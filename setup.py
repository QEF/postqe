#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

import os
import re
import pathlib
import platform

from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from distutils.command.clean import clean  # type: ignore[attribute]

VERSION_NUMBER_PATTERN = re.compile(r"version_number\s*=\s*(\'[^\']*\'|\"[^\"]*\")")

QE_SOURCE_URL = "https://github.com/QEF/q-e/archive/refs/tags/qe-6.8.zip"
QE_SOURCE_MD5SUM = "787c4aad3b203f7dd61d03a849a1c4e9"


def find_pyqe_module():
    """
    Returns the absolute pathname of the pyqe module built for the running platform.

    :return: A pathname string or `None` if no suitable module is found.
    """
    python_version = ''.join(platform.python_version_tuple()[:2])

    for filename in map(str, pathlib.Path(__file__).parent.glob('postqe/_pyqe.*.so')):
        if f'-{python_version}-' not in filename and f'{python_version}m' not in filename:
            continue
        elif platform.python_implementation().lower() not in filename:
            continue
        elif platform.system().lower() not in filename:
            continue
        elif platform.machine() not in filename:
            continue
        return filename


class BuildExtCommand(build_ext):

    def run(self):
        build_dir = pathlib.Path(__file__).parent.joinpath('postqe/fortran')
        qe_topdir = os.environ.get('QE_TOPDIR') or next(build_dir.rglob('q-e-*'), None)

        if qe_topdir is not None:
            qe_topdir = pathlib.Path(qe_topdir)
        else:
            import zipfile

            qe_source_path = str(build_dir.joinpath('q-e.zip'))
            if not os.path.isfile(qe_source_path):
                print("Download QE source code ...")
                os.system(f'curl -SL {QE_SOURCE_URL} -o "{qe_source_path}"')

            print("Checks MD5SUM of source archive ... ")
            if os.system(f'echo `md5sum {qe_source_path}` | grep --quiet "^{QE_SOURCE_MD5SUM}[[:blank:]]"'):
                raise ValueError("Checksum of source code archive doesn't match!")

            print("Extract source files from archive ...")
            with zipfile.ZipFile(qe_source_path, 'r') as zfp:
                zfp.extractall(build_dir)

            qe_topdir = pathlib.Path(next(build_dir.rglob('q-e-*')))

        print("QE top directory {}".format(str(qe_topdir)))

        if 'QE_TOPDIR' not in os.environ:
            os.environ['QE_TOPDIR'] = str(qe_topdir)

        # Check QE installation
        for version_file in qe_topdir.glob('include/*version.h'):
            with version_file.open() as fp:
                version_number = VERSION_NUMBER_PATTERN.search(fp.read())

            if version_number is None:
                raise ValueError("Missing version label in {}".format(str(version_file)))

            print("QE {}".format(version_number.group(0)))
            if version_number.group(1).strip('\'"') < '6.8':
                raise ValueError("Build aborted: QE version too old!")

            break
        else:
            raise FileNotFoundError("Missing QE version file info!")

        if not qe_topdir.joinpath('make.inc').is_file():
            print("Configure Quantum Espresso ...")
            for configure_file in qe_topdir.rglob('**/configure'):
                os.chmod(configure_file, 0o755)
            os.system(str(qe_topdir.joinpath('configure')))

            with qe_topdir.joinpath('make.inc').open() as fp:
                make_inc_lines = fp.readlines()

            changed = False
            for k in range(len(make_inc_lines)):
                line = make_inc_lines[k]
                if line.startswith("CFLAGS ") or line.startswith("FFLAGS "):
                    if '-fPIC' not in line:
                        make_inc_lines[k] = line.replace(' = ', ' = -fPIC ')
                        changed = True

            if changed:
                with qe_topdir.joinpath('make.inc').open(mode='w') as fp:
                    fp.writelines(make_inc_lines)

        print("Build pyqe module ...")
        os.system('make -C {} all'.format(str(build_dir)))

        # Modify python wrapper module, fixing import in postqe package.
        with build_dir.parent.joinpath('pyqe.py').open() as fp:
            python_wrapper_lines = fp.readlines()

        with build_dir.parent.joinpath('pyqe.py').open(mode='w') as fp:
            python_wrapper_lines[0] = "# Altered wrapper for postqe\n"
            python_wrapper_lines[1] = "from . import _pyqe\n"
            fp.writelines(python_wrapper_lines)

        build_ext.run(self)


class CleanCommand(clean):

    def run(self):
        build_dir = pathlib.Path(__file__).parent.absolute().joinpath('postqe/fortran')
        qe_topdir = os.environ.get('QE_TOPDIR') or next(build_dir.rglob('q-e-*'), None)
        if qe_topdir is None:
            qe_topdir = pathlib.Path(next(build_dir.rglob('q-e-*')))
        else:
            qe_topdir = pathlib.Path(qe_topdir)

        if 'QE_TOPDIR' not in os.environ:
            os.environ['QE_TOPDIR'] = str(qe_topdir.absolute())

        print("Clean QE and pyqe module building files ...")
        os.system('make -C {} clean'.format(str(build_dir.absolute())))
        super().run()


class InstallCommand(install):

    def run(self):
        if find_pyqe_module() is None:
            print("A suitable pyqe module not found, invoke build_ext ...")
            self.run_command('build_ext')
        install.run(self)


setup(
    name='postqe',
    version='1.0.0',
    packages=['postqe'],
    package_data={'postqe': ['_pyqe.*.so']},
    install_requires=[
        'numpy>=1.17.0', 'ase~=3.20.0', 'qeschema~=1.1', 'scipy',
        'h5py', 'matplotlib', 'colormath', 'natsort', 'moviepy',
    ],
    entry_points={
        'console_scripts': [
            'postqe=postqe.cli:main'
        ]
    },
    cmdclass={
        'build_ext': BuildExtCommand,
        'clean': CleanCommand,
        'install': InstallCommand,
    },
    author='Mauro Palumbo, Pietro Delugas, Davide Brunato',
    author_email='pdelugas@sissa.it, brunato@sissa.it',
    license='LGPL-2.1',
    long_description='Post processing tools for Quantum Espresso',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Physics'
    ]
)
