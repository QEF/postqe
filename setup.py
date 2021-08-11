#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

import os
import re
import glob
import pathlib
import platform
import shutil
import subprocess
import itertools

from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install


VERSION_NUMBER_PATTERN = re.compile(r"version_number\s*=\s*(\'[^\']*\'|\"[^\"]*\")")


def find_qe_installation():
    """
    Find a Quantum Espresso installation to be linked to postqe. The best match
    installation found is used for building Fortran 90 wrappers.
    Priority is for subdir installations and then for the most recent version.

    :return: QE installation dynamic library path or `None` if no installation is found.
    """
    installations_found = []

    pw_path = shutil.which('pw.x', path=os.environ.get('QE_BIN_DIR'))
    if pw_path is not None:
        cmd = f'readelf -Wwi {pw_path} | grep DW_AT_comp_dir'
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
        dir_paths = {line.decode('utf-8').split('): ')[-1] for line in output.splitlines()}
        print(dir_paths)

        for path in map(lambda x: pathlib.Path(x), dir_paths):
            print(list(path.glob('*.f90')))


def find_pyqe_module():
    """
    Returns the absolute pathname of the pyqe module build for the running platform.

    :return: A pathname string or `None` if no suitable module is found.
    """
    project_dir = os.path.dirname(__file__)
    python_version = ''.join(platform.python_version_tuple()[:2])
    pyqe_modules = glob.iglob(os.path.join(project_dir, 'postqe/pyqe.*.so'))

    for filename in pyqe_modules:
        if '{}m'.format(python_version) not in filename:
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
        qe_topdir = os.environ.get('QE_TOPDIR')
        if qe_topdir is None:
            # TODO: discovery in certain paths?
            raise KeyError("QE_TOPDIR environment variable not found")

        qe_topdir = pathlib.Path(qe_topdir)
        print("QE top directory {}".format(str(qe_topdir)))

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
        os.system('make -C postqe/fortran all')

        build_ext.run(self)


class InstallCommand(install):

    def run(self):
        if find_pyqe_module() is None:
            print("A suitable pyqe module not found, invoke build_ext ...")
            self.run_command('build_ext')
        return
        install.run(self)


setup(
    name='postqe',
    version='1.0.0',
    packages=['postqe'],
    package_data={'postqe': ['pyqe.*.so']},
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
