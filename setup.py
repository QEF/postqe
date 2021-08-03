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

    LIB_DIRS = ['UtilXlib']  #, 'Modules'  #'upflib', 'LAXlib']
    EXCLUDED_FILES = {'mp_base_gpu.f90', 'mp_bands_util.f90'}

    def run(self):
        qe_topdir = os.environ.get('QE_TOPDIR')
        if qe_topdir is None:
            # TODO: discovery in certain paths?
            raise KeyError("QE_TOPDIR environment variable not found")

        qe_topdir = pathlib.Path(qe_topdir)

        """
        pw_path = shutil.which('pw.x', path=os.environ.get('QE_BIN_DIR'))

        if pw_path is None:
            if 'QE_BIN_DIR' in os.environ:
                msg = "Build aborted: pw.x executable not found in QE_BIN_DIR directory path!"
            else:
                msg = "Build aborted: pw.x executable not found in executable search paths!"
            raise FileNotFoundError(msg)

        cmd = f'readelf -Wwi {pw_path} | grep DW_AT_comp_dir'
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)

        dir_paths = {line.decode('utf-8').split('): ')[-1] for line in output.splitlines()}
        if not dir_paths:
            raise FileNotFoundError("Missing QE build directory!")

        # for d in dir_paths:
        #    print(d)

        # qe_build_dir = pathlib.Path(os.path.commonpath(dir_paths))
        """
        print("QE top directory {}".format(str(qe_topdir)))

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

        # invocare il configure aggiungendo l'opzione -fPIC
        # ./configure FFLAGS = "-O3 -g -fPIC -fallow-argument-mismatch"

        # ci pensa lui a fare la verifica.

        fortran_files = []
        for path in map(lambda x: qe_topdir.joinpath(x), self.LIB_DIRS):
            if not path.is_dir():
                raise FileNotFoundError("Missing directory {}".format(str(path)))

            fortran_files.append(path.glob('mp_base.f90'))

        fortran_files = [str(x) for x in itertools.chain.from_iterable(fortran_files)
                         if x.name not in self.EXCLUDED_FILES]

        if not fortran_files:
            raise FileNotFoundError("Missing QE Fortran files!")
        print("Found {} Fortran source files ...".format(len(fortran_files)))

        print("Remove old wrappers ...")
        for filepath in glob.iglob('postqe/wrappers/f90wrap_*.f90'):
            os.unlink(filepath)

        print("Create wrappers using f90wrap ...\n")
        os.chdir('postqe/wrappers')
        os.system('f90wrap -k kind_map -m wrapper_module {}'.format(' '.join(fortran_files)))

        breakpoint()
        print("Compile wrappers with f2py ...\n")
        os.system('f2py-f90wrap -c -m _wrapper_module f90wrap_*.f90 -L')
        return

        print("Build f2py extension module ...")
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
