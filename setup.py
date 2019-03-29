#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

import os
import glob
import platform
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install


def find_pyqe_module():
    """
    Returns the absolute pathname of the pyqe module build for the running platform.

    :return: A pathname string or `None` if no suitable module is found.
    """
    project_dir = os.path.dirname(__file__)
    python_version = ''.join(platform.python_version_tuple()[:2])
    platform.python_implementation()
    pyqe_modules = glob.glob(os.path.join(project_dir, 'postqe/pyqe.*.so'))

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
        print("Build f2py extension module ...")
        os.system('make -C postqe/fortran all')
        build_ext.run(self)


class InstallCommand(install):

    def run(self):
        if find_pyqe_module() is None:
            print("A suitable pyqe module not found, invoke build_ext ...")
            self.run_command('build_ext')
        install.run(self)


setup(
    name='postqe',
    version='0.3',
    packages=['postqe', 'postqe/ase'],
    package_data={'postqe': ['schemas/*.xsd', 'pyqe.*.so']},
    install_requires=[
        'numpy>=1.10.1', 'ase>=3.10', 'scipy', 'h5py', 'matplotlib',
        'xmlschema>=0.9.10', 'colormath', 'natsort', 'moviepy'
    ],
    data_files=[
        # ('/usr/share/doc/postqe/example1', glob.glob('examples/example1/*')),
        # ('/usr/share/doc/postqe/example2', glob.glob('examples/example2/*')),
        # ('/usr/share/doc/postqe/example3', glob.glob('examples/example3/*')),
        #('/usr/share/doc/postqe/RGB', [fn for fn in glob.glob('examples/RGB/*') if os.path.isfile(fn)]),
        #('/usr/share/doc/postqe/RGB/EIG', glob.glob('examples/RGB/EIG/*')),
        #('/usr/share/doc/postqe/RGB/plot', glob.glob('examples/RGB/plot/*')),
        #('/usr/share/doc/postqe/RGB/spectra', glob.glob('examples/RGB/spectra/*')),
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
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Fortran',
        'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Physics'
    ]
)
