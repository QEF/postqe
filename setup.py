#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

import os
import glob
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.sdist import sdist
from setuptools.command.install import install


class MyBuildExt(build_ext):

    def run(self):
        os.system('make -C postqe/fortran all')
        build_ext.run(self)


class MySDist(sdist):

    def run(self):
        sdist.run(self)


class MyInstall(install):

    def run(self):
        install.run(self)

setup(
    name='postqe',
    version='0.2',
    packages=['postqe', 'postqe.ase'],
    package_data={'postqe': [
        'schemas/*.xsd'
    ]},
    install_requires=[
        'numpy>=1.10.1', 'ase>=3.10', 'scipy', 'h5py', 'matplotlib',
        'xmlschema>=0.9.10', 'colormath', 'natsort', 'moviepy'
    ],
    data_files=[
        # ('/usr/share/doc/postqe/example1', glob.glob('examples/example1/*')),
        # ('/usr/share/doc/postqe/example2', glob.glob('examples/example2/*')),
        # ('/usr/share/doc/postqe/example3', glob.glob('examples/example3/*')),
        ('/usr/share/doc/postqe/RGB', [fn for fn in glob.glob('examples/RGB/*') if os.path.isfile(fn)]),
        ('/usr/share/doc/postqe/RGB/EIG', glob.glob('examples/RGB/EIG/*')),
        ('/usr/share/doc/postqe/RGB/plot', glob.glob('examples/RGB/plot/*')),
        ('/usr/share/doc/postqe/RGB/spectra', glob.glob('examples/RGB/spectra/*')),
    ],
    entry_points={
        'console_scripts': [
            'postqe=postqe.cli:main'
        ]
    },
    cmdclass={
        'build_ext': MyBuildExt,
        'sdist': MySDist,
        'install': MyInstall
    },

    author='Mauro Palumbo',
    author_email='mpalumbo@sissa.it',
    license='LGPL-2.1',
    long_description='Post processing tools for Quantum Espresso',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Fortran',
        'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Physics'
    ]
)
