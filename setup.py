import os
import glob
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.sdist import sdist
from setuptools.command.install import install


class MyBuildExt(build_ext):

    def run(self):
        os.system('make -C postqe/fortranmodules2/flib all')
        os.system('make -C postqe/fortranmodules2 all')
        build_ext.run(self)


class MySDist(sdist):

    def run(self):
        os.system('make -C postqe/fortranmodules2/flib all')
        os.system('make -C postqe/fortranmodules2 all')
        sdist.run(self)


class MyInstall(install):

    def run(self):
        os.system('make -C postqe/fortranmodules2/flib all')
        os.system('make -C postqe/fortranmodules2 all')
        install.run(self)

setup(
    name='postqe',
    version='0.1',
    packages=['postqe'],
    package_data={'postqe': ['schemas/*.xsd']},
    install_requires=[
        'xmlschema', 'numpy', 'ase', 'scipy', 'h5py', 'colormath', 'natsort', 'moviepy', 'matplotlib'
    ],
    author='Mauro Palumbo',
    author_email='mpalumbo@sissa.it',

    # summary = 'Just another Python package for the cheese shop',
    license='MIT',
    long_description='Post processing tools for Quantum Espresso',
    data_files=[
        #('/usr/share/doc/postqe/example1', glob.glob('examples/example1/*')),
        #('/usr/share/doc/postqe/example2', glob.glob('examples/example2/*')),
        #('/usr/share/doc/postqe/example3', glob.glob('examples/example3/*')),
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
    requires=['python (>=3.3)'],
    cmdclass={
        'build_ext': MyBuildExt,
        'sdist': MySDist,
        'install': MyInstall
    },

    # could also include long_description, download_url, classifiers, etc.
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Fortran',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Physics'
    ]
)
