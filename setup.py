from setuptools import setup, find_packages

setup (
    name='postqe',
    version='0.1',
    packages=find_packages(),

    # Declare your packages' dependencies here, for eg:
    install_requires=[
        'numpy', 'scipy', 'h5py', 'colormath', 'natsort', 'moviepy', 'wx'
    ],

    # Fill in these to make your Egg ready for upload to
    # PyPI
    author='Mauro Palumbo',
    author_email='mpalumbo@sissa.it',

    #summary = 'Just another Python package for the cheese shop',
    url='',
    license='MIT',
    long_description='Post processing tools for Quantum Espresso',

    # could also include long_description, download_url, classifiers, etc.
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Fortran',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Physics'
    ]
)
