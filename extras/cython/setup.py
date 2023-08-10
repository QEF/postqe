from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = 'cython functions',
    ext_modules = cythonize("compute_vs_cython.pyx"),
)

