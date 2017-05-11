from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'cython functions',
  ext_modules = cythonize("cextensions/cythonfun.pyx"),
)

