from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'cython plot',
  ext_modules = cythonize("cythonplot.pyx"),
)

