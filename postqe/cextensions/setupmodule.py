#!/usr/bin/env python3

from distutils.core import setup, Extension

gvect = Extension('gvect', sources = ['gvect.c'])

setup(name = 'gvect',
      version = '1.0',
      description = 'Python Package for g vectors generation C Extension',
      ext_modules = [gvect],
      author='Mauro Palumbo',
      author_email=''
     )
