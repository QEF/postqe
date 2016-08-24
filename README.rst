=========
postqe
=========

This is a first experimental implementation of post processing tools for QE
developed in Python 3 (tested with Python 3.4.3).

A first version of a GUI has also been created as "postqeGUI.py". It is based on Python version of wxWidgets, i.e. WxPython. Since the postqe code is written in Python 3, we are using here the Phoenix version of WxPython. Please note that the porting of WxPython to Python 3 is still ongoing. More details on https://wiki.wxpython.org/ProjectPhoenix.

Features
--------
- Read charge and wavefunctions files from Quantum Espresso output (both binary and HDF5). 
- Print the charge data as a text file for further postprocessing
- Calculate the bare potential of the atoms.
- Calculate the bare plus the Hartree potential
- Calculate the bare + Hartree + exchange-correlation potential 
- 1D and 2D plots of the charge along certain directions/planes


