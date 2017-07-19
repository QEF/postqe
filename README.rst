======
postqe
======

This is a first experimental implementation of post processing tools for QE developed in for Python 2.7
and Python 3.3+. It relies on numpy, which must be available in your platform.
Some functionalities, for example RGB calculations, need additional Python modules (moviepy, colormath,
natsort) which can be easily downloaded and installed from the web or with "pip".

Some examples are provided in the directory examples, running the command line version of the code.
RGA factors are not yet available in the GUI.
 

Features
--------
- Read charge and wavefunctions files from Quantum Espresso output (both binary and HDF5). 
- Print the charge data as a text file for further postprocessing
- Calculate the bare potential of the atoms.
- Calculate the bare plus the Hartree potential
- Calculate the bare + Hartree + exchange-correlation potential 
- Print the potentials as a text file for further postprocessing
- 1D and 2D plots of the charge along certain directions/planes using Fourier interpolation
- Calculate RGB factors


