.. _introduction:


****************
Introduction
****************

:py:mod:`postqe` is a Python package for postprocessing of results obtained with the Quantum Espresso (*QE*) code [#QE]_. The package provides Python API functions for the most common tasks, such as plotting the charge density or fitting the total energy with an equation of state (EOS). It also makes available in Python some QE functionalities using the F2PY code [#F2PY]_ and wrappers to generate Python modules from QE dynamically linked libraries. Finally, it also includes an interface with the popular Atomic Simulation Environment (*ASE*) [#ASE]_, from which in fact leverages for some functionalities.

The package implents some Python classes for handling the most important quantities, such as the charge or a potential. Some classes are derived from *ASE* and adapted to the postprocessing needs of QE.

It is meant to be imported into your own Python scripts or used from the command line interface (see the Tutorial for some examples). It is also meant for people who want to tinker with the code and adapt it to their own needs. The package is based on *numpy*, *scipy*, *matplotlib* and *ASE* libraries.

Finally it is a software framework where Quantum Espresso developers or advanced users may implement new functionalities which are needed by the community. In this respect, it offers the possibility to develop code in different languages (Python, C/C++, Fortran) and then use Python to "glue" everything together.

Current features of the package include: 

* Fit the total energy :math:`E_{tot}(V)` with an equation of state (Murnaghan, Vinet, Birch, etc.)
* Calculate and plot the electronic band structure
* Calculate and plot the electronic density of states (DOS)
* Plot 1D, 2D or 3D sections of the charge density 
* Plot 1D, 2D or 3D sections of different potentials (Hartree, exchange-correlation, etc.)

.. [#QE] http://www.quantum-espresso.org/
.. [#F2PY]  https://docs.scipy.org/doc/numpy-dev/f2py/
.. [#ASE] https://wiki.fysik.dtu.dk/ase/


================
Installation
================

You can download all package files from GitHub  and then install it with the command:

.. code-block:: bash 

   sudo python setup.py install


The most useful functions for the user are directly accessible. You can import all of them as:

.. code-block:: python 

   from postqe import *

or you can import only the ones you need. The above command also makes available a number of useful constants that you can use for unit conversions.

More functions are available as submodules. See the related documentation for more details. Note, however, that most of these functions are less well documented and are meant for advanced users or if you want to tinker with the code.

================
General notes
================

----------------------------
Version
----------------------------

The package is still under development. Hence, it may contain bugs and some features may not be fully implemented. Use at your own risk.

----------------------------
Plotting
----------------------------

:py:mod:`postqe` uses the *matplotlib* library for Plotting. Some functions in the package are simply useful wrappers for *matplotlib* functionalities of common uses. They return a *matplotlib* object which can be further adapted to specific needs and personal taste. Alternatively, you can of course manipulate and plot the post-procecessed data with any other Python tool of your choice.

It is also possible to export the charge (and various potentials) into text files according to different available formats (XSF XCrySDen format, cube Gaussian format, Gnuplot formats, contour.x and plotrho.x formats).  
