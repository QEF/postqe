======
postqe
======

:py:mod:`postqe` is a Python package to perform postprocessing calculations for results
obtained with the Quantum Espresso (QE) code [#QE]_. The package provides Python functions
to post-process the results, such as plotting the charge density (or the bare/Hartree/total
potentials) on 1D or 2D sections, fitting the total energy with Murnaghan Equation of State
(EOS), etc. The package also exposes some QE functionalities in Python using the F2PY code
[#F2PY]_ and wrappers to generate Python modules from QE dynamically linked libraries.

It is meant to be imported in your own Python code or used from the command line (see the
Tutorial part of this documentation). It is also meant for people who want to tinker with
the code and adapt it to their own needs. The package is based on numpy, scipy and matplotlib libraries.


Current features of the package include:

* Fit the total energy :math:`E_{tot}(V)` with Murnaghan's equation of state
* Calculate and plot the electronic band structure
* Calculate and plot the electronic density of states (DOS)
* Plot 1D or 2D sections of the charge density
* Plot 1D or 2D sections of different potentials (Hartree, exchange-correlation, etc.)


.. [#QE] http://www.quantum-espresso.org/
.. [#F2PY]  https://docs.scipy.org/doc/numpy-dev/f2py/


Installation
------------

You can download all package files from GitHub  and then install it with the command:

.. code-block:: bash

   sudo python setup.py install


Usage
-----

The most useful functions for the user are directly accessible. You can import all of them as:

.. code-block:: python

   from postqe import *

or you can import only the ones you need. The above command also makes available a number of
useful constants that you can use for unit conversions.

More functions are available as submodules. See the related documentation for more details.
Note, however, that most of these functions are less well documented and are meant for advanced
users or if you want to tinker with the code.


Status
------

Development(Alpha)
