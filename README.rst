======
postqe
======

*Postqe* is a Python package to perform postprocessing calculations for results
obtained with the Quantum Espresso (QE) code [#QE]_. The package provides Python functions
to post-process the results, such as plotting the charge density (or the bare/Hartree/total
potentials) on 1D or 2D sections, fitting the total energy with Murnaghan Equation of State
(EOS), etc. The package also exposes some QE functionalities in Python using the F2PY code
[#F2PY]_ and wrappers to generate Python modules from QE dynamically linked libraries.

It is meant to be imported in your own Python code or used from the command line (see the
Tutorial part of this documentation). It is also meant for people who want to tinker with
the code and adapt it to their own needs. The package is based on numpy, scipy and
matplotlib libraries.


Current features of the package include:

* Fit the total energy :math:`E_{tot}(V)` with Murnaghan's equation of state
* Calculate and plot the electronic band structure
* Calculate and plot the electronic density of states (DOS)
* Plot 1D or 2D sections of the charge density
* Plot 1D or 2D sections of different potentials (Hartree, exchange-correlation, etc.)


.. [#QE] http://www.quantum-espresso.org/
.. [#F2PY]  https://docs.scipy.org/doc/numpy-dev/f2py/


Package requirements
--------------------
There are some non-Python packages required to be installed in your system before installing *postqe*:

Fortran compiler:
    Maybe *gfortran* on a Gnu Linux platform for free or another commercial package, like PGI or Intel Fortran.

Lapack development libraries:
    the basename of this package should be *lapack-devel* on a RHEL-based system or *liblapack-dev*
    on a Debian-based system).


Installation
------------

Download all package's files from GitHub (unpack the downloaded archive
or clone the repository).
If you are in a virtual environment you can install it with the command:

.. code-block:: bash

   pip install -r requirements.txt
   python setup.py install

Otherwise avoid system wide installations with root privileges but install it in user-space
with the command:

.. code-block:: bash

   pip install -r requirements.txt --user
   python setup.py install --user

The setup script takes care to download, configure and compile a fresh installation of
Quantum ESPRESSO v6.8 into the setup's directory `build/`.

Alternatively you can download Quantum ESPRESSO v6.8 sources by your own and then unpack
the archive in another directory. Then provide the variable `QE_TOPDIR` to the install
command, setting it with the path of the directory containing the Quantum ESPRESSO sources,
for example:

.. code-block:: bash

   QE_TOPDIR=<path-to-QE-sources-base-dir> python setup.py install

The advantage of the alternate method is that you don't need to recompile Quantum ESPRESSO
after a clean of the build directories (`python setup.py clean`).


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


Authors
-------
Mauro Palumbo
Davide Brunato
Pietro Delugas


License
-------
This software is distributed under the terms of the LGPL-2.1 license. See
the file 'LICENSE' in the root directory of the present distribution, or
https://opensource.org/licenses/LGPL-2.1 .

