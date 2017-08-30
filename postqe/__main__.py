#!/usr/bin/env python
#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
"""
Execute postqe module as a script (see PEP-338).
"""
if not locals()['__package__']:
    # When this module is loaded before the package then __package__ is None or ''
    # and the relative imports are disabled. In this case import the package and
    # set __package__.
    #
    # $ python postqe --> __package__ == ''
    # $ python postqe/__main__.py --> __package__ is None
    #
    # Ref: https://www.python.org/dev/peps/pep-0366/ for details.
    import os
    import sys
    pkg_search_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if sys.path[0] != pkg_search_path:
        sys.path.insert(0, pkg_search_path)
    import postqe
    __package__ = postqe.__name__

from .cli import main
main()
