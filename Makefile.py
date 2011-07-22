#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for fortran90

from make import *

# We clean undesirable files. Indeed, we will compile everything everytime.
clean(["o", "mod"])

sources_filename = lister("*.f90")

# We create the binaries
make_binaries(sources_filename, ["mercury.f90", "element.f90", "close.f90"], debug=True, gdb=False, profiling=False)
