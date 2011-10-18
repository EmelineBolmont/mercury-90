#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for fortran90

from make import *

# We clean undesirable files. Indeed, we will compile everything everytime.
clean(["o", "mod"])

#~ sourceFile.setCompilator("ifort")
#~ sourceFile.setCompilingOptions("-check all")

sourceFile.setCompilator("gfortran")
sourceFile.setCompilingOptions("-O3 -march=native -pipe -finit-real=snan")
#~ sourceFile.setCompilingOptions("")

# pour tester les bornes des tableaux : -fbounds-check (il faut ensuite faire tourner le programme, des tests sont effectués au cours de l'exécution)

sources_filename = lister("*.f90")

# We create the binaries
make_binaries(sources_filename, ["mercury.f90", "element.f90", "close.f90"], debug=False, gdb=False, profiling=False)

