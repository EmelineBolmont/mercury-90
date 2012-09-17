#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for fortran90

from make import *
import sys # to be able to retrieve arguments of the script


###############################################
## Beginning of the program
###############################################

filename = "test_disk.f90"

isProblem = False
problem_message = "Given in parameters is the name of the source code you want to compile." + "\n" + \
"The default filename is : "+filename + "\n" + \
"For instance :" + "\n" + \
"maketest.py resultant_torque.f90"

# We get arguments from the script
for arg in sys.argv[1:]:
	if (arg == 'help'):
		isProblem = True
	else:
		filename = arg

if isProblem:
	print(problem_message)
	exit()

# We clean undesirable files. Indeed, we will compile everything everytime.
clean(["o", "mod"])

sourceFile.setCompilator("gfortran")

sourceFile.setCompilingOptions("-O0 -finit-real=nan")

sources_filename = lister("*.f90")

# We create the binaries
make_binaries(sources_filename, [filename], debug=True, gdb=False, profiling=False)
