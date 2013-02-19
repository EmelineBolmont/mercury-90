#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for fortran90

from make import *
import sys # to be able to retrieve arguments of the script


###############################################
## Beginning of the program
###############################################

filename = "test_disk.f90"
force = False # To force the compilation of every module

isProblem = False
problem_message = "Given in parameters is the name of the source code you want to compile." + "\n" + \
"The default filename is : "+filename + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * force : To force the compilation of every module even those not modified" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * name : (default='%s')The filename of the sourcecode you want to compile" % filename + "\n" + \
"Example :" + "\n" + \
"maketest.py resultant_torque.f90"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'help'):
    isProblem = True
  elif (key == 'force'):
    force = True
  elif (key == 'name'):
    filename = value
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

# We clean undesirable files. Indeed, we will compile everything everytime.
if force:
  clean(["o", "mod"])

sourceFile.setCompilator("gfortran")

sourceFile.setCompilingOptions("")

sources_filename = lister("*.f90")

# We create the binaries
make_binaries(sources_filename, [filename], debug=True, gdb=False, profiling=False)
