#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for fortran90

from make import *
import sys # to be able to retrieve arguments of the script


###############################################
## Beginning of the program
###############################################

LOG_NAME = 'compilation.log'

filename = "test_disk.f90"
force = False # To force the compilation of every module
isOptimize = False

isProblem = False
problem_message = "Given in parameters is the name of the source code you want to compile." + "\n" + \
"The default filename is : "+filename + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * force : To force the compilation of every module even those not modified" + "\n" + \
" * no_debug : will skip all debug options and use optimisation options instead" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * name : (default='%s')The filename of the sourcecode you want to compile" % filename + "\n" + \
"Example :" + "\n" + \
"maketest.py name=resultant_torque.f90"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'help'):
    isProblem = True
  elif (key == 'no_debug'):
    isOptimize = True
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

sources_filename = lister("*.f90")

# Before compiling, we delete the previous compilation log. Indeed, we need to append the several warnings in the same file
# But we do not want to have infos of the previous compilation in it.
if os.path.isfile(LOG_NAME):
  os.remove(LOG_NAME)

if isOptimize:
  sourceFile.setCompilingOptions("-O3 -march=native -pipe -finit-real=nan")
  # We create the binaries
  make_binaries(sources_filename, [filename], debug=False, gdb=False, profiling=False)
else:  
  sourceFile.setCompilingOptions("")
  
  # We create the binaries
  make_binaries(sources_filename, [filename], debug=True, gdb=False, profiling=False)
