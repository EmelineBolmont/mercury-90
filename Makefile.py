#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for fortran90

from make import *

import git_infos

# Parameters
debug = False
gdb = False
profiling = False

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * debug : [%s] activate debug options" % debug + "\n" + \
" * gdb : [%s] activate options for gdb" % gdb + "\n" + \
" * profiling : [%s] activate options for profiling" % profiling + "\n" + \
" Example : " + "\n" + \
" Makefile.py gdb"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'debug'):
    debug = True
  elif (key == 'gdb'):
    gdb = True
  elif (key == 'profiling'):
    profiling = True
  elif (key == 'help'):
    isProblem = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

git_infos.write_infos_in_f90_file()

isModifs = git_infos.is_non_committed_modifs()

# We clean undesirable files. Indeed, we will compile everything everytime.
clean(["o", "mod"])

#~ sourceFile.setCompilator("ifort")
#~ sourceFile.setCompilingOptions("-check all")

sourceFile.setCompilator("gfortran")
sourceFile.setCompilingOptions("-O3 -march=native -pipe -finit-real=nan")
#~ sourceFile.setCompilingOptions("")

# pour tester les bornes des tableaux : -fbounds-check (il faut ensuite faire tourner le programme, des tests sont effectués au cours de l'exécution)

sources_filename = lister("*.f90")

# We create the binaries
make_binaries(sources_filename, ["mercury.f90", "element.f90", "close.f90"], debug=debug, gdb=gdb, profiling=profiling)
#~ make_binaries(sources_filename, {"mercury.f90":"mercury2", "element.f90":"element2", "close.f90":"close2"}, debug=False, gdb=False, profiling=False)

if (isModifs):
	print("Warning: There is non committed modifs!")