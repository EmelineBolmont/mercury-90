#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# Of course, everything is not tested, but it is planed to test as many things as possible

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "07 juin 2011"
__version__ = "$Revision: 2.0 $"
__credits__ = """We run a test simulation and erase all the files created after the tests. The simulations files are thought to be 
in a "simu_test" subdirectory of the directory were are the sources (and binaries) of mercury (and this script)"""

from make import clean
import os
import subprocess # To launch various process, get outputs et errors, returnCode and so on.
import difflib # To compare two strings
import pdb # To debug

# We import all variables linked to original output files.
from original_files import *

def run(command):
    """run a command that will be a string.
    The function return a tuple with the output, 
    the stderr and the return code. 
    It is fitted to run only command placed in the parent directory (of type "../command")"""

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
    (process_stdout, process_stderr) = process.communicate()
    returnCode = process.poll()
    
    if (returnCode == 0):
        print(command+" was runned successfully")
        
        # We delete the last \n introduced by subprocess and that currupt the diff managed to see differences between outputs.
        return (process_stdout[0:-1], process_stderr[0:-1])
    else:
        print(command+" has terminated with an error code: "+str(returnCode))
        
        temp = command.lstrip("../")
        
        errlogname = temp+".err"
        print("Writing the error of the compilation in '"+errlogname+"'")
        f = open(errlogname,'w')
        f.write(process_stderr)
        f.close()
        
        outlogname = temp+".out"
        print("Writing the output of the compilation in '"+outlogname+"'")
        f = open(outlogname,'w')
        f.write(process_stdout)
        f.close()
        
        return 0

def compare(original, new):
    """function that compare and print differences between to strings that are compared line by line."""
    
    modified = ["+ ", "- ", "? "]
    
    # We want lists, because difflib.Differ only accept list of lines.
    original = original.split("\n")
    new = new.split("\n")
    
    d = difflib.Differ()

    result = list(d.compare(original, new))
    
    differences = []
    for (number, line) in enumerate(result):
      if (line[0:2] in modified):
        differences.append("("+line[0]+") l"+str(number+1)+" :"+line[2:])

    # We print separately because it is more convenient if we want to store in a file instead.
    if (differences != []):
      for line in differences:
        print(line)
    else:
      print "OK : No differences"
    return 0

def compare2file(original,new_file):
  """Function that will use compare to see differences between 'original' 
  that is thought to be a variable and 'new_file' that is the name of a 
  file to read then use as input
  """
  
  f = open(new_file, 'r')
  new_lines = f.readlines()
  f.close()
  
  compare(original, ''.join(new_lines))
  
  return 0

FOLDER = "simu_test"

##################
# Outputs of various binaries and tests to compare with the actual ones. 
# Theses outputs are those of the original version of mercury, that is, mercury6_2.for
##################



os.chdir(FOLDER)
# We clean undesirable files. 
clean(["out", "clo", "aei", "dmp", "tmp"])

print("""###################
# Test of mercury #
###################""")

(process_stdout, process_stderr) = run("../mercury")

print("for the Output of mercury")
compare(MERCURY_OUTPUT_ORIGINAL, process_stdout)

for (new_file, oldfile_output) in zip(MERCURY_NAMES, MERCURY_FILES):
  print("\nfor "+new_file+" :")
  compare2file(oldfile_output, new_file)
  
print("""#################
# Test of close #
#################""")

(process_stdout, process_stderr) = run("../close")
for (new_file, oldfile_output) in zip(CLOSE_NAMES, CLOSE_FILES):
  print("\nfor "+new_file+" :")
  compare2file(oldfile_output, new_file)

print("""###################
# Test of element #
###################""")

(process_stdout, process_stderr) = run("../element")
for (new_file, oldfile_output) in zip(ELEMENT_NAMES, ELEMENT_FILES):
  print("\nfor "+new_file+" :")
  compare2file(oldfile_output, new_file)





