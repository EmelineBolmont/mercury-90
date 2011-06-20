#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# Of course, everything is not tested, but it is planed to test as many things as possible

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "09 juin 2011"
__version__ = "$Revision: 2.12 $"
__credits__ = """We run a test simulation and erase all the files created after the tests. The simulations files are thought to be 
in a "simu_test" subdirectory of the directory were are the sources (and binaries) of mercury (and this script)"""

from make import clean
import os
import subprocess # To launch various process, get outputs et errors, returnCode and so on.
import difflib # To compare two strings
import pdb # To debug
from time import time

# Here are all the name of the output files, for each binary
MERCURY_FILENAMES = ["info.out", "big.dmp", "small.dmp", "param.dmp", "restart.dmp", "big.tmp", "small.tmp", "param.tmp", "restart.tmp"]
ELEMENT_FILENAMES = ["APOLLO.aei", "JUPITER.aei", "MERCURY.aei", "ORPHEUS.aei", "TOUTATIS.aei", "EARTHMOO.aei", "KHUFU.aei", "MINOS.aei", 
  "PLUTO.aei", "URANUS.aei", "JASON.aei", "MARS.aei", "NEPTUNE.aei", "SATURN.aei", "VENUS.aei"]
CLOSE_FILENAMES = ["APOLLO.clo", "JUPITER.clo", "MERCURY.clo", "ORPHEUS.clo", "TOUTATIS.clo", "EARTHMOO.clo", "KHUFU.clo", "MINOS.clo", 
  "PLUTO.clo", "URANUS.clo", "JASON.clo", "MARS.clo", "NEPTUNE.clo", "SATURN.clo", "VENUS.clo"]

EXTENTION_ORIGINAL = ".ori"

def run(command):
  """run a command that will be a string.
  The function return a tuple with the output, 
  the stderr and the return code. 
  It is fitted to run only command placed in the parent directory (of type "../command")"""

  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  
  (process_stdout, process_stderr) = process.communicate()
  returnCode = process.poll()
  
  if (returnCode == 0):
    print(command+" was runned successfully\n")
    
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
    
    return (process_stdout, process_stderr)

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
    return "\n".join(differences)
  else:
    return None

def compare2file(ori_files,new_files):
  """Function that will use compare to see differences between 'original' 
  that is thought to be a variable and 'new_file' that is the name of a 
  file to read then use as input
  """
  no_diff = []
  diff = []
  
  for (original, new) in zip(ori_files, new_files):
    f_old = open(original, 'r')
    old_lines = f_old.readlines()
    f_old.close()
    
    f_new = open(new, 'r')
    new_lines = f_new.readlines()
    f_new.close()
    
    
    difference = compare(''.join(old_lines), ''.join(new_lines))
    if (difference == None):
      no_diff.append(new)
    else:
      diff.append([new, difference])
  
  # Now we output results
  if (no_diff != []):
    print "No differences seen on :",', '.join(no_diff)
  
  if (diff != []):
    for (file, comp) in diff:
      print("\nFor "+file)
      print(comp)
  
  return 0

FOLDER = "simu_test"

##################
# Outputs of various binaries and tests to compare with the actual ones. 
# Theses outputs are those of the original version of mercury, that is, mercury6_2.for
##################



os.chdir(FOLDER)
# We clean undesirable files. 
clean(["out", "clo", "aei", "dmp", "tmp"])

# We clean original files from mercury, element and close
clean(["ori"])

print("""###################
# Running original binaries #
###################""")
start = time()
(merc_or_stdout, merc_or_stderr) = run("../mercury_original/mercury")
t_merc_ori = time() - start

start = time()
(clo_or_stdout, clo_or_stderr) = run("../mercury_original/close")
t_clo_ori = time() - start

start = time()
(ele_or_stdout, ele_or_stderr) = run("../mercury_original/element")
t_ele_ori = time() - start

for file in MERCURY_FILENAMES:
  os.rename(file, file+EXTENTION_ORIGINAL)

for file in CLOSE_FILENAMES:
  os.rename(file, file+EXTENTION_ORIGINAL)
  
for file in ELEMENT_FILENAMES:
  os.rename(file, file+EXTENTION_ORIGINAL)
  
MERCURY_FILENAMES_OLD = []
for file in MERCURY_FILENAMES:
  MERCURY_FILENAMES_OLD.append(file+EXTENTION_ORIGINAL)

ELEMENT_FILENAMES_OLD = []
for file in ELEMENT_FILENAMES:
  ELEMENT_FILENAMES_OLD.append(file+EXTENTION_ORIGINAL)
  
CLOSE_FILENAMES_OLD = []
for file in CLOSE_FILENAMES:
  CLOSE_FILENAMES_OLD.append(file+EXTENTION_ORIGINAL)


print("""###################
# Running new binaries #
###################""")

clean(["out"])

start = time()
(merc_new_stdout, merc_new_stderr) = run("../mercury")
t_merc_new = time() - start

start = time()
(clo_new_stdout, clo_new_stderr) = run("../close")
t_clo_new = time() - start

start = time()
(ele_new_stdout, ele_new_stderr) = run("../element")
t_ele_new = time() - start

  
print("""###################
# Test of mercury #
###################""")


print("For the Output of mercury")
diff = compare(merc_or_stdout, merc_new_stdout)
if (diff != None):
  print diff
else:
  print("OK : No differences\n")

compare2file(MERCURY_FILENAMES_OLD, MERCURY_FILENAMES)
  
print("""
#################
# Test of close #
#################""")


print("For the Output of close")
diff = compare(clo_or_stdout, clo_new_stdout)
if (diff != None):
  print diff
else:
  print("OK : No differences\n")

compare2file(CLOSE_FILENAMES_OLD, CLOSE_FILENAMES)

print("""
###################
# Test of element #
###################""")


print("For the Output of element")
diff = compare(ele_or_stdout, ele_new_stdout)
if (diff != None):
  print diff
else:
  print("OK : No differences\n")

compare2file(ELEMENT_FILENAMES_OLD, ELEMENT_FILENAMES)

# Not representative, it depend really on the execution. We should do a mean value on several runs, or a longer run.
#~ print("""
#~ #############################
#~ # Execution time comparison #
#~ #############################""")
#~ 
#~ # We determine the pourcentage of difference from original and new binaries
#~ pourcent = (t_merc_new - t_merc_ori) / t_merc_ori
#~ merc_pourcent = " ("+str(round(pourcent*100))+")"
#~ 
#~ pourcent = (t_ele_new - t_ele_ori) / t_ele_ori
#~ ele_pourcent = " ("+str(round(pourcent*100))+")"
#~ 
#~ pourcent = (t_clo_new - t_clo_ori) / t_clo_ori
#~ clo_pourcent = " ("+str(round(pourcent*100))+")"
#~ 
#~ print("binary\toriginal\tnew")
#~ print("mercury\t"+str(t_merc_ori)+"\t"+str(t_merc_new)+merc_pourcent)
#~ print("element\t"+str(t_ele_ori)+"\t"+str(t_ele_new)+ele_pourcent)
#~ print("close\t"+str(t_clo_ori)+"\t"+str(t_clo_new)+clo_pourcent)

