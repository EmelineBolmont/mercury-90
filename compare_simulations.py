#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# The goal is to compare with a given version of the code, either an old one, the previous or current one.

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "3 august 2014"
__version__ = "$Revision: 3.0 $"

import sys
import os
import difflib # To compare two strings
import subprocess # To launch various process, get outputs et errors, returnCode and so on.
import pdb # To debug
import glob # to get list of file through a given pattern

NEW_TEST = "example_simulation"
PREVIOUS_TEST = "old_simulation"
REVISION = "HEAD"
PROGRAM_NAME = "mercury"

# Parameters
force_source = False # To force the compilation of every module
force_simulation = False # To force generation of simulation outputs for the "old" version of the code

isProblem = False
problem_message = """Script that run a mercury simulation and test if the outputs and binaries have correct behaviour.
The goal is to compare with a given version of the code, either an old one, the previous or current one.
By default, we take the current version of the code (the first time) and compare with the non committed modifications. 
but one can force the change of reference version by using the "rev=xxx" option. Without changing the code, 
one can force the comparison with a different simulation (the one in "example_simulation" will be copied in the other folder)
by using the "force" option.

The script can take various arguments :
(no spaces between the key and the values, only separated by '=')
 * help : display a little help message on HOW to use various options
 * force : To force generation of outputs for the 'old' program (after copying simulation files from the example)
 * actual : To force copying HEAD simulation, compyling it, then generating simulation outputs
 * faq : Display possible problems that might occurs during comparison
 * rev=%s : (previous, actual, current) are possible. Else, every Git ID syntax is OK.
the reference revision for the comparison with actual code

 Example : 
(examples are ordered. From the more common, to the more drastic.)
 compare_simulations.py #only generate new outputs, do nothing for the old one
 compare_simulations.py force # copy inputs in old folder, then generate outputs for both binaries
 compare_simulations.py actual # compile HEAD, copy inputs in old folder, then generate outputs for for both binaries
 compare_simulations.py rev=cdabb998 # compile the given revision, copy input in old folder then generate outputs for both binaries""" % REVISION

isFAQ = False
faq_message = """* If you have differences, ensure that all 
your modules have been compiled with 'test' options. 
To make sure of that:
> Makefile.py test force

Once this is done, make:
> Makefile.py force
to recompile all module with speed options.
"""

value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'force'):
    force_simulation = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'actual'):
    force_source = True
    force_simulation = True
    REVISION = "HEAD"
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'help'):
    isProblem = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'faq'):
    isFAQ = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'rev'):
    # If a revision is specified, we force the actualisation of source code, compilation and simulation.
    force_source = True
    force_simulation = True
    if (value in ['actual', 'current']):
      REVISION = "HEAD"
    elif (value in ['previous']):
      REVISION = "HEAD^"
    else:
      REVISION = value
  else:
    print("the key '%s' does not match" % key)
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

if isFAQ:
  print(faq_message)
  exit()


def run(commande):
  """lance une commande qui sera typiquement soit une liste, soit une 
  commande seule. La fonction renvoit un tuple avec la sortie, 
  l'erreur et le code de retour"""
  if (type(commande)==list):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  elif (type(commande)==str):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  else:
    raise TypeError("The command is neither a string nor a list.")
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  # there is .poll() or .wait() but I don't remember the difference. For some kind of things, one of the two was not working
  return (process_stdout, process_stderr, returncode)

def clean():
  """Delete all outputs files for a given nautilus simulation"""
  run("rm *.tmp")
  run("rm *.dmp")
  run("rm *.out")

def ASCIICompare(original, new):
  """function that compare and print differences between to strings that are compared line by line."""
  
  modified = ["+ ", "- ", "? "]
  
  # We want lists, because difflib.Differ only accept list of lines.
  original = original.split("\n")
  new = new.split("\n")
  
  d = difflib.Differ()

  result = list(d.compare(original, new))
  
  differences = []
  line_number_original = 0
  line_number_new = 0
  for line in result:
    if (line[0:2] == "  "):
      line_number_original += 1
      line_number_new += 1
    elif (line[0:2] == "- "):
      line_number_original += 1
      differences.append("[ori] l"+str(line_number_original)+" :"+line[2:])
    elif (line[0:2] == "+ "):
      line_number_new += 1
      differences.append("[new] l"+str(line_number_new)+" :"+line[2:])
    elif (line[0:2] == "? "):
      differences.append("      l"+str(max(line_number_new, line_number_original))+" :"+line[2:])

  # We print separately because it is more convenient if we want to store in a file instead.
  if (differences != []):
    return "\n".join(differences)
  else:
    return None

def compare2files(ori_files,new_files):
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
    
    
    difference = ASCIICompare(''.join(old_lines), ''.join(new_lines))
    if (difference == None):
      no_diff.append(new)
    else:
      diff.append([new, difference])
  
  # Now we output results
  if (diff != []):
    for (file, comp) in diff:
      print("\nFor "+file)
      print(comp)
      
    if (no_diff != []):
      print "No differences seen on :",', '.join(no_diff)
  else:
    print("Everything OK")  
  
  return 0

def compare2Binaries(ori_files, new_files):
  """Function that will use compare to see differences between 'original' 
  that is thought to be a variable and 'new_file' that is the name of a 
  file to read then use as input
  """
  no_diff = []
  diff = []
  
  for (original, new) in zip(ori_files, new_files):
    (stdout, stderr, returnCode) = run("md5sum %s" % original)
    md5_ori = stdout.split()[0]
    
    (stdout, stderr, returnCode) = run("md5sum %s" % new)
    md5_new = stdout.split()[0]
    
    if (md5_new != md5_ori):
      diff.append(original)
    else:
      no_diff.append(original)
  
  # Now we output results
  if (diff != []):
    for filename in diff:
      print("\ndifferences with binary  %s" % filename)
      
    if (no_diff != []):
      print("No differences seen on :%s" % ', '.join(no_diff))
  else:
    print("Everything OK")
  
  return 0

##################
# Outputs of various binaries and tests to compare with the actual ones. 
# Theses outputs are those of the original version of mercury, that is, mercury6_2.for
##################


os.chdir(NEW_TEST)
# We clean undesirable files beforehand, because we will copy if necessary the input simulation files to the other folder. 
clean()
os.chdir("..")

# We create folder and force old simulation generation if this is the first time we run the script
if not(os.path.isdir(PREVIOUS_TEST)):
  os.mkdir(PREVIOUS_TEST)
  force_source = True
  force_simulation = True

# We delete old files, get the desired revision of the code, and the corresponding simulation files, compile it and so on.
if force_source:
  print("Preparing old binaries ...")
  run("rm %s/*" % PREVIOUS_TEST) # Delete all files in the test directory
  
  # Copy the given revision REVISION in the required sub-folder (PREVIOUS_TEST)
  ## REVISION can either be a given commit, or HEAD, or READ^ and so on. Any commit ID available in Git.
  #> @Warning Problems for comparison can occurs if output files changes format between the two revisions.
  ## This script is intended to compare recent version of the code, when outputs changes are not expected. When outputs changes,
  ## one must be carefull with implementation and do the test themselves.
  get_revision = "git archive %s --format=tar --prefix=%s/ | tar xf -" % (REVISION, PREVIOUS_TEST)
  
  print(get_revision)
  run(get_revision)
  
  # We retrieve the commit ID from the possible alias stored in 'REVISION'
  (REVISION_ID, dummy, returnCode) = run("git rev-parse %s" % REVISION)
  (HEAD_ID, dummy, returnCode) = run("git rev-parse HEAD")
  
  os.chdir(PREVIOUS_TEST)
  
  # We store information about old commit used for comparison, in the corresponding folder.
  revision_file = open("revision.in", 'w')
  revision_file.write("Old revision: %s\n" % REVISION)
  revision_file.write("Old revision ID: %s\n" % REVISION_ID)
  revision_file.write("Current revision ID (HEAD): %s\n(but uncommitted changes might exists)\n" % HEAD_ID)
  revision_file.close()
  
  # Compilation of previous code
  previous_compilation = "Makefile.py test"
  print(previous_compilation)
  (stdout, stderr, returnCode) = run(previous_compilation)
  
  if (returnCode != 0):
    print(stdout)
    print(stderr)
  
  os.chdir("..")

if force_simulation:
  os.chdir(PREVIOUS_TEST)
  # Copy of simulation files
  copy_files = "cp ../%s/* ." % NEW_TEST
  print(copy_files)
  run(copy_files)
  
  # We run the simulation
  print("##########################################")
  sys.stdout.write("Running old binaries ...\r")
  sys.stdout.flush()
  (naut_or_stdout, naut_or_stderr, returnCode) = run("./%s" % PROGRAM_NAME)
  print("Running old binaries ...ok")
  print("##########################################")

  # Go back in parent directory (containing the current code and test script)
  os.chdir("..")
else:
  print("Skipping running original Nautilus, output already exists")

print("##########################################")
sys.stdout.write("Running new binaries ...\r")
sys.stdout.flush()
os.chdir(NEW_TEST)

(naut_new_stdout, naut_new_stderr, returnCode) = run("../%s" % PROGRAM_NAME)

OUTPUT_FILENAMES = ["xv.out"]
CLOSE_FILENAMES = ["ce.out"]

# list are sorted to ensure we compare the right files between actual and original outputs
OUTPUT_FILENAMES.sort()
CLOSE_FILENAMES.sort()

os.chdir("..")
print("Running new binaries ...ok")
print("##########################################")

if not(force_simulation):
  # To prevent finding differences in the standard output and error when the old simulation is not re-generated here
  naut_or_stdout = naut_new_stdout
  naut_or_stderr = naut_new_stderr  

os.chdir(PREVIOUS_TEST)
OUTPUT_FILENAMES_OLD = ["xv.out"]
CLOSE_FILENAMES_OLD = ["ce.out"]

# list are sorted to ensure we compare the right files between actual and original outputs
OUTPUT_FILENAMES_OLD.sort()
CLOSE_FILENAMES_OLD.sort()

os.chdir("..")

print("##########################################")

# We make the comparison

if (len(OUTPUT_FILENAMES_OLD) != len(OUTPUT_FILENAMES)):
  print("Error: number of abundances files is different")
  print("Old: %s" % OUTPUT_FILENAMES_OLD)
  print("Actual: %s" % OUTPUT_FILENAMES)

if (len(CLOSE_FILENAMES_OLD) != len(CLOSE_FILENAMES)):
  print("Error: number of rates files is different")
  print("Old: %s" % CLOSE_FILENAMES_OLD)
  print("Actual: %s" % CLOSE_FILENAMES)

diff = ASCIICompare(naut_or_stdout, naut_new_stdout)
if (diff != None):
  print("\nTest of nautilus")
  print("\tFor the Output of nautilus")
  print diff

# We create names including the folder in which they are
CLOSE_FILENAMES = [os.path.join(NEW_TEST, filename) for filename in CLOSE_FILENAMES]
OUTPUT_FILENAMES = [os.path.join(NEW_TEST, filename) for filename in OUTPUT_FILENAMES]

# We create names including the folder in which they are
CLOSE_FILENAMES_OLD = [os.path.join(PREVIOUS_TEST, filename) for filename in CLOSE_FILENAMES_OLD]
OUTPUT_FILENAMES_OLD = [os.path.join(PREVIOUS_TEST, filename) for filename in OUTPUT_FILENAMES_OLD]

print("comparing outputs:")
compare2Binaries(OUTPUT_FILENAMES_OLD, OUTPUT_FILENAMES)

print("comparing close encounters:")
compare2Binaries(CLOSE_FILENAMES_OLD, CLOSE_FILENAMES)

ASCII_FILES = ['info.out']

# We include the folder name because we are in the parent folder.
ASCII_OLD = [os.path.join(PREVIOUS_TEST, filename) for filename in ASCII_FILES]
ASCII_NEW = [os.path.join(NEW_TEST, filename) for filename in ASCII_FILES]

#~ pdb.set_trace()
print("comparing info.out")
compare2files(ASCII_OLD, ASCII_NEW)
