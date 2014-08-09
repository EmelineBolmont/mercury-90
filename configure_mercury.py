#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that add commands in the .bash_profile

import pdb
import os
import sys
import subprocess

BASH_PROFILE = os.path.join(os.path.expanduser("~"),".bash_profile")

BEGIN_MARK = "@mercury_begin"
END_MARK = "@mercury_end"

MERCURY_PATH = os.path.dirname(os.path.realpath(sys.argv[0]))

MERCURY_ADDON = """# %s DO NOT EDIT !
#########################################

# To access mercury folder via 'cd $mercury'
mercury="%s"

# import bash commands from a Git file
source $mercury/.mercury_profile

#########################################
# %s DO NOT EDIT !""" % (BEGIN_MARK, MERCURY_PATH, END_MARK)

is_uninstall = False

isProblem = False
problem_message = """AIM : Add commands to the .bash_profile
"The script can take various arguments :
"(no spaces between the key and the values, only separated by '=')
" * uninstall: To suppress any command that were previously added to the .bash_profile
" * help : display a little help message on HOW to use various options."""


value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

def run(command):
  """Run a command that can be a list or a string of list. 
  The function return a tuple with standard output, 
  standard error and return code."""
  if (type(command)==list):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  elif (type(command)==str):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  else:
    raise TypeError("The command is neither a string nor a list.")
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  # there is .poll() or .wait() but I don't remember the difference. For some kind of things, one of the two was not working
  return (process_stdout, process_stderr, returncode)

def delete_mercury_addon():
  """
  """
  object_file = open(BASH_PROFILE, 'r')
  
  lines = []
  skip = False
  for line in object_file:
    # If we find the begin marker, we flag a region where we start to ignore lines
    if BEGIN_MARK in line:
      skip = True
    
    # Before testing for ignore mark, we skip or store lines (because we want to ignore the end marker also)
    if skip:
      pass
    else:
      lines.append(line)
    
    if END_MARK in line:
      skip = False
  
  object_file.close()
  
  object_file = open(BASH_PROFILE, 'w')
  
  for line in lines:
    object_file.write(line)
  
  object_file.close()
  

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'uninstall'):
    is_uninstall = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'help'):
    isProblem = True
    if (value != None):
      print(value_message % (key, key, value))
  else:
    print("the key '%s' does not match" % key)
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

# Test if the bash profile exist. Test also if one installation was already done.
if (os.path.isfile(BASH_PROFILE)):
  (process_stdout, process_stderr, returnCode) = run("grep '%s' %s|wc -l" % (BEGIN_MARK, BASH_PROFILE))
  
  test = float(process_stdout)
  
  if (test == 0):
    already_installed = False
  elif (test == 1):
    already_installed = True
  else:
    raise ValueError("Error: '%s' exist more than once in %s" % (BEGIN_MARK, BASH_PROFILE))
    exit()
else:
  already_installed = False



if (already_installed):
  delete_mercury_addon()

if (is_uninstall):
  if (not(already_installed)):
    print("Nothing to uninstall in %s." % BASH_PROFILE)
  exit()

object_file = open(BASH_PROFILE, 'a')
  
object_file.write(MERCURY_ADDON)

object_file.close()

print("Modifs will be available at next start.")
print("If you do not want to wait, run:")
print("source %s" % BASH_PROFILE)

