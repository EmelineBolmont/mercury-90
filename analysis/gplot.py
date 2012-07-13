#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that search for every gnuplot script available
from __future__ import print_function

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "13 juillet 2012"
__version__ = "$Revision: 1.1 $"
__credits__ = """Display each gnuplot script available in the current working directory and allow us to choose which one we want to run"""

import os
import subprocess
import pdb
import sys # to allow retrieving of parameters of the script

GNUPLOT_EXTENSION = "gnuplot"

VISUALISEUR = {"svg":"gthumb", "png":"gthumb", "pdf":"gv", "jpg":"gthumb"}

def lister(scheme):
  """list all the files corresponding to the given expression (might 
  be regular, in fact the same rules that goes for bash).  The function
  return the list of files, or the return code if their was an error."""
  
  if (type(scheme) != str):
    print("The argument must be a string")
    return None
  
  commande = "ls "+scheme
  
  process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

  (process_stdout, process_stderr) = process.communicate()
  returnCode = process.poll()
  
  # If returnCode is not 0, then there was a problem
  if (returnCode==0):
    files = process_stdout.split("\n")
    files.remove('')
    return files
  else:
    return returnCode

def run(command):
  """run a command that will be a string.
  The function return a tuple with the output, 
  the stderr and the return code. 
  It is fitted to run only command placed in the parent directory (of type "../command")"""
  
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  
  (process_stdout, process_stderr) = process.communicate()
  returnCode = process.poll()
  
  if (returnCode == 0):
    #print(command+" was runned successfully\n")
    
    # We delete the last \n introduced by subprocess and that corrupt the diff managed to see differences between outputs.
    return (process_stdout, process_stderr)
  else:
    print(command+" has terminated with an error code: "+str(returnCode))
    
    print("\nstderr:")
    if (type(process_stderr) == list):
      for line in process_stderr:
        print(line)
    else:
      print(process_stderr)
    
    return (process_stdout, process_stderr)
    
def enumeration(liste):
  """liste, un élément par ligne, les éléments de l'argument
  avec les indices associés"""
  output_list = []
  for (indice,element) in enumerate(liste):
    txt = "%2i : %-34s" % (indice, os.path.splitext(element)[0])
    output_list.append(txt)
  
  print(" ".join(output_list))

def getOutput(scriptname):
  """function that return the name of the output file of the given gnuplot script, if there is any"""
  
  f=open(scriptname,'r')
  lines = f.readlines()
  f.close()
  
  for line in lines:
    words = line.split()
    if (words[0] == "set") and (words[1] == "output"):
      output = words[2][1:-1]
      
      # in case of output file printed not in the directory of the .dat file
      temp = output.split("/")
      output = temp[-1]
      
      return output
  
  return None


###############################################
## Beginning of the program
###############################################

#>>>>>>>>>>>>> We read parameters of the script <<<<<<<<<<<<<<<<

indice_script = None # to allow default actions if nothing is given in parameters
generate_all = False

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * script (the script we want to launch, 'all' if we want to execute them all)" + "\n" + \
" * help (display a little help message on HOW to use various options"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'script'):
    if (value == 'all'):
      generate_all = True
      indice_script = 0
    else:
      generate_all = False
      indice_script = int(value)
  elif (key == 'help'):
    isProblem = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()
  
#>>>>>>>>>>>>> Get input if no parameters specified the gnuplot script to launch <<<<<<<<<<<<<<<<
scripts = lister("*."+GNUPLOT_EXTENSION)
nb_scripts = len(scripts)

# If the action is not given in parameters, we ask the user through standard input
if (indice_script == None):
  enumeration(scripts)

  while not(0<=indice_script<nb_scripts):
    try:
      txt_input = raw_input("Quel graphique voulez-vous afficher? (0-"+str(nb_scripts-1)+" ; 'all' pour tous les traiter ; 'l' pour la liste)\n")
      
      indice_script = int(txt_input)
    except ValueError:
      if (txt_input == 'l'):
        enumeration(scripts)
      elif (txt_input == 'all'):
        generate_all = True
        indice_script = 0
      else:
        print("The parameter must be between 0 and %i" % (nb_scripts-1))

if not(generate_all):
  indexes = [indice_script]
else:
  indexes = range(nb_scripts)

# We run gnuplot for every script needed. 
for index in indexes:
  script = scripts[index]

  (stdout, stderr) = run("gnuplot "+script)

# We only display the output file for the last script (in order to display only one graphic if we choose to run them all
output_file = getOutput(script)
if (output_file != None):
  output_extension = os.path.splitext(output_file)[1][1:] # [1:] is here to get rid of the prefixed dot  "." of the extension

  run(VISUALISEUR[output_extension]+" "+output_file)
else:
  print("Warning: The current gnuplot file has no output file")




