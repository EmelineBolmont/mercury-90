#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.1
# To check if there is NaN in the orbits, or if the simulation did not have time to finish itself before the allowed time by the server.

import os
import pdb
import autiwa
import sys
import subprocess
import simulations_utilities
import mercury_utilities

# Get current working directory
rep_exec = os.getcwd()

# Get the machine hostname
hostname = simulations_utilities.getHostname()

scriptFolder = os.path.dirname(os.path.realpath(__file__)) # the folder in which the module is. 
binaryPath = os.path.join(scriptFolder, os.path.pardir)

#    .-.     .-.     .-.     .-.     .-.     .-.     .-.     .-.     .-. 
#  .'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `.
# (    .     .-.     .-.     .-.     .-.     .-.     .-.     .-.     .    )
#  `.   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   .'
#    )    )                                                       (    (
#  ,'   ,'                                                         `.   `.
# (    (                     DEBUT DU PROGRAMME                     )    )
#  `.   `.                                                         .'   .' 
#    )    )                                                       (    (
#  ,'   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   `.
# (    '  _  `-'  _  `-'  _  `-'  _  `-'  _  `-'  _  `-'  _  `-'  _  `    )
#  `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .'
#    `-'     `-'     `-'     `-'     `-'     `-'     `-'     `-'     `-'
isOK = True # If no problems are found, this boolean help display a message that says that everything went fine.

isRestart = False # Say if we want to re-run the simulation or not.
isForcedStart = False # will erase output file if they exists and run all the simulations.
isForcedContinue = False # will continue the simulation no matter what
isContinue = False # Do we want to continue simulations that did not have time to finish?
WALLTIME = None

isProblem = False
problem_message = "AIM : Check if the status of the mercury simulations in the CWD (either sub or sub-sub folder --- meta option)" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * continue : Say that we want to continue the simulations that did not have enough time to finish" + "\n" + \
" * restart : Say if we want to re-run the simulation in case of NaN problems" + "\n" + \
" * force-start : will erase output file if they exists and run all the simulations." + "\n" + \
" * force-continue : will erase output file if they exists and continue all the simulations" + "\n" + \
" * walltime : (in hours) the estimated time for the job. Only used for avakas" + "\n" + \
"" + "\n" + \
"Example : \n" + \
"> mercury-check.py restart continue\n" + \
"will continue non finished simulation or restart it if NaN are present\n" + \
"> mercury-check.py walltime=48 continue\n"

value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'restart'):
    isRestart = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'force-start'):
    isForcedStart = True
    isRestart = False
    isContinue = False
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'force-continue'):
    isForcedContinue = True
    isForcedStart = False
    isRestart = False
    isContinue = False
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'continue'):
    isContinue = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'walltime'):
    WALLTIME = int(value)
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

if (('avakas' in hostname) and WALLTIME == None and (isContinue or isRestart or isForcedStart)):
  print("Walltime option must be set. type 'help' for a description of the options")
  exit()

absolute_parent_path = rep_exec

if not(os.path.isfile("param.in")):
  print("Doesn't look like a regular simulation folder")
  print("\t 'param.in' does not exist, script termianted")
  
  exit()


# We check informations on the current simulation
simulation_status = 0 # if 0, then the simulation ended correctly

# We check which folders contain NaN
isNaN = False
if os.path.isfile("big.dmp"):
  (stdout, stderr, returnCode) = autiwa.lancer_commande('grep -l "NaN" big.dmp')
  if (returnCode == 0):
    isNaN = True

# We check if the simulation had time to finish
command = 'tail info.out|grep "Integration complete"|wc -l'
(stdout, stderr, returnCode) = autiwa.lancer_commande(command)
if (returnCode != 0):
  print("The command '%s' did not end correctly" % command)
  print(stderr)
  pdb.set_trace()
is_finished = int(stdout.split("\n")[0]) # We get the number of times "Integration complete" is present in the end of the 'info.out' file

# If there is Nan, we do not want to continue the simulation, but restart it, or check manually, so theses two kinds of problems are separated.
if not(isNaN):
  if (is_finished == 0):
    simulation_status = 1
    
    isOK = False
else:
  simulation_status = 2
  
  isOK = False

if (simulation_status == 0):
  log_message = "The simulation is finished"
elif (simulation_status == 1):
  log_message = "The simulation is not finished"
elif (simulation_status == 2):
  log_message = "NaN are present"
else:
  log_message = None

if (log_message != None):
  print(log_message)

if (((simulation_status != 0) and (isContinue or isRestart)) or (isForcedStart or isForcedContinue)):
  mercury_utilities.prepareSubmission(BinaryPath=binaryPath, walltime=WALLTIME)

# If the option 'start' is given, we force the run of the simulation, whatever there is an old simulation or not in the folder.
if (isForcedStart):
  mercury_utilities.mercury_restart()

if (((simulation_status == 1) and isContinue) or isForcedContinue):
  mercury_utilities.mercury_continue()

if ((simulation_status == 2) and isRestart):
  mercury_utilities.mercury_restart()

# TODO Check in a folder if a simulation is currently running (don't know how to test that)
