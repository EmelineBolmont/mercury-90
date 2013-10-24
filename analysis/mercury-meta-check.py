#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.1
# To check if there is NaN in the orbits, or if the simulation did not have time to finish itself before the allowed time by the server.

import os
import pdb
import autiwa
import sys
import subprocess
import mercury
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
isMeta = False # If we consider the current folder as a folder that list sub-meta-simulations where the simulations really are
isContinue = False # Do we want to continue simulations that did not have time to finish?
isVerbose = False # Force display of individual information when the 'meta' option is active
showFinished = False # By default, We only show problems. 
WALLTIME = None

isProblem = False
problem_message = "AIM : Check if the status of the mercury simulations in the CWD (either sub or sub-sub folder --- meta option)" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * meta : option that will consider the current folder as a folder that list meta simulation instead of simple simulations" + "\n" + \
" * continue : Say that we want to continue the simulations that did not have enough time to finish" + "\n" + \
" * restart : Say if we want to re-run the simulation in case of NaN problems" + "\n" + \
" * force-start : will erase output file if they exists and run all the simulations." + "\n" + \
" * force-continue : will erase output file if they exists and continue all the simulations" + "\n" + \
" * walltime : (in hours) the estimated time for the job. Only used for avakas" + "\n" + \
" * finished : Instead of show problems, will only show finished simulations" + "\n" + \
" * verbose : Froce show of individual informations when the 'meta' option is active" + "\n" + \
"" + "\n" + \
"Example : \n" + \
"> mercury-check-simulation.py meta restart continue\n" + \
"will continue non finished simulation and restart thoses with NaN considering\n" + \
" that each subfolder of the current folder contain a simulation in \n" + \
"each sub-folder (subsubfolder of the PWD).\n" + \
">mercury-check-simulation.py finished meta walltime=48\n" + \
"will show the simulations that are finished. "

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
  elif (key == 'verbose'):
    isVerbose = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'meta'):
    isMeta = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'walltime'):
    WALLTIME = int(value)
  elif (key == 'finished'):
    showFinished = True
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

if (('avakas' in hostname) and WALLTIME == None and (isContinue or isRestart or isForcedStart)):
  print("Walltime option must be set. type 'help' for a description of the options")
  exit()

# We go in each sub folder of the current working directory

# If sub folders are meta simulation folders instead of folders, we list the meta simulation folder to run the test in each sub folder.
if (isMeta):
  meta_list = [dir for dir in os.listdir(".") if (os.path.isdir(dir))]
  meta_list.sort()
else:
  meta_list = ["."]

logs = {}
nb_simulations = {}
nb_finished = {}
meta_OK = {}

for meta in meta_list:
  os.chdir(meta)
  
  if (meta == '.'):
    absolute_parent_path = rep_exec
  else:
    absolute_parent_path = os.path.join(rep_exec, meta)
  
  # We get the list of simulations
  simu_list = [dir for dir in os.listdir(".") if (os.path.isdir(dir))]
  #autiwa.suppr_dossier(liste_simu,dossier_suppr)
  simu_list.sort()
  
  logs[meta] = []
  nb_simulations[meta] = len(simu_list)
  finished = 0 # We initialize the number of simulations finished for this meta simulation
  
  # We check which folders contain NaN
  (stdout, stderr, returnCode) = autiwa.lancer_commande('grep -l "NaN" */big.dmp')
  if (returnCode == 0):
    NaN_folder = stdout.split("/big.dmp\n")
    NaN_folder.remove('') # we remove an extra element that doesn't mean anything
  else:
    NaN_folder = []

  for simu in simu_list:
    os.chdir(simu)
    
    simulation_status = 0 # if 0, then the simulation ended correctly
    
    if not(os.path.isfile("param.in")):
      print("%s/%s : doesn't look like a regular simulation folder" % (absolute_parent_path, simu))
      print("\t 'param.in' does not exist, folder skipped")
      os.chdir("..")
      break
    
    

    # We check if the simulation had time to finish
    # We get the integration time
    paramin = mercury.Param(algorithme="HYBRID", start_time=0, stop_time=0, output_interval=0,h=0)
    paramin.read(filename="param.dmp")
    stop_time = paramin.get_stop_time()

    bigdmp = open("big.dmp", 'r')
    for i in range(5):
        tmp = bigdmp.readline()
    # We round the result because the result can be slightly lower to the expected time, and the simulation still be over.
    big_time = float(tmp.split()[-1]) + 1.
    
    if (big_time > stop_time):
        is_finished = 1
    else:
        is_finished = 0

    
    # If there is Nan, we do not want to continue the simulation, but restart it, or check manually, so theses two kinds of problems are separated.
    if simu not in NaN_folder:
      if (is_finished == 0):
        simulation_status = 1
        
        isOK = False

    else:
      simulation_status = 2
      
      isOK = False
    
    if (simulation_status == 0 and showFinished):
      log_message = "%s/%s : The simulation is finished" % (absolute_parent_path, simu)
    
    elif (simulation_status == 1 and not(showFinished)):
      log_message = "%s/%s : The simulation is not finished" % (absolute_parent_path, simu)
    
    elif (simulation_status == 2 and not(showFinished)):
      log_message = "%s/%s : NaN are present" % (absolute_parent_path, simu)
    else:
		log_message = None
		
    if (log_message != None):
      if not(isMeta):
        print(log_message)
      else:
        logs[meta].append(log_message)
    
    if (((simulation_status != 0) and (isContinue or isRestart)) or (isForcedStart or isForcedContinue)):
      mercury_utilities.prepareSubmission(BinaryPath=binaryPath, walltime=WALLTIME)
    
    # If the option 'start' is given, we force the run of the simulation, whatever there is an old simulation or not in the folder.
    if (isForcedStart):
      mercury_utilities.mercury_restart()
    
    if (((simulation_status == 1) and isContinue) or isForcedContinue):
      mercury_utilities.mercury_continue()
    
    if ((simulation_status == 2) and isRestart):
      mercury_utilities.mercury_restart()
    
    
    
    if (simulation_status == 0):
      finished += 1
    
    # We get back in the parent directory
    os.chdir("..")
  
  nb_finished[meta] = finished
  meta_OK[meta] = isOK
  os.chdir(rep_exec)

if isMeta:
  if isVerbose:
    for meta in meta_list:
      for item in logs[meta]:
        print(item)
  
  for meta in meta_list:
    print("%s : %d/%d simulations are finished" % (meta, nb_finished[meta], nb_simulations[meta]))
else:
  if (isOK and not(isForcedStart)):
    print("All the simulations finished correctly and without NaN")
  
# TODO Check in a folder if a simulation is currently running (don't know how to test that)
