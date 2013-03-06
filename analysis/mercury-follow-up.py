#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.1
# To display the evolution of a running simulation

import os
import sys
import autiwa
import pdb
import subprocess

# Get the machine hostname
#~ hostname = simulations_utilities.getHostname()

def getJobInfos(jobID):
  """
  With the jobID (only integers) in parameter, will return in a tuple : 
  
  ellapsed time as an object 'Temps', and Curent workind directory (absolute path) of the job as a string. 
  
  Hence : (ellapsed_time, cwd)
  
  """
  
  (process_stdout, process_stderr, return_code) = autiwa.lancer_commande("qstat -j %d" % jobID)
  
  remaining_time = "???"

  # cpu time
  index = process_stdout.index("cpu=")
  length = process_stdout[index:].index(',')
  cpu_info = process_stdout[index+4:index+length]
  
  ellapsed_time = autiwa.Temps(*map(float,cpu_info.split(":")))
  
  
  index = process_stdout.index("cwd:")
  length = process_stdout[index:].index('\n')
  cwd = process_stdout[index+4:index+length].lstrip()
    
  return (ellapsed_time, cwd)

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

isAll = False # to have info on all the running simulations
isProblem = False
isVerbose = True

problem_message = " This script will show information about a running simulation" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * all : to have info on all the running simulations"  + "\n" + \
" * help : display a little help message on HOW to use various options" 

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'help'):
    isProblem = True
  elif (key == 'all'):
    isAll = True
    isVerbose = False
  elif (key == 'verbose'):
    isVerbose = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

if isAll:
  logname = os.getlogin()
  
  (stdout, stderr, returnCode) = autiwa.lancer_commande('qstat -u %s' % logname)
  if (returnCode == 0):
    lines = stdout.split("\n")
  
  # We delete the two first lines
  for i in range(2):
    del(lines[0])
  lines.remove('')
    
  #~ pdb.set_trace()
  jobIDs = [int(line.split()[0]) for line in lines]
else:
  (process_stdout, process_stderr, return_code) = autiwa.lancer_commande("ls simulation.sh.o*")
  if (return_code != 0):
    print("the command return an error "+str(return_code))
    print(process_stderr)
    exit()
    
  list_jobs = process_stdout.split("\n")
  list_jobs.remove('') # we remove an extra element that doesn't mean anything
  list_jobs.sort()

  current_job = list_jobs[-1]

  jobID = int(current_job.split("simulation.sh.o")[1])
  jobIDs = [jobID]




####################
# We get the execution time, if any (must be a job, and not a manual execution)
####################


for jobID in jobIDs:
  (ellapsed_time, cwd) = getJobInfos(jobID)
  
  os.chdir(cwd)
  
  # We count the initial number of bodies
  (stdout, stderr, returnCode) = autiwa.lancer_commande('grep "m=" big.in|wc -l')
  if (returnCode == 0):
    init_nb_bodies = int(stdout)
    #~ init_nb_bodies.remove('') # we remove an extra element that doesn't mean anything
  else:
    init_nb_bodies = 0

  # We check which folders contain NaN
  (stdout, stderr, returnCode) = autiwa.lancer_commande('grep "m=" big.dmp|wc -l')
  if (returnCode == 0):
    current_nb_bodies = int(stdout)
    #~ init_nb_bodies.remove('') # we remove an extra element that doesn't mean anything
  else:
    current_nb_bodies = 0

  # Integration Time
  parameters = open("param.in", 'r')
  lines = parameters.readlines()
  parameters.close()

  # in reversed order to be sure to test all the line (because once we delete a line, the index are shifted)
  for line in reversed(lines):
    if (line[0] == ')'):
      lines.remove(line)
  t_start = float(lines[1].split("=")[1])
  t_stop = float(lines[2].split("=")[1])

  integration_time = (t_stop - t_start) / 365.25

  # Current time
  big = open("big.dmp", 'r')
  for i in range(5):
    line = big.readline()
  big.close()
  current_time = float(line.split("=")[1]) / 365.25 # In years

  percentage = current_time / integration_time * 100.
  

  remaining_time = autiwa.Temps(ellapsed_time.temps / (percentage/100.))
    
    
  # Print infos
  print("jobID %d : %s" % (jobID, cwd))
  print("  Integration time : %g / %g years (%.1f%%)" % (current_time, integration_time, percentage))

  if isVerbose:
    print("  Number of bodies : Initial=%d ; Current=%d" % (init_nb_bodies, current_nb_bodies))
    print("  Remaining time < %s" % remaining_time)
    