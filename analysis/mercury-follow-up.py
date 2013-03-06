#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.1
# To display the evolution of a running simulation

import sys
import autiwa
import pdb
import subprocess

# Get the machine hostname
#~ hostname = simulations_utilities.getHostname()


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

isProblem = False
problem_message = " This script will show information about a running simulation" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" 

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'help'):
    isProblem = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()


# We check which folders contain NaN
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

####################
# We get the execution time, if any (must be a job, and not a manual execution)
####################
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

(process_stdout, process_stderr, return_code) = autiwa.lancer_commande("qstat -j %d" % jobID)
cpu_info = ''
remaining_time = "???"
if (return_code == 0):
  index = process_stdout.index("cpu")
  length = process_stdout[index:].index(',')
  cpu_info = process_stdout[index+4:index+length]
  
  ellapsed_time = autiwa.Temps(*map(float,cpu_info.split(":")))

  remaining_time = autiwa.Temps(ellapsed_time.temps / (percentage/100.))
  
  
# Print infos
print("Number of bodies : Initial=%d ; Current=%d" % (init_nb_bodies, current_nb_bodies))
print("Integration time : %7.2g / %7.2g years" % (current_time, integration_time))
print("Completed : %.1f%%" % percentage)
print("Estimated remaining time < %s" % remaining_time)
