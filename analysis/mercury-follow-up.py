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

class Temps(object):
  """Classe qui dÃ©finit un temps. 
  Ceci permet d'additionner deux objets temps, les afficher Ã  l'aide 
  de print, et ainsi de les manipuler en toute transparence
  
  *values : (an, jour, heure, minute, seconde) liste de valeurs, 
  au maxi 5, au mini 1. Si une seule valeur est donnÃ©e, elle sera prise
   comme Ã©tant le nombre de seconde. Si 2 sont donnÃ©es, la premiÃ¨re 
   sera le nombre de minutes et la 2e le nombre de secondes, et ainsi 
   de suite au fur et Ã  mesure que les nombres de valeurs augmentent.
  """
  
  NB_DIGITS = 2
  
  # Valeur en seconde de diffÃ©rentes durÃ©es
  MINUTE = 60
  HEURE = 60 * MINUTE
  JOUR = 24 * HEURE
  AN = 365.25 * JOUR
  
  def __init__(self, *values):
    nb_values = len(values)
    values = list(values) # On converti le tuple pour pouvoir ajouter des valeurs
    if (nb_values > 5):
      raise ValueError("The number of parameters must be of 5 maximum")
    for i in range(5-nb_values):
      values.insert(0,0)
    
    self.temps = float(values[0] * Temps.AN + values[1] * Temps.JOUR + values[2] * Temps.HEURE + values[3] * Temps.MINUTE + values[4])
    self.__update()
  
  def __update(self):
    """mÃ©thode privÃ©e qui met Ã  jour toutes les variables internes
    
    met Ã  jour self.nb_an, self.nb_jour, self.nb_heure, self.nb_minute et self.nb_seconde
    """
    reste = self.temps
  
    self.nb_an = int(reste / Temps.AN)
    reste -= self.nb_an * Temps.AN
    
    self.nb_jour = int(reste / Temps.JOUR)
    reste -= self.nb_jour * Temps.JOUR
    
    self.nb_heure = int(reste / Temps.HEURE)
    reste = reste - self.nb_heure * Temps.HEURE
    
    self.nb_minute = int(reste / Temps.MINUTE)
    reste = reste - self.nb_minute * Temps.MINUTE
    
    self.nb_seconde = round(reste, Temps.NB_DIGITS)
  
  def __str__(self):
    """Permet, via print, d'afficher le temps
    """
    # On prÃ©pare la chaÃ®ne de caractÃ¨re
    str_temps = ''
    if (self.nb_an == 1):
      str_temps += str(self.nb_an)+"an "
    elif (self.nb_an > 1):
      str_temps += str(self.nb_an)+"ans "
      
    if (self.nb_jour == 1):
      str_temps += str(self.nb_jour)+"jour "
    elif (self.nb_jour > 1):
      str_temps += str(self.nb_jour)+"jours "
      
    if (self.nb_heure != 0):
      str_temps += str(self.nb_heure)+"h "
      
    if (self.nb_minute != 0):
      str_temps += str(self.nb_minute)+"min "
      
    str_temps += str(self.nb_seconde)+"s"
    return str_temps
  
  def __add__(self, other):
    """Surcharge de +"""
    return Temps(self.temps + other.temps)
  
  def __sub__(self, other):
    """Surcharge de -"""
    return Temps(self.temps - other.temps)
  
  def __eq__(self, other):
    """surcharge de =="""
    if (self.temps == other.temps):
      return True
    else:
      return False
  
  def __ne__(self, other):
    """surcharge de !="""
    return not(self.__eq__(other))
    
def lancer_commande(commande):
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

def getJobInfos(jobID):
  """
  With the jobID (only integers) in parameter, will return in a tuple : 
  
  ellapsed time as an object 'Temps', and Curent workind directory (absolute path) of the job as a string. 
  
  Hence : (ellapsed_time, cwd)
  
  """
  
  (process_stdout, process_stderr, return_code) = lancer_commande("qstat -j %d" % jobID)
  
  remaining_time = "???"

  # cpu time
  try:
    index = process_stdout.index("cpu=")
    length = process_stdout[index:].index(',')
    cpu_info = process_stdout[index+4:index+length]
    
    ellapsed_time = Temps(*map(float,cpu_info.split(":")))
  except:
    ellapsed_time = Temps(0)
  
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
" * verbose : to display more infos about the simulations (must be defined after 'all' option)"  + "\n" + \
" * help : display a little help message on HOW to use various options"   + "\n" + \
"\nExample:"  + "\n" + \
"> mercury-follow-up.py"  + "\n" + \
"> mercury-follow-up.py all"  + "\n" + \
"> mercury-follow-up.py all verbose"

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
  
  (stdout, stderr, returnCode) = lancer_commande('qstat -u %s' % logname)
  if (returnCode == 0):
    lines = stdout.split("\n")
  
  # We delete the two first lines
  for i in range(2):
    del(lines[0])
  lines.remove('')
    
  
  jobIDs = [int(line.split()[0]) for line in lines]
else:
  (process_stdout, process_stderr, return_code) = lancer_commande("ls *.o[0-9]*")
  if (return_code != 0):
    print("the command return an error "+str(return_code))
    print(process_stderr)
    exit()
    
  list_jobs = process_stdout.split("\n")
  list_jobs.remove('') # we remove an extra element that doesn't mean anything
  list_jobs.sort()

  current_job = list_jobs[-1]

  jobID = int(current_job.split(".o")[1])
  jobIDs = [jobID]




####################
# We get the execution time, if any (must be a job, and not a manual execution)
####################

infoAll = [] # The array where to store display infos as strings (in tuple to allow sorting)
waiting_list = []
for jobID in jobIDs:
  (ellapsed_time, cwd) = getJobInfos(jobID)
  
  os.chdir(cwd)
  
  # In case it exists several old jobs for the same simulation :
  (process_stdout, process_stderr, return_code) = lancer_commande("ls *.o[0-9]*")
  if (return_code != 0):
    waiting_list.append("%d" % jobID)
    continue
    
  list_jobs = process_stdout.split("\n")
  list_jobs.remove('') # we remove an extra element that doesn't mean anything
  list_jobs.sort()
  
  if (len(list_jobs) > 1):
    ellapsed_time = Temps(0)
    for job in list_jobs:
      current_ID = int(job.split(".o")[1])
      (tmp_time, dummy) = getJobInfos(jobID)
      ellapsed_time = ellapsed_time + tmp_time
  
  # NORMAL CALCULATION
  
  # We count the initial number of bodies
  (stdout, stderr, returnCode) = lancer_commande('grep "m=" big.in|wc -l')
  if (returnCode == 0):
    init_nb_bodies = int(stdout)
    #~ init_nb_bodies.remove('') # we remove an extra element that doesn't mean anything
  else:
    init_nb_bodies = 0

  # We check which folders contain NaN
  (stdout, stderr, returnCode) = lancer_commande('grep "m=" big.dmp|wc -l')
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
  
  remaining_time = Temps(ellapsed_time.temps * (100. / percentage - 1.))
    
    
  # Print infos
  infos = "jobID %d\n" % jobID
  infos += "    %s\n" % cwd

  if isVerbose:
    infos += "    Integration time : %g / %g years (%.1f%%)\n" % (current_time, integration_time, percentage)
    infos += "    Number of bodies : Initial=%d ; Current=%d\n" % (init_nb_bodies, current_nb_bodies)
    infos += "    Remaining time < %s\n" % remaining_time
  else:
    infos += "    [end] %s (%.1f%%)\n" % (remaining_time, percentage)
  
  infoAll.append((remaining_time.temps, infos))

# To display infos in some order
infoAll.sort(reverse=True)

for (ti, info) in infoAll:
  print(info)

nb_wait = len(waiting_list)
print("%d jobs on waiting list : %s" % (nb_wait, " ".join(waiting_list)))