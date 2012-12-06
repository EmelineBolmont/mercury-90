#!/usr/bin/env python 
# script to continue a mercury simulation assuming that a new folder with all the files of the previous 
#simulations has been created and that we are in this directory 
# Version 1.0

import os
import subprocess
import mercury

def lancer_commande(commande):
  """lance une commande qui sera typiquement soit une liste, soit une 
  commande seule. La fonction renvoit un tuple avec la sortie, 
  l'erreur et le code de retour"""
  if (type(commande)==list):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  elif (type(commande)==str):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  else:
    raise TypeError("La commande n'est ni une liste, ni une chaine de caractere.")
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  return (process_stdout, process_stderr, returncode)
    
# First, we need to get the dump files 
files_to_move = ["param", "big"]

if os.path.exists("small.dmp"):
  files_to_move.append("small")

for filename in files_to_move:
  lancer_commande("mv "+filename+".dmp "+filename+".in")

# We delete output files
lancer_commande("rm *.aei *.clo *.out")

# We delete temporary files (since the ones we were interested in were moved to *.in)
lancer_commande("rm *.tmp *.dmp")

paramin = mercury.Param('', '', '', '', '')
paramin.read()

integration_time = raw_input("For how many time do you want to launch the integration (in years)?\n")# In years

paramin.stop_time = paramin.start_time + float(integration_time) * 365.25 # in days

paramin.write()
