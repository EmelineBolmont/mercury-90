#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.2
# To display the manual torque profile, if any, of a given mercury simulation. The file "torque_profile.dat" is needed

import pylab as pl
import sys # to be able to retrieve arguments of the script


###############################################
## Beginning of the program
###############################################

isLog = False # We set the false option before. Because if not, we will erase the 'true' with other option that are not log, and 
# thus will lead to be in the else and put log to false.
OUTPUT_EXTENSION = 'pdf' # default value

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help (display a little help message on HOW to use various options" + "\n" + \
" * ext=png (The extension for the output files)"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'help'):
    isProblem = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

####################
# On lit, pour chaque planète, le contenu du fichier et on stocke les variables qui nous intéressent.
####################
distance = [] # demi-grand axe en ua
torque = [] # torque in units of Gamma_0

# On récupère les données orbitales


tableau = open("torque_profile.dat", 'r')

# On passe la 1e ligne d'entête.
for indice in range(1):
  tableau.readline()

for ligne in tableau:
  # Pour chaque ligne du tableau, on découpe suivant les 
  # espaces. (par défaut, s'il y a plusieurs espaces à la 
  # suite, il ne va pas créer 50 000 colonnes, ce qui ne 
  # serait pas le cas avec un 'split(" ")')
  colonne = ligne.split()
  
  # On essaye de rajouter les éléments. Si un seul d'entre eux
  # à un soucis, on élimine toute la ligne (en gros, lors de 
  # l'éjection d'une planète, certains paramètres peuvent 
  # devenir ****** au lieu d'un float, et ça génère une erreur
  # lors de la conversion en float)
  try:
    ai = float(colonne[0])
    Ti = float(colonne[1])

  except:
    pass
  distance.append(ai)     # demi-grand axe en ua
  torque.append(Ti)       # torque in units of Gamma_0

tableau.close()

# on trace les plots

pl.figure(1)
pl.clf()
pl.plot(distance, torque)

pl.xlabel("Distance [AU]")
pl.ylabel("Torque [Gamma/Gamma_0]")
pl.autoscale(axis='x', tight=True)
#~ pl.legend()
pl.grid(True)



#~ dossier_output = "output"
#~ system("mkdir dossier_output")
#~ system("cd dossier_output")

pl.figure(1)
nom_fichier_plot = "torque_profile"
#~ pl.savefig(nom_fichier_plot+'.svg', format='svg')
pl.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

