#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.1

from math import *
import pylab as pl
import os, pdb, autiwa
import numpy as np
from constants import MT, MS
import sys # to be able to retrieve arguments of the script

###############################################
## Beginning of the program
###############################################
OUTPUT_EXTENSION = "pdf"

####################
# We read OPTIONS
####################
isProblem = False
problem_message = "AIM : Plot the evolution of mass with time, of all planets in the simulation" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * tmax=1e6 : the end of the output [years]" + "\n" + \
" * tmin=1e3 : the beginning of the output [years]" + "\n" + \
" * ext=%s : The extension for the output files" % OUTPUT_EXTENSION + "\n" + \
" * help : Display a little help message on HOW to use various options"

value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'tmin'):
    t_min = float(value)
  elif (key == 'tmax'):
    t_max = float(value)
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

####################
# On récupère la liste des fichiers planètes.aei
####################
(process_stdout, process_stderr, return_code) = autiwa.lancer_commande("ls *.aei")
if (return_code != 0):
  print("the command return an error "+str(return_code))
  print(process_stderr)
  exit()
  
liste_aei = process_stdout.split("\n")
liste_aei.remove('') # we remove an extra element that doesn't mean anything
nb_planete = len(liste_aei)


####################
# On lit, pour chaque planète, le contenu du fichier et on stocke les variables qui nous intéressent.
####################
t = [] # temps en année
m = [] # masse de la planète en masse terrestre



# On récupère les données orbitales
for planete in range(nb_planete):
  
  fichier_source = liste_aei[planete]
  tableau = open(fichier_source, 'r')
  
  t.append([]) # temps en année
  m.append([]) # masse de la planète en masse solaire

  # On passe les 3 premières lignes d'entête.
  for indice in range(3):
    tableau.readline()

  entete = tableau.readline()
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
      ti = float(colonne[0])
      mi = float(colonne[7])
    except:
      pass
    t[-1].append(ti)
    m[-1].append(mi * MS / MT)

    
  tableau.close()

dt = t[0][1] - t[0][0]
if ('t_max' in locals()):
  id_max = int((t_max - t[0][0]) / dt)
else:
  #~ t_max = max([max(ti) for ti in i])
  id_max = max([len(ti) for ti in t]) - 1
  t_max = t[0][0] + id_max * dt

if ('t_min' in locals()):
  id_min = int((t_min - t[0][0]) / dt)
else:
  t_min = t[0][0]
  id_min = 0
  
# on trace les plots

fig = pl.figure()
plot_m = fig.add_subplot(1, 1, 1)

for planet in range(nb_planete):
  plot_m.plot(t[planet][id_min:id_max], m[planet][id_min:id_max])
plot_m.set_xlim([t_min, t_max])
plot_m.set_xlabel("Time [years]")
plot_m.set_ylabel("Mass [Earths]")
plot_m.grid(True)

nom_fichier_plot = "mass"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

