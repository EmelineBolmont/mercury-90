#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.1
# Pour lire le fichier temperature_profile et 
# afficher le profil de profondeur optique et de diffusivité thermique

import pylab as pl
import numpy as np
import sys

AU = 1.4959787e13 # AU = astronomical unit in [cm]
day = 86400 # s

OUTPUT_EXTENSION = "pdf"

####################
# We read OPTIONS
####################
isProblem = False
problem_message = "AIM : Plot the optical depth profile of the simulation." + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
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
# On lit, pour chaque planète, le contenu du fichier et on stocke les variables qui nous intéressent.
####################
tau = [] # optical depth
chi = [] # thermal diffusivity
index = [] # temperature index 
a = [] # demi-grand axe en ua
T = [] # temperature in K

# On récupère les données orbitales


tableau = open("temperature_profile.dat", 'r')

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
    indexi = float(colonne[2])
    chi_i = float(colonne[3])
    tau_i = float(colonne[4])
  except:
    pass
  tau.append(tau_i) # optical depth
  chi.append(chi_i) # thermal diffusivity
  index.append(indexi) # temperature index 
  a.append(ai)       # demi-grand axe en ua
  T.append(Ti)       # temperature in K

tau = np.array(tau)
chi = np.array(chi)
index = np.array(index)
a = np.array(a)
T = np.array(T)

chi = chi * AU**2 / day # we change the units to put in cgs, i.e in cm^2/s

tableau.close()

# on trace les plots

fig = pl.figure()
pl.clf()
# We create subplots. add_subplot(2, 3, 1) means we have 2 lines, 3 columns, 
# and that the active plot is the first, starting from top left (for 6 plots in total)
plot_tau = fig.add_subplot(2, 1, 1)
plot_tau.semilogy(a, tau)

plot_tau.set_xlabel("Semi-major axis [AU]")
plot_tau.set_ylabel("Optical depth")
#~ pl.legend()
plot_tau.grid(True)

plot_chi = fig.add_subplot(2, 1, 2)

plot_chi.plot(a, chi)

plot_chi.set_xlabel("Semi-major axis [AU]")
plot_chi.set_ylabel("Thermal diffusivity [cm^2/s]")
plot_chi.grid(True)



#~ dossier_output = "output"
#~ system("mkdir dossier_output")
#~ system("cd dossier_output")

nom_fichier_plot = "optical_depth_profile"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

