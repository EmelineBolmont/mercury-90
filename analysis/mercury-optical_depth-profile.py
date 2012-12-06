#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.1
# Pour lire le fichier temperature_profile et 
# afficher le profil de profondeur optique et de diffusivité thermique

import pylab as pl
import numpy as np

AU = 1.4959787e13 # AU = astronomical unit in [cm]
day = 86400 # s

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
    tau_i = float(colonne[3])
    chi_i = float(colonne[4])
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

pl.figure(1)
pl.clf()
# On crée des sous plots. Pour subplot(231), ça signifie qu'on a 2 lignes, 3 colonnes, et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
pl.subplot(211)
pl.plot(a, tau)

pl.xlabel(unicode("a [UA]",'utf-8'))
pl.ylabel(unicode("optical depth",'utf-8'))
#~ pl.legend()
pl.grid(True)

pl.subplot(212)
pl.plot(a, chi)

pl.xlabel(unicode("a [UA]",'utf-8'))
pl.ylabel(unicode("thermal diffusivity [cm^2/s]",'utf-8'))
pl.grid(True)



#~ dossier_output = "output"
#~ system("mkdir dossier_output")
#~ system("cd dossier_output")

pl.figure(1)
nom_fichier_plot = "optical_depth_profile"
#~ pl.savefig(nom_fichier_plot+'.svg', format='svg')
pl.savefig(nom_fichier_plot+'.pdf', format='pdf')

pl.show()

