#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0.1.1
# To study the initial properties of the planets generated. Must be in a meta-simulation folder.


# TODO
# _ écrire un fichier qui stocke les paramètres orbitaux des planètes qu'il reste.
# _ sélectionner les planètes qui sont plus proche qu'une certaine valeur de leur étoile.

from math import *
import pylab as pl
import os, pdb, autiwa
from simu_constantes import *
import random
import numpy as np

# Mass threshold
MASS_THRESHOLD = 5 # Earth mass

CONVERGENCE_ZONE = 6. # in AU, the location of the convergence zone (as a reference for plotting period ratios)

DELTA_RATIO = 0.005
DELTA_M = 0.1

#######################
# On prépare le fichier log
#######################
# Get current working directory
rep_exec = os.getcwd()

resonances = ["3:2", "4:3", "5:4", "6:5", "7:6", "8:7", "9:8", "10:9", "11:10"]

#######################
# On prépare les plots
#######################
OUTPUT_EXTENSION = 'pdf' # default value



#######################
# On définit et prépare les différents dossiers 
# où on doit lire les résultats des simulations
#######################
# On définit une liste de nom de dossier que, quoi qu'il arrive
# on ne souhaite pas traiter et donc qu'on vire de la liste.
dossier_suppr = ["output"]

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

# We go in each sub folder of the current working directory

# On récupère la liste des sous-dossiers
liste_simu = [dir for dir in os.listdir(".") if os.path.isdir(dir)]
autiwa.suppr_dossier(liste_simu,dossier_suppr)
liste_simu.sort()


nb_simu = len(liste_simu)

# We initialize the variables where to store datas

# List of orbital elements of ALL the simulations
a = [] # in AU
e = [] # eccentricity
I = [] # inclination in degre??
m = [] # mass in earth mass
m_relat = [] # mass expressed in function of the mass of the most massive planet of the system

# List of orbital elements of the closest planet from the convergence zone
a_clo = [] # in AU
e_clo = [] # eccentricity
I_clo = [] # inclination in degre??
m_clo = [] # mass in earth mass

# List of period ratios for each planet with the reference being, in each system, the closest planet from the resonance.
period_ratio = []

most_massive = [] # the most massive planet of the system
second_massive = [] # the second most massive planet of the system

# Lists for masses of the first and second planets of resonances
m_first = []  # All the coorbitals
m_second = [] # All the coorbitals
m_clo_first = []  # The coorbitals with the closest planet from the convergence zone
m_clo_second = [] # The coorbitals with the closest planet from the convergence zone

final_nb_planets = [] # the final number of planets in the system

print("\t Reading datas")

for simu in liste_simu:
  os.chdir(os.path.join(rep_exec, simu))
  
  # We do nothing of there is no 'big.in' file
  if not(os.path.isfile("big.in")):
    continue
  
  dataFile = open("big.in", 'r')
  

  
  # List initialized for each new simulation that contains values of the current simulation
  a_system = [] # in AU
  e_system = [] # eccentricity
  I_system = [] # inclination in degre??
  m_system = [] # mass in earth mass
  
  # We get rid of the header
  for i in range(6):
    dataFile.readline()
  
  lines = dataFile.readlines()
  dataFile.close()
  
  nb_planets = len(lines) / 2
  final_nb_planets.append(nb_planets)
  
  for line in lines[1::2]: # only odd lines
    datas = line.split()
    a_system.append(float(datas[0]))
    #~ e_system.append(float(datas[2]))
    #~ I_system.append(float(datas[3]))
    #~ m_system.append(float(datas[4]))
  #~ 
  # We add all the elements of the system into the list of values corresponding to all the simulations
  a.extend(a_system)
  #~ e.extend(e_system)
  #~ I.extend(I_system)
  #~ m.extend(m_system)
  
  #~ m_max = max(m_system)
  #~ m_relat.extend([min((mi + random.uniform(-0.5,0.5))/m_max,1.) for mi in m_system])
  #~ 
  #~ ###########
  #~ # Statistical studies
  #~ 
  #~ # We search for the closest planet from the convergence zone
  #~ a_system = np.array(a_system)
  #~ a_temp = abs(a_system - CONVERGENCE_ZONE)
  #~ idx_clo = a_temp.argmin()
  #~ a_ref = a_system[idx_clo]
  #~ 
  #~ a_clo.append(a_system[idx_clo])
  #~ e_clo.append(e_system[idx_clo])
  #~ I_clo.append(I_system[idx_clo])
  #~ m_clo.append(m_system[idx_clo])
  #~ 
  
  #~ # We search for all the previous planets that are in coorbit with the reference one
  #~ idx_before = idx_clo # 
  #~ while (idx_before > 0 and (a_system[idx_before-1]/a_ref)**1.5>(1-DELTA_RATIO)):
    #~ idx_before -= 1
  #~ 
  #~ # We search for all the outer planets that are in coorbit with the reference one
  #~ idx_after = idx_clo # 
  #~ while (idx_after < nb_planets-1 and (a_system[idx_after+1]/a_ref)**1.5<(1+DELTA_RATIO)):
    #~ idx_after += 1
    #~ 
  #~ 
  #~ if (idx_before > 0):
    #~ idx_before -= 1
  #~ 
  #~ if (idx_after < nb_planets-1):
    #~ idx_after += 1
  #~ 
  #~ # We search for the two most massives planets of a system
  #~ tmp = list(m_system)
  #~ tmp.sort()
  #~ most_massive.append(tmp[-1])
  #~ second_massive.append(tmp[-2])
  #~ 
  #~ 
  #~ 
  #~ period_ratio.extend((a_system[idx_before:idx_clo]/a_ref)**(1.5))# We do not add 1 to the last index to exclude the central planet
  #~ period_ratio.extend((a_system[idx_clo+1:idx_after+1]/a_ref)**(1.5))# We need to add 1 the the outer index because the last index is excluded
  #~ 
  #~ 
  #~ 
  #~ i = 0
  #~ while (i<nb_planets):
    #~ a_ref = a_system[i]
    #~ m_ref = m_system[i]
    #~ i += 1
    #~ while (i<nb_planets):# We search for all the coorbitals with the current reference planet
      #~ tmp = (a_system[i] / a_ref)**1.5
      #~ res = (tmp if tmp > 1 else 1/tmp)
      #~ 
      #~ if (res<(1+DELTA_RATIO)):
        #~ if (m_ref>=m_system[i]):
          #~ m_min = m_system[i]
          #~ m_max = m_ref
        #~ else:
          #~ m_max = m_system[i]
          #~ m_min = m_ref
        #~ 
        #~ # We do not store in the same array if the coorbitals are close or not from the convergence zone
        #~ if (1+abs(1.-a_ref/a_system[idx_clo]) < (1+DELTA_RATIO)):
          #~ m_clo_first.append(m_max)
          #~ m_clo_second.append(m_min)
          #~ 
        #~ else:
          #~ m_first.append(m_max)
          #~ m_second.append(m_min)
        #~ i += 1
      #~ else:
        #~ break
      #~ 
  #~ ##########
  #~ 


#######################
#   Tracé des plots   #
#######################
print("\t Computing Plots")
os.chdir(rep_exec)
nb_bins = 50

fig = pl.figure()
hist_a = fig.add_subplot(1, 1, 1)
hist_a.set_xlabel("Semi-major axis [AU]")
hist_a.set_ylabel("Distribution")
hist_a.hist(a, 50, histtype='step')




print("Storing plots")
nom_fichier_plot = "histogram_initial_a"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)  
  
pl.show()

# TODO tracer la zone de convergence et les résonnances sur les plots
