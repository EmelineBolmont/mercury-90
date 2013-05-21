#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0.1.1
# Pour générer des histogrammes sur plusieurs simulations afin de 
# regarder les caractéristiques statistiques des planètes qui restent.


# TODO
# _ écrire un fichier qui stocke les paramètres orbitaux des planètes qu'il reste.
# _ sélectionner les planètes qui sont plus proche qu'une certaine valeur de leur étoile.

from math import *
#~ import pylab as pl
import matplotlib.pyplot  as pl
import os, pdb, autiwa
from simu_constantes import *
import random
import numpy as np
import sys
import shutil

DELTA_RATIO = 0.005
DELTA_M = 0.1

PREFIX = "simu" # the prefix for the directories we want to take into the statistic.
t_max = None # the integration time of the simulations. Assumed to be the final time of the first simulation read

isDisk = True

# We get the path toward the binaries
scriptFolder = os.path.dirname(os.path.realpath(__file__)) # the folder in which the module is. 
binaryPath = os.path.join(scriptFolder, os.path.pardir)

#######################
# On prépare le fichier log
#######################
# Get current working directory
rep_exec = os.getcwd()

resonances = ["2:1", "3:2", "4:3", "5:4", "6:5", "7:6", "8:7", "9:8", "10:9", "11:10"]

#######################
# We prepare the plots
#######################

# Extensions voulues pour les fichiers de sortie
OUTPUT_EXTENSION = 'png' # default value in bitmap, because vectoriel can take time and space if there is a lot of data

nom_fichier_plot = [] # list of names for each plot
figures = [] # list of figures
plots = [] # list of plots (the axes of the fig objects, assuming there is no subplots in there)

nom_fichier_plot.append("miscellaneous")
figures.append(pl.figure())
fig = figures[-1].add_subplot(2, 3, 1)
plots.append([]) # Because it is a subplots figure
plots[-1].append(fig)
fig.set_xlabel("mass [Earths]")
fig.set_ylabel("Distribution")

fig = figures[-1].add_subplot(2, 3, 2)
plots[-1].append(fig)
fig.set_xlabel("mass [Earths]")
fig.set_ylabel("eccentricity")

fig = figures[-1].add_subplot(2, 3, 3)
plots[-1].append(fig)
fig.set_xlabel("distance (in AU)")
fig.set_ylabel("eccentricity")

fig = figures[-1].add_subplot(2, 3, 4)
plots[-1].append(fig)
fig.set_xlabel("I (in degrees)")
fig.set_ylabel("Distribution")

fig = figures[-1].add_subplot(2, 3, 5)
plots[-1].append(fig)
fig.set_xlabel("nb_final")
fig.set_ylabel("Distribution")

fig = figures[-1].add_subplot(2, 3, 6)
plots[-1].append(fig)
fig.set_xlabel("Orbital Distance [AU]")
fig.set_ylabel("Distribution")
fig.set_xscale("log")

nom_fichier_plot.append("m_fct_a")
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("a (in AU)")
fig.set_ylabel("mass [Earths]")
fig.grid()

nom_fichier_plot.append("histogrammes_res")
figures.append(pl.figure())
hist_res = figures[-1].add_subplot(1, 1, 1)
plots.append(hist_res)
hist_res.set_xlabel("Period ratio relative to the most massive planet")
hist_res.set_ylabel("Distribution")


nom_fichier_plot.append("")
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("most massive [Earths]")
fig.set_ylabel("second most massive [Earths]")


nom_fichier_plot.append("coorbital_pl_mass")
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_title("Mass of the coorbital planets")
fig.set_xlabel("mass of the first planet (m_earth)")
fig.set_ylabel("mass of the second planet (m_earth)")

nom_fichier_plot.append("most_massive_position")
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("Position [AU]")
fig.set_ylabel("Mass of most massive [Earths]")
fig.grid()

#######################
# On définit et prépare les différents dossiers 
# où on doit lire les résultats des simulations
#######################
# We define a list of forbidden names that will not be considered as simulations
dossier_suppr = ["output", "test", "reference_simu"]

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

isLog = False # We set the false option before. Because if not, we will erase the 'true' with other option that are not log, and 
# thus will lead to be in the else and put log to false.
isStat = False # Will create a global element.out file that contains all the surviving planets

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * nodisk : to avoid torque diagram display" + "\n" + \
" * stat :  to create a global 'element.out' file that contain all the surviving planets" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * ext=pdf : The extension for the output files"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'nodisk'):
    isDisk = False
  elif (key == 'stat'):
    isStat = True
  elif (key == 'help'):
    isProblem = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()


# We go in each sub folder of the current working directory

# We get the list of simulations
liste_simu = [dir for dir in os.listdir(".") if os.path.isdir(dir)]
autiwa.suppr_dossier(liste_simu,dossier_suppr)
liste_simu.sort()

nb_simu = len(liste_simu)

# If there is a disk, we get the torque diagram
if isDisk:
  shutil.copy2("%s/disk.in" % liste_simu[0], ".")
  shutil.copy2("%s/param.in" % liste_simu[0], ".")
  (process_stdout, process_stderr, returncode) = autiwa.lancer_commande(os.path.join(binaryPath, "torque_diagram"))
  if (returncode != 0):
    print(process_stdout)
    print(process_stderr)
    print("Process terminated with error %d" % returncode)
    pdb.set_trace()
    
  
  # We read the contour of the zero torque zones
  contour_filename = "contour_total_torque.dat"
  
  contours_a = [[]]
  contours_m = [[]]
  contours_aapp = contours_a[-1].append
  contours_mapp = contours_m[-1].append

  contour_file = open(contour_filename, 'r')
  for line in contour_file:
    if (line == ' \n'):
      contours_a.append([])
      contours_m.append([])
      contours_aapp = contours_a[-1].append
      contours_mapp = contours_m[-1].append
    else:
      (a, m, dummy) = line.split()
      contours_aapp(float(a))
      contours_mapp(float(m))
  
  # We delete all the empty arrays created in the list of contours (one between each contour)
  while True:
    try:
      contours_a.remove([])
      contours_m.remove([])
    except:
      break
  
      
# We initialize the variables where to store datas

# List of orbital elements of ALL the simulations
a = [] # in AU
e = [] # eccentricity
I = [] # inclination in degre
m = [] # mass in earth mass
m_relat = [] # mass expressed in function of the mass of the most massive planet of the system

# List of orbital elements of the closest planet from the convergence zone
a_clo = [] # in AU
e_clo = [] # eccentricity
I_clo = [] # inclination in degre
m_clo = [] # mass in earth mass

# List of period ratios for each planet with the reference being, in each system, the closest planet from the resonance.
period_ratio = []

most_massive = [] # the most massive planet of the system
most_massive_a = [] # the position (AU) of the most massive planet of the system
second_massive = [] # the second most massive planet of the system

# Lists for masses of the first and second planets of resonances
m_first = []  # All the coorbitals
m_second = [] # All the coorbitals
m_clo_first = []  # The coorbitals with the closest planet from the convergence zone
m_clo_second = [] # The coorbitals with the closest planet from the convergence zone

final_nb_planets = [] # the final number of planets in the system
final_time = [] # Final time of the integration. Used to check if there was problems with the simulations

print("\t Reading datas")

if (isStat):
  globalElementOut = open("element.out", 'w')

for simu in liste_simu:
  os.chdir(os.path.join(rep_exec, simu))
  try:
    dataFile = open("element.out", 'r')
  except:
    print("folder "+simu+" does not contains any 'element.out' file")
    continue
  
  # List initialized for each new simulation that contains values of the current simulation
  a_system = [] # in AU
  e_system = [] # eccentricity
  I_system = [] # inclination in degre
  m_system = [] # mass in earth mass
  
  header = []
  
  # We get rid of the header
  for i in range(5):
    header.append(dataFile.readline())
  tmp = header[1]
  tmp = tmp.split(":")[1]
  t_max0 = float(tmp)
  
  if (t_max == None):
    # we set the t_max with the first simulation
    t_max = t_max0
    if isStat:
      globalElementOut.write("".join(header))
  
  if (t_max0 != t_max):
    print("folder "+simu+" output time is "+tmp+" instead of "+str(t_max))
    continue
  
  lines = dataFile.readlines()
  dataFile.close()
  
  if (isStat):
    globalElementOut.write("".join(lines))
  
  nb_planets = len(lines)
  final_nb_planets.append(nb_planets)
  
  for line in lines:
    datas = line.split()
    ai = float(datas[1])
    ei = float(datas[2])
    Ii = float(datas[3])
    mi = float(datas[4])
    a_system.append(ai)
    e_system.append(ei)
    I_system.append(Ii)
    m_system.append(mi)
    #~ if (ai < 23 and ai > 19 and mi > 16 and mi < 20):
      #~ print("in %s a planet has a=%f AU and m=%f mt" % (simu, ai, mi))
    #~ if (ai>2.5) and (mi > 10.):
      #~ print("in %s a planet has a=%f AU and m=%f mt" % (simu, ai, mi))
    #~ if (ai > 1.3):
      #~ print("in %s %s has a=%f and m=%f mt" % (simu, datas[0], ai, mi))
  
  # We add all the elements of the system into the list of values corresponding to all the simulations
  a.extend(a_system)
  e.extend(e_system)
  I.extend(I_system)
  m.extend(m_system)
  
  ###########
  # Statistical studies
  
  if (final_nb_planets[-1]>1):
    # We search for the most massive planet of the system
    a_system = np.array(a_system)
    m_system = np.array(m_system)
    idx_clo = m_system.argmax()
    a_ref = a_system[idx_clo]
    
    a_clo.append(a_system[idx_clo])
    e_clo.append(e_system[idx_clo])
    I_clo.append(I_system[idx_clo])
    m_clo.append(m_system[idx_clo])
  
    
    # We search for all the previous planets that are in coorbit with the reference one
    idx_before = idx_clo
    while (idx_before > 0 and (a_system[idx_before-1]/a_ref)**1.5>(1-DELTA_RATIO)):
      idx_before -= 1
    
    # We search for all the outer planets that are in coorbit with the reference one
    idx_after = idx_clo
    while (idx_after < nb_planets-1 and (a_system[idx_after+1]/a_ref)**1.5<(1+DELTA_RATIO)):
      idx_after += 1
      
    #~ idx_before:idx_after is the range of the coorbitals of the reference planet
    if (idx_before > 0):
      idx_before -= 1
    
    if (idx_after < nb_planets-1):
      idx_after += 1
    
    # idx_before:idx_after is the range of planets between the first non coorbital inner 
    # and outer the position of the reference planet (if they exists)
    period_ratio.extend((a_system[idx_before:idx_clo]/a_ref)**(1.5))# We do not add 1 to the last index to exclude the central planet
    period_ratio.extend((a_system[idx_clo+1:idx_after+1]/a_ref)**(1.5))# We need to add 1 the the outer index because the last index is excluded
    
    # We search for the two most massives planets of a system
    tmp = [(mi, ai) for (ai, mi) in zip(a_system, m_system)]
    tmp.sort()
    most_massive.append(tmp[-1][0])
    most_massive_a.append(tmp[-1][1])
    second_massive.append(tmp[-2][0])

    if (tmp[-1] > 6. and tmp[-1] < 9.):
      print("most massive = %.1f in %s:" % (tmp[-1],simu))
    #~ if (final_nb_planets[-1] == 7):
      #~ print("nb_planets :",final_nb_planets[-1],simu)
    #~ if (max(e_system) > 0.8):
      #~ print("max eccentricity :",max(e_system),simu)
  else:
    print("in %s there is only %i planet left" % (simu, final_nb_planets[-1]))
  
  i = 0
  while (i<nb_planets):
    a_ref = a_system[i]
    m_ref = m_system[i]
    i += 1
    while (i<nb_planets):# We search for all the coorbitals with the current reference planet
      tmp = (a_system[i] / a_ref)**1.5
      res = (tmp if tmp > 1 else 1/tmp)
      
      if (res<(1+DELTA_RATIO)):
        if (m_ref>=m_system[i]):
          m_min = m_system[i]
          m_max = m_ref
        else:
          m_max = m_system[i]
          m_min = m_ref
        
        # We do not store in the same array if the coorbitals are close or not from the convergence zone
        if (1+abs(1.-a_ref/a_system[idx_clo]) < (1+DELTA_RATIO)):
          m_clo_first.append(m_max)
          m_clo_second.append(m_min)
          
        else:
          m_first.append(m_max)
          m_second.append(m_min)
        i += 1
      else:
        break
      
  ##########

if isStat:
  globalElementOut.close()

print("less massive in the list of most massive : %f" % min(most_massive))
print("most massive planet formed in all simulations : %f" % max(most_massive))

#######################
#   Tracé des plots   #
#######################
print("\t Computing Plots")
os.chdir(rep_exec)
nb_bins = 50

plots[0][0].hist(m_clo, bins=range(25), normed=True, histtype='step')

plots[0][1].plot(m, e, 'o', markersize=5)

plots[0][2].plot(a, e, 'o', markersize=5)

plots[0][3].hist(I, bins=[0.002*i for i in range(25)], normed=True, histtype='step')

plots[0][4].hist(final_nb_planets, bins=range(25), histtype='step')

a_max = max(a)
a_min = min(a)
plots[0][5].hist(a, bins=np.logspace(log10(a_min), log10(a_max), 100))


plots[1].semilogx(a, m, 'o', markersize=5)
if (isDisk and (len(contours_a) > 0)):
  plots[1].fill(contours_a[0], contours_m[0], facecolor="#ff0000", alpha=0.3, edgecolor='none', label="Outward migration")
  for (c_a, c_m) in zip(contours_a[1:], contours_m[1:]):
    plots[1].fill(c_a, c_m, facecolor="#ff0000", alpha=0.3, edgecolor='#000000')


plots[2].hist(period_ratio, bins=[0.5+0.0025*i for i in range(400)], normed=True, histtype='step')

plots[3].plot(most_massive, second_massive, 'o', markersize=5)


plots[4].plot(m_clo_first, m_clo_second, 'ro', markersize=5, label="closest from biggest")
plots[4].plot(m_first, m_second, 'bo', markersize=5, label="all others")
ylims = list(plots[4].get_ylim())
xlims = list(plots[4].get_xlim())
limits = [max(xlims[0], ylims[0]), min(xlims[1], ylims[1])]
plots[4].plot(limits, limits, 'k--', label="equal mass")

plots[5].semilogx(most_massive_a, most_massive, 'bo', markersize=5)
if (isDisk and (len(contours_a) > 0)):
  plots[5].fill(contours_a[0], contours_m[0], facecolor="#ff0000", alpha=0.3, edgecolor='none', label="Outward migration")
  for (c_a, c_m) in zip(contours_a[1:], contours_m[1:]):
    plots[5].fill(c_a, c_m, facecolor="#ff0000", alpha=0.3, edgecolor='#000000')

# We add the resonances
ylims = list(hist_res.get_ylim())
xlims = list(hist_res.set_xlim([0.6, 1.4]))
for res in resonances:
  nb_period = map(float, res.split(":")) # We get the two integers value of the resonance.
  ratio = nb_period[0] / nb_period[1]
  hist_res.plot([ratio, ratio], ylims, 'k--')
  hist_res.plot([1./ratio, 1./ratio], ylims, 'k--')
  hist_res.text(ratio, ylims[1], " "+res, horizontalalignment='center', verticalalignment='bottom', rotation='vertical', size=7)
  hist_res.text(1./ratio, ylims[1], " "+res, horizontalalignment='center', verticalalignment='bottom', rotation='vertical', size=7)

print("Storing plots")
for (name, fig) in zip(nom_fichier_plot, figures):
  fig.savefig("%s.%s" % (name,OUTPUT_EXTENSION))

pl.show()