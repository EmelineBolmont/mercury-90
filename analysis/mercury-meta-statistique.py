#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0.1.1
# Goal : Generate several histograms over several folders that each represent meta-simulations that contains several simulation following a common rule. 
# The name of the meta-simulation will be assumed to be the name of the folder.

from math import *
import matplotlib.pyplot as pl
import os, pdb, autiwa
import sys # To handle options of the script, for instance
from simu_constantes import *
import random
import numpy as np
from matplotlib.ticker import FormatStrFormatter

# Mass threshold
MASS_THRESHOLD = 5 # Earth mass

DELTA_RATIO = 0.005
DELTA_M = 0.1

# Get current working directory
rep_exec = os.getcwd()

resonances = ["3:2", "4:3", "5:4", "6:5", "7:6", "8:7", "9:8", "10:9", "11:10"]

#######################
# We prepare the plots
#######################

# Extensions voulues pour les fichiers de sortie
OUTPUT_EXTENSION = 'png' # default value in bitmap, because vectoriel can take time and space if there is a lot of data

nom_fichier_plot = [] # list of names for each plot
figures = [] # list of figures
plots = [] # list of plots (the axes of the fig objects, assuming there is no subplots in there)

distance_format = FormatStrFormatter("%.3g")

nom_fichier_plot.append("histogrammes_m")
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("Mass [Earths]")
fig.set_ylabel("Distribution")

nom_fichier_plot.append("e_fct_m")
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("Mass [Earths]")
fig.set_ylabel("Eccentricity")

nom_fichier_plot.append("e_fct_a")
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("Distance [AU]")
fig.set_ylabel("Eccentricity")
fig.xaxis.set_major_formatter(distance_format)

nom_fichier_plot.append('histogrammes_I')
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("Inclination (in degrees)")
fig.set_ylabel("Distribution")

nom_fichier_plot.append('histogrammes_nb_pl')
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("Final Number of planets")
fig.set_ylabel("Distribution")

nom_fichier_plot.append("m_fct_a")
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("Distance [AU]")
fig.set_ylabel("Mass [Earths]")
fig.set_xscale("log")
fig.xaxis.set_major_formatter(distance_format)

nom_fichier_plot.append("histogrammes_res")
figures.append(pl.figure())
hist_res = figures[-1].add_subplot(1, 1, 1)
plots.append(hist_res)
hist_res.set_xlabel("Period ratio relative to the most massive planet")
hist_res.set_ylabel("Distribution")

nom_fichier_plot.append("most_massives")
figures.append(pl.figure())
fig = figures[-1].add_subplot(1, 1, 1)
plots.append(fig)
fig.set_xlabel("most massive [Earths]")
fig.set_ylabel("second most massive [Earths]")

#######################
# We prepare the various folders we we will read the simulations
#######################
# We define a list of folder names we don't want to read. 
# Theses names will be suppressed from the final list, in case they exist.
dossier_suppr = ["output", "indiv_simu_01"]

#######################
# We list all the simulations of the current working directory each one of them being a 'meta simulation'
#######################
liste_meta_simu = [dir for dir in os.listdir(".") if os.path.isdir(dir)]
autiwa.suppr_dossier(liste_meta_simu,dossier_suppr)
liste_meta_simu.sort()

nb_meta_simu = len(liste_meta_simu)

# We generate a list of colors
colors = [ '#'+li for li in autiwa.colorList(nb_meta_simu, exclude=['000000'])]

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

# We read options, if any

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * ext=png : [%s] The extension for the output files" % OUTPUT_EXTENSION

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 't_min'):
    t_min = float(value)
  elif (key == 't_max'):
    t_max = float(value)
  elif (key == 'log'):
    isLog = True
  elif (key == 'alog'):
    isaLog = True
  elif (key == 'xs'):
    isXS = True
  elif (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'help'):
    isProblem = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

# We go in each sub folder of the current working directory
for (meta_index, meta_simu) in enumerate(liste_meta_simu):
  print("Traitement de %s"%meta_simu)
  os.chdir(os.path.join(rep_exec,meta_simu))
  
  meta_prefix = meta_simu # it is thought to have the initial number of planets in here
  

  # For each meta-simulation, we get the list of simulation that it contains
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
  
  #~ dataFile = open("element.out", 'w')
  #~ dataFile.write("""
  #~ Time (years):   10000000.0001293
#~ 
              #~ a        e       i      mass    Rot/day  Obl
#~ 
  #~ """)
  #~ dataFile.close()
  
  print("\t Reading datas")
  for simu in liste_simu:
    os.chdir(os.path.join(rep_exec, meta_simu, simu))
    dataFile = open("element.out", 'r')
    
    # List initialized for each new simulation that contains values of the current simulation
    a_system = [] # in AU
    e_system = [] # eccentricity
    I_system = [] # inclination in degre
    m_system = [] # mass in earth mass
    
    # We get rid of the header
    for i in range(5):
      dataFile.readline()
    
    lines = dataFile.readlines()
    dataFile.close()
    
    nb_planets = len(lines)
    final_nb_planets.append(nb_planets)
    
    for line in lines:
      datas = line.split()
      a_system.append(float(datas[1]))
      e_system.append(float(datas[2]))
      I_system.append(float(datas[3]))
      m_system.append(float(datas[4]))
    
    # We add all the elements of the system into the list of values corresponding to all the simulations
    a.extend(a_system)
    e.extend(e_system)
    I.extend(I_system)
    m.extend(m_system)
    
    #~ dataFile = open("../element.out", 'a')
    #~ dataFile.write("".join(lines))
    #~ dataFile.close()
   
    if (nb_planets > 0):
      m_max = max(m_system)
      m_relat.extend([min((mi + random.uniform(-0.5,0.5))/m_max,1.) for mi in m_system])
      
      ###########
      # Statistical studies
      
      # We search for the closest planet from the convergence zone
      a_system = np.array(a_system)
      idx_clo = np.argmax(m_system)
      a_temp = abs(a_system - a_system[idx_clo])
      a_ref = a_system[idx_clo]
      
      a_clo.append(a_system[idx_clo])
      e_clo.append(e_system[idx_clo])
      I_clo.append(I_system[idx_clo])
      m_clo.append(m_system[idx_clo])
      
      
      
      # We search for all the previous planets that are in coorbit with the reference one
      idx_before = idx_clo # 
      while (idx_before > 0 and (a_system[idx_before-1]/a_ref)**1.5>(1-DELTA_RATIO)):
        idx_before -= 1
      
      # We search for all the outer planets that are in coorbit with the reference one
      idx_after = idx_clo # 
      while (idx_after < nb_planets-1 and (a_system[idx_after+1]/a_ref)**1.5<(1+DELTA_RATIO)):
        idx_after += 1
        
      #~ idx_before:idx_after is the range of the coorbitals of the reference planet
      
      if (idx_before > 0):
        idx_before -= 1
      
      if (idx_after < nb_planets-1):
        idx_after += 1
      
      # We search for the two most massives planets of a system
      tmp = list(m_system)
      tmp.sort()
      try:
        mass_first = tmp[-1]
        mass_second = tmp[-2]
        most_massive.append(mass_first)
        second_massive.append(mass_second)
        #~ if (abs(mass_first - 16) < DELTA_M):
          #~ print("meta_simu :"+meta_prefix+"\t simu :"+simu)
      except:
        pass
      
      #~ idx_before:idx_after is the range of planets between the first non coorbital inner and outer the position of the reference planet (if they exists)
      
      period_ratio.extend((a_system[idx_before:idx_clo]/a_ref)**(1.5))# We do not add 1 to the last index to exclude the central planet
      period_ratio.extend((a_system[idx_clo+1:idx_after+1]/a_ref)**(1.5))# We need to add 1 the the outer index because the last index is excluded
      
      
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
    else:
      print("%d planets in %s/%s" % (nb_planets, meta_simu, simu))
  
  
  #######################
  #   Plots
  #######################
  print("\t Computing Plots")
  os.chdir(rep_exec)
  nb_bins = 50

  plots[0].hist(m_clo, bins=range(25), normed=True, histtype='step', label=meta_prefix)

  plots[1].plot(m, e, 'o', markersize=5, label=meta_prefix)
  
  plots[2].plot(a, e, 'o', markersize=5, label=meta_prefix)

  plots[3].hist(I, bins=[0.002*i for i in range(25)], normed=True, histtype='step', label=meta_prefix)
  
  plots[4].hist(final_nb_planets, bins=range(25), histtype='step', label=meta_prefix)
  
  plots[5].plot(a, m, 'o', markersize=5, color=colors[meta_index], label=meta_prefix)
  
  plots[6].hist(period_ratio, bins=[0.5+0.0025*i for i in range(400)], normed=True, histtype='step', label=meta_prefix)
  
  plots[7].plot(most_massive, second_massive, 'o', markersize=5, label=meta_prefix)
  

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
for (name, fig, plot) in zip(nom_fichier_plot, figures, plots):
  plot.legend()
  fig.savefig("%s.%s" % (name,OUTPUT_EXTENSION))
  


pl.show()

