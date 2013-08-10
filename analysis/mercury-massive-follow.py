#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version 1.2

# The script will display in a a-m plot the evolution of the most massive planets, in surimpression, the migration map.

import pdb # Pour le debug
import numpy as np
#~ from fractions import Fraction
import pylab as pl
import autiwa
import sys # to get access to arguments of the script
#~ import shutil, os
import os
from constants import MT, MS
import mercury_utilities

################
## Parameters ##
################
# We get the path toward the binaries
scriptFolder = os.path.dirname(os.path.realpath(__file__)) # the folder in which the module is. 
binaryPath = os.path.join(scriptFolder, os.path.pardir)

NOM_FICHIER_PLOT = "massive_follow"
OUTPUT_EXTENSION = "pdf"

# Maximum number of planets (the most massives) that will be colored
MAX_COLORED = 3
###############################################
## Beginning of the program
###############################################

# We get arguments from the script

#~ pdb.set_trace()
isProblem = False
isDisk = True
isLogX = False
problem_message = "AIM : Display in a m = f(a) diagram, the most massive planets and their evolution in m=f(a)" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * nodisk : avoid torque diagram display" + "\n" + \
" * log : display distances in log" + "\n" + \
" * massive : (%d) the number of most massive planets to be tracked" % MAX_COLORED + "\n" + \
" * ext=png (The extension for the output files)" + "\n" + \
" * help : display this current message"

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'nodisk'):
    isDisk = False
  elif (key == 'log'):
    isLogX = True
  elif (key == 'massive'):
    MAX_COLORED = int(value)
  elif (key == 'help'):
    print(problem_message)
    exit()
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)

if isDisk:
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

# We prepare the plots
fig = pl.figure()
# On crée des sous plots. Pour subplot(231), ça signifie qu'on a 2 lignes, 3 colonnes, et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
plot_AM = fig.add_subplot(1, 1, 1)
plot_AM.set_xlabel("Semi-major axis [AU]")
plot_AM.set_ylabel("mass [Earths]")
plot_AM.xaxis.grid(True, which='major', color='#cccccc', linestyle='-')
plot_AM.yaxis.grid(True, which='major', color='#cccccc', linestyle='-')


# To each planet named "planet01" exist a "planet01.aei" file that contains its coordinates with time.
name = [] # the name of the planets in the current simulation. 
mass = [] # mass in earth mass

table = open("element.out", 'r')

# We get rid of the header
for i in range(5):
  table.readline()

lines = table.readlines()
table.close()

nb_planets = len(lines)
for line in lines:
  datas = line.split()
  name.append(datas[0])
  mass.append(float(datas[4]))

mass = np.array(mass)
name = np.array(name)
sorted_by_mass = np.argsort(mass)
name = name[sorted_by_mass][::-1] # most massive first

# We only keep the most massive ones
name = ["%s.aei" % ni for ni in name[0:MAX_COLORED]]

nb_planets = len(name)

a = [] # Semi-major axis in AU
m = [] # planet mass in earth mass

aappend = a.append
mappend = m.append

for aei_file in name:
  (ai, mi) = np.loadtxt(aei_file, skiprows=4, usecols = (1,7), dtype=float, unpack=True)
  mi = (MS / MT) * mi
  
  aappend(ai)
  mappend(mi)


# We generate a list of colors
colors = [ '#'+li for li in autiwa.colorList(nb_planets)]

if isLogX:
  plot = plot_AM.semilogx
  
  plot_AM.xaxis.grid(True, which='minor', color='#cccccc', linestyle='--')

else:
  plot = plot_AM.plot

for planet in range(nb_planets):
  # We display a circle for the planet
  plot(a[planet], m[planet], color=colors[planet])
  plot(a[planet][-1], m[planet][-1], 'o', color=colors[planet])

ylims = list(plot_AM.get_ylim())
xlims = list(plot_AM.get_xlim())

if (isDisk and (len(contours_a) > 0)):
  plot_AM.fill(contours_a[0], contours_m[0], facecolor="#ffb2b2", edgecolor='#000000', label="Outward migration")
  for (c_a, c_m) in zip(contours_a[1:], contours_m[1:]):
    plot_AM.fill(c_a, c_m, facecolor="#ffb2b2", edgecolor='#000000')

#~ plot_AM.set_xlim(xlims)
#~ plot_AM.set_ylim(ylims)

plot_AM.legend()
fig.savefig("%s.%s" % (NOM_FICHIER_PLOT, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

# TODO
# add an option to plot resonances at a different time of the evolution.


