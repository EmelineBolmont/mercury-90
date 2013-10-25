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
from matplotlib.ticker import FormatStrFormatter

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
COLOR_FINAL = '000000'
COLOR_EJECTED = '909090'
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
" * massive=%d : the number of most massive planets to be tracked" % MAX_COLORED + "\n" + \
" * ext=%s : The extension for the output files" % OUTPUT_EXTENSION + "\n" + \
" * help : display a little help message on HOW to use various options"


value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'nodisk'):
    isDisk = False
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'log'):
    isLogX = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'massive'):
    MAX_COLORED = int(value)
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

if isDisk:
  (process_stdout, process_stderr, returncode) = autiwa.lancer_commande(os.path.join(binaryPath, "migration_map"))
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
# We create subplots. add_subplot(2, 3, 1) means we have 2 lines, 3 columns, 
# and that the active plot is the first, starting from top left (for 6 plots in total)
plot_AM = fig.add_subplot(1, 1, 1)
plot_AM.set_xlabel("Semi-major axis [AU]")
plot_AM.set_ylabel("mass [Earths]")
plot_AM.xaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_AM.yaxis.grid(True, which='major', color='#222222', linestyle='-')


######################################
####################
# On rÃ©cupÃ¨re la liste des fichiers planÃ¨tes.aei
####################
(process_stdout, process_stderr, return_code) = autiwa.lancer_commande("ls *.aei")
if (return_code != 0):
  print("the command return an error "+str(return_code))
  print(process_stderr)
  exit()
  
liste_aei = process_stdout.split("\n")
liste_aei.remove('') # we remove an extra element that doesn't mean anything
nb_planets = len(liste_aei)


####################
# On lit, pour chaque planete, le contenu du fichier et on stocke les variables qui nous intÃÂ©ressent.
####################
t = [] # time in years
a = [] # Semi-major axis in AU
m = [] # planet mass in earth mass

# We retrieve the orbital data
for planete in range(nb_planets):
  
  fichier_source = liste_aei[planete]
  (ti, ai, mi) = np.loadtxt(fichier_source, skiprows=4, usecols = (0,1,7), dtype=float, unpack=True)
  mi = (MS / MT) * mi
  
  if (type(ti) == np.ndarray):
    t.append(ti/1e6)
    a.append(ai)
    m.append(mi)
  else:
    # In case the is only one point, we force to have a list, to avoid plotting problems
    t.append(np.array([ti/1e6]))
    a.append(np.array([ai]))
    m.append(np.array([mi]))


# We get the array of reference time, i.e, one of the longuest list of time available in the list of planets. 
len_t = [len(ti) for ti in t]
ref_len = max(len_t)
ref_id = len_t.index(ref_len) # The ID of the longuest time array
ref_time = t[ref_id] # The longuest time array

delta_t = ref_time[1] - ref_time[0]

# We get the index for the t_max value
if ('t_max' in locals()):
  id_max = int((t_max - ref_time[0]) / delta_t)
  t_max = ref_time[id_max]
else:
  id_max = ref_len - 1
  t_max = ref_time[-1]

# We get the index for the t_min value
if ('t_min' in locals()):
  id_min = int((t_min - ref_time[0]) / delta_t)
  t_min = ref_time[id_min]
else:
  id_min = 0
  t_min = ref_time[0]

# We get the list of masses at the 't_max' time. By default, it will be the final time of the simulation.
final_name = []
final_mass = []
for (aei_file,mi) in zip(liste_aei,m):
  if (len(mi) > id_max):
    final_name.append(aei_file.rstrip(".aei"))
    final_mass.append(mi[id_max])

nb_final = len(final_name)
decorated = [(final_mass[i], i, final_name[i]) for i in range(nb_final)]
decorated.sort(reverse=True) # We get the final planets, the bigger the firsts

for (i , (mass, dumb,  name)) in enumerate(decorated):
  final_name[i] = name
  final_mass[i] = mass


# We get all the association of name/index of planets
ind_of_planet = {}
for (ind, filename) in enumerate(liste_aei):
  basename = os.path.splitext(filename)[0]
  ind_of_planet[basename] = ind

# We generate a list of colors
tmp = autiwa.colorList(MAX_COLORED, exclude=['ffffff', COLOR_FINAL, COLOR_EJECTED])
tmp.extend([COLOR_FINAL]*(nb_final-MAX_COLORED))

# We initialize the total color list. 
# By default, all the planets have the same colors, that will illustrate ejection and collisions with the sun
colors = ['#'+COLOR_EJECTED] * nb_planets

# For the planets that remains in the simulation at the end, we change the color (either a flashy one for the most massives, or a default one)
for (name, color) in zip(final_name, tmp):
  colors[ind_of_planet[name]] = '#'+color



tableau = open('info.out', 'r')
# We must read the file backward because at the beginning, we only have the color of a few planets.
lines = tableau.readlines()
tableau.close()

# We retrieve the history of collision to give the same colors to all the bodies that collided with the remaining bodies
# We start from the last collisions, because the color of the lost body will be the same as the color of the remaining one. 
#  If only one body is left, we want to see each of the bodies that collides with him to get his color.
lost_in_collisions = [] 
for line in reversed(lines):
  if (line.count('was hit by') > 0):
    words = line.split()
    collision_time = float(words[-2]) / 1e6
    if ((collision_time > t_max) or (collision_time < t_min)):
      continue
    remaining_planet = words[0]
    lost_planet = words[4]
    colors[ind_of_planet[lost_planet]] = colors[ind_of_planet[remaining_planet]]
    lost_in_collisions.append(ind_of_planet[lost_planet])
#######################################

if isLogX:
  plot = plot_AM.semilogx
  
  plot_AM.xaxis.grid(True, which='minor', color='#888888', linestyle='-')

else:
  plot = plot_AM.plot

for planet in range(nb_planets):
  # We display a circle for the planet
  plot(a[planet], m[planet], color=colors[planet])


# We display surviving planets
for name in final_name:
  planet = ind_of_planet[name]
  plot(a[planet][-1], m[planet][-1], 'o', color=colors[planet])

for planet in lost_in_collisions:
  plot(a[planet][-1], m[planet][-1], 'o', markerfacecolor='None', markeredgewidth=2, markeredgecolor=colors[planet])

ylims = list(plot_AM.get_ylim())
xlims = list(plot_AM.get_xlim())

if (isDisk and (len(contours_a) > 0)):
  plot_AM.fill(contours_a[0], contours_m[0], facecolor="#ffb2b2", edgecolor='#000000', label="Outward migration")
  for (c_a, c_m) in zip(contours_a[1:], contours_m[1:]):
    plot_AM.fill(c_a, c_m, facecolor="#ffb2b2", edgecolor='#000000')

#~ plot_AM.set_xlim(xlims)
#~ plot_AM.set_ylim(ylims)

distance_format = FormatStrFormatter("%.3g")
plot_AM.xaxis.set_major_formatter(distance_format)

plot_AM.legend()
fig.savefig("%s.%s" % (NOM_FICHIER_PLOT, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

# TODO
# add an option to plot resonances at a different time of the evolution.


