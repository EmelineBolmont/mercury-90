#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.1
# Script that will display the evolution of the Semi-major axis, mass, 
# eccentricity and inclination for all the planets of the current simulation 
# (You launch the script in the folder of the mercury simulation)

from math import *
import os, pdb, autiwa
import numpy as np
from constants import MT, MS, RADTODEG
import sys # to be able to retrieve arguments of the script
import pylab as pl
#~ import matplotlib.pyplot as pl
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from analysis import get_x_s
import mercury_utilities


###############################################
## Beginning of the program
###############################################
isaLog = False # If true, the distance axis will be in log
isLog = False # We set the false option before. Because if not, we will erase the 'true' with other option that are not log, and 
# thus will lead to be in the else and put log to false.
isXS = False # If true, will display the semiwitdh of the horseshoe region in the eccentricity plot
OUTPUT_EXTENSION = 'png' # default value in bitmap, because vectoriel can take time and space if there is a lot of data
isAM = False # If activated, only display semi major axis and masses

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * t_max : the end of the output (in years)" + "\n" + \
" * t_min : the beginning of the output (in years)" + "\n" + \
" * log : time will be displayed in log" + "\n" + \
" * alog : distance will be displayed in log" + "\n" + \
" * xs : will display the semiwitdh of the horseshoe region in the eccentricity plot" + "\n" + \
" * am : will only display semi-major axis and masses" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * ext=png : [%s] The extension for the output files" % OUTPUT_EXTENSION

isValue = False # Set to True if someone define 'key=value' when key doesn't need a parameter
value_message = "/!\ Warning: You defined a value ('key=value') to some parameters that to do need them ('key')"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 't_min'):
    t_min = float(value) / 1e6
  elif (key == 't_max'):
    t_max = float(value) / 1e6
  elif (key == 'log'):
    isLog = True
    if (value != None):
      print("/!\ Warning: %s doesn't need any value.")
  elif (key == 'alog'):
    isaLog = True
    if (value != None):
      print("/!\ Warning: %s doesn't need any value.")
  elif (key == 'xs'):
    isXS = True
    if (value != None):
      print("/!\ Warning: %s doesn't need any value.")
  elif (key == 'am'):
    isAM = True
    if (value != None):
      print("/!\ Warning: %s doesn't need any value.")
  elif (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'help'):
    isProblem = True
    if (value != None):
      print("/!\ Warning: %s doesn't need any value.")
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

if isValue:
  print(value_message)

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
# On lit, pour chaque planete, le contenu du fichier et on stocke les variables qui nous intÃ©ressent.
####################
t = [] # time in years
a = [] # Semi-major axis in AU
e = [] # eccentricity
q = [] # perihelion
Q = [] # aphelion
I = [] # inclinaison (degrees)
m = [] # planet mass in earth mass

i = 0
error = True
positions = []
while (len(positions) == 0):
  try:
    positions = mercury_utilities.get_column_position(liste_aei[i])
  except:
    i += 1

# We retrieve the orbital data
for planete in range(nb_planete):
  
  ti = [] # time in years
  ai = [] # Semi-major axis in AU
  ei = [] # eccentricity
  Ii = [] # inclinaison (degrees)
  mi = [] # planet mass in earth mass
  
  fichier_source = liste_aei[planete]
  object_file = open(fichier_source, 'r')
  
  for i in range(4):
    object_file.readline()
  
  tiapp = ti.append
  aiapp = ai.append
  eiapp = ei.append
  Iiapp = Ii.append
  miapp = mi.append

  for line in object_file:
    # When using [a:b], the actual range will be from a to b-1 included.
    tiapp(float(line[positions[0]:positions[1]]) / 1e6)
    aiapp(float(line[positions[1]:positions[2]]))
    eiapp(float(line[positions[2]:positions[3]]))
    Iiapp(float(line[positions[3]:positions[4]]))
    miapp(float(line[positions[7]:positions[8]]))

  object_file.close()
  
  ti = np.array(ti)
  ai = np.array(ai)
  ei = np.array(ei)
  Ii = np.array(Ii)
  mi = np.array(mi)
  
  qi = ai * (1 - ei)
  Qi = ai * (1 + ei)
  mi = (MS / MT) * mi
  
  t.append(ti)
  a.append(ai)
  e.append(ei)
  q.append(qi)
  Q.append(Qi)
  I.append(Ii)
  m.append(mi)
  

if (isXS):
  x_s = []
  for mi in m:
    x_s.append(get_x_s((MT/MS) * mi))

# We get the array of reference time, i.e, one of the longuest list of time available in the list of planets. 
len_t = [ti.size for ti in t]
ref_id = np.argmax(len_t) # The ID of the longuest time array
ref_time = t[ref_id] # The longuest time array
ref_len = ref_time.size

#~ pdb.set_trace()

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

# We generate a list of colors
colors = [ '#'+li for li in autiwa.colorList(nb_planete)]

# We display plots
if isAM:
  nb_rows = 1
else:
  nb_rows = 2

fig = pl.figure()
pl.clf()
fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)
# We create subplots. add_subplot(2, 3, 1) means we have 2 lines, 3 columns, 
# and that the active plot is the first, starting from top left (for 6 plots in total)
plot_a = fig.add_subplot(2, nb_rows, 1)

if (isaLog and isLog):
  plot = plot_a.loglog
elif (isaLog and not(isLog)):
  plot = plot_a.semilogy
elif (not(isaLog) and isLog):
  plot = plot_a.semilogx
else:
  plot = plot_a.plot

for planet in range(nb_planete):
  plot(t[planet][id_min:id_max+1], a[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
  plot(t[planet][id_min:id_max+1], q[planet][id_min:id_max+1], color=colors[planet])
  plot(t[planet][id_min:id_max+1], Q[planet][id_min:id_max+1], color=colors[planet])

plot_a.set_xlabel("Time [million years]")
plot_a.set_ylabel("Semi-major axis [AU]")
plot_a.grid(True)

plot_m = fig.add_subplot(2, nb_rows, nb_rows + 1, sharex=plot_a)
if isLog:
  plot = plot_m.semilogx
else:
  plot = plot_m.plot
  
for planet in range(nb_planete):
  plot(t[planet][id_min:id_max+1], m[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
plot_m.set_xlabel("Time [million years]")
plot_m.set_ylabel("Mass [Earths]")
plot_m.grid(True)

if not(isAM):
  plot_e = fig.add_subplot(2, 2, 2, sharex=plot_a)
  if isLog:
    plot = plot_e.loglog
  else:
    plot = plot_e.semilogy

  for planet in range(nb_planete):
    plot(t[planet][id_min:id_max+1], e[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))

  if isXS:
    for planet in range(nb_planete):
      print("At the end, %s has x_s=%f" % (liste_aei[planet], x_s[planet][id_max]))
      plot(t[planet][id_min:id_max+1], x_s[planet][id_min:id_max+1], color=colors[planet])


  plot_e.set_xlabel("Time [million years]")
  plot_e.set_ylabel("Eccentricity")
  plot_e.grid(True)
  
  plot_I = fig.add_subplot(2, 2, 4, sharex=plot_a)
  if isLog:
    plot = plot_I.semilogx
  else:
    plot = plot_I.plot

  for planet in range(nb_planete):
    plot(t[planet][id_min:id_max+1], I[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
  plot_I.set_xlabel("Time [million years]")
  plot_I.set_ylabel("Inclination [degrees]")
  plot_I.grid(True)

  myyfmt = ScalarFormatter(useOffset=False)
  myyfmt.set_scientific(True)
  myyfmt.set_powerlimits((-2, 3)) 
  plot_I.yaxis.set_major_formatter(myyfmt)

#~ myxfmt = ScalarFormatter(useOffset=True)
#~ myxfmt._set_offset(1e5)
#~ myxfmt.set_scientific(True)
#~ myxfmt.set_powerlimits((-3, 3)) 
myxfmt = FormatStrFormatter("%.3g")

afmt = FormatStrFormatter("%.3g")

plot_a.xaxis.set_major_formatter(myxfmt)
plot_a.yaxis.set_major_formatter(afmt)

plot_a.autoscale(axis='x', tight=True)

nom_fichier_plot = "evolution_planete"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

