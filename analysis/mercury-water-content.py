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
import mercury
import mercury_utilities
from matplotlib.ticker import FormatStrFormatter
import glob

################
## Parameters ##
################
# We get the path toward the binaries
scriptFolder = os.path.dirname(os.path.realpath(__file__)) # the folder in which the module is. 
binaryPath = os.path.join(scriptFolder, os.path.pardir)

NOM_FICHIER_PLOT = "water_content"
OUTPUT_EXTENSION = "pdf"

###############################################
## Beginning of the program
###############################################

# We get arguments from the script

#~ pdb.set_trace()
isProblem = False
problem_message = "AIM : Display in a m = f(a) diagram, the most massive planets and their evolution in m=f(a)" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * tmax=1e6 : the end of the output [years]" + "\n" + \
" * tmin=1e3 : the beginning of the output [years]" + "\n" + \
" * ext=%s : The extension for the output files" % OUTPUT_EXTENSION + "\n" + \
" * help : display a little help message on HOW to use various options"


value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'tmin'):
    t_min = float(value)
  elif (key == 'tmax'):
    t_max = float(value)
  elif (key == 'ext'):
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

def get_water_color(mass_ratio):
  """From the mass ratio between water mass and total mass, return a color, 
  gradient from totally dry color and totally wet color
  """
  
  #~ DRY_COLOR = (255., 0., 0.) # red
  #~ WET_COLOR = (0., 0., 255.) # blue
  #~ 
  #~ RGB_COLOR = [int(DRY_COLOR[j] + mass_ratio*(WET_COLOR[j]-DRY_COLOR[j])) for j in range(3)]
  #~ HEX_COLOR = "#"+"".join(["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in RGB_COLOR])
  
  if (mass_ratio < 0.05):
    HEX_COLOR = '#000000' # black
  elif (mass_ratio < 0.1):
    HEX_COLOR = '#ff0000' # red
  elif (mass_ratio < 0.2):
    HEX_COLOR = '#ffa500' # orange
  elif (mass_ratio < 0.3):
    HEX_COLOR = '#000000' # yellow
  else:
    HEX_COLOR = '#0000ff' # blue
  
  return HEX_COLOR

# We prepare the plots
fig = pl.figure()
# We create subplots. add_subplot(2, 3, 1) means we have 2 lines, 3 columns, 
# and that the active plot is the first, starting from top left (for 6 plots in total)
plot_AM = fig.add_subplot(1, 1, 1)
plot_AM.set_xlabel("Semi-major axis [AU]")
plot_AM.set_ylabel("Mass [Earths]")
plot_AM.xaxis.grid(True, which='major', color='#cccccc', linestyle='-')
plot_AM.yaxis.grid(True, which='major', color='#cccccc', linestyle='-')


# We list the ASCII output files
aei_files = glob.glob("*.aei")
aei_files.sort()
nb_planets = len(aei_files)

t = []
a = []
m = []
planet_name = []


for filename in aei_files:
  name = os.path.splitext(filename)[0]
  
  (ti, ai, mi) = np.loadtxt(filename, skiprows=4, usecols = (0,1,7), dtype=float, unpack=True)
  mi = (MS / MT) * mi
  
  if (type(ti) == np.ndarray):
    t.append(ti)
    a.append(ai)
    m.append(mi)
  else:
    # In case the is only one point, we force to have a list, to avoid plotting problems
    t.append(np.array([ti]))
    a.append(np.array([ai]))
    m.append(np.array([mi]))
  planet_name.append(name)
  

# We set the initial amount of water., null below 2.5 AU, 5% between 2.5 and 5 AU and 50% above 5 AU
## Initially, we fill all the array with the same amount of water. Only through info.out will we modify the desired items.
water_content = []
for (ai,mi) in zip(a,m):
  if (ai[0] < 2.5):
    water_mass = 0.
  elif (ai[0] < 5.):
    water_mass = 0.05 * mi
  else:
    water_mass = 0.5 * mi
  water_content.append(water_mass)

# For each planet name, we make a dictionnary to know the corresponding planet index
planet_index = {}
for (index, name) in enumerate(planet_name):
  planet_index[name] = index

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

if (t_max > 1e6):
  unit_time = "Million years"
  time_conversion = 1e6
elif (t_max > 1e3):
  unit_time = " 000 years"
  time_conversion = 1e3
else:
  unit_time = "years"
  time_conversion = 1.

tableau = open('info.out', 'r')
# We must read the file backward because at the beginning, we only have the color of a few planets.
lines = tableau.readlines()
tableau.close()

# We retrieve the history of collision to give the same colors to all the bodies that collided with the remaining bodies
# We start from the last collisions, because the color of the lost body will be the same as the color of the remaining one. 
#  If only one body is left, we want to see each of the bodies that collides with him to get his color.
lost_in_collisions = [] 
for line in lines:
  if (line.count('was hit by') > 0):
    words = line.split()
    
    collision_time = float(words[-2]) / 1e6
    if ((collision_time > t_max) or (collision_time < t_min)):
      continue
    
    remaining_planet = words[0]
    remaining_index = planet_index[remaining_planet]
    
    lost_planet = words[4]
    lost_index = planet_index[lost_planet]
    
    # Starting index to modify the remaining planet water content
    tstart_index = len(t[lost_index])
    new_water_mass = water_content[remaining_index][tstart_index -1] + water_content[lost_index][tstart_index - 1]
    
    # We set the new water content right after the collision time for the remaining planet
    water_content[remaining_index][tstart_index:] = new_water_mass

plot = plot_AM.plot

id_time = id_max

min_mass = min([min(mi) for mi in m]) # earth mass
markersize_prefactor = 4 / (min_mass**0.33)
for planet in range(nb_planets):
  try:
    water_mass_ratio = water_content[planet][id_time] / m[planet][id_time]
    plot(a[planet][id_time], m[planet][id_time], 'o', color='#000000', 
         markersize=max(int(markersize_prefactor * (m[planet][id_time])**0.33),markersize_prefactor), 
         markeredgecolor=get_water_color(water_mass_ratio), markeredgewidth=4)
  except:  
    pass

pl.show()

pdb.set_trace()
