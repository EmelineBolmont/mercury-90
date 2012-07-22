#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version 1.0
# The script will calculate the resonances between each planet through time

import pdb # Pour le debug
import numpy as np
from fractions import Fraction
import pylab as pl
import autiwa
import sys # to get access to arguments of the script
import mercury_utilities
import os

################
## Parameters ##
################
NOM_FICHIER_PLOT = "timed_resonances"
OUTPUT_EXTENSION = "png"

DENOMINATOR_LIMIT = 15 # Maximum value allowed of the denominator when we want to get a fraction from a decimal value
NUMBER_OF_VALUES = 100 # sampling for period ratio around the given value
UNCERTAINTY = 5 # In percentage
NB_LAST_POINTS = 15 # Number of points we want to test the libration of angles.

TIME_DELAY = 5 # Interval in index before two tests of the resonance between planets. 

# the threshold of the standard deviation of a resonant angle, 
# below which we consider there is libration and thus, a resonance.
STD_THRESHOLD = 40 

# Extreme angles for angle folding. We will apply a mod 2pi to all angles, 
# but theses values determines between wich values we will constrains the angles.
ANGLE_MIN = - 180.
ANGLE_MAX = 180.

def get_possible_resonances(periodRatio):
  periodMin = periodRatio * (1 - uncertainty)
  periodMax = periodRatio * (1 + uncertainty)
  
  # We do not want period ratios less than 1 (this only happens for coorbitals I think)
  if (periodMin < 1.):
    periodMin = 1.
  
  periodWidth = periodMax - periodMin
  deltaPeriod = periodWidth/NUMBER_OF_VALUES
  
  
  periods = [periodMin + deltaPeriod * i for i in range(NUMBER_OF_VALUES)]

  resonances = []
  for period_i in periods:
    fraction = Fraction(period_i).limit_denominator(DENOMINATOR_LIMIT)
    #print(fraction)
    resonances.append(fraction)

  # We exclude all values that appears several time to only keep one occurence of each value
  resonances = list(set(resonances))

  # We sort the resonances to get the more interesting first (3:2 before 32:27 for instance)
  tmp = [(res.numerator, res) for res in resonances]
  tmp.sort()
  resonances = [element[1] for element in tmp]
  
  return resonances

def isResonance(res, g_inner, n_inner, M_inner, g_outer, n_outer, M_outer):
  """Given a resonance as a Fraction object, and g, n M for inner and
  outer planet, the function return if there is the resonance between
  the two planets"""
  outer_period_nb = res.denominator
  inner_period_nb = res.numerator
  
  # Resonances are usually displayed as (p+q):p where q is the order of
  # the resonance. We retreive thoses parameters
  p = outer_period_nb
  q = inner_period_nb - outer_period_nb

  # We calculate the resonant angles
  long_of_peri_inner = g_inner + n_inner
  mean_longitude_inner = M_inner + long_of_peri_inner

  long_of_peri_outer = g_outer + n_outer
  mean_longitude_outer = M_outer + long_of_peri_outer

  phi = np.empty((q+1, NB_LAST_POINTS)) # we create an array, empty for the moment, that will contain all the resonant angles associated with the supposed resonance.

  temp_value = inner_period_nb * mean_longitude_outer - outer_period_nb * mean_longitude_inner

  for i in range(q+1):
    phi[i] = temp_value - i * long_of_peri_inner - (q - i) * long_of_peri_outer

  # We take modulo 2*pi of the resonant angle
  phi = phi%(360.)
  too_low = phi < ANGLE_MIN
  too_high = phi > ANGLE_MAX
  phi[too_low] = phi[too_low] + 360.
  phi[too_high] = phi[too_high] - 360.

  delta_longitude = long_of_peri_outer - long_of_peri_inner
  delta_longitude = delta_longitude%(360.)
  too_low = delta_longitude < ANGLE_MIN
  too_high = delta_longitude > ANGLE_MAX
  delta_longitude[too_low] = delta_longitude[too_low] + 360.
  delta_longitude[too_high] = delta_longitude[too_high] - 360.
  
  # If one of the std's is small (Typically, around 25, but 
  # I had once a 80 that was not a resonance, so I think a threshold around 40 is a good one)
  standard_deviation = min(phi.std(1))
    
  if (standard_deviation < STD_THRESHOLD):
    return True
  else:
    return False
###############################################
## Beginning of the program
###############################################

# We get arguments from the script
isLog = False
isProblem = False
problem_message = "AIM : Display in a m = f(a) diagram, all the planets of the current mercury simulation" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * log (time will be displayed in log)" + "\n" + \
" * ext=png (The extension for the output files)" + "\n" + \
" * help : display this current message"

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'log'):
    isLog = True
  elif (key == 'help'):
    print(problem_message)
    exit()
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)

####################
# we get the name of the files for all the planets in the system
####################
liste_aei = mercury_utilities.get_aei_files()
nb_planets = len(liste_aei)

# We also initialize the list for the names of the planets
planet_names = []
for filename in liste_aei:
  planet_names.append(os.path.splitext(filename)[0])

####################
# Initialization of some variables and tests
####################
if (TIME_DELAY > NB_LAST_POINTS):
  print("Warning: The interval between two checks of resonance is greater than the number of points used to test it. As the consequence, all the integration time will not be tested")

uncertainty = 0.01 * float(UNCERTAINTY)

# We initialize the arrays
t = [] # time in years
a = [] # demi-grand axe en ua
g = [] # argument of pericentre (degrees)
n = [] # longitude of ascending node (degrees)
M = [] # Mean anomaly (degrees)

# to store the successive resonances for each planet. If a resonance is written at the index i, that mean that the planet i is the inner planet of the resonance given in the array. 
resonance_index_range = [[] for planet in range(nb_planets)] # For each planet, each tuple is a range of time during which the planet is in resonance with the planet just after it.
resonance_type = [[] for planet in range(nb_planets)] # For each planet, store the type of resonance 3:2, 4:3 and so on. The index i of the sublist (of a given planet) correspond to the range in time of the array resonance_time

planet_order = [[] for planet in range(nb_planets)] # a list of the order of the planets, in order to display the resonances nicely
time_order = [[] for planet in range(nb_planets)] # the corresponding time in years
resonance_with = [[] for planet in range(nb_planets)] # the index of the outer planet for the resonance (the inner planet index being the index of the sublist considered)
resonance_inner_rank = [[] for planet in range(nb_planets)] # The rank, in distance, of the inner planet in the resonance (useful to place the resonance in the second plot)

####################
# We read the datas for all the planets
####################
for planet_datafile in liste_aei:
  (ti, ai, gi, ni, Mi) = np.loadtxt(planet_datafile, skiprows=4, usecols=(0,1,4,5,6), dtype=float, unpack=True)
  t.append(ti)
  a.append(ai)
  g.append(gi)
  n.append(ni)
  M.append(Mi)

# We calculate the period
period = []
for ai in a:
  period.append(ai**1.5)

### We create the first point for some arrays
# We sort the planet in the order of distance from the host star
distance = [ai[0] for ai in a]
planet_index_sorted_by_distance = np.argsort(distance) # the i-th closest planet in distance is the planet planet_index_sorted_by_distance[i]
ordering_planets = np.argsort(planet_index_sorted_by_distance) # the i-th planet is the ordering_planets[i] in distance

# we append the ordering of the current planets in order to display resonances later
for (planet_idx, order) in enumerate(ordering_planets):
  planet_order[planet_idx].append(order+1)
  time_order[planet_idx].append(t[planet_idx][0])
###

# We retrieve the longuest array (the planet that stay the longuest in the simulation
lengths = []
for ai in a:
  lengths.append(ai.size)
max_lengths = max(lengths)

for instant_index in range(NB_LAST_POINTS,max_lengths,TIME_DELAY):
  range_start = instant_index - (NB_LAST_POINTS-1)
  range_stop  = instant_index

  # We make sublist only with the planets that still exist at the present time
  still_here_planets = []
  a_temp = []
  g_temp = []
  n_temp = []
  M_temp = []

  for planet in range(nb_planets):
    if (t[planet].size > instant_index):
      still_here_planets.append(planet)
      a_temp.append(a[planet][range_start:range_stop+1])
      g_temp.append(g[planet][range_start:range_stop+1])
      n_temp.append(n[planet][range_start:range_stop+1])
      M_temp.append(M[planet][range_start:range_stop+1])
      
  # We sort the planet in the order of distance from the host star
  distance = [ai[-1] for ai in a_temp]
  planet_index_sorted_by_distance = np.argsort(distance) # the i-th closest planet in distance is the planet planet_index_sorted_by_distance[i]
  ordering_planets = 1 + np.argsort(planet_index_sorted_by_distance) # the i-th planet is the ordering_planets[i] in distance (starting at 1)

  # we append the ordering of the current planets in order to display resonances later
  for (planet_idx, order) in zip(still_here_planets, ordering_planets):
    planet_order[planet_idx].append(order)
    time_order[planet_idx].append(t[planet_idx][instant_index])

  for (inner, outer) in zip(planet_index_sorted_by_distance[0:-1], planet_index_sorted_by_distance[1:]):
    g_inner = g_temp[inner]
    n_inner = n_temp[inner]
    M_inner = M_temp[inner]

    g_outer = g_temp[outer]
    n_outer = n_temp[outer]
    M_outer = M_temp[outer]

    # We get the various possible resonance between the two planets
    periodRatio = (distance[outer] / distance[inner])**1.5
    resonances = get_possible_resonances(periodRatio)

    # For each resonance we check if this one exist between the two considered planets
    index = 0
    while (index < len(resonances)):
      res = resonances[index]
      is_resonance = isResonance(res, g_inner, n_inner, M_inner, g_outer, n_outer, M_outer)
      if (is_resonance):
        isExtend = False # boolean that say if the current resonance is the extension of the last resonance listed for the inner planet
        if (len(resonance_type[inner]) != 0):
          last_type = resonance_type[inner][-1]
          last_index_range = resonance_index_range[inner][-1]
          
          # if the two index ranges overlap
          if ((last_type == res) and (last_index_range[1] >= range_start-1)):
            isExtend = True
        
        # if the current resonance already existed before, we only extend the index range of validity for the last resonance of the inner planet index
        if isExtend:
          resonance_index_range[inner][-1][1] = instant_index
        else:
          resonance_type[inner].append(res)
          resonance_index_range[inner].append([range_start, range_stop])
          resonance_with[inner].append(outer)
          resonance_inner_rank[inner].append(ordering_planets[inner])

        # If we find a resonance, we do not search for another one
        break
      index += 1

####################
# We now want to display in a fashion way the resonances
####################

# We generate a list of colors
colors = [ '#'+li for li in autiwa.colorList(nb_planets)]

fig = pl.figure(1)
pl.clf()
fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)
# On crée des sous plots. Pour subplot(311), ça signifie qu'on a 2 lignes, 3 colonnes, et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
plot_a = fig.add_subplot(211)
for planet in range(nb_planets):
  if isLog:
    plot_a.semilogx(t[planet], a[planet], color=colors[planet], label=planet_names[planet])
  else:
    plot_a.plot(t[planet], a[planet], color=colors[planet], label=planet_names[planet])

ylims = list(pl.ylim())
for planet in range(nb_planets):
  for (outer, (id_begin, id_end)) in zip(resonance_with[planet], resonance_index_range[planet]):
    pl.plot([t[planet][id_begin], t[planet][id_begin]], [a[planet][id_begin], a[outer][id_begin]], 'k--')
    pl.plot([t[planet][id_end], t[planet][id_end]], [a[planet][id_end], a[outer][id_end]], 'k--')


pl.xlabel("time [years]")
pl.ylabel("a [AU]")
pl.grid(True)

plot_res = fig.add_subplot(212, sharex=plot_a)
for planet in range(nb_planets):
  for (res, inner_rank, (id_begin, id_end)) in zip(resonance_type[planet], resonance_inner_rank[planet], resonance_index_range[planet]):
    x_position = (t[planet][id_begin] + t[planet][id_end]) / 2.
    y_position = inner_rank + .5
    
    x_left = t[planet][id_begin]
    x_right = t[planet][id_end]
    y_bottom = inner_rank + 0.1
    y_top = inner_rank + 0.9
    
    text = "%i:%i" % (res.numerator, res.denominator)
    pl.text(x_position, y_position, text, horizontalalignment='center', verticalalignment='center', rotation='vertical')
    plot_res.fill([x_left, x_right, x_right, x_left],[y_bottom, y_bottom, y_top, y_top], fill=True, color='#cccccc')

for planet in range(nb_planets):
  if isLog:
    plot_res.semilogx(time_order[planet], planet_order[planet], color=colors[planet], label=planet_names[planet])
  else:
    plot_res.plot(time_order[planet], planet_order[planet], color=colors[planet], label=planet_names[planet])
pl.xlabel("time [years]")
pl.ylabel("planet order")
pl.ylim((0, nb_planets+1))
pl.grid(True)
pl.show()


# TODO
# trier les resonances pour commencer par les résonances de nombre les plus faibles (vu qu'on s'arrête de chercher quand on en trouve une.
# how to represent the resonances in a graph?
# To display resonances, one plot of a un function of time, and another plot with orbital period in function of time ?

