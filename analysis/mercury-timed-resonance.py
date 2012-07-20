#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version 1.2

# The script will calculate the resonances between each planet through time

import pdb # Pour le debug
import numpy as np
from fractions import Fraction
import pylab as pl
import autiwa
import sys # to get access to arguments of the script
import mercury_utilities

################
## Parameters ##
################
NOM_FICHIER_PLOT = "system_resonances"
OUTPUT_EXTENSION = "pdf"

DENOMINATOR_LIMIT = 30 # Maximum value allowed of the denominator when we want to get a fraction from a decimal value
NUMBER_OF_VALUES = 100 # sampling for period ratio around the given value
UNCERTAINTY = 5 # In percentage
NB_LAST_POINTS = 50 # Number of points we want to test the libration of angles.

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

#~ pdb.set_trace()
isProblem = False
problem_message = "AIM : Display in a m = f(a) diagram, all the planets of the current mercury simulation" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * cz=1 (The position of a convergence zone in AU)" + "\n" + \
"   cz=[[1,60],[4,30]] (the list of mass (earth mass) and zero torque position (in AU) successively)" + "\n" + \
" * ext=png (The extension for the output files)" + "\n" + \
" * help : display this current message"

for arg in sys.argv[1:]:
	try:
		(key, value) = arg.split("=")
	except:
		key = arg
	if (key == 'ext'):
		OUTPUT_EXTENSION = value
	elif (key == 'cz'):
		CZ = eval(value)
	elif (key == 'help'):
		print(problem_message)
		exit()
	else:
		print("the key '"+key+"' does not match")
		isProblem = True

if isProblem:
	print(problem_message)


####################
# Initialization of some variables
####################
uncertainty = 0.01 * float(UNCERTAINTY)

# We initialize the arrays
t = [] # time in years
a = [] # demi-grand axe en ua
g = [] # argument of pericentre (degrees)
n = [] # longitude of ascending node (degrees)
M = [] # Mean anomaly (degrees)

####################
# we get the name of the files for all the planets in the system
####################
liste_aei = get_aei_files()
nb_planets = len(liste_aei)

# We also initialize the list for the names of the planets
planet_names = []
for filename in liste_aei:
  planet_names.append(os.path.splitext(filename)[0]


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

# we calculate the period
period = []
for ai in a:
  period.append(ai**1.5)

# We retrieve the longuest array (the planet that stay the longuest in the simulation
lengths = []
for ai in a:
  length.append(ai.size)
longuest_index = np.argmax(length)

#~ for instant_index in t[longuest_index][NB_LAST_POINTS::NB_LAST_POINTS]: TODO
instant_index = 51 # the index for the current time where we check the resonances

# We make sublist only with the planets that still exist at the present time
a_temp = []
g_temp = []
n_temp = []
M_temp = []

for planet in range(nb_planets):
  if (t[planet].size > instant_index):
    a_temp.append([planet][instant_index - NB_LAST_POINTS:instant_index])
    g_temp.append([planet][instant_index - NB_LAST_POINTS:instant_index])
    n_temp.append([planet][instant_index - NB_LAST_POINTS:instant_index])
    M_temp.append([planet][instant_index - NB_LAST_POINTS:instant_index])
    

# We sort the planet in the order of distance from the host star
distance = [ai[-1] for ai in a_temp]
ordering_planets = np.argsort(distance)

for (inner, outer) in zip(ordering_planets[0:-1], ordering_planets[1:]):
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
  res = resonances[i]
  is_resonance = isResonance(res, g_inner, n_inner, M_inner, g_outer, n_outer, M_outer)
  if (is_resonance):
    TODO add resonance to a list. But don't know how to display it yet
  index += 1



# TODO
# trier les resonances pour commencer par les résonances de nombre les plus faibles (vu qu'on s'arrête de chercher quand on en trouve une.
# how to represent the resonances in a graph?

