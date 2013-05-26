#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version 1.4
# 09-08-12
# The script will calculate the resonances between each planet through time
#
# /!\ Warning : It seems to be problems in some python version in the line 62 (with the module Fraction. 
#               But it works perfectly on another server, so I don't know what to think about this.

import pdb # Pour le debug
import numpy as np
import analysis
import pylab as pl
import autiwa
import sys # to get access to arguments of the script
import mercury_utilities
import os
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter


# README
# There a several parameters to deal with to obtain an optimised calculation. 
# First, you can control the number of instants where you will check for resonances
# Second, you can control the number of points in the past you will use, for a given instant, to check if there were resonance. Beware
#  here, because too many points will make you loose some resonances (because the point of stability can vary in time, even with libration
#  and if there is too few points, you will find resonances that might not exist in reality
# And finally, you can modify the tolerance over the standard deviation for the resonant angle. 
# But you should not modify that unless you know what you are doing. The tolerance must not be too small, because weak resonances with
# secular change in the resonance angle libration can exist. But if there is not enough outputs, the risk is to find a resonance
# where there is none, only because the huge time between outputs can make the resonant angle to be in libration just by a sample effect.
#
# /!\ If a resonance is found at only one instant, she will not be displayed, becarefull of the number of instant where to test them


################
## Parameters ##
################
NOM_FICHIER_PLOT = "timed_resonances"
OUTPUT_EXTENSION = "png"

## The running time of the script is very sensitive to this parameter that determine the number of time, for each resonance, 
# that we search for a corresping fraction to the orbital period ratio
NUMBER_OF_VALUES = 10 # sampling for period ratio around the given value 
DENOMINATOR_LIMIT = 12 # Maximum value allowed of the denominator when we want to get a fraction from a decimal value
NUMERATOR_LIMIT = 20 # maximum value allowed for the numerator
UNCERTAINTY = 5 # In percentage
NB_LAST_POINTS = 15 # Number of points we want to test the libration of angles.

NB_MEASUREMENTS = 500 # The number of times we test the resonances between planets (because the total number of output can vary from one simulation to another)

# the threshold of the standard deviation of a resonant angle, 
# below which we consider there is libration and thus, a resonance.
STD_THRESHOLD = 70 

# Extreme angles for angle folding. We will apply a mod 2pi to all angles, 
# but theses values determines between wich values we will constrains the angles.
ANGLE_CENTER_VALUE = 90.
 
###############################################
## Beginning of the program
###############################################

# We get arguments from the script
isLog = False
isProblem = False
problem_message = "AIM : Display all the resonances of all the planet during the simulation" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * t_max : the end of the output (in years)" + "\n" + \
" * t_min : the beginning of the output (in years)" + "\n" + \
" * instants : (%d) the number of points in time where to search for resonances between planets" % NB_MEASUREMENTS + "\n" + \
" * sample : (%d) The number of successive points used to test a resonance" % NB_LAST_POINTS + "\n" + \
" * log : (%s) time will be displayed in log" % isLog + "\n" + \
" * ext=png : (%s) The extension for the output files" % OUTPUT_EXTENSION + "\n" + \
" * help : display this current message"

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 't_min'):
    t_min = float(value)
  elif (key == 't_max'):
    t_max = float(value)
  elif (key == 'instants'):
    NB_MEASUREMENTS = int(value)
  elif (key == 'sample'):
    NB_LAST_POINTS = int(value)
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
uncertainty = 0.01 * float(UNCERTAINTY)

# We initialize the arrays
t = [] # time in years
a = [] # demi-grand axe en ua
e = [] # eccentricity
g = [] # argument of pericentre (degrees)
n = [] # longitude of ascending node (degrees)
M = [] # Mean anomaly (degrees)
q = [] # periastron (AU)
Q = [] # apoastron (AU)


####################
# We read the datas for all the planets
####################

i = 0
error = True
positions = []
while (len(positions) == 0):
  try:
    positions = mercury_utilities.get_column_position(liste_aei[i])
  except:
    i += 1

# We retrieve the orbital data
for (planete, planet_datafile) in enumerate(liste_aei):
  sys.stdout.write("Reading data files %5.1f %% : %s        \r" % ((planete+1.) * 100. / float(nb_planets), planet_datafile))
  sys.stdout.flush()
  
  ti = [] # time in years
  ai = [] # Semi-major axis in AU
  ei = [] # eccentricity
  gi = [] # (in degrees)
  ni = [] # (in degrees)
  Mi = [] # (in degrees)

  
  fichier_source = liste_aei[planete]
  object_file = open(fichier_source, 'r')
  
  for i in range(4):
    object_file.readline()
  
  tiapp = ti.append
  aiapp = ai.append
  eiapp = ei.append
  giapp = gi.append
  niapp = ni.append
  Miapp = Mi.append
  

  for line in object_file:
    # When using [a:b], the actual range will be from a to b-1 included.
    tiapp(float(line[positions[0]:positions[1]]))
    aiapp(float(line[positions[1]:positions[2]]))
    eiapp(float(line[positions[2]:positions[3]]))
    giapp(float(line[positions[4]:positions[5]]))
    niapp(float(line[positions[5]:positions[6]]))
    Miapp(float(line[positions[6]:positions[7]]))

  object_file.close()
  
  ti = np.array(ti)
  ai = np.array(ai)
  ei = np.array(ei)
  gi = np.array(gi)
  ni = np.array(ni)
  Mi = np.array(Mi)
  
  qi = ai * (1 - ei)
  Qi = ai * (1 + ei)
  
  t.append(ti)
  a.append(ai)
  e.append(ei)
  g.append(gi)
  n.append(ni)
  M.append(Mi)
  q.append(qi)
  Q.append(Qi)

####################
# We change the range in time, if needed and do some pre-calculations
####################
# We get the array of reference time, i.e, one of the longuest list of time available in the list of planets. 
lengths = [ai.size for ai in a]

ref_id = np.argmax(lengths)# The ID of the longuest time array
ref_time = t[ref_id] # The longuest time array

delta_t = ref_time[1] - ref_time[0]

# We get the index for the t_max value
if ('t_max' in locals()):
  id_max = int((t_max - ref_time[0]) / delta_t)
else:
  id_max = len(ref_time)

# We get the index for the t_min value
if ('t_min' in locals()):
  id_min = int((t_min - ref_time[0]) / delta_t)
else:
  id_min = 0

for planet in range(nb_planets):
  t[planet] = t[planet][id_min:id_max]
  a[planet] = a[planet][id_min:id_max]
  e[planet] = e[planet][id_min:id_max]
  g[planet] = g[planet][id_min:id_max]
  n[planet] = n[planet][id_min:id_max]
  M[planet] = M[planet][id_min:id_max]

# We go in reverse order because once we delete the index 'i', all the indexes after are shifted and 4 become 3 if we delete the index 3
for planet in range(nb_planets)[::-1]: 
  if (len(t[planet]) == 0):
    del(t[planet])
    del(a[planet])
    del(e[planet])
    del(g[planet])
    del(n[planet])
    del(M[planet])
    del(liste_aei[planet])
    del(planet_names[planet])

# Once we eventually have deleted the planet that have dissapeared in the range considered, we calculate once again the total number of planets
nb_planets = len(t)

# We retrieve the longuest array (the planet that stay the longuest in the simulation) We do that again in case of a range selection in time
lengths = [ai.size for ai in a]
max_lengths = max(lengths)

# If the user require more measurement than timestep available, we force to have nb_measurements equal the number of timestep
if (NB_MEASUREMENTS < max_lengths):
  time_delay = max_lengths / NB_MEASUREMENTS
else:
  time_delay = 1

# We want at least to test every single point once. 
# Thus, we make sure that "NB_LAST_POINTS" is at least equal to the space between two tests of the resonances
if (time_delay > NB_LAST_POINTS):
  print("Warning: The interval between two checks of resonance is greater than the number of points used to test it.")
  print("         'NB_LAST_POINTS' extended to 'time_delay'=%d" % time_delay)
  NB_LAST_POINTS = time_delay

# We calculate q and Q
q = [ai * (1 - ei) for (ai, ei) in zip(a,e)]
Q = [ai * (1 + ei) for (ai, ei) in zip(a,e)]

####################
# We declare the arrays needed for the resonances
####################
# to store the successive resonances for each planet. If a resonance is written at the index i, that mean that the planet i is the inner planet of the resonance given in the array. 
resonance_index_range = [[] for planet in range(nb_planets)] # For each planet, each tuple is a range of time during which the planet is in resonance with the planet just after it.
resonance_type = [[] for planet in range(nb_planets)] # For each planet, store the type of resonance 3:2, 4:3 and so on. The index i of the sublist (of a given planet) correspond to the range in time of the array resonance_time

dynamic_order = [[] for planet in range(nb_planets)] # a list, for each planet, of its position through time. If 2, this means that the plant is the 2nd closest from its host star.
time_order = [[] for planet in range(nb_planets)] # the corresponding time in years
resonance_with = [[] for planet in range(nb_planets)] # the index of the outer planet for the resonance (the inner planet index being the index of the sublist considered)
resonance_inner_rank = [[] for planet in range(nb_planets)] # The rank, in distance, of the inner planet in the resonance (useful to place the resonance in the second plot)


### We create the first point for some arrays
# We sort the planet in the order of distance from the host star
distance = [ai[0] for ai in a]
planet_index_sorted_by_distance = np.argsort(distance) # the i-th closest planet in distance is the planet planet_index_sorted_by_distance[i]
ordering_planets = 1 + np.argsort(planet_index_sorted_by_distance) # the i-th planet is the ordering_planets[i] in distance starting at 1

# we append the ordering of the current planets in order to display resonances later
for (planet_idx, order) in enumerate(ordering_planets):
  dynamic_order[planet_idx].append(order)
  time_order[planet_idx].append(t[planet_idx][0])
###
for instant_index in range(NB_LAST_POINTS,max_lengths,time_delay):
  
  # We display a progress bar of the computation
  # The extra spaces are to make sure that no old character from the previous line will appear
  sys.stdout.write("Progression %6.2f %% : %i / %i                   \r" % ((instant_index * 100. / float(max_lengths)), instant_index, max_lengths))
  sys.stdout.flush()
  
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
  distance_end = [ai[-1] for ai in a_temp] # The distance of each planet from the host star, in AU, at the end of the current sub-range (of a_temp)
  distance_begin = [ai[0] for ai in a_temp] # The distance of each planet from the host star, in AU, at the beginning of the current sub-range (of a_temp)
  planet_index_sorted_by_distance = np.argsort(distance_end) # the i-th closest planet in distance is the planet planet_index_sorted_by_distance[i]
  ordering_planets = 1 + np.argsort(planet_index_sorted_by_distance) # the i-th planet is the ordering_planets[i] in distance (starting at 1)

  # we append the ordering of the current planets in order to display resonances later
  for (planet_idx, order) in zip(still_here_planets, ordering_planets):
    dynamic_order[planet_idx].append(order)
    time_order[planet_idx].append(t[planet_idx][instant_index])

  for (inner, outer) in zip(planet_index_sorted_by_distance[0:-1], planet_index_sorted_by_distance[1:]):
    g_inner = g_temp[inner]
    n_inner = n_temp[inner]
    M_inner = M_temp[inner]

    g_outer = g_temp[outer]
    n_outer = n_temp[outer]
    M_outer = M_temp[outer]

    # We get the various possible resonance between the two planets
    periodRatio_begin = (distance_begin[outer] / distance_begin[inner])**1.5
    periodRatio_end = (distance_end[outer] / distance_end[inner])**1.5
    
    # If the period ratio is too different between the beginning and the end of the range, we do not calculate possible resonances to gain time.
    if (abs(periodRatio_begin - periodRatio_end) > 0.02):
      continue
    
    resonances = analysis.get_possible_resonances(periodRatio_end, uncertainty=uncertainty, 
                 denominator_limit=DENOMINATOR_LIMIT, numerator_limit=NUMERATOR_LIMIT, sampling=NUMBER_OF_VALUES)

    # For each resonance we check if this one exist between the two considered planets
    index = 0
    while (index < len(resonances)):
      res = resonances[index]
      is_resonance = analysis.isResonance(res, g_inner, n_inner, M_inner, g_outer, n_outer, M_outer, 
                    nb_points=NB_LAST_POINTS, angle_center_value=ANGLE_CENTER_VALUE, std_threshold=STD_THRESHOLD)
      if (is_resonance):
        isExtend = False # boolean that say if the current resonance is the extension of the last resonance listed for the inner planet
        if (len(resonance_type[still_here_planets[inner]]) != 0):
          last_type = resonance_type[still_here_planets[inner]][-1]
          last_index_range = resonance_index_range[still_here_planets[inner]][-1]
          
          # if the two index ranges overlap
          if ((last_type == res) and (last_index_range[1] >= range_start-1)):
            isExtend = True
        
        # if the current resonance already existed before, we only extend the index range of validity for the last resonance of the inner planet index
        if isExtend:
          resonance_index_range[still_here_planets[inner]][-1][1] = instant_index
        else:
          # We test here the previous resonance. If she was only at one instant, we delete it
          if (len(resonance_index_range[still_here_planets[inner]]) != 0):
            tmp = resonance_index_range[still_here_planets[inner]][-1]
            if (tmp[0] == tmp[1]):
              del(resonance_type[still_here_planets[inner]][-1])
              del(resonance_index_range[still_here_planets[inner]][-1])
              del(resonance_with[still_here_planets[inner]][-1])
              del(resonance_inner_rank[still_here_planets[inner]][-1])
            
          resonance_type[still_here_planets[inner]].append(res)
          # To avoid overlap, we define the resonance at the position of the range_stop, without associating, by default, any length.
          resonance_index_range[still_here_planets[inner]].append([range_stop, range_stop])
          resonance_with[still_here_planets[inner]].append(still_here_planets[outer])
          resonance_inner_rank[still_here_planets[inner]].append(ordering_planets[inner])

        # If we find a resonance, we do not search for another one
        break
      index += 1

####################
# We now want to display in a fashion way the resonances
####################
sys.stdout.write("Generating graphics                          \r")
sys.stdout.flush()

# We generate a list of colors
colors = [ '#'+li for li in autiwa.colorList(nb_planets)]

fig = pl.figure(1)
pl.clf()
fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)
# On crée des sous plots. Pour subplot(311), ça signifie qu'on a 2 lignes, 3 colonnes, et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
plot_a = fig.add_subplot(311)
for planet in range(nb_planets):
  sys.stdout.write("Generating graphics  %5.1f %%                          \r" % ((planet+1) * 25. / float(nb_planets)))
  sys.stdout.flush()
  if isLog:
    plot_a.semilogx(t[planet], a[planet], color=colors[planet], label=planet_names[planet])
    plot_a.semilogx(t[planet], q[planet], color=colors[planet])
    plot_a.semilogx(t[planet], Q[planet], color=colors[planet])
  else:
    plot_a.plot(t[planet], a[planet], color=colors[planet], label=planet_names[planet])
    plot_a.plot(t[planet], q[planet], color=colors[planet])
    plot_a.plot(t[planet], Q[planet], color=colors[planet])

ylims = list(pl.ylim())
for planet in range(nb_planets):
  sys.stdout.write("Generating graphics  %5.1f %%                          \r" % (20.+(planet+1) * 20. / float(nb_planets)))
  sys.stdout.flush()
  for (outer, (id_begin, id_end)) in zip(resonance_with[planet], resonance_index_range[planet]):
    plot_a.plot([t[planet][id_begin], t[planet][id_begin]], [a[planet][id_begin], a[outer][id_begin]], 'k--')
    plot_a.plot([t[planet][id_end], t[planet][id_end]], [a[planet][id_end], a[outer][id_end]], 'k--')

# For a huge number of planet, the legend will be horrible, so we skip this line in that case
if (nb_planets<=10):
  plot_a.legend()

plot_a.set_xlabel("time [years]")
plot_a.set_ylabel("a [AU]")
plot_a.grid(True)

plot_res = fig.add_subplot(312, sharex=plot_a)
for planet in range(nb_planets):
  sys.stdout.write("Generating graphics  %5.1f %%                          \r" % (40.+(planet+1) * 20. / float(nb_planets)))
  sys.stdout.flush()
  for (res, inner_rank, (id_begin, id_end)) in zip(resonance_type[planet], resonance_inner_rank[planet], resonance_index_range[planet]):
    x_position = (t[planet][id_begin] + t[planet][id_end]) / 2.
    y_position = inner_rank + .5
    
    x_left = t[planet][id_begin]
    x_right = t[planet][id_end]
    y_bottom = inner_rank + 0.1
    y_top = inner_rank + 0.9
    
    text = "%i:%i" % (res.numerator, res.denominator)
    plot_res.text(x_position, y_position, text, horizontalalignment='center', verticalalignment='center', rotation='vertical')
    plot_res.fill([x_left, x_right, x_right, x_left],[y_bottom, y_bottom, y_top, y_top], fill=True, color='#cccccc')
    plot_res.plot([x_left, x_left], [y_bottom, y_top], 'k-')
    plot_res.plot([x_right, x_right], [y_bottom, y_top], 'k-')

for planet in range(nb_planets):
  sys.stdout.write("Generating graphics  %5.1f %%                          \r" % (60.+(planet+1) * 20. / float(nb_planets)))
  sys.stdout.flush()
  if isLog:
    plot_res.semilogx(time_order[planet], dynamic_order[planet], color=colors[planet], label=planet_names[planet])
  else:
    plot_res.plot(time_order[planet], dynamic_order[planet], color=colors[planet], label=planet_names[planet])
plot_res.set_xlabel("time [years]")
plot_res.set_ylabel("planet order")
plot_res.set_ylim((0, nb_planets+1))
plot_res.yaxis.set_ticks(range(1, nb_planets+1))
plot_res.yaxis.set_ticklabels(range(1, nb_planets+1))
plot_res.grid(True)

plot_e = fig.add_subplot(313, sharex=plot_a)
for planet in range(nb_planets):
  sys.stdout.write("Generating graphics  %5.1f %%                          \r" % ((planet+1) * 25. / float(nb_planets)))
  sys.stdout.flush()
  if isLog:
    plot_e.loglog(t[planet], e[planet], color=colors[planet], label=planet_names[planet])
  else:
    plot_e.semilogy(t[planet], e[planet], color=colors[planet], label=planet_names[planet])
plot_e.set_xlabel("time [years]")
plot_e.set_ylabel("eccentricity")
plot_e.grid(True)

myxfmt = ScalarFormatter(useOffset=True)
myxfmt._set_offset(1e5)
myxfmt.set_scientific(True)
myxfmt.set_powerlimits((-3, 3)) 
plot_a.xaxis.set_major_formatter(myxfmt)
plot_res.yaxis.set_major_formatter(FormatStrFormatter('%i'))

sys.stdout.write("Saving graphics                          \r")
sys.stdout.flush()
pl.savefig('%s.%s' % (NOM_FICHIER_PLOT, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

sys.stdout.write("Displaying graphics                          \n")
sys.stdout.flush()
pl.show()

## TODO
# Store the list of fraction for a given period ratio. If a period ratio is less than 
#  delta ratio/2 of an existing one, we do not calculate again the fraction, but instead retrieve the existing list

## Tricks
# One thing to understand is the fact that when checking resonances at t=ti, we actually search for a resonance between t=ti-dt and t=ti. 
# If, for a particular reason, the current resonance can be extended, the resonance validity is extended to ti-dt to ti'=ti+step 
# where step is the separation in time between each check of resonances. 
# Though, if step is smaller than dt, then a new resonance can be declare in a range that overlap the previous one.
