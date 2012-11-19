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
#	here, because too many points will make you loose some resonances (because the point of stability can vary in time, even with libration
#	and if there is too few points, you will find resonances that might not exist in reality
# And finally, you can modify the tolerance over the standard deviation for the resonant angle. 
# But you should not modify that unless you know what you are doing. The tolerance must not be too small, because weak resonances with
# secular change in the resonance angle libration can exist. But if there is not enough outputs, the risk is to find a resonance
# where there is none, only because the huge time between outputs can make the resonant angle to be in libration just by a sample effect.
#
# /!\ If a resonance is found at only one instant, she will not be displayed, becarefull of the number of instant where to test them


################
## Parameters ##
################
NOM_FICHIER_PLOT = "plot-periods"
OUTPUT_EXTENSION = "png"

 
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

# We initialize the arrays
t = [] # time in years
a = [] # demi-grand axe en ua
e = [] # eccentricity


t_temp = [] # time in years
a_temp = [] # demi-grand axe en ua
e_temp = [] # eccentricity

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
	ai = [] # semi major axis in AU
	ei = [] # eccentricity
	
	fichier_source = liste_aei[planete]
	object_file = open(fichier_source, 'r')
	
	for i in range(4):
		object_file.readline()
	
	tiapp = ti.append
	aiapp = ai.append
	eiapp = ei.append
	
	for line in object_file:
		# When using [a:b], the actual range will be from a to b-1 included.
		tiapp(float(line[positions[0]:positions[1]]))
		aiapp(float(line[positions[1]:positions[2]]))
		eiapp(float(line[positions[2]:positions[3]]))
	object_file.close()
	
	if (ti != []):
		t_temp.append(np.array(ti))
		a_temp.append(np.array(ai))
		e_temp.append(np.array(ei))

a_init = [ai[0] for ai in a_temp]
initial_order = np.argsort(a_init)

for order in initial_order:
	t.append(t_temp[order])
	a.append(a_temp[order])
	e.append(e_temp[order])

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


# Once we eventually have deleted the planet that have dissapeared in the range considered, we calculate once again the total number of planets
nb_planets = len(t)

for planet in range(nb_planets):
  t[planet] = t[planet][id_min:id_max]
  a[planet] = a[planet][id_min:id_max]
  e[planet] = e[planet][id_min:id_max]

# We retrieve the longuest array (the planet that stay the longuest in the simulation) We do that again in case of a range selection in time
lengths = [ai.size for ai in a]
max_lengths = max(lengths)

pl_min = [] # The extremum to plot datas
pl_max = [] # The extremum to plot datas
for planet in range(nb_planets):
	pl_min.append(max(id_min, 0))
	pl_max.append(min(id_max, lengths[planet]))

OP = [ai**1.5 for ai in a] # The orbital period in years

#~ ####################
#~ # We declare the arrays needed for the resonances
#~ ####################


orbital_periods = np.empty((nb_planets, max_lengths))
orbital_periods.fill(np.NaN)
for i in range(nb_planets):
	orbital_periods[i, 0:lengths[i]] = OP[i]


planet_index_sorted_by_distance = orbital_periods.argsort(axis=0)
planet_rank = 1 + planet_index_sorted_by_distance.argsort(axis=0) # for each planet (one line), the various rank in orbital distance through time


orbital_periods.sort(axis=0)
period_ratio = np.empty((nb_planets-1, max_lengths))
for planet in range(nb_planets-1):
	period_ratio[planet] = orbital_periods[planet+1] / orbital_periods[planet]


####################
# We now want to display in a fashion way the resonances
####################

# We generate a list of colors
colors = [ '#'+li for li in autiwa.colorList(nb_planets)]

fig = pl.figure(1)
pl.clf()

plot_a = fig.add_subplot(311)
for planet in range(nb_planets):
	plot_a.plot(t[planet], a[planet], color=colors[planet], label=planet_names[planet])

plot_a.set_xlabel("time [years]")
plot_a.set_ylabel("a [AU]")
plot_a.grid(True)
plot_a.legend()

plot_order = fig.add_subplot(312, sharex=plot_a)
for planet in range(nb_planets):
	plot_order.plot(t[planet], planet_rank[planet][pl_min[planet]:pl_max[planet]], color=colors[planet])

plot_order.set_xlabel("time [years]")
plot_order.set_ylabel("order")
plot_order.grid(True)

ylims = plot_order.get_ylim()
plot_order.set_ylim([ylims[0]-1, ylims[1]+1])


plot_PR = fig.add_subplot(313, sharex=plot_a)
for planet in range(nb_planets-1):
	plot_PR.plot(ref_time, period_ratio[planet], color=colors[planet], label="period ratio %i/%i" % (planet+2, planet+1))

plot_PR.set_xlabel("time [years]")
plot_PR.set_ylabel("period ratio")
plot_PR.grid(True)
plot_PR.legend()
plot_PR.set_ylim(ymin=0.95)

print("                                         ")

pl.show()
