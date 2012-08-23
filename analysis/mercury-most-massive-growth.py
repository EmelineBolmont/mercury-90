#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.1
# Pour lire des fichiers de simulations, récupérer les caractéristiques
# des planètes qu'il reste en fin de simulation et les écrire dans un 
# seul fichier que le script "analyse_simu" va lire

from math import *
import pylab as pl
import os, pdb, autiwa
import numpy as np
from constants import MT, MS
import sys # to be able to retrieve arguments of the script

# Maximum number of planets (the most massives) that will be colored
MAX_COLORED = 3
COLOR_FINAL = '000000'
COLOR_EJECTED = '909090'

###############################################
## Beginning of the program
###############################################
# We set the false option before. Because if not, we will erase the 'true' with other option that are not log, and 
# thus will lead to be in the else and put log to false.
isLog = False 

isaLog = False # Put a semilog in 'a' if true

OUTPUT_EXTENSION = "pdf" # default extension for outputs

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * t_max : the end of the output, in years)" + "\n" + \
" * t_min : the beginning of the output (in years)" + "\n" + \
" * massive : the number of most massive planets to be tracked" + "\n" + \
" * log : time (x-axis) will be displayed in log" + "\n" + \
" * alog : semi major axis (y-axis) will be displayed in log" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * ext=pdf : The extension for the output files"

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
	elif (key == 'massive'):
		MAX_COLORED = int(value)
	elif (key == 'log'):
		isLog = True
	elif (key == 'alog'):
		isaLog = True
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
a = [] # semi major axis in AU
#~ e = [] # eccentricity
#~ q = [] # perihelion
#~ Q = [] # aphelion
#~ I = [] # inclinaison (degrees)
m = [] # planet mass in earth mass

# We retrieve the orbital data
for planete in range(nb_planete):
	
	fichier_source = liste_aei[planete]
	(ti, ai, mi) = np.loadtxt(fichier_source, skiprows=4, usecols = (0,1,7), dtype=float, unpack=True)
	#~ qi = ai * (1 - ei)
	#~ Qi = ai * (1 + ei)
	mi = (MS / MT) * mi
	
	if (type(ti) == np.ndarray):
		t.append(ti)
		a.append(ai)
		#~ e.append(ei)
		#~ q.append(qi)
		#~ Q.append(Qi)
		#~ I.append(Ii)
		m.append(mi)
	else:
		# In case the is only one point, we force to have a list, to avoid plotting problems
		t.append(np.array([ti]))
		a.append(np.array([ai]))
		#~ e.append(np.array([ei]))
		#~ q.append(np.array([qi]))
		#~ Q.append(np.array([Qi]))
		#~ I.append(np.array([Ii]))
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
colors = ['#'+COLOR_EJECTED] * nb_planete

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
		if (float(words[-2]) > t_max):
			continue
		remaining_planet = words[0]
		lost_planet = words[4]
		colors[ind_of_planet[lost_planet]] = colors[ind_of_planet[remaining_planet]]
		lost_in_collisions.append(ind_of_planet[lost_planet])

#~ # We generate a list of colors
#~ tmp = autiwa.colorList(nb_planete)
#~ colors = [ '#'+li for li in autiwa.colorList(nb_planete)]

# on trace les plots

fig = pl.figure(1)
pl.clf()
fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)

# On crée des sous plots. Pour subplot(321), ça signifie qu'on a 2 lignes, 3 colonnes, et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
plot_a = fig.add_subplot(211)
for planet in range(nb_planete):
	if isaLog:
		if isLog:
			plot_a.loglog(t[planet][id_min:id_max+1], a[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
		else:
			plot_a.semilogy(t[planet][id_min:id_max+1], a[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
	else:
		if isLog:
			plot_a.semilogx(t[planet][id_min:id_max+1], a[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
		else:
			plot_a.plot(t[planet][id_min:id_max+1], a[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))


for planet in lost_in_collisions:
	if isaLog:
		if isLog:
			plot_a.loglog(t[planet][-1], a[planet][-1], 'o', markerfacecolor='None', markeredgewidth=2, markeredgecolor=colors[planet])
		else:
			plot_a.semilogy(t[planet][-1], a[planet][-1], 'o', markerfacecolor='None', markeredgewidth=2, markeredgecolor=colors[planet])
	else:
		if isLog:
			plot_a.semilogx(t[planet][-1], a[planet][-1], 'o', markerfacecolor='None', markeredgewidth=2, markeredgecolor=colors[planet])
		else:
			plot_a.plot(t[planet][-1], a[planet][-1], 'o', markerfacecolor='None', markeredgewidth=2, markeredgecolor=colors[planet])

plot_a.set_xlabel("time [years]")
plot_a.set_ylabel("a [UA]")
plot_a.grid(True)

plot_mass = fig.add_subplot(212, sharex=plot_a)
for planet in range(nb_planete):
	if isLog:
		plot_mass.semilogx(t[planet][id_min:id_max+1], m[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
	else:
		plot_mass.plot(t[planet][id_min:id_max+1], m[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
plot_mass.set_xlabel("time [years]")
plot_mass.set_ylabel("mass [Earths]")
plot_mass.grid(True)


nom_fichier_plot = "plot_am"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

