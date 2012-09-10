#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0.1
# Script that will display the evolution of the semi major axis, mass, 
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



###############################################
## Beginning of the program
###############################################

isLog = False # We set the false option before. Because if not, we will erase the 'true' with other option that are not log, and 
# thus will lead to be in the else and put log to false.
OUTPUT_EXTENSION = 'png' # default value in bitmap, because vectoriel can take time and space if there is a lot of data

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * t_max (the end of the output, in years)" + "\n" + \
" * t_min (the beginning of the output (in years)" + "\n" + \
" * log (time will be displayed in log)" + "\n" + \
" * help (display a little help message on HOW to use various options" + "\n" + \
" * ext=png (The extension for the output files)"

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
e = [] # eccentricity
q = [] # perihelion
Q = [] # aphelion
I = [] # inclinaison (degrees)
m = [] # planet mass in earth mass

# We retrieve the orbital data
for planete in range(nb_planete):
	
	fichier_source = liste_aei[planete]
	(ti, ai, ei, Ii, mi) = np.loadtxt(fichier_source, skiprows=4, usecols = (0,1,2,3,7), dtype=float, unpack=True)
	qi = ai * (1 - ei)
	Qi = ai * (1 + ei)
	mi = (MS / MT) * mi
	
	if (type(ti) == np.ndarray):
		t.append(ti)
		a.append(ai)
		e.append(ei)
		q.append(qi)
		Q.append(Qi)
		I.append(Ii)
		m.append(mi)
	else:
		# In case the is only one point, we force to have a list, to avoid plotting problems
		t.append(np.array([ti]))
		a.append(np.array([ai]))
		e.append(np.array([ei]))
		q.append(np.array([qi]))
		Q.append(np.array([Qi]))
		I.append(np.array([Ii]))
		m.append(np.array([mi]))


# We get the array of reference time, i.e, one of the longuest list of time available in the list of planets. 
len_t = [ti.size for ti in t]
ref_id = np.argmax(len_t) # The ID of the longuest time array
ref_time = t[ref_id] # The longuest time array
ref_len = ref_time.size

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

# on trace les plots

fig = pl.figure(1)
pl.clf()
fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)
# On crée des sous plots. Pour subplot(311), ça signifie qu'on a 2 lignes, 3 colonnes, et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
plot_a = fig.add_subplot(221)
for planet in range(nb_planete):
	if isLog:
		plot_a.semilogx(t[planet][id_min:id_max+1], a[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
		plot_a.semilogx(t[planet][id_min:id_max+1], q[planet][id_min:id_max+1], color=colors[planet])
		plot_a.semilogx(t[planet][id_min:id_max+1], Q[planet][id_min:id_max+1], color=colors[planet])
	else:
		try:
			plot_a.plot(t[planet][id_min:id_max+1], a[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
			plot_a.plot(t[planet][id_min:id_max+1], q[planet][id_min:id_max+1], color=colors[planet])
			plot_a.plot(t[planet][id_min:id_max+1], Q[planet][id_min:id_max+1], color=colors[planet])
		except:
			pdb.set_trace()

pl.xlabel("time [years]")
pl.ylabel("a [AU]")
#~ myyfmt = ScalarFormatter(useOffset=True)
#~ myyfmt._set_offset(1e9)
#~ myxfmt = ticker.ScalarFormatter()
#~ myxfmt.set_scientific(True)
#~ myxfmt.set_powerlimits((-3, 3)) 

#~ pl.legend()
pl.grid(True)

plot_e = fig.add_subplot(2, 2, 2, sharex=plot_a)
for planet in range(nb_planete):
	if isLog:
		plot_e.semilogx(t[planet][id_min:id_max+1], e[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
	else:
		plot_e.plot(t[planet][id_min:id_max+1], e[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
pl.xlabel("time [years]")
pl.ylabel("eccentricity")
pl.grid(True)

plot_m = fig.add_subplot(2, 2, 3, sharex=plot_a)
for planet in range(nb_planete):
	if isLog:
		plot_m.semilogx(t[planet][id_min:id_max+1], m[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
	else:
		plot_m.plot(t[planet][id_min:id_max+1], m[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
pl.xlabel("time [years]")
pl.ylabel("mass [Earths]")
pl.grid(True)

plot_I = fig.add_subplot(2, 2, 4, sharex=plot_a)
for planet in range(nb_planete):
	if isLog:
		plot_I.semilogx(t[planet][id_min:id_max+1], I[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
	else:
		plot_I.plot(t[planet][id_min:id_max+1], I[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
pl.xlabel("time [years]")
pl.ylabel("Inclination [degrees]")
pl.grid(True)

myyfmt = ScalarFormatter(useOffset=False)
#~ myyfmt._set_offset(1e5)
myyfmt.set_scientific(True)
myyfmt.set_powerlimits((-2, 3)) 
#~ myyfmt._set_offset(1e9)
myxfmt = ScalarFormatter(useOffset=True)
myxfmt._set_offset(1e5)
myxfmt.set_scientific(True)
myxfmt.set_powerlimits((-3, 3)) 
#~ myxfmt = FormatStrFormatter('%0.0e')
#~ pdb.set_trace()
#~ plot_a.xaxis.set_major_formatter(FormatStrFormatter('%0.0e'))
plot_a.xaxis.set_major_formatter(myxfmt)
plot_I.yaxis.set_major_formatter(myyfmt)
plot_e.yaxis.set_major_formatter(myyfmt)

#~ dossier_output = "output"
#~ system("mkdir dossier_output")
#~ system("cd dossier_output")


nom_fichier_plot = "evolution_planete"
#~ pl.savefig(nom_fichier_plot+'.svg', format='svg')
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

