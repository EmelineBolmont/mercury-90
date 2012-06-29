#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0
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
OUTPUT_EXTENSION = 'pdf' # default value

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
# On lit, pour chaque planète, le contenu du fichier et on stocke les variables qui nous intéressent.
####################
t = [] # temps en année
a = [] # demi-grand axe en ua
e = [] # excentricité
q = [] # perihelion
Q = [] # aphelion
I = [] # inclinaison (degrees)
n = [] # longitude of ascending node (degrees)
m = [] # masse de la planète en masse solaire



# On récupère les données orbitales
for planete in range(nb_planete):
	
	fichier_source = liste_aei[planete]
	tableau = open(fichier_source, 'r')
	
	t.append([]) # temps en année
	a.append([]) # demi-grand axe en ua
	e.append([]) # excentricité
	q.append([]) # perihelion
	Q.append([]) # aphelion
	I.append([]) # inclinaison (degrees)
	n.append([]) # longitude of ascending node (degrees)
	m.append([]) # masse de la planète en masse solaire

	# On passe les 3 premières lignes d'entête.
	for indice in range(3):
		tableau.readline()

	entete = tableau.readline()
	for ligne in tableau:
		# Pour chaque ligne du tableau, on découpe suivant les 
		# espaces. (par défaut, s'il y a plusieurs espaces à la 
		# suite, il ne va pas créer 50 000 colonnes, ce qui ne 
		# serait pas le cas avec un 'split(" ")')
		colonne = ligne.split()
		
		# On essaye de rajouter les éléments. Si un seul d'entre eux
		# à un soucis, on élimine toute la ligne (en gros, lors de 
		# l'éjection d'une planète, certains paramètres peuvent 
		# devenir ****** au lieu d'un float, et ça génère une erreur
		# lors de la conversion en float)
		try:
			ti = float(colonne[0])
			ai = float(colonne[1])
			ei = float(colonne[2])
			Ii = float(colonne[3])
			ni = float(colonne[5])
			mi = float(colonne[7])
		except:
			pass
		t[-1].append(ti)
		if (ai < 0):
			print("t="+str(ti)+" years a="+str(ai)+" AU for "+fichier_source)
		a[-1].append(ai)
		e[-1].append(ei)
		q[-1].append(ai * (1 - ei))
		Q[-1].append(ai * (1 + ei))
		I[-1].append(Ii)
		n[-1].append(ni)
		m[-1].append(mi * MS / MT)

		
	tableau.close()

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

# We generate a list of colors
tmp = autiwa.colorList(nb_planete)
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
		plot_a.plot(t[planet][id_min:id_max+1], a[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
		plot_a.plot(t[planet][id_min:id_max+1], q[planet][id_min:id_max+1], color=colors[planet])
		plot_a.plot(t[planet][id_min:id_max+1], Q[planet][id_min:id_max+1], color=colors[planet])

pl.xlabel(unicode("time [years]",'utf-8'))
pl.ylabel(unicode("a [AU]",'utf-8'))
#~ myyfmt = ScalarFormatter(useOffset=True)
#~ myyfmt._set_offset(1e9)
#~ myxfmt = ticker.ScalarFormatter()
#~ myxfmt.set_scientific(True)
#~ myxfmt.set_powerlimits((-3, 3)) 

#~ pl.legend()
pl.grid(True)

plot_e = fig.add_subplot(222, sharex=plot_a)
for planet in range(nb_planete):
	if isLog:
		plot_e.semilogx(t[planet][id_min:id_max+1], e[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
	else:
		plot_e.plot(t[planet][id_min:id_max+1], e[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
pl.xlabel(unicode("time [years]",'utf-8'))
pl.ylabel(unicode("eccentricity",'utf-8'))
pl.grid(True)

plot_m = fig.add_subplot(223, sharex=plot_a)
for planet in range(nb_planete):
	if isLog:
		plot_m.semilogx(t[planet][id_min:id_max+1], m[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
	else:
		plot_m.plot(t[planet][id_min:id_max+1], m[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
pl.xlabel(unicode("time [years]",'utf-8'))
pl.ylabel(unicode("mass [Earths]",'utf-8'))
pl.grid(True)

plot_I = fig.add_subplot(224, sharex=plot_a)
for planet in range(nb_planete):
	if isLog:
		plot_I.semilogx(t[planet][id_min:id_max+1], I[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
	else:
		plot_I.plot(t[planet][id_min:id_max+1], I[planet][id_min:id_max+1], color=colors[planet], label='PLANETE'+str(planet))
pl.xlabel(unicode("time [years]",'utf-8'))
pl.ylabel(unicode("Inclination [degrees]",'utf-8'))
pl.grid(True)

myyfmt = ScalarFormatter(useOffset=True, useMathText=True)
#~ myyfmt._set_offset(1e5)
myyfmt.set_scientific(True)
myyfmt.set_powerlimits((-2, 3)) 
#~ myyfmt._set_offset(1e9)
myxfmt = ScalarFormatter(useOffset=True, useMathText=True)
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

pl.figure(1)
nom_fichier_plot = "evolution_planete"
#~ pl.savefig(nom_fichier_plot+'.svg', format='svg')
pl.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

