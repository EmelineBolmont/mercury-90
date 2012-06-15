#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0.1.1
# Pour générer des histogrammes sur plusieurs simulations afin de 
# regarder les caractéristiques statistiques des planètes qui restent.


# TODO
# _ écrire un fichier qui stocke les paramètres orbitaux des planètes qu'il reste.
# _ sélectionner les planètes qui sont plus proche qu'une certaine valeur de leur étoile.

from math import *
import pylab as pl
import os, pdb, autiwa
from simu_constantes import *
import random
import numpy as np
import sys

# Mass threshold
MASS_THRESHOLD = 5 # Earth mass

CONVERGENCE_ZONE = 2.5 # in AU, the location of the convergence zone (as a reference for plotting period ratios)

DELTA_RATIO = 0.005
DELTA_M = 0.1

PREFIX = "simu" # the prefix for the directories we want to take into the statistic.
t_max = 1e7 # the integration time of the simulations

#######################
# On prépare le fichier log
#######################
# Get current working directory
rep_exec = os.getcwd()

resonances = ["4:3", "5:4", "6:5", "7:6", "8:7", "9:8", "10:9", "11:10"]

#######################
# On prépare les plots
#######################

OUTPUT_EXTENSION = "pdf" # default extension for outputs

dossier_plot = "output"

#######################
# On définit et prépare les différents dossiers 
# où on doit lire les résultats des simulations
#######################
# On définit une liste de nom de dossier que, quoi qu'il arrive
# on ne souhaite pas traiter et donc qu'on vire de la liste.
dossier_suppr = ["output", "test", "reference_simu"]

#    .-.     .-.     .-.     .-.     .-.     .-.     .-.     .-.     .-. 
#  .'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `.
# (    .     .-.     .-.     .-.     .-.     .-.     .-.     .-.     .    )
#  `.   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   .'
#    )    )                                                       (    (
#  ,'   ,'                                                         `.   `.
# (    (                     DEBUT DU PROGRAMME                     )    )
#  `.   `.                                                         .'   .' 
#    )    )                                                       (    (
#  ,'   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   `.
# (    '  _  `-'  _  `-'  _  `-'  _  `-'  _  `-'  _  `-'  _  `-'  _  `    )
#  `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .'
#    `-'     `-'     `-'     `-'     `-'     `-'     `-'     `-'     `-'

isLog = False # We set the false option before. Because if not, we will erase the 'true' with other option that are not log, and 
# thus will lead to be in the else and put log to false.

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help (display a little help message on HOW to use various options" + "\n" + \
" * ext=pdf (The extension for the output files)"

# We get arguments from the script
for arg in sys.argv[1:]:
	try:
		(key, value) = arg.split("=")
	except:
		key = arg
	if (key == 'ext'):
		OUTPUT_EXTENSION = value
	elif (key == 'help'):
		isProblem = True
	else:
		print("the key '"+key+"' does not match")
		isProblem = True

if isProblem:
	print(problem_message)
	exit()


# We go in each sub folder of the current working directory

# On récupère la liste des sous-dossiers
liste_simu = [dir for dir in os.listdir(".") if (os.path.isdir(dir) and dir.startswith(PREFIX))]
autiwa.suppr_dossier(liste_simu,dossier_suppr)
liste_simu.sort()


nb_simu = len(liste_simu)

# We initialize the variables where to store datas

# List of orbital elements of ALL the simulations
a = [] # in AU
e = [] # eccentricity
I = [] # inclination in degre
m = [] # mass in earth mass
m_relat = [] # mass expressed in function of the mass of the most massive planet of the system

# List of orbital elements of the closest planet from the convergence zone
a_clo = [] # in AU
e_clo = [] # eccentricity
I_clo = [] # inclination in degre
m_clo = [] # mass in earth mass

# List of period ratios for each planet with the reference being, in each system, the closest planet from the resonance.
period_ratio = []

most_massive = [] # the most massive planet of the system
second_massive = [] # the second most massive planet of the system

# Lists for masses of the first and second planets of resonances
m_first = []  # All the coorbitals
m_second = [] # All the coorbitals
m_clo_first = []  # The coorbitals with the closest planet from the convergence zone
m_clo_second = [] # The coorbitals with the closest planet from the convergence zone

final_nb_planets = [] # the final number of planets in the system
final_time = [] # Final time of the integration. Used to check if there was problems with the simulations

print("\t Reading datas")

for simu in liste_simu:
	os.chdir(os.path.join(rep_exec, simu))
	try:
		dataFile = open("element.out", 'r')
	except:
		print("folder "+simu+" does not contains any 'element.out' file")
		continue
	
	# List initialized for each new simulation that contains values of the current simulation
	a_system = [] # in AU
	e_system = [] # eccentricity
	I_system = [] # inclination in degre
	m_system = [] # mass in earth mass
	
	# We get rid of the header
	dataFile.readline()
	tmp = dataFile.readline()
	tmp = tmp.split(":")[1]
	t_max0 = float(tmp)
	if (t_max0 != t_max):
		print("folder "+simu+" output time is "+tmp+" instead of "+str(t_max))
		continue
	#~ final_time.append(t_max0)
	
	for i in range(3):
		dataFile.readline()
	
	lines = dataFile.readlines()
	dataFile.close()
	
	nb_planets = len(lines)
	final_nb_planets.append(nb_planets)
	
	for line in lines:
		datas = line.split()
		ai = float(datas[1])
		ei = float(datas[2])
		Ii = float(datas[3])
		mi = float(datas[4])
		a_system.append(ai)
		e_system.append(ei)
		I_system.append(Ii)
		m_system.append(mi)
		if (ai < 23 and ai > 19 and mi > 16 and mi < 20):
			print("in "+simu+" a planet has a="+str(ai)+" AU and m="+str(mi)+" mt")
		#~ if (a_system[-1]> 20.):
			#~ print("in "+simu+" "+datas[0]+" has a="+datas[1])
	
	# We add all the elements of the system into the list of values corresponding to all the simulations
	a.extend(a_system)
	e.extend(e_system)
	I.extend(I_system)
	m.extend(m_system)
	
	###########
	# Statistical studies
	
	if (final_nb_planets[-1]>1):
		# We search for the most massive planet of the system
		a_system = np.array(a_system)
		#~ a_temp = abs(a_system - CONVERGENCE_ZONE)
		#~ idx_clo = a_temp.argmin()
		m_system = np.array(m_system)
		idx_clo = m_system.argmax()
		a_ref = a_system[idx_clo]
		
		a_clo.append(a_system[idx_clo])
		e_clo.append(e_system[idx_clo])
		I_clo.append(I_system[idx_clo])
		m_clo.append(m_system[idx_clo])
	
	
		#~ # We delete the element corresponding to the reference. Hence, each planet at '1' will be a coorbital
		#~ a_system = np.delete(a_system, idx_clo)
		#~ period_ratio.extend((a_system / a_ref)**(1.5))
		
		# We search for all the previous planets that are in coorbit with the reference one
		idx_before = idx_clo # 
		while (idx_before > 0 and (a_system[idx_before-1]/a_ref)**1.5>(1-DELTA_RATIO)):
			idx_before -= 1
		
		# We search for all the outer planets that are in coorbit with the reference one
		idx_after = idx_clo # 
		while (idx_after < nb_planets-1 and (a_system[idx_after+1]/a_ref)**1.5<(1+DELTA_RATIO)):
			idx_after += 1
			
		#~ idx_before:idx_after is the range of the coorbitals of the reference planet
		
		if (idx_before > 0):
			idx_before -= 1
		
		if (idx_after < nb_planets-1):
			idx_after += 1
		
		# We search for the two most massives planets of a system
		tmp = list(m_system)
		tmp.sort()
		most_massive.append(tmp[-1])
		second_massive.append(tmp[-2])

		if (tmp[-1] > 6. and tmp[-1] < 9.):
			print("most massive = %.1f in %s:" % (tmp[-1],simu))
		#~ if (final_nb_planets[-1] == 7):
			#~ print("nb_planets :",final_nb_planets[-1],simu)
		#~ if (max(e_system) > 0.8):
			#~ print("max eccentricity :",max(e_system),simu)
	else:
		print("in "+simu+" there is only "+str(final_nb_planets[-1])+" planet left")
	
	#~ idx_before:idx_after is the range of planets between the first non coorbital inner and outer the position of the reference planet (if they exists)
	
	period_ratio.extend((a_system[idx_before:idx_clo]/a_ref)**(1.5))# We do not add 1 to the last index to exclude the central planet
	period_ratio.extend((a_system[idx_clo+1:idx_after+1]/a_ref)**(1.5))# We need to add 1 the the outer index because the last index is excluded
	
	
	
	i = 0
	while (i<nb_planets):
		a_ref = a_system[i]
		m_ref = m_system[i]
		i += 1
		while (i<nb_planets):# We search for all the coorbitals with the current reference planet
			tmp = (a_system[i] / a_ref)**1.5
			res = (tmp if tmp > 1 else 1/tmp)
			
			if (res<(1+DELTA_RATIO)):
				if (m_ref>=m_system[i]):
					m_min = m_system[i]
					m_max = m_ref
				else:
					m_max = m_system[i]
					m_min = m_ref
				
				# We do not store in the same array if the coorbitals are close or not from the convergence zone
				if (1+abs(1.-a_ref/a_system[idx_clo]) < (1+DELTA_RATIO)):
					m_clo_first.append(m_max)
					m_clo_second.append(m_min)
					
				else:
					m_first.append(m_max)
					m_second.append(m_min)
				i += 1
			else:
				break
			
	##########
print("less massive in the list of most massive :",min(most_massive))
print("most massive planet formed in all simulations :",max(most_massive))

#######################
#   Tracé des plots   #
#######################
print("\t Computing Plots")
os.chdir(rep_exec)
nb_bins = 50

nom_fichier_plot = [] # list of names for each plot

nom_fichier_plot1 = "miscellaneous"
pl.figure(1)
pl.subplot(231)
pl.xlabel(unicode("masse (en mj)",'utf-8'))
pl.ylabel("density of probability")
pl.hist(m_clo, bins=range(25), normed=True, histtype='step')

m2 = [mi + random.uniform(-0.5,0.5) for mi in m]
pl.subplot(232)
pl.xlabel("mass (in earth mass)")
pl.ylabel("eccentricity")
pl.plot(m2, e, 'o', markersize=5)

pl.subplot(233)
pl.xlabel("distance (in AU)")
pl.ylabel("eccentricity")
pl.plot(a, e, 'o', markersize=5)

pl.subplot(234)
pl.xlabel(unicode("I (in degrees)",'utf-8'))
pl.ylabel("density of probability")
pl.hist(I, bins=[0.002*i for i in range(25)], normed=True, histtype='step')

pl.subplot(235)
pl.xlabel(unicode("nb_final",'utf-8'))
pl.ylabel("density of probability")
pl.hist(final_nb_planets, bins=range(25), histtype='step')

nom_fichier_plot2 = "m_fct_a"
m2 = [mi + random.uniform(-0.5,0.5) for mi in m]
pl.figure(2)
pl.xlabel(unicode("a (in AU)",'utf-8'))
pl.ylabel("mass (in m_earth)")
pl.plot(a, m2, 'o', markersize=5)


nom_fichier_plot3 = "histogrammes_res"
pl.figure(3)
#~ pl.clf()
pl.xlabel("Period ratio relative to the most massive planet")
pl.ylabel("density of probability")

pl.hist(period_ratio, bins=[0.5+0.0025*i for i in range(400)], normed=True, histtype='step')

nom_fichier_plot4 = "most_massives"
most_massive = [mi + random.uniform(-0.5, 0.5) for mi in most_massive]
second_massive = [mi + random.uniform(-0.5, 0.5) for mi in second_massive]

pl.figure(4)
#~ pl.clf()
pl.xlabel("most massive (m_earth)")
pl.ylabel("second most massive (m_earth)")
pl.plot(most_massive, second_massive, 'o', markersize=5)


pl.figure(5)
m_first = [mi + random.uniform(-0.25,0.25) for mi in m_first]
m_second = [mi + random.uniform(-0.25,0.25) for mi in m_second]
m_clo_first = [mi + random.uniform(-0.25,0.25) for mi in m_clo_first]
m_clo_second = [mi + random.uniform(-0.25,0.25) for mi in m_clo_second]
nom_fichier_plot5 = "coorbital_pl_mass"
pl.title("Mass of the coorbital planets")
pl.xlabel("mass of the first planet (m_earth)")
pl.ylabel("mass of the second planet (m_earth)")

pl.plot(m_clo_first, m_clo_second, 'ro', markersize=5, label="closest from CZ")
pl.plot(m_first, m_second, 'bo', markersize=5, label="all others")
ylims = list(pl.ylim())
xlims = list(pl.xlim())
limits = [max(xlims[0], ylims[0]), min(xlims[1], ylims[1])]
pl.plot(limits, limits, 'k--', label="equal mass")



print("Storing plots")
pl.figure(1)
pl.savefig(nom_fichier_plot1+"."+OUTPUT_EXTENSION)

pl.figure(2)
pl.savefig(nom_fichier_plot2+"."+OUTPUT_EXTENSION)

pl.figure(3)
#~ pl.axis('tight')
ylims = list(pl.ylim())
xlims = list(pl.xlim([0.6, 1.4]))
for res in resonances:
	#~ pdb.set_trace()
	nb_period = map(float, res.split(":")) # We get the two integers value of the resonance.
	ratio = nb_period[0] / nb_period[1]
	pl.plot([ratio, ratio], ylims, 'k--')
	pl.plot([1./ratio, 1./ratio], ylims, 'k--')
	pl.text(ratio, ylims[1], " "+res, horizontalalignment='center', verticalalignment='bottom', rotation='vertical', size=7)
	pl.text(1./ratio, ylims[1], " "+res, horizontalalignment='center', verticalalignment='bottom', rotation='vertical', size=7)
pl.savefig(nom_fichier_plot3+"."+OUTPUT_EXTENSION)

pl.figure(4)
pl.savefig(nom_fichier_plot4+"."+OUTPUT_EXTENSION)

pl.figure(5)
pl.legend()
pl.savefig(nom_fichier_plot5+"."+OUTPUT_EXTENSION)
	
pl.show()

# TODO tracer la zone de convergence et les résonnances sur les plots
