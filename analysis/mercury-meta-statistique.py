#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0.1.1
# Pour gÃ©nÃ©rer des histogrammes sur plusieurs simulations afin de 
# regarder les caractÃ©ristiques statistiques des planÃ¨tes qui restent.


# TODO
# _ Ã©crire un fichier qui stocke les paramÃ¨tres orbitaux des planÃ¨tes qu'il reste.
# _ sÃ©lectionner les planÃ¨tes qui sont plus proche qu'une certaine valeur de leur Ã©toile.

from math import *
import pylab as pl
import os, pdb, autiwa
from simu_constantes import *
import random
import numpy as np

# Mass threshold
MASS_THRESHOLD = 5 # Earth mass

CONVERGENCE_ZONE = 3. # in AU, the location of the convergence zone (as a reference for plotting period ratios)

DELTA_RATIO = 0.005
DELTA_M = 0.1

# Get current working directory
rep_exec = os.getcwd()

resonances = ["3:2", "4:3", "5:4", "6:5", "7:6", "8:7", "9:8", "10:9", "11:10"]

#######################
# On prÃ©pare les plots
#######################

# Extensions voulues pour les fichiers de sortie
extensions = ['pdf', 'png', 'svg']

dossier_plot = "output"

#######################
# We prepare the various folders we we will read the simulations
#######################
# We define a list of folder names we don't want to read. 
# Theses names will be suppressed from the final list, in case they exist.
dossier_suppr = ["output", "indiv_simu_01"]

#######################
# On se place dans le dossier de simulation souhaitÃ©
#######################
liste_meta_simu = [dir for dir in os.listdir(".") if os.path.isdir(dir)]
autiwa.suppr_dossier(liste_meta_simu,dossier_suppr)
liste_meta_simu.sort()

nb_meta_simu = len(liste_meta_simu)

# We generate a list of colors
colors = [ '#'+li for li in autiwa.colorList(nb_meta_simu, exclude=['000000'])]

# We chose the number of plot in x and y axis for the p.multi
# environment in order to plot ALL the resonant angles and have x and
# y numbers as close on from another as possible. There are q+1
# resonant angles, the period ratio and w1 - w2, i.e q+3 plots
nb_plots_x = 1
nb_plots_y = 1
while (nb_plots_x * nb_plots_y < nb_meta_simu):
  if (nb_plots_x == nb_plots_y):
    nb_plots_y += 1
  else:
    nb_plots_x += 1
subplot_index = nb_plots_x * 100 + nb_plots_y * 10

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

# We go in each sub folder of the current working directory

for (meta_index, meta_simu) in enumerate(liste_meta_simu):
	print("Traitement de %s"%meta_simu)
	os.chdir(os.path.join(rep_exec,meta_simu))
	
	meta_prefix = meta_simu # it is thought to have the initial number of planets in here
	

	# On rÃ©cupÃ¨re la liste des sous-dossiers
	liste_simu = [dir for dir in os.listdir(".") if os.path.isdir(dir)]
	autiwa.suppr_dossier(liste_simu,dossier_suppr)
	liste_simu.sort()
	

	nb_simu = len(liste_simu)
	
	# We initialize the variables where to store datas
	
	intial_nb_planets = int(meta_prefix)
	
	# List of orbital elements of ALL the simulations
	a = [] # in AU
	e = [] # eccentricity
	I = [] # inclination in degre??
	m = [] # mass in earth mass
	m_relat = [] # mass expressed in function of the mass of the most massive planet of the system
	
	# List of orbital elements of the closest planet from the convergence zone
	a_clo = [] # in AU
	e_clo = [] # eccentricity
	I_clo = [] # inclination in degre??
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
	
	print("\t Reading datas")
	
	for simu in liste_simu:
		os.chdir(os.path.join(rep_exec, meta_simu, simu))
		dataFile = open("element.out", 'r')
		
		# List initialized for each new simulation that contains values of the current simulation
		a_system = [] # in AU
		e_system = [] # eccentricity
		I_system = [] # inclination in degre??
		m_system = [] # mass in earth mass
		
		# We get rid of the header
		for i in range(5):
			dataFile.readline()
		
		lines = dataFile.readlines()
		dataFile.close()
		
		nb_planets = len(lines)
		final_nb_planets.append(nb_planets)
		
		for line in lines:
			datas = line.split()
			a_system.append(float(datas[1]))
			e_system.append(float(datas[2]))
			I_system.append(float(datas[3]))
			m_system.append(float(datas[4]))
		
		# We add all the elements of the system into the list of values corresponding to all the simulations
		a.extend(a_system)
		e.extend(e_system)
		I.extend(I_system)
		m.extend(m_system)
		
		m_max = max(m_system)
		m_relat.extend([min((mi + random.uniform(-0.5,0.5))/m_max,1.) for mi in m_system])
		
		###########
		# Statistical studies
		
		# We search for the closest planet from the convergence zone
		a_system = np.array(a_system)
		a_temp = abs(a_system - CONVERGENCE_ZONE)
		idx_clo = a_temp.argmin()
		a_ref = a_system[idx_clo]
		
		a_clo.append(a_system[idx_clo])
		e_clo.append(e_system[idx_clo])
		I_clo.append(I_system[idx_clo])
		m_clo.append(m_system[idx_clo])
		
		
		
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
		try:
			mass_first = tmp[-1]
			mass_second = tmp[-2]
			most_massive.append(mass_first)
			second_massive.append(mass_second)
			#~ if (abs(mass_first - 16) < DELTA_M):
				#~ print("meta_simu :"+meta_prefix+"\t simu :"+simu)
		except:
			pass
		
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
				
	
	
	#######################
	#   TracÃ© des plots   #
	#######################
	print("\t Computing Plots")
	os.chdir(rep_exec)
	nb_bins = 50
	
	nom_fichier_plot = [] # list of names for each plot

	nom_fichier_plot1 = "histogrammes_m"
	pl.figure(1)
	pl.xlabel(unicode("masse (en mj)",'utf-8'))
	pl.ylabel("density of probability")
	pl.hist(m_clo, bins=range(25), normed=True, histtype='step', label='nb='+meta_prefix)

	nom_fichier_plot2 = "e_fct_m"
	#~ m2 = [mi + random.uniform(-0.5,0.5) for mi in m]
	pl.figure(2)
	pl.xlabel("mass [Earths]")
	pl.ylabel("eccentricity")
	pl.plot(m, e, 'o', markersize=5, label='nb='+meta_prefix)
	
	nom_fichier_plot3 = "e_fct_a"
	#~ dist = [ai - CONVERGENCE_ZONE for ai in a]
	pl.figure(3)
	pl.xlabel("distance [AU]")
	pl.ylabel("eccentricity")
	pl.plot(a, e, 'o', markersize=5, label='nb='+meta_prefix)

	nom_fichier_plot4 = 'histogrammes_I'
	pl.figure(4)
	pl.xlabel(unicode("I (in degrees)",'utf-8'))
	pl.ylabel("density of probability")
	pl.hist(I, bins=[0.002*i for i in range(25)], normed=True, histtype='step', label='nb='+meta_prefix)
	
	nom_fichier_plot5 = 'histogrammes_nb_pl'
	pl.figure(5)
	pl.xlabel(unicode("nb_final",'utf-8'))
	pl.ylabel("density of probability")
	pl.hist(final_nb_planets, bins=range(25), histtype='step', label='nb='+meta_prefix)
	
	
	nom_fichier_plot6 = "m_fct_a"
	m2 = [mi + random.uniform(-0.5,0.5) for mi in m]
	pl.figure(6)
	subplot_index += 1
	pl.subplot(subplot_index)
	pl.xlabel(unicode("a [AU]",'utf-8'))
	pl.ylabel("mass [Earths]")
	pl.plot(a, m2, 'o', markersize=2, color=colors[meta_index], label='nb='+meta_prefix)
	pl.ylim(0, 12)
	#~ pl.xlim(1, 10)
	pl.legend()
	
	
	nom_fichier_plot7 = "histogrammes_res"
	pl.figure(7)
	#~ pl.clf()
	pl.xlabel("Period ratio relative to the closest planet of the system from the convergence zone")
	pl.ylabel("density of probability")
	
	pl.hist(period_ratio, bins=[0.5+0.0025*i for i in range(400)], normed=True, histtype='step', label='nb='+meta_prefix)
	
	nom_fichier_plot8 = "most_massives"
	most_massive = [mi + random.uniform(-0.5, 0.5) for mi in most_massive]
	second_massive = [mi + random.uniform(-0.5, 0.5) for mi in second_massive]
	
	pl.figure(8)
	#~ pl.clf()
	pl.xlabel("most massive [Earths]")
	pl.ylabel("second most massive [Earths]")
	pl.plot(most_massive, second_massive, 'o', markersize=5, label='nb='+meta_prefix)
	
	#~ pl.figure(9)
	#~ if (meta_prefix == '50'):
		#~ m_first = [mi + random.uniform(-0.25,0.25) for mi in m_first]
		#~ m_second = [mi + random.uniform(-0.25,0.25) for mi in m_second]
		#~ m_clo_first = [mi + random.uniform(-0.25,0.25) for mi in m_clo_first]
		#~ m_clo_second = [mi + random.uniform(-0.25,0.25) for mi in m_clo_second]
		#~ nom_fichier_plot9 = "coorbital_pl_mass"
		#~ pl.title("Mass of the coorbital planets")
		#~ pl.xlabel("mass of the first planet [Earths]")
		#~ pl.ylabel("mass of the second planet [Earths]")
		#~ 
		#~ pl.plot(m_clo_first, m_clo_second, 'ro', markersize=5, label="closest from CZ")
		#~ pl.plot(m_first, m_second, 'bo', markersize=5, label="all others")
		#~ ylims = list(pl.ylim())
		#~ xlims = list(pl.xlim())
		#~ limits = [max(xlims[0], ylims[0]), min(xlims[1], ylims[1])]
		#~ pl.plot(limits, limits, 'k--', label="equal mass")

fig = pl.figure(6)
ax1 = fig.add_subplot(221)
#~ ax1.set_xticklabels([])

ax1 = fig.add_subplot(222)
#~ ax1.set_xticklabels([])
#~ ax1.set_yticklabels([])

ax1 = fig.add_subplot(223)


ax1 = fig.add_subplot(224)
#~ ax1.set_yticklabels([])

print("Storing plots")
for ext in [".png"]:
	pl.figure(1)
	pl.legend()
	pl.savefig(nom_fichier_plot1+ext)
	
	pl.figure(2)
	pl.legend()
	pl.savefig(nom_fichier_plot2+ext)
	
	pl.figure(3)
	pl.legend()
	pl.savefig(nom_fichier_plot3+ext)
	
	pl.figure(4)
	pl.legend()
	pl.savefig(nom_fichier_plot4+ext)
	
	pl.figure(5)
	pl.legend()
	pl.savefig(nom_fichier_plot5+ext)
	
	pl.figure(6)
	#~ pl.subplots_adjust(hspace = 0, wspace = 0)
	pl.legend()
	pl.savefig(nom_fichier_plot6+ext)
	
	pl.figure(7)
	#~ pl.axis('tight')
	ylims = list(pl.ylim())
	xlims = list(pl.xlim([0.6, 1.4]))
	#~ pl.text(xlims[0], 0.9*ylims[1], " For "+meta_prefix+" planets", horizontalalignment='left', verticalalignment='top', size = 15)
	for res in resonances:
		#~ pdb.set_trace()
		nb_period = map(float, res.split(":")) # We get the two integers value of the resonance.
		ratio = nb_period[0] / nb_period[1]
		pl.plot([ratio, ratio], ylims, 'k--')
		pl.plot([1./ratio, 1./ratio], ylims, 'k--')
		pl.text(ratio, ylims[1], " "+res, horizontalalignment='center', verticalalignment='bottom', rotation='vertical', size=7)
		pl.text(1./ratio, ylims[1], " "+res, horizontalalignment='center', verticalalignment='bottom', rotation='vertical', size=7)
	pl.legend()
	pl.savefig(nom_fichier_plot7+".pdf")
	
	pl.figure(8)
	pl.legend()
	pl.savefig(nom_fichier_plot8+ext)
	
	#~ pl.figure(9)
	#~ pl.legend()
	#~ pl.savefig(nom_fichier_plot9+ext)


pl.show()

# TODO tracer la zone de convergence et les rÃ©sonnances sur les plots
