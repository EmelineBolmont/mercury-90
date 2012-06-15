#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.2
# Pour lire des fichiers de simulations, récupérer les caractéristiques
# des planètes qu'il reste en fin de simulation et les écrire dans un 
# seul fichier que le script "analyse_simu" va lire

import pylab as pl

####################
# On lit, pour chaque planète, le contenu du fichier et on stocke les variables qui nous intéressent.
####################
tau = [] # optical depth
chi = [] # thermal diffusivity
index = [] # temperature index 
a = [] # demi-grand axe en ua
T = [] # temperature in K

# On récupère les données orbitales


tableau = open("temperature_profile.dat", 'r')

# On passe la 1e ligne d'entête.
for indice in range(1):
	tableau.readline()

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
		ai = float(colonne[0])
		Ti = float(colonne[1])
		indexi = float(colonne[2])
		tau_i = float(colonne[3])
		chi_i = float(colonne[4])
	except:
		pass
	tau.append(tau_i) # optical depth
	chi.append(chi_i) # thermal diffusivity
	index.append(indexi) # temperature index 
	a.append(ai)       # demi-grand axe en ua
	T.append(Ti)       # temperature in K

tableau.close()

# on trace les plots

pl.figure(1)
pl.clf()
# On crée des sous plots. Pour subplot(231), ça signifie qu'on a 2 lignes, 3 colonnes, et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
ax1 = pl.subplot(211)
pl.plot(a, T)

pl.xlabel(unicode("a [UA]",'utf-8'))
pl.ylabel(unicode("Temperature [K]",'utf-8'))
#~ pl.legend()
pl.grid(True)

pl.subplot(212, sharex=ax1)
pl.plot(a, index)

pl.xlabel(unicode("a [UA]",'utf-8'))
pl.ylabel(unicode("Temperature power law",'utf-8'))
pl.grid(True)



#~ dossier_output = "output"
#~ system("mkdir dossier_output")
#~ system("cd dossier_output")

pl.figure(1)
nom_fichier_plot = "temperature_profile"
#~ pl.savefig(nom_fichier_plot+'.svg', format='svg')
pl.savefig(nom_fichier_plot+'.pdf', format='pdf')

pl.show()

