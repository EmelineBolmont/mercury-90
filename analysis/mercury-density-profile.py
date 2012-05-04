#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.1
# Pour lire des fichiers de simulations, récupérer les caractéristiques
# des planètes qu'il reste en fin de simulation et les écrire dans un 
# seul fichier que le script "analyse_simu" va lire

import pylab as pl
import numpy as np

####################
# On lit, pour chaque planète, le contenu du fichier et on stocke les variables qui nous intéressent.
####################
a = [] # demi-grand axe en ua
Sigma = [] # surface density (in g/cm^2)
index = [] # surface density index


# On récupère les données orbitales


tableau = open("density_profile.dat", 'r')

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
		Sigmai = float(colonne[1])
		#~ indexi = float(colonne[2])
	except:
		pass

	a.append(ai)       # demi-grand axe en ua
	Sigma.append(Sigmai) # surface density (in g/cm^2)
	#~ index.append(indexi) # surface density index
tableau.close()

a = np.array(a)
Sigma = np.array(Sigma)
#~ index = np.array(index)

# on trace les plots


pl.figure(1)
pl.clf()
pl.plot(a, Sigma)
pl.axis('tight')
pl.xlim(xmin=0)
pl.ylim(ymin=0)

pl.xlabel(unicode("a [UA]",'utf-8'))
pl.ylabel(unicode("surface density [g/cm^2]",'utf-8'))
#~ pl.legend()
pl.grid(True)



#~ dossier_output = "output"
#~ system("mkdir dossier_output")
#~ system("cd dossier_output")

pl.figure(1)
nom_fichier_plot = "density_profile"
#~ pl.savefig(nom_fichier_plot+'.svg', format='svg')
pl.savefig(nom_fichier_plot+'.pdf', format='pdf')

pl.show()

