#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.1
# Pour lire des fichiers de simulations, récupérer les caractéristiques
# des planètes qu'il reste en fin de simulation et les écrire dans un 
# seul fichier que le script "analyse_simu" va lire

from math import *
import pylab as pl
import os, pdb, autiwa
import numpy as np
from constants import MT, MS
import sys # to be able to retrieve arguments of the script

###############################################
## Beginning of the program
###############################################

# We get arguments from the script
try:
	#~ pdb.set_trace()
	if (len(sys.argv)-1) == 1:
		t_max = float(sys.argv[1])
	elif (len(sys.argv)-1) == 2:
		t_min = float(sys.argv[1])
		t_max = float(sys.argv[2])
except:
	print("If one argument, then this must be the max time for output. \nIf two arguments, the first is t_min and the second is t_max")
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
m = [] # masse de la planète en masse terrestre



# On récupère les données orbitales
for planete in range(nb_planete):
	
	fichier_source = liste_aei[planete]
	tableau = open(fichier_source, 'r')
	
	t.append([]) # temps en année
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
			mi = float(colonne[7])
		except:
			pass
		t[-1].append(ti)
		m[-1].append(mi * MS / MT)

		
	tableau.close()

dt = t[0][1] - t[0][0]
if ('t_max' in locals()):
	id_max = int((t_max - t[0][0]) / dt)
else:
	#~ t_max = max([max(ti) for ti in i])
	id_max = max([len(ti) for ti in t]) - 1
	t_max = t[0][0] + id_max * dt

if ('t_min' in locals()):
	id_min = int((t_min - t[0][0]) / dt)
else:
	t_min = t[0][0]
	id_min = 0
	
# on trace les plots

pl.figure(1)
pl.clf()

for planet in range(nb_planete):
	pl.plot(t[planet][id_min:id_max], m[planet][id_min:id_max], label='PLANETE'+str(planet))
pl.xlim([t_min, t_max])
pl.xlabel(unicode("time [years]",'utf-8'))
pl.ylabel(unicode("mass [mt]",'utf-8'))
#~ pl.legend()
pl.grid(True)


#~ dossier_output = "output"
#~ system("mkdir dossier_output")
#~ system("cd dossier_output")

pl.figure(1)
nom_fichier_plot = "mass"
#~ pl.savefig(nom_fichier_plot+'.svg', format='svg')
pl.savefig(nom_fichier_plot+'.pdf', format='pdf')

pl.show()

