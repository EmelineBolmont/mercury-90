#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0
# Pour lire des fichiers de simulations, récupérer les caractéristiques
# des planètes qu'il reste en fin de simulation et les écrire dans un 
# seul fichier que le script "analyse_simu" va lire

import pylab as pl
import numpy as np

OUTPUT_EXTENSION = 'pdf'

####################
# On lit, pour chaque planète, le contenu du fichier et on stocke les variables qui nous intéressent.
####################
# a : demi-grand axe en ua
# H : scaleheight (in AU)
# h : aspect ratio h=H/R

# On récupère les données orbitales
(a, H, h) = np.loadtxt("scaleheight_profile.dat", dtype=float, unpack=True)

# on trace les plots
fig = pl.figure(1)
# On crée des sous plots. Pour subplot(231), ça signifie qu'on a 2 lignes, 3 colonnes, 
# et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
plot_H = fig.add_subplot(211)
plot_H.plot(a, H, 'k-')
plot_H.plot(a,-H, 'k-')

plot_H.set_xlabel("a [UA]")
plot_H.set_ylabel("Scaleheight H [AU]")
#~ pl.legend()
plot_H.grid(True)

plot_h = fig.add_subplot(212, sharex=plot_H)
plot_h.plot(a, h, 'k-')
plot_h.plot(a,-h, 'k-')

plot_h.set_xlabel("a [UA]")
plot_h.set_ylabel("aspect ratio h")
plot_h.grid(True)

nom_fichier_plot = "scaleheight_profile"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)


pl.show()

