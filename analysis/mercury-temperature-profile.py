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
# tau : optical depth
# chi : thermal diffusivity
# index : temperature index 
# a : demi-grand axe en ua
# T : temperature in K

# On récupère les données orbitales
(a, T, index, tau, chi) = np.loadtxt("temperature_profile.dat", skiprows=1, usecols = (0,1,2,3,4), dtype=float, unpack=True)

# on trace les plots

fig = pl.figure(1)
# On crée des sous plots. Pour subplot(231), ça signifie qu'on a 2 lignes, 3 colonnes, et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
plot_T = fig.add_subplot(211)
plot_T.loglog(a, T)

plot_T.set_xlabel("a [UA]")
plot_T.set_ylabel("Temperature [K]")
plot_T.grid(True)

plot_idx = fig.add_subplot(212, sharex=plot_T)
plot_idx.semilogx(a, index)

plot_idx.set_xlabel("a [UA]")
plot_idx.set_ylabel("Temperature power law")
plot_idx.grid(True)

# We generate the output file
nom_fichier_plot = "temperature_profile"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

