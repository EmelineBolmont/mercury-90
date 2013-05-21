#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.1
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
# sigma : surface density (in g/cm^2)
# index : surface density index

# We retrieve the datas
(a, sigma, index) = np.loadtxt("density_profile.dat", skiprows=1, dtype=float, unpack=True)

# on trace les plots

fig = pl.figure(1)
plot_density = fig.add_subplot(211)
plot_density.plot(a, sigma)

plot_density.set_xlabel("a [UA]")
plot_density.set_ylabel("surface density [g/cm^2]")
plot_density.axis('tight')
plot_density.set_xlim(xmin=0)
plot_density.set_ylim(ymin=0)
plot_density.grid(True)

plot_idx = fig.add_subplot(212, sharex=plot_density)
plot_idx.plot(a, index)

plot_idx.set_xlabel("a [UA]")
plot_idx.set_ylabel("density power law")
plot_idx.set_ylim((-2,2))
plot_idx.grid(True)

# We generate the output file
nom_fichier_plot = "density_profile"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)


pl.show()

