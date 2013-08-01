#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.1
# Pour lire des fichiers de simulations, récupérer les caractéristiques
# des planètes qu'il reste en fin de simulation et les écrire dans un 
# seul fichier que le script "analyse_simu" va lire

import pylab as pl
import numpy as np
from matplotlib.ticker import FormatStrFormatter


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
plot_density.loglog(a, sigma)

plot_density.set_xlabel("Semi-major axis [AU]")
plot_density.set_ylabel("Surface density [g/cm^2]")
#~ plot_density.axis('tight')
#~ plot_density.set_xlim(xmin=0)
#~ plot_density.set_ylim(ymin=0)
plot_density.xaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_density.yaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_density.xaxis.grid(True, which='minor', color='#888888', linestyle='-')
plot_density.yaxis.grid(True, which='minor', color='#888888', linestyle='-')

plot_idx = fig.add_subplot(212, sharex=plot_density)
plot_idx.semilogx(a, index)

plot_idx.set_xlabel("Semi-major axis [AU]")
plot_idx.set_ylabel("density power law")
plot_idx.set_ylim((-2,2))
plot_idx.xaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_idx.yaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_idx.xaxis.grid(True, which='minor', color='#888888', linestyle='-')

myxfmt = FormatStrFormatter("%.3g")
plot_density.xaxis.set_major_formatter(myxfmt)

# We generate the output file
nom_fichier_plot = "density_profile"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)


pl.show()

