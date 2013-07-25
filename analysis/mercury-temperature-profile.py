#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0
# Pour lire des fichiers de simulations, récupérer les caractéristiques
# des planètes qu'il reste en fin de simulation et les écrire dans un 
# seul fichier que le script "analyse_simu" va lire

import pylab as pl
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter


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
(a, T, index, chi, tau) = np.loadtxt("temperature_profile.dat", skiprows=1, usecols = (0,1,2,3,4), dtype=float, unpack=True)

# on trace les plots
fig = pl.figure(1)
# On crée des sous plots. Pour subplot(231), ça signifie qu'on a 2 lignes, 3 colonnes, et que le subplot courant est le 1e. (on a donc 2*3=6 plots en tout)
plot_T = fig.add_subplot(211)
plot_T.loglog(a, T, label="my profile")

(inner_edge, outer_edge) = plot_T.get_xlim()
rect_width = outer_edge - inner_edge
snowline_T_min = 145. # K
snowline_T_max = 170. # K
snowline_T_mid = (snowline_T_max + snowline_T_min)/2.
snowline = Rectangle((inner_edge, snowline_T_min), rect_width, (snowline_T_max - snowline_T_min), color='#dddddd')

plot_T.add_patch(snowline)
plot_T.text(outer_edge, snowline_T_mid, "Snow line", horizontalalignment='right', verticalalignment='bottom', size=10)
plot_T.plot([inner_edge, outer_edge], [snowline_T_mid, snowline_T_mid], color='#000000')

plot_T.set_xlabel("Semi-major axis [AU]")
plot_T.set_ylabel("Temperature [K]")
plot_T.xaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_T.yaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_T.xaxis.grid(True, which='minor', color='#888888', linestyle='-')
plot_T.yaxis.grid(True, which='minor', color='#888888', linestyle='-')

plot_idx = fig.add_subplot(212, sharex=plot_T)
plot_idx.semilogx(a, index)

plot_idx.set_xlabel("Semi-major axis [AU]")
plot_idx.set_ylabel("Temperature power law")
plot_idx.xaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_idx.yaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_idx.xaxis.grid(True, which='minor', color='#888888', linestyle='-')


myxfmt = FormatStrFormatter("%.3g")
plot_T.xaxis.set_major_formatter(myxfmt)

# We generate the output file
nom_fichier_plot = "temperature_profile"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show()

