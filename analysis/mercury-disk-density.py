#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v0.1
# Pour lire des fichiers de simulations, récupérer les caractéristiques
# des planètes qu'il reste en fin de simulation et les écrire dans un 
# seul fichier que le script "analyse_simu" va lire

import pylab as pl
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import sys

OUTPUT_EXTENSION = 'pdf'

####################
# We read OPTIONS
####################
isProblem = False
problem_message = "AIM : Plot the density profile of the simulation." + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * ext=%s : The extension for the output files" % OUTPUT_EXTENSION + "\n" + \
" * help : Display a little help message on HOW to use various options"

value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'help'):
    isProblem = True
    if (value != None):
      print(value_message % (key, key, value))
  else:
    print("the key '%s' does not match" % key)
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

####################
# On lit, pour chaque planète, le contenu du fichier et on stocke les variables qui nous intéressent.
####################
# a : demi-grand axe en ua
# sigma : surface density (in g/cm^2)
# index : surface density index

# We retrieve the datas
(a, sigma, index) = np.loadtxt("density_profile.dat", skiprows=1, dtype=float, unpack=True)

# on trace les plots

fig = pl.figure()
plot_density = fig.add_subplot(2, 1, 1)
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

plot_idx = fig.add_subplot(2, 1, 2, sharex=plot_density)
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

