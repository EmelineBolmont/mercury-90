#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0
# Pour lire des fichiers de simulations, récupérer les caractéristiques
# des planètes qu'il reste en fin de simulation et les écrire dans un 
# seul fichier que le script "analyse_simu" va lire

import pylab as pl
import numpy as np
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
import sys

OUTPUT_EXTENSION = 'pdf'

####################
# We read OPTIONS
####################
isProblem = False
problem_message = "AIM : Plot the scaleheight profile of the simulation." + "\n" + \
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
# H : scaleheight (in AU)
# h : aspect ratio h=H/R

# On récupère les données orbitales
(a, H, h) = np.loadtxt("scaleheight_profile.dat", dtype=float, unpack=True)

# on trace les plots
fig = pl.figure()
# We create subplots. add_subplot(2, 3, 1) means we have 2 lines, 3 columns, 
# and that the active plot is the first, starting from top left (for 6 plots in total)
plot_H = fig.add_subplot(2, 1, 1)
plot_H.loglog(a, H, 'k-')
#~ plot_H.plot(a,-H, 'k-')

plot_H.set_xlabel("Semi-major axis [AU]")
plot_H.set_ylabel("Scaleheight H [AU]")
plot_H.xaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_H.yaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_H.xaxis.grid(True, which='minor', color='#888888', linestyle='-')
plot_H.yaxis.grid(True, which='minor', color='#888888', linestyle='-')

plot_h = fig.add_subplot(2, 1, 2, sharex=plot_H)
plot_h.semilogx(a, h, 'k-')
#~ plot_h.plot(a,-h, 'k-')

plot_h.set_xlabel("Semi-major axis [AU]")
plot_h.set_ylabel("aspect ratio h")
plot_h.xaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_h.yaxis.grid(True, which='major', color='#222222', linestyle='-')
plot_h.xaxis.grid(True, which='minor', color='#888888', linestyle='-')
#~ plot_h.yaxis.grid(True, which='minor', color='#888888', linestyle='-')

myxfmt = FormatStrFormatter("%.3g")
plot_H.xaxis.set_major_formatter(myxfmt)

nom_fichier_plot = "scaleheight_profile"
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)


pl.show()

