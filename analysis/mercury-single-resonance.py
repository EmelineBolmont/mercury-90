#!/usr/bin/env python
# -*- coding: utf-8 -*-
# script to plot resonance between planets
# Version 1.3

#  this routine plots the resonant angles for a inner_period_nb:outer_period_nb resonance
#  (inner_period_nb>outer_period_nb)

# Parameters :
# inner_planet: name of the .aei file associated with the inner planet
# outer_planet: name of the .aei file associated with the outer planet
# inner_period_nb: number of periods for the inner planet
# outer_period_nb: number of periods for the outer planet

#~ inner_planet    = "PLANETE1"
#~ outer_planet    = "PLANETE2"
#~ inner_period_nb = 3
#~ outer_period_nb = 2

# HOW TO # There is a resonance if w1 - w2 is in libration.
# You have found the right resonance if at least one phi angle of the
# current resonance is in libration

# EXAMPLE
# To search for a 4:3 resonance between inner PLANET66 and outer PLANET68 :
# mercury-single-resonance.py PLANET66.aei PLANET68.aei 4 3


import pdb        # to debug via pdb.set_trace()
import os         # to create folder, change directory and so on
import subprocess # to launch 'runjob'
import autiwa
from math import *
import numpy as np
import pylab as pl
import sys # to use in particuliar the sys.argv list to retrieve parameters of the script
from matplotlib.ticker import ScalarFormatter

ANGLE_MIN = - 90.
ANGLE_MAX = 270.

OUTPUT_EXTENSION = 'png' # default value in bitmap, because vectoriel can take time and space if there is a lot of data



def lancer_commande(commande):
  """lance une commande qui sera typiquement soit une liste, soit une
  commande seule. La fonction renvoit un tuple avec la sortie,
  l'erreur et le code de retour"""
  if (type(commande)==list):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  elif (type(commande)==str):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  else:
    raise TypeError("This command is neither a list nor a string.")
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  return (process_stdout, process_stderr, returncode)


###############################################
## Beginning of the program
###############################################

# We get arguments from the script
if (len(sys.argv)-1) != 4:
  print("Not enough arguments. There must be 4\n" + \
        "1: name of the .aei file associated with the inner planet.\n" + \
        "2: name of the .aei file associated with the outer planet.\n" + \
        "3: the number of inner period.\n" + \
        "4: the number of outer period.\n\n" + \
        "Example : mercury-resonances.py PLANETE1.aei PLANETE2.aei 3 2\n" + \
        "will search for a 3:2 resonance between planet 1 and 2")
  exit()

inner_planet = sys.argv[1].rstrip(".aei")
outer_planet = sys.argv[2].rstrip(".aei")
inner_period_nb = int(sys.argv[3])
outer_period_nb = int(sys.argv[4])



####################
# On récupère la liste des fichiers planètes.aei
####################
(process_stdout, process_stderr, return_code) = lancer_commande("ls *.aei")
if (return_code != 0):
  print("the command return an error "+str(return_code))
  print(process_stderr)
  exit()

liste_aei = process_stdout.split("\n")
liste_aei.remove('') # we remove an extra element that doesn't mean anything
nb_planete = len(liste_aei)

if not(liste_aei.count(inner_planet+".aei")):
  print("The planet "+inner_planet+" does not exist")
  print("Available planets are :")
  print(liste_aei)
  exit()

if not(liste_aei.count(outer_planet+".aei")):
  print("The planet "+outer_planet+" does not exist")
  print("Available planets are :")
  print(liste_aei)
  exit()



#~ pdb.set_trace()

####################
# On lit, pour chaque planète, le contenu du fichier et on stocke les variables qui nous intéressent.
####################
t_inner = [] # time in years
a_inner = [] # demi-grand axe en ua
g_inner = [] # argument of pericentre (degrees)
n_inner = [] # longitude of ascending node (degrees)
M_inner = [] # Mean anomaly (degrees)

t_outer = [] # time in years
a_outer = [] # demi-grand axe en ua
g_outer = [] # argument of pericentre (degrees)
n_outer = [] # longitude of ascending node (degrees)
M_outer = [] # Mean anomaly (degrees)


values = [[t_inner, a_inner, g_inner, n_inner, M_inner], [t_outer, a_outer, g_outer, n_outer, M_outer]]
filenames = [inner_planet, outer_planet]
for ind in range(2):

  source_file = filenames[ind]+".aei"

  sys.stdout.write("Reading data files %5.1f %% : %s                \r" % ((ind) * 25., source_file))
  sys.stdout.flush()

  tableau = open(source_file, 'r')

  # we skip the four first lines
  for indice in range(4):
      tableau.readline()

  tapp = values[ind][0].append
  aapp = values[ind][1].append
  gapp = values[ind][2].append
  napp = values[ind][3].append
  Mapp = values[ind][4].append

  for ligne in tableau:
    colonne = ligne.split()

    tapp(float(colonne[0]))
    aapp(float(colonne[1]))
    gapp(float(colonne[4]))
    napp(float(colonne[5]))
    Mapp(float(colonne[6]))
  tableau.close()


sys.stdout.write("Calculating angles %5.1f %%                          \r" % (50.))
sys.stdout.flush()

# We search for an ejection and the the low numbered list between inner and outer planet
nb_times = min(len(t_inner), len(t_outer))

t = np.array(t_inner[:nb_times]) # no matter inner or outer since we take the minimun of the two

a_inner = np.array(a_inner[:nb_times])
g_inner = np.array(g_inner[:nb_times])
n_inner = np.array(n_inner[:nb_times])
M_inner = np.array(M_inner[:nb_times])

a_outer = np.array(a_outer[:nb_times])
g_outer = np.array(g_outer[:nb_times])
n_outer = np.array(n_outer[:nb_times])
M_outer = np.array(M_outer[:nb_times])

# Resonances are usually displayed as (p+q):p where q is the order of
# the resonance. We retreive thoses parameters
p = outer_period_nb
q = inner_period_nb - outer_period_nb

# We calculate the resonant angles
long_of_peri_inner = g_inner + n_inner
mean_longitude_inner = M_inner + long_of_peri_inner

long_of_peri_outer = g_outer + n_outer
mean_longitude_outer = M_outer + long_of_peri_outer

phi = np.empty((q+1, nb_times)) # we create an array, empty for the moment, that will contain all the resonant angles associated with the supposed resonance.

temp_value = inner_period_nb * mean_longitude_outer - outer_period_nb * mean_longitude_inner

for i in range(q+1):
  phi[i] = temp_value - i * long_of_peri_inner - (q - i) * long_of_peri_outer

# We take modulo 2*pi of the resonant angle
phi = phi%(360.)
too_low = phi < ANGLE_MIN
too_high = phi > ANGLE_MAX
phi[too_low] = phi[too_low] + 360.
phi[too_high] = phi[too_high] - 360.

delta_longitude = long_of_peri_outer - long_of_peri_inner
delta_longitude = delta_longitude%(360.)
too_low = delta_longitude < ANGLE_MIN
too_high = delta_longitude > ANGLE_MAX
delta_longitude[too_low] = delta_longitude[too_low] + 360.
delta_longitude[too_high] = delta_longitude[too_high] - 360.


# We chose the number of plot in x and y axis for the p.multi
# environment in order to plot ALL the resonant angles and have x and
# y numbers as close on from another as possible. There are q+1
# resonant angles, the period ratio and w1 - w2, i.e q+3 plots
(nb_lines, nb_rows) = autiwa.get_subplot_shape(q+3)

# on trace les plots
myxfmt = ScalarFormatter(useOffset=True)
myxfmt._set_offset(1e5)
myxfmt.set_scientific(True)
myxfmt.set_powerlimits((-3, 3))

sys.stdout.write("Displaying %5.1f %%                          \r" % (75.))
sys.stdout.flush()

fig = pl.figure(1)
pl.clf()
fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)

# We create a variable to store the index of the subplot
subplot_index = 0
fig.suptitle("resonance "+str(inner_period_nb)+":"+str(outer_period_nb)+" between "+inner_planet+" and "+outer_planet)

subplot_index += 1
plot_period = fig.add_subplot(nb_lines, nb_rows, subplot_index)
plot_period.plot(t, (a_outer / a_inner)**1.5)
plot_period.set_xlabel("time [years]")
plot_period.set_ylabel("period ratio")
#~ pl.legend()
plot_period.grid(True)

# For each resonant angle, it is important to show only the points and not lines,
# so that we do not mask interesting features by meaningless lines.
subplot_index += 1
plot_dl = fig.add_subplot(nb_lines, nb_rows, subplot_index, sharex=plot_period)
plot_dl.plot(t, delta_longitude, '.')
plot_dl.set_xlabel("time [years]")
plot_dl.set_ylabel("w2 - w1")
plot_dl.grid(True)

for i in range(q+1):
  sys.stdout.write("Displaying %5.1f %%                          \r" % (75. + float((i+3) / (q+3)) * 25.))
  sys.stdout.flush()
  subplot_index += 1
  plot_phi = fig.add_subplot(nb_lines, nb_rows, subplot_index, sharex=plot_period)
  plot_phi.plot(t, phi[i], '.')
  plot_phi.set_xlabel("time [years]")
  plot_phi.set_ylabel(unicode("φ%i" % i, 'utf8'))
  plot_phi.grid(True)

plot_period.xaxis.set_major_formatter(myxfmt)

nom_fichier_plot = "resonance_%i_%i_between_%s_and_%s" % (inner_period_nb, outer_period_nb, inner_planet, outer_planet)
fig.savefig('%s.%s' % (nom_fichier_plot, OUTPUT_EXTENSION), format=OUTPUT_EXTENSION)

pl.show() # We display the plot
