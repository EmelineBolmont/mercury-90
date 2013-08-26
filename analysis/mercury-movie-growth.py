#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.1
# To display the orbits of the planets in the system in the (x,y) plane

from math import *
import pylab as pl
import os, pdb, autiwa
import numpy as np
import sys # to be able to retrieve arguments of the script
import mercury
from matplotlib.lines import Line2D

# We get the path toward the binaries
scriptFolder = os.path.dirname(os.path.realpath(__file__)) # the folder in which the module is. 
BINARY_FOLDER = os.path.join(scriptFolder, os.path.pardir)

FRAME_PREFIX = "frame_growth_"
OUTPUT_EXTENSION = 'png'
OUTPUT_FOLDER = "movie"


NB_POINTS = 50 # Number of points for the display or circles
NB_FRAMES = 2

###############################################
## Beginning of the program
###############################################


# We get arguments from the script
isDisk = True

#~ pdb.set_trace()
isProblem = False
problem_message = "AIM : Movie in a m = f(a) diagram, to see how the planet grow" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * tmax=1e6 : the end of the output [years]" + "\n" + \
" * tmin=1e3 : the beginning of the output [years]" + "\n" + \
" * range=1. : the farthest location in the disk that will be displayed (in AU)" + "\n" + \
" * nodisk : to avoid torque diagram display" + "\n" + \
" * frames=1 : the number of frames you want" + "\n" + \
" * ext=%s : The extension for the output files" % OUTPUT_EXTENSION + "\n" + \
" * help : display a little help message on HOW to use various options"

value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'tmin'):
    t_min = float(value)
  elif (key == 'tmax'):
    t_max = float(value)
  elif (key == 'nodisk'):
    isDisk = False
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'frames'):
    NB_FRAMES = int(value)
  elif (key == 'range'):
    plot_range = float(value)
  elif (key == 'ext'):
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

if not(os.path.exists(OUTPUT_FOLDER)):
    os.mkdir(OUTPUT_FOLDER)

if isDisk:
  (process_stdout, process_stderr, returncode) = autiwa.lancer_commande(os.path.join(BINARY_FOLDER, "torque_diagram"))
  if (returncode != 0):
    print(process_stdout)
    print(process_stderr)
    print("Process terminated with error %d" % returncode)
    pdb.set_trace()
    
  
  # We read the contour of the zero torque zones
  contour_filename = "contour_total_torque.dat"
  
  contours_a = [[]]
  contours_m = [[]]
  contours_aapp = contours_a[-1].append
  contours_mapp = contours_m[-1].append

  contour_file = open(contour_filename, 'r')
  for line in contour_file:
    if (line == ' \n'):
      contours_a.append([])
      contours_m.append([])
      contours_aapp = contours_a[-1].append
      contours_mapp = contours_m[-1].append
    else:
      (a, m, dummy) = line.split()
      contours_aapp(float(a))
      contours_mapp(float(m))
  
  # We delete all the empty arrays created in the list of contours (one between each contour)
  while True:
    try:
      contours_a.remove([])
      contours_m.remove([])
    except:
      break

####################
# On recupere la liste des fichiers planetes.aei
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
# On lit, pour chaque planete, le contenu du fichier et on stocke les variables qui nous interessent.
####################
t = [] # temps en annee
a = [] # the smei major axis in AU
m = [] # mass in earth mass


# On recupere les donnees orbitales
for planete in range(nb_planete):
  
  fichier_source = liste_aei[planete]
  tableau = open(fichier_source, 'r')
  
  t.append([]) # time in year
  a.append([]) # the Semi-major axis in AU
  m.append([]) # mass in earth mass

  # On passe les 3 premieres lignes d'entete.
  for indice in range(3):
    tableau.readline()

  entete = tableau.readline()
  for ligne in tableau:
    # Pour chaque ligne du tableau, on decoupe suivant les 
    # espaces. (par defaut, s'il y a plusieurs espaces Ã  la 
    # suite, il ne va pas creer 50 000 colonnes, ce qui ne 
    # serait pas le cas avec un 'split(" ")')
    colonne = ligne.split()
    
    # On essaye de rajouter les elements. Si un seul d'entre eux
    # Ã  un soucis, on elimine toute la ligne (en gros, lors de 
    # l'ejection d'une planete, certains parametres peuvent 
    # devenir ****** au lieu d'un float, et ca genere une erreur
    # lors de la conversion en float)
    try:
      ti = float(colonne[0])
      ai = float(colonne[1])
      mi = float(colonne[7]) / 3.00374072e-6 # in earth mass
      
      # We must append inside the 'try' to avoid appending a value twice. The problem should occurs before anyway, 
      # when trying to convert into float
      t[-1].append(ti)
      a[-1].append(ai)
      m[-1].append(mi)
    except:
      pass

    
  tableau.close()

a1 = [ai[0] for ai in a]
a2 = [ai[-1] for ai in a]
a1.extend(a2)
a_max = max(a1) # We get the biggest Semi-major axis of the simulation (either at the beginning or the end of the simulation)
m_max = max([mi[-1] for mi in m]) # We get the biggest mass of the simulation
m_min = 0
a_min = 0
#~ a_max = 1.5 * CZ_LOCATION
delta_t = t[0][1] - t[0][0]

# If the timestep between two outputs is to big, we do not display a tail, because the planet will have time to do more than one 
# orbit between two values 
if (delta_t > 10.):
  isTail = False

# We get the array of reference time, i.e, one of the longuest list of time available in the list of planets. 
ref_len = 0
ref_id = 0
for planet in range(nb_planete):
  len_i = len(t[planet])
  if (len_i > ref_len):
    ref_len = len_i
    ref_id = planet
ref_time = t[ref_id]

# We get the index for the t_max value
if ('t_max' in locals()):
  id_max = int((t_max - ref_time[0]) / delta_t)
  t_max = ref_time[id_max]
else:
  id_max = ref_len - 1
  t_max = ref_time[-1]

# We get the index for the t_min value
if ('t_min' in locals()):
  id_min = int((t_min - ref_time[0]) / delta_t)
  t_min = ref_time[id_min]
else:
  id_min = 0
  t_min = ref_time[0]

if (t_max > 1e6):
  unit_time = "Million years"
  time_convertion = 1e6
elif (t_max > 1e3):
  unit_time = "Thousand years"
  time_convertion = 1e3
else:
  unit_time = "years"
  time_convertion = 1.
  


# on trace les plots
autiwa.lancer_commande("rm %s/%s*" % (OUTPUT_FOLDER, FRAME_PREFIX)) # We delete the previous frames

delta_t_min = (t_max - t_min) / (float(NB_FRAMES -1.))
# Number of timestep between each frame
# real number to be as close as possible from the real value, and do not encounter rounding problems. 
# The conversion to an integer is done at the very end.
ts_per_frame = delta_t_min / delta_t 

if (ts_per_frame < 1):
  ts_per_frame = 1
  NB_FRAMES = id_max - id_min +1

# We generate a list of colors
colors = [ '#'+li for li in autiwa.colorList(nb_planete)]

fig = pl.figure()
plot_orbits = fig.add_subplot(1, 1, 1)
plot = plot_orbits.plot
MAX_LENGTH = len(str(NB_FRAMES)) # The maximum number of characters needed to display

# We prepare the timeline
# Might work only if a_min and m_min are equal to 0
timeline_width = 0.7 # total width of the plot is "1"
timeline_height = 1.05 # total height of the plot is "1"
timetick_length = 0.02 # The semi-length of the extremal ticks of the timeline, in units of the total height of the plot

timeline_start = Line2D([0., 0.], [timeline_height - timetick_length, timeline_height + timetick_length], 
                        clip_on=False, color="#000000", linewidth=3, transform=plot_orbits.transAxes)
timeline_stop = Line2D([timeline_width, timeline_width], [timeline_height - timetick_length, timeline_height + timetick_length], 
                        clip_on=False, color="#000000", linewidth=3, transform=plot_orbits.transAxes)


t_frame = -1.
for frame_i in range(NB_FRAMES):
  id_time = id_min + int(frame_i * ts_per_frame)
  t_frame = t_min + int(frame_i * ts_per_frame) * delta_t
  
  percentage = (frame_i) / float(NB_FRAMES - 1)
  sys.stdout.write("%3.0f%% frame %*d : T = %#.2e years\r" % (percentage * 100., MAX_LENGTH, frame_i, t_frame))
  sys.stdout.flush()
  
  plot_orbits.clear()
  
  plot_orbits.fill([0, 0.004, 0.004, 0, 0], [0, 0, m_max, m_max, 0], color='yellow')

  if (isDisk and (len(contours_a) > 0)):
    plot_orbits.fill(contours_a[0], contours_m[0], facecolor="#ff0000", alpha=0.3, edgecolor='none', label="Outward migration")
    for (c_a, c_m) in zip(contours_a[1:], contours_m[1:]):
      plot_orbits.fill(c_a, c_m, facecolor="#ff0000", alpha=0.3, edgecolor='#000000')
  
  min_mass = 0.2 # earth mass
  markersize_prefactor = 4 / (min_mass**0.33)
  for planet in range(nb_planete):
    try:
      
      plot(a[planet][id_time], m[planet][id_time], 'o', color=colors[planet], 
           markersize=max(int(markersize_prefactor * (m[planet][id_time])**0.33),markersize_prefactor))
    except:
      #~ if (frame_i == 9):
        #~ pdb.set_trace()
      pass
      #~ # The planet has been ejected  

  plot_orbits.text(timeline_width, timeline_height, " %.1f %s" % (t_frame / time_convertion, unit_time), 
                   horizontalalignment='left', verticalalignment='center', 
                   size=15, transform=plot_orbits.transAxes)
  
  plot_orbits.add_line(timeline_start)
  plot_orbits.add_line(timeline_stop)
  
  timeline = Line2D([0., percentage * (timeline_width-0.01)], [timeline_height, timeline_height], 
                    marker=">", markevery=(1,1), color="#000000", linewidth=3, markersize=10, 
                    clip_on=False, transform=plot_orbits.transAxes)
  plot_orbits.add_line(timeline)
  
  plot_orbits.set_xlabel("Distance [AU]")
  plot_orbits.set_ylabel("Mass [Earths]")
  
  plot_orbits.axis('tight')
  plot_orbits.set_ylim(m_min, m_max)
  plot_orbits.set_xlim(a_min, a_max)
  plot_orbits.grid(True)

  nom_fichier_plot = "%s%0*d" % (FRAME_PREFIX, MAX_LENGTH, frame_i)  
  fig.savefig(os.path.join(OUTPUT_FOLDER, "%s.%s" % (nom_fichier_plot, OUTPUT_EXTENSION)), format=OUTPUT_EXTENSION)

sys.stdout.write("Movie Completed. Total number of frames : %d\n" % NB_FRAMES)
sys.stdout.flush()

if (NB_FRAMES<3):
  pl.show()

