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

BINARY_FOLDER = '$HOME/bin/mercury'
OUTPUT_EXTENSION = "png"

NB_POINTS = 50 # Number of points for the display or circles
NB_FRAMES = 2
CZ = None # POSITION OF THE CONVERGENCE ZONE IN AU

###############################################
## Beginning of the program
###############################################


# We get arguments from the script

#~ pdb.set_trace()
isProblem = False
problem_message = "AIM : Display in a m = f(a) diagram, all the planets of the current mercury simulation" + "\n" + \
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * t_max (the end of the output, in years)" + "\n" + \
" * t_min (the beginning of the output (in years)" + "\n" + \
" * frames=1 (the number of frames you want)" + "\n" + \
" * cz=1 (The position of a convergence zone in AU)" + "\n" + \
"   cz=[[1,60],[4,30]] (the list of mass (earth mass) and zero torque position (in AU) successively)" + "\n" + \
" * ext=png (The extension for the output files)" + "\n" + \
" * help : display this current message"

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 't_min'):
    t_min = float(value)
  elif (key == 't_max'):
    t_max = float(value)
  elif (key == 'frames'):
    NB_FRAMES = int(value)
  elif (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'cz'):
    CZ = eval(value)
  elif (key == 'help'):
    print(problem_message)
    exit()
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)

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
  a.append([]) # the semi major axis in AU
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
a_max = max(a1) # We get the biggest semi major axis of the simulation (either at the beginning or the end of the simulation)
m_max = max([mi[-1] for mi in m]) # We get the biggest mass of the simulation
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



# on trace les plots

delta_t_min = (t_max - t_min) / (float(NB_FRAMES -1.))
# Number of timestep between each frame
# real number to be as close as possible from the real value, and do not encounter rounding problems. 
# The conversion to an integer is done at the very end.
ts_per_frame = delta_t_min / delta_t 

if (ts_per_frame < 1):
  ts_per_frame = 1
  NB_FRAMES = id_max - id_min +1



# We generate a list of colors
tmp = autiwa.colorList(nb_planete)
colors = [ '#'+li for li in autiwa.colorList(nb_planete)]

angles = [2 * pi / NB_POINTS * i for i in range(NB_POINTS)]
angles.append(angles[0]) # we want to have a full circle, perfectly closed
angles = np.array(angles)


t_frame = -1.
for frame_i in range(NB_FRAMES):
  id_time = id_min + int(frame_i * ts_per_frame)
  t_frame = t_min + int(frame_i * ts_per_frame) * delta_t
  print(t_frame)
  
  pl.figure(1)
  pl.clf()
  # We put a yellow star to display the central body
  #~ pl.plot(0, 0, '*', color='yellow', markersize=20) 
  pl.fill([0, 0.004, 0.004, 0, 0], [0, 0, m_max, m_max, 0], color='yellow')

  # We display the convergence zone
  if (type(CZ) == list):
    pl.plot(CZ[1], CZ[0], '-.', color="#000000")
  elif (type(CZ) in [float, int]):
    pl.plot([CZ, CZ], [0, m_max], '-.', color="#000000")
  
  # We draw circles for each orbit. This part might be used if a delta_t is more than one orbit of a planet.
  for planet in range(nb_planete):
    try:
      pl.plot(a[planet][id_time], m[planet][id_time], 'o', color=colors[planet], markersize=int(5* (m[planet][id_time])**0.33))
    except:
      #~ if (frame_i == 9):
        #~ pdb.set_trace()
      pass
      #~ # The planet has been ejected
  #~ pl.figtext(0.40, 0.95, "T = "+str(t_frame)+" years", horizontalalignment='left', verticalalignment="top")
  pl.title("T = %#.2e years" % t_frame)
  pl.xlabel("a [AU]")
  pl.ylabel("mass [Earths]")
  
  pl.axis('tight')
  pl.ylim(0, m_max)
  pl.xlim(0, a_max)
  #~ pl.legend()
  pl.grid(True)
  #~ pdb.set_trace()
  nom_fichier_plot = "frame_growth_"+autiwa.number_fill(frame_i,len(str(NB_FRAMES)))
  
  pl.savefig(nom_fichier_plot+'.'+OUTPUT_EXTENSION, format=OUTPUT_EXTENSION)

#~ dossier_output = "output"
#~ system("mkdir dossier_output")
#~ system("cd dossier_output")


pl.show()

