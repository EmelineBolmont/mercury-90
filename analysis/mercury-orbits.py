#!/usr/bin/env python
# -*- coding: utf-8 -*-
# v1.0
# To display the orbits of the planets in the system in the (x,y) plane

from math import *
import pylab as pl
import os, pdb, autiwa
import numpy as np
import sys # to be able to retrieve arguments of the script
import mercury

BINARY_FOLDER = '$HOME/bin/mercury'
OUTPUT_EXTENSION = 'png'

NB_POINTS = 50 # Number of points for the display of circles
NB_FRAMES = 2
#~ CZ_LOCATION = 6. # POSITION OF THE CONVERGENCE ZONE IN ai

isTail = True # There is a part where this boolean is changed automatically if the timestep between two output is to huge.

###############################################
## Beginning of the program
###############################################


# We get arguments from the script

#~ pdb.set_trace()
isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * t_max (the end of the output, in years)" + "\n" + \
" * t_min (the beginning of the output (in years)" + "\n" + \
" * forceCircle True or False, to display circle instead of real orbits if the output interval is huge" + "\n" + \
" * zoom (the farthest location in the disk that will be displayed (in AU)" + "\n" + \
" * frames=1 (the number of frames you want)" + "\n" + \
" * ext=png (The extension for the output files)"

for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 't_min'):
    t_min = float(value)
  elif (key == 't_max'):
    t_max = float(value)
  elif (key == 'forceCircle'):
    isTail = True
  elif (key == 'frames'):
    NB_FRAMES = int(value)
  elif (key == 'ext'):
    OUTPUT_EXTENSION = value
  elif (key == 'zoom'):
    plot_range = float(value)
  elif (key == 'help'):
    print(problem_message)
    exit()
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)

if (NB_FRAMES <= 1):
  print("The number of frames cannot be lower than 2")
  NB_FRAMES = 2
  
####################
# We erase old output files, add 'x, y and z' to the outputs values, and generate new output files
####################

autiwa.lancer_commande("rm *.aei element.out")

elementin = mercury.Element()
elementin.read()
elementin.set_format_sortie(" a21e e21e i8.4 g8.4 n8.4 l8.4 m21e x21e y21e z21e")
elementin.write()

(stdout, stderr, returnCode) = autiwa.lancer_commande(os.path.join(BINARY_FOLDER, "element"))
if (returnCode != 0):
  print("Unable to Launch 'element'")

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
x = [] # x cartesian coordinate in AU
y = [] # y cartesian coordinate in AU
z = [] # z cartesian coordinate in AU


# On recupere les donnees orbitales
for planete in range(nb_planete):
  
  fichier_source = liste_aei[planete]
  tableau = open(fichier_source, 'r')
  
  t.append([]) # time in year
  a.append([]) # the semi major axis in AU
  m.append([]) # mass in earth mass
  x.append([]) # demi-grand axe en ua
  y.append([]) #
  z.append([]) #

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
      mi = float(colonne[7]) / 3.00374072e-6 # in solar mass
      xi = float(colonne[8])
      yi = float(colonne[9])
      zi = float(colonne[10])
      
      # We must append inside the 'try' to avoid appending a value twice. The problem should occurs before anyway, when trying to convert into float
      t[-1].append(ti)
      a[-1].append(ai)
      m[-1].append(mi)
      x[-1].append(xi)
      y[-1].append(yi)
      z[-1].append(zi)
    except:
      pass
    

    
  tableau.close()

#~ a_max = max([ai[0] for ai in a]) # We get the biggest semi major axis of the simulation initially
#~ a_max = 1.5 * CZ_LOCATION
# The separation between outputs is not always the same because the real output interval in mercury and element might be slightly different.
delta_t = (t[0][-1] - t[0][0]) / float(len(t[0]))

# If the timestep between two outputs is to big, we do not display a tail, because the planet will have time to do more than one 
# orbit between two values 
if (delta_t > 10.):
  print("/!\ time between output will have the effect to display ugly orbits. \
  Try the option 'isTail=False', but orbits will be circle and not representative of reality")

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

# We get the index for the t_max value
if ('t_min' in locals()):
  id_min = int((t_min - ref_time[0]) / delta_t)
  t_min = ref_time[id_min]
else:
  id_min = 0
  t_min = ref_time[0]




#~ idx_tail = []
#~ for planet in range(nb_planete):
  #~ tmp = a[planet][id_max]**1.5 # the period of the planet in years
  #~ tmp = int(id_max - tmp / delta_t + 2)
  #~ # Negative index is not possible. So if the planet did not have the time to do one orbit, the tail will be artitificially put to 0
  #~ if (tmp < 0):
    #~ tmp = 0
  #~ idx_tail.append(tmp) 

# on trace les plots

delta_t_min = (t_max - t_min) / (float(NB_FRAMES -1.))
# Number of timestep between each frame
# real number to be as close as possible from the real value, and do not encounter rounding problems. 
# The conversion to an integer is done at the very end.
ts_per_frame = delta_t_min / delta_t 

# If there is too many frames for the outputs availables, we impose 1 output between each frames and reduce the total number of frames
if (ts_per_frame < 1):
  ts_per_frame = 1
  NB_FRAMES = id_max - id_min +1


# We generate a list of colors
tmp = autiwa.colorList(nb_planete)
colors = [ '#'+li for li in autiwa.colorList(nb_planete)]

angles = [2 * pi / NB_POINTS * i for i in range(NB_POINTS)]
angles.append(angles[0]) # we want to have a full circle, perfectly closed
angles = np.array(angles)



for frame_i in range(NB_FRAMES):
  id_time = id_min + int(frame_i * ts_per_frame)
  t_frame = t_min + int(frame_i * ts_per_frame) * delta_t
  
  if (frame_i == NB_FRAMES - 1):
    id_time = id_max
    t_frame = t_max
  
  print("frame %d : T = %#.2e years" % (frame_i, t_frame))
  
  pl.figure(1)
  pl.clf()
  # We put a yellow star to display the central body
  pl.plot(0, 0, '*', color='yellow', markersize=20) 

  # We display a circle for the convergence zone
  #~ CZ_x = CZ_LOCATION * np.cos(angles)
  #~ CZ_y = CZ_LOCATION * np.sin(angles)
  #~ pl.plot(CZ_x, CZ_y, '-.', color="#000000")
  
  if isTail:
    idx_tail = [None] * nb_planete
    for planet in range(nb_planete):
      # If the planet is still in the system at that time, we display it, else, the 'except' make us pass to the next planet.
      try:
        tmp = a[planet][id_time]**1.5 # the period of the planet in years
        tmp = int(id_time - tmp / delta_t + 2)
        # Negative index is not possible. So if the planet did not have the time to do one orbit, the tail will be artitificially put to 0
        if (tmp < 0):
          tmp = 0
        idx_tail[planet] = tmp
        pl.plot(x[planet][idx_tail[planet]:id_time+1], y[planet][idx_tail[planet]:id_time+1], color=colors[planet], label='PLANETE'+str(planet))
        pl.plot(x[planet][id_time], y[planet][id_time], 'o', color=colors[planet], markersize=int(5* (m[planet][id_time])**0.33))
      except:
        pass
        # The planet has been ejected
  else:
    # We draw circles for each orbit. This part might be used if a delta_t is more than one orbit of a planet.
    for planet in range(nb_planete):
      try:
        r = sqrt(x[planet][id_time]**2 + y[planet][id_time]**2)
        x_circ = r * np.cos(angles)
        y_circ = r * np.sin(angles)
        pl.plot(x_circ, y_circ, color=colors[planet], label='PLANETE'+str(planet))
        pl.plot(x[planet][id_time], y[planet][id_time], 'o', color=colors[planet], markersize=int(5* (m[planet][id_time])**0.33))
      except:
        pass
        #~ # The planet has been ejected
  pl.title("T = %#.2e years" % t_frame)
  pl.xlabel("x (in AU)")
  pl.ylabel("y (in AU)")
  
  pl.axis('equal')
  if ('plot_range' in vars()):
    pl.ylim(-plot_range, plot_range)
    pl.xlim(-plot_range, plot_range)
  #~ pl.legend()
  pl.grid(True)
  #~ pdb.set_trace()
  nom_fichier_plot = "frame_"+autiwa.number_fill(frame_i,len(str(NB_FRAMES)))
  
  pl.savefig(nom_fichier_plot+'.'+OUTPUT_EXTENSION, format=OUTPUT_EXTENSION)

#~ dossier_output = "output"
#~ system("mkdir dossier_output")
#~ system("cd dossier_output")

pl.show()

