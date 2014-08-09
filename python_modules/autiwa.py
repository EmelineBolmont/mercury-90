#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Module librement utilisable
"""module autiwa qui contient les fonctions générales que je trouve utile"""
from __future__ import print_function

__author__ = "Autiwa <autiwa@gmail.com>"
__date__ = "24 juin 2011"
__version__ = "$Revision: 2.1 $"
__credits__ = """Thanks to Bastien for his help when I was learning Python and his usefull tutorial"""



import sys
import os
import time
import commands
import subprocess
import math
import pdb


# TODO
# _ Modifier findRepExec() pour qu'elle renvoit vraiment le répertoire
#   d'exécution du tout premier script ou alors changer ma façon de m'en servir
# _ tester la méthode doc() de la classe AutiwaObject


SEPARATEUR_PARAMETRES = ';'

class AutiwaObject(object):
  """class that define general methods usefull for inheritance for other class. To sum up, this class doesn't do anything and is defined to be used for inheritance only"""
  
  def __init__(self):
    """ for the moment, do nothing. Nothing to initialize"""
  
  def doc(self, method="class"):
    """method to display the documentation of a class. 
    
    method="class" : if nothing is given, by default, it prints the documentation of the class. Else, the string must be the name of a method
    
    return : nothing. It prints the documentation of the class, by default, or of the method whose name is given in parameter
    """

    if (method == "class"):
      print("Documentation of the Class of this object :",type(self))
      print(self.__doc__)
    elif (method in dir(self)):  # test if a method with the name contained in the var 'method' exist within the current class
      print("Documentation of '"+method+"'")
      print(eval("self."+method+".__doc__"))
    else:
      raise Warning("There is no method named '"+method+"' in this class")
  
# Test unitaires
if (__name__ == '__main__'):
  test = AutiwaObject()
  test.doc()
  test.doc("__init__")
  try:
    test.doc("toto")
  except Warning:
    print("doc() renvoie bien une erreur quand test n'existe pas")
  except:
    raise Exception
  

class Temps(object):
  """Classe qui définit un temps. 
  Ceci permet d'additionner deux objets temps, les afficher à l'aide 
  de print, et ainsi de les manipuler en toute transparence
  
  *values : (an, jour, heure, minute, seconde) liste de valeurs, 
  au maxi 5, au mini 1. Si une seule valeur est donnée, elle sera prise
   comme étant le nombre de seconde. Si 2 sont données, la première 
   sera le nombre de minutes et la 2e le nombre de secondes, et ainsi 
   de suite au fur et à mesure que les nombres de valeurs augmentent.
  """
  
  NB_DIGITS = 2
  
  # Valeur en seconde de différentes durées
  MINUTE = 60
  HEURE = 60 * MINUTE
  JOUR = 24 * HEURE
  AN = 365.25 * JOUR
  
  def __init__(self, *values):
    nb_values = len(values)
    values = list(values) # On converti le tuple pour pouvoir ajouter des valeurs
    if (nb_values > 5):
      raise ValueError("The number of parameters must be of 5 maximum")
    for i in range(5-nb_values):
      values.insert(0,0)
    
    self.temps = float(values[0] * Temps.AN + values[1] * Temps.JOUR + values[2] * Temps.HEURE + values[3] * Temps.MINUTE + values[4])
    self.__update()
  
  def __update(self):
    """méthode privée qui met à jour toutes les variables internes
    
    met à jour self.nb_an, self.nb_jour, self.nb_heure, self.nb_minute et self.nb_seconde
    """
    reste = self.temps
  
    self.nb_an = int(reste / Temps.AN)
    reste -= self.nb_an * Temps.AN
    
    self.nb_jour = int(reste / Temps.JOUR)
    reste -= self.nb_jour * Temps.JOUR
    
    self.nb_heure = int(reste / Temps.HEURE)
    reste = reste - self.nb_heure * Temps.HEURE
    
    self.nb_minute = int(reste / Temps.MINUTE)
    reste = reste - self.nb_minute * Temps.MINUTE
    
    self.nb_seconde = round(reste, Temps.NB_DIGITS)
  
  def __str__(self):
    """Permet, via print, d'afficher le temps
    """
    # On prépare la chaîne de caractère
    str_temps = ''
    if (self.nb_an == 1):
      str_temps += str(self.nb_an)+"an "
    elif (self.nb_an > 1):
      str_temps += str(self.nb_an)+"ans "
      
    if (self.nb_jour == 1):
      str_temps += str(self.nb_jour)+"jour "
    elif (self.nb_jour > 1):
      str_temps += str(self.nb_jour)+"jours "
      
    if (self.nb_heure != 0):
      str_temps += str(self.nb_heure)+"h "
      
    if (self.nb_minute != 0):
      str_temps += str(self.nb_minute)+"min "
      
    str_temps += str(self.nb_seconde)+"s"
    return str_temps
  
  def __add__(self, other):
    """Surcharge de +"""
    return Temps(self.temps + other.temps)
  
  def __sub__(self, other):
    """Surcharge de -"""
    return Temps(self.temps - other.temps)
  
  def __eq__(self, other):
    """surcharge de =="""
    if (self.temps == other.temps):
      return True
    else:
      return False
  
  def __ne__(self, other):
    """surcharge de !="""
    return not(self.__eq__(other))

if (__name__ == '__main__'):
  t1 = Temps(100000)
  print(t1)
  t2 = Temps(50000)
  print(t2)
  
  # test de __eq__
  assert not(t1 == t2)
  
  # test de __ne__
  assert t1 != t2
  
  #test of __add__ method
  t3 = Temps(150000)
  assert t3 == t1 + t2
  
  #test of __sub__ method
  assert t2 == t1 - t2


######################
# définition des fonctions
######################
def find_rep_exec():
  """fonction qui renvoie le chemin absolu où se trouve le script, 
  c'est à dire le répertoire d'exécution. Si on appelle des scripts 
  dans le script principal, le dossier renvoyé sera celui du dernier 
  script lancé. Je sais pas comment faire pour afficher tout le temps 
  le répertoire du script principal.
  """
  # Version 1.0
  temp = sys.argv[0]
  path = os.path.dirname(temp)
  full_path = os.path.abspath(path)
  return full_path

def write_log(nom_fichier_log, texte_log, rep_exec=find_rep_exec()):
  """se place dans le répertoire d'exécution du script, écrit une ligne
  dans le fichier log, puis retourne dans le répertoire où on était
  avant l'exécution du script. 
  
  De plus, affiche à l'écran la ligne de texte que l'on écrit dans le fichier log."""
  # Version 1.21
  
  print(texte_log)
  
  rep_initial = os.getcwd()
  os.chdir(rep_exec)
  log = open(nom_fichier_log,'a')
  log.write("["+str(time.strftime('%d/%m/%Y %H:%M:%S'))+"] "+texte_log+"\n")
  log.close()
  os.chdir(rep_initial)

def fusionner_fichier(fusion, liste):
  """On passe en argument le nom du fichier résultant, et la liste des
  fichiers que l'on souhaite fusionner (des textes donc)"""
  fichier_fusion = open(fusion, 'w')
  for file in liste:
    fichier = open(file, 'r')
    for ligne in fichier:
      fichier_fusion.write(ligne)
    fichier.close()
  fichier_fusion.close()

def copier_fichier(liste, rep_arrive):
  """On déplace les fichiers de sortie dans le répertoire prévu
   à cet effet, afin d'avoir une trace."""
  for nom in liste:
    os.system("cp "+nom+" "+os.path.join(rep_arrive,nom))

def mise_en_forme_temps(temps):
  """Prend en argument un temps en seconde et renvoie une chaîne de 
  caractères qui mettra en forme ce temps pour le rendre plus lisible 
  si celui-ci dépasse plusieurs heures"""

  temps = float(temps)
  
  nb_heure = int(temps / 3600)
  temps = temps - nb_heure * 3600
  
  nb_minute = int(temps / 60)
  temps = temps - nb_minute * 60
  
  nb_seconde = round(temps, 2)
  
  # On prépare la chaîne de caractère
  str_temps = ''
  if (nb_heure!=0):
    str_temps += str(nb_heure)+"h "
  if (nb_minute!=0):
    str_temps += str(nb_minute)+"min "
  str_temps += str(nb_seconde)+"s"
  return str_temps

def nettoyage(fichiers_effacer):
  """Essaye d'effacer la liste de fichier passée en arguments"""
  for nom in fichiers_effacer:
    if os.path.exists(nom):
      os.remove(nom)

def enumeration(liste):
  """liste, un élément par ligne, les éléments de l'argument 
  avec les indices associés"""
  for (indice,element) in enumerate(liste):
    print(indice, ":",element)

def lancer_commande(commande):
  """lance une commande qui sera typiquement soit une liste, soit une 
  commande seule. La fonction renvoit un tuple avec la sortie, 
  l'erreur et le code de retour"""
  if (type(commande)==list):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  elif (type(commande)==str):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  else:
    raise TypeError("The command is neither a string nor a list.")
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  # there is .poll() or .wait() but I don't remember the difference. For some kind of things, one of the two was not working
  return (process_stdout, process_stderr, returncode)
    
def suppr_dossier(liste,liste_suppr):
  """Supprime de la 1e liste les éléments de la 2e liste"""
  for dir in liste_suppr:
    try:
      liste.remove(dir)
    except ValueError:
      pass

def nom_fichier(chemin):
  """À partir du chemin absolu ou relatif d'un fichier, récupère le 
  nom du fichier sans l'extension. Pour l'instant, version simplifiée 
  qui suppose que le nom est de la forme './nom.py'"""  
  
  temp = chemin
  temp = temp.rstrip(".py")
  temp = temp.lstrip("./")
  
  return temp
  
def number_fill(number, fill):
  number = str(number)
  # tant qu'on n'a pas la bonne taille, on rajoute des "0" à gauche.
  while (len(number) < fill):
    number = "0"+number
  
  return number

def ecriture_variable(nom, entete, *variables):
  """ 'nom' est le nom du fichier dans lequel seront écrites les 
  variables. L'extension sera '.dat' et 'nom' ne doit donc pas avoir 
  d'extension. 'entete' sera l'entete du fichier, permettant notamment
   de donner les noms des différentes variables que l'on s'apprete à 
   stocker. 'variables' contiendra la liste de toutes les variables 
   qui nous intéressent. Le nombre de variable n'est pas déterminé à l'avance."""
  if not(type(nom)==str):
    raise TypeError("Le nom de la variable n'est pas une chaine de caracteres")
  
  extension = "dat"
  
  fichier = open(nom+"."+extension, 'w')  
  # on écrit l'entete même si elle est vide pour ne pas avoir de problèmes lors de la lecture
  fichier.write(entete+"\n") 

  nb_var = len(variables)
  nb_ligne = len(variables[0])
  
  # On boucle d'abord sur l'indice représentant chaque élément
  for ligne in range(nb_ligne):
    temp = []
    
    ## On boucle sur les différentes variables de la liste et on 
    # récupère les éléments correspondant à la ligne courante.
    for var in range(nb_var):
      temp.append(str(variables[var][ligne]))
    # On fusionne les éléments en les séparant par le séparateur déterminé à l'avance
    temp = SEPARATEUR_PARAMETRES.join(temp)
    
    # On écrit tout ça dans le fichier.
    fichier.write(temp+"\n")
  fichier.close()
  
def lecture_variable(nom, colonnes=''):
  """On va lire le fichier s'appellant 'nom' (l'extension étant donnée
  dans le corps de la fonction) et renverra les colonnes que l'on veut.
   Ces colonnes sont données via la variable 'colonnes'. 
   Si celle ci est un entier, alors cela signifie simplement que l'on 
   veut ce nombre de colonnes (et ce même si le fichier en contient 
   plus.). Si rien n'est donné, alors on détermine le nombre de 
   colonnes du fichier et on renvoit la totalité des données"""
  if not(type(nom)==str):
    raise TypeError("Le premier argument, nom, doit etre une chaine de caracteres")
  
  extension = 'dat'
  nom_fichier = nom+"."+extension
  
  # si c'est un entier, c'est simplement le nombre de colonnes
  if type(colonnes)==int:
    liste_colonnes = range(colonnes)
  # Si c'est une liste ou un tuple, ce sont les indices des colonnes qu'on veut.
  elif type(colonnes) in (list,tuple):
    liste_colonnes = colonnes
  # Si rien ne correspond, on détermine le nombre de colonnes et on lit tout.
  else:
    fichier = open(nom_fichier, 'r')
    entete = fichier.readline()
    line = fichier.readline()
    fichier.close()
    line = line.split(SEPARATEUR_PARAMETRES)
    liste_colonnes = range(len(line))
  fichier = open(nom_fichier, 'r')
  
  entete = fichier.readline()
  
  tableau = [[] for i in liste_colonnes]
  for line in fichier:
    line = line.split(SEPARATEUR_PARAMETRES)
    ## On essaye de lire et stocker la ligne. 
    # Si ça marche pas, on passe à la suivante.
    try:
      # Pour chaque indice voulu, on stocke la valeur.
      for (indice, colonne) in enumerate(liste_colonnes):
        tableau[indice].append(float(line[colonne]))
    except:
      pass
  fichier.close()
  
  ## S'il n'y a qu'une seule variable, alors on ne renvoit pas un 
  # tuple mais la variable elle même.
  if (len(tableau)==1):
    sortie = tableau[0]
  else:
    sortie = tuple(tableau)
  
  return sortie
  
def printCR(text):
  """function that allow us to display on the screen some text in order to rewrite it afterwards"""
  
  sys.stdout.write(text+chr(13))
  sys.stdout.flush()

#http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
def rotation_matrix(angle, direction):
  """
  Return the rotation matrix for the given angle (in radian) and vector
  angle : the angle of the rotation, this angle is algebric
  direction : the vector axis of the rotation
  
  return : return the matrix given by the angle and the axis rotation
  """
  import numpy
  
  s = math.sin(angle)
  c = math.cos(angle)
  
  ## Transformation du vecteur directeur en array numpy (à la base, ça 
  # peut être un tuple, un array ou une liste
  direction = numpy.array(direction)
  
  #normalisation du vecteur
  (ux, uy, uz) = direction / math.sqrt(sum(direction**2))
  
  mat = numpy.array([[ux**2 + (1 - ux**2) * c, ux * uy * (1 - c) - uz * s, ux * uz * (1 - c) + uy * s], 
                    [ux * uy * (1 - c) + uz * s, uy**2 + (1 - uy**2) * c, uy * uz * (1 - c) - ux * s],
                    [ux * uz * (1 - c) - uy * s, uy * uz * (1 - c) + ux * s, uz**2 + (1 - uz**2) * c]], dtype=numpy.float64)
  
  return mat

def global_rotation(*matrices):
  """effectue le produit des matrices données en argument, dans l'ordre
  donné. Si on veut faire plusieurs matrices de rotation dans l'ordre 
  1, 2 puis 3, il faudra utiliser 
  mat = global_rotation(mat3, mat2, mat1). 
  Marche pour des matrices 3x3"""
  import numpy
  
  temp = numpy.identity(3)
  for i in matrices:
    temp = numpy.dot(temp, i)
  return temp

def affiche_cadre_texte(text, char="*", size=40, shell_size=80):
  """
  function that print a text surrounded by a character that make a frame arround it.
  text : the text that will be printed in the middle of the frame
  char = "*" : the character that will form the frame
  size = 40 : total size of the frame. If the text is longer than that, it will be printed in several lines.
  shell_size = 80 : size of the shell. By default 80. Allow the function to center the frame in comparison to the shell.
  
  return : return nothing, only 0 if everything is ok, but print the text in a console
  """
  from math import ceil
  
  if (size > shell_size):
    size = shell_size
  
  horiz_edge = char * size
  
  # each edge has a width of 1 character
  max_size = size - 2
  # at least, we want to have one space between the text and each vertical edge.
  line_size = max_size - 2
  
  text_size = len(text)
  
  # we slice the text in enough lines to fit the frame width.
  # ceil give the superior integer of the argument that must be a floating number. So we mutiply the argument by 1. to have a real. After this, we convert the given superior integer into an integer because ceil return a floating number. 
  nb_lignes = int(ceil(1. * text_size/line_size))
  lines = []
  for line in range(nb_lignes-1):
    lines.append(text[line * line_size:(line + 1) * line_size])
  # On rajoute la dernière ligne, qui probablement ne sera pas entière, après la boucle
  lines.append(text[(nb_lignes - 1) * line_size:])  
  
  lines_centered = []
  for line in lines:
    lines_centered.append(line.center(max_size))
  
  # We define the frame :
  frame_lines = []
  frame_lines.append(horiz_edge)
  for line in lines_centered:
    frame_lines.append(char+line+char)
  frame_lines.append(horiz_edge)
  
  # We display the centered frame
  for line in frame_lines:
    print(line.center(shell_size))
  
  return 0
  
def significativeRound(real,roundNumber):
  """function that return a truncative floating number of the number 
  in parameter, given the roundNumber we want
  Version 2.0
  
  Parameters
  real : The floating point number we want to truncate
  roundNumber : the number of significative figures we want
  """
  #print "Warning:pour test, ne fait rien!!"
  #return real
  import math
  
  if (roundNumber == 0):
    print("'roundNumber' is set to 0")
    exit
  elif (type(roundNumber) != int):
    print("'roundNumber' is not an integer")
    exit
  
  # If an integer, we'll have problem for divisions.
  if (type(real) != float):
    real = float(real)
  
  if (real == 0.):
    return 0.
  
  # in case we have negative number
  if (real < 0.):
    negative = True
    real = - real
  else:
    negative = False
  
  significativeNumbers = []
  
  # First we initialise various parameters for the loop. The principle
  # is to extract the integer part of the logarithm. This way, we 
  # retrieve the first significative number. By substracting this 
  # number to the original number, we can successively retrieve all 
  # the significative figures.
  
  # We first get the order of magnitude of the number to be able to 
  # retrieve all the significative numbers in the loop after that.
  logNumber = int(math.log10(real))

  significativeNumbers = [0] # we force the 'while' loop to do at least one turn
  # For numbers under 1, we search for the first significative number
  while (significativeNumbers[-1] == 0):
    ordMagn = 10**logNumber
    significativeNumbers = [int(real/ordMagn)]
    decimalPart = real/ordMagn - significativeNumbers[-1]
    logNumber = logNumber - 1 # if this is the last turn, this line will do nothing. This avoid to substract one in the first turn (if we had done it at the beginning)

  if (significativeNumbers[-1] == 0):
    ordMagn = ordMagn / 10
    try:
      significativeNumbers = [int(real/ordMagn)]
    except:
      print(real)
      pdb.set_trace()
    decimalPart = real/ordMagn - significativeNumbers[-1]
    
  for i in range(roundNumber-2):
    significativeNumbers.append(int(decimalPart * 10))
    decimalPart = decimalPart * 10 - significativeNumbers[-1]
  
  # For the last number we must round.
  significativeNumbers.append(int(round(decimalPart * 10 + 0.5)))
    
  numberWithoutMagn = 0.
  for (i, si) in enumerate(significativeNumbers):
    numberWithoutMagn += si * 10**(-i)
  
  truncatedNumber = ordMagn * numberWithoutMagn
  
  if not(negative):
    return truncatedNumber
  else:
    return -truncatedNumber

def colorList(nb_colors, exclude=['ffffff']):
  """Function that return a list of colors given the number of different colors we want. 
  It is still a simple version of the function that may contain bugs, especially for large values of colors. 
  The color 'white' is prohibited here, but not 'black'
  
  Parameter
  nb_colors : The number of different colors we want.
  exclude=['ffffff'] : The colors we do not want in our list
  
  Return 
  A list of colors, with the desired number of elements
  """
  from math import ceil
  from random import shuffle
  
  excludeColors = list(exclude)
  
  # We add the length of the excludedColors array. Maybe some of them will not be generated in our case with the individual 
  # colors taken into account, but we want, at least, to generate enough colors. 
  individual_colors = int(ceil((nb_colors + len(exclude))**0.33)) 
  
  colors = [0, 255, 127]
  iteration = 0
  # We search for all individual values (for R, G, B) that we need to define at least the number of colors we want.
  while (len(colors) < individual_colors):
    step = 255 / (2**(iteration + 1))
    
    iteration += 1

    #For the current iteration, the number of subdivision will be a power of 2
    nb_el = 2**iteration

    # We start at half the step
    temp = step/2
    for i in range(nb_el):
      colors.append(temp)
      temp += step
  
  # If the total number of colors is low, no need to have 127 that is in fact only needed for the algorithm.
  if (nb_colors <=6):
    colors.remove(127)
  
  
  
  HEXcolors = []
  nb = 0
  
  if (nb_colors == 3):
    for colorTemp in ['ff0000', '00ff00', '0000ff']:
      if (colorTemp not in excludeColors):
        HEXcolors.append(colorTemp)
        nb += 1
    if (nb == 3):
      return HEXcolors
  
  
  HEXcolors = []
  nb = 0
  
  for i in range(len(colors)):
    subcolors = colors[0:i+1]
    for R in subcolors:
      for G in subcolors:
        for B in subcolors:
          colorTemp = hexColor(R)+hexColor(G)+hexColor(B)
          if (colorTemp not in excludeColors):
            HEXcolors.append(colorTemp)
            excludeColors.append(colorTemp)
            
            nb += 1
          
          if (nb == nb_colors):
            return HEXcolors
  
  pdb.set_trace()

def hexColor(integer):
  """given a number between 0 and 255, return the hexadecimal value as a two character string
  """
  
  if (type(integer) != int):
    raise TypeError("The parameter must be an integer")
  
  if (integer <0 or integer > 255):
    raise ValueError("The parameter must be between 0 and 255")

  value = hex(integer)
  value = value.split("0x")[1]

  if (len(value) < 2):
    value = "0"+value

  return value

def contrastColor(ref_color):
  """
  Return '000000' by default, but if 'color' is quite dark, then return 'ffffff'

  The principle is simple, we calculate the brightness of the color and regarding a tolerance parameter, 
  we put white or black in order to have a visible color on the given one.

  Example :
  contrastColor('ffffff') = '000000'
  contrastColor('000000') = 'ffffff'
  """
  
  tolerance = 130

  R = int('0x'+ref_color[0:2],16)
  G = int('0x'+ref_color[2:4],16)
  B = int('0x'+ref_color[4:6],16)
  
  # We calculate the square of the brightness in order to avoid the use of sqrt.
  brightness_square = .241 * R**2 + .691 * G**2 + 0.068 * B**2 # source : http://alienryderflex.com/hsp.html
  

  if (brightness_square > tolerance**2):
    color = "000000"
  else:
    color = "ffffff"

  return color

def get_subplot_shape(number_of_plots):
  """We chose the number of plot in x and y axis for the p.multi
  environment in order to plot ALL the resonant angles and have x and
  y numbers as close on from another as possible.

  HOW TO :
  fig = pl.figure(1)
  subplot_i = 0
  (nb_line, nb_row) = get_subplot_shape(5)
  subplot_i += 1
  plot_dof = fig.add_subplot(nb_line, nb_row, subplot_i)"""
  nb_plots_x = 1
  nb_plots_y = 1
  while (nb_plots_x * nb_plots_y < number_of_plots):
    if (nb_plots_x == nb_plots_y):
      nb_plots_y += 1
    else:
      nb_plots_x += 1

  return (nb_plots_x, nb_plots_y)

# Si le module est lancé directement, on teste certaines méthodes.
if __name__=='__main__':



  # Test de afficheCadreTexte
  text = (4) * '012345678901234567890123456789123456'
  affiche_cadre_texte(text)
  affiche_cadre_texte(text+'rab')
  affiche_cadre_texte("Partie Principale")

  # Test de significativeRound
  liste = [1.23875486757865, 1.2487598986e18, 1.8975896457e-30, -1.76458764, 0.00003487586, 0.983297103893]
  rounding=4
  print("numbers rounded with ",rounding," significative numbers : ")
  for ele in liste:
    print(ele," rounded give : ",significativeRound(ele,4))
