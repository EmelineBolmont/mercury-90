#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Librairies contenant des constantes utiles
# pour importer les variables, faire "from constants import *", 
# sinon, les variables s'appeleront constants.nom_variable
import math

__author__ = "Autiwa <autiwa@gmail.com>"
__date__ = "24 août 2011"
__version__ = "$Revision: 1.0 $"
__credits__ = """Script that define several constants"""

"""The module is not retro-compatible with previous version where the constants were defined in small characters instead of CAPSLOCK"""

######################
# Définition des constantes
######################
MS = 1.9891e30 # kg masse du soleil
MJ = 1.8986e27 # kg masse de jupiter
MT = 5.9736e24 # kg masse de la terre
MSAT = 5.6846e26 # kg masse de saturne

DJ = 1.326 # g/cm³ densité moyenne d'une planète géante (ici, densité moyenne de jupiter)
DT = 5.515 # g/cm³ densité moyenne d'une planète de type terrestre (ici, densité moyenne de la terre)

# en mètre valeur d'une unité astronomique
AU = 1.495979e11

# nombre de jours dans un an, c'est plus simple ensuite pour calculer T
YEAR = 365.25

# nombre de secondes dans une journée
DAY = 86400

# Constante de gravitation universelle en unité SI
G = 6.6726e-11

# Constante de la gravitation avec distance en ua, 
# temp en jours, et masse en masse solaire
G0 = 2.959122082855911e-4 # G * jour**2 * ms / ua**3 

######################
# Facteurs de conversion
######################
DEGTORAD = math.pi / 180.
RADTODEG = 180. / math.pi
