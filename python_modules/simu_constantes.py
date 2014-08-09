#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Librairies contenant des constantes utiles
# Version 0.3
# pour importer les variables, faire "from simu_constantes import *", 
# sinon, les variables s'appeleront simu_constantes.nom_variable
import math

######################
# Définition des constantes
######################
MS = 1.9891e30# kg masse du soleil
MJ = 1.8986e27# kg masse de jupiter
MT = 5.9736e24# kg masse de la terre
MSAT = 5.6846e26# kg masse de saturne

DJ = 1.326 #g/cm³ densité moyenne d'une planète géante (ici, densité moyenne de jupiter)
DT = 5.515 #g/cm³ densité moyenne d'une planète de type terrestre (ici, densité moyenne de la terre)

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
G0 = 2.959122082855911e-4#G * jour**2 * ms / ua**3 

######################
# Facteurs de conversion
######################
DEGTORAD = math.pi/180
RADTODEG = 180/math.pi

#######################
# On prépare le stockage des fichiers
#######################
# liste des fichiers à supprimer avant de lancer mercury. 
# !!! Par défaut, ces valeurs ne sont pas utilisées car elles sont 
# définies dans les fonctions de déplacement et de suppression.
fichier_sortie = ['xv.out', 'ce.out', 'info.out', "mercury_stdout.txt", "mercury_stderr.txt"]
fichier_dmp = ['big.dmp', 'small.dmp', 'param.dmp', 'restart.dmp']
fichier_tmp = ['big.tmp', 'small.tmp', 'param.tmp', 'restart.tmp']
fichier_entree = ['param.in', 'big.in', 'element.in']
fichiers_effacer = fichier_dmp + fichier_tmp + fichier_sortie
fichiers_deplacer = fichier_dmp + fichier_tmp + fichier_sortie + fichier_entree

#######################
# Variables pour la syntaxe des fichiers .simu
####################### 
prefixe_commentaire = "#"
prefixe_nom_simulation = "&"
# à droite se trouve deux membres séparés par un égal. le membre de 
# gauche doit avoir un nom rigoureux que l'on peut tester ensuite pour 
# stocker le membre de droite dans un nom de variable. (pour l'instant ne fonctionne pas)
prefixe_variable = "@"
__liste_variable__ = ['temps_integration', 'nb_planete_total']
operateur_assignation = "="

separateur_parametres = ";"

fichier_global_simu = "fin_simu.simu"

#######################
# variables diverses
#######################
near_zero = 1e-15
