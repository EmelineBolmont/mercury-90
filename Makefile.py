#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for fortran90

from make import *

# We clean undesirable files. Indeed, we will compile everything everytime.
clean(["o", "mod"])

sourceFile.setCompilator("gfortran")

sourceFile.setCompilingOptions("")#-O3 -march=native -pipe")

sources_filename = lister("*.f90")

# We create the binaries
make_binaries(sources_filename, ["mercury.f90", "element.f90", "close.f90"], debug=False, gdb=False, profiling=False)

#-O2 déclenche toutes les options d’optimisation de code garantissant un gain de performance sans
#altérer les réactions du programme. -O1 n’active que les options d’optimisation n’entraînant pas une
#augmentation de la taille de l’exécutable, tandis que -O3 active toutes les options d’optimisation, y
#compris celles qui ne garantissent pas nécessairement un gain dans les performances ou qui peuvent
#modiﬁer les réactions du programme.

#-fomit-frame-pointer indique au compilateur de ne pas conserver de pointeur de cadre dans un re-
#gistre pour les fonctions qui n’en ont pas besoin — ce pointeur ne sert, globalement, qu’au déboguage.
#Ceci évite les instructions pour sauver, mettre à jour et restaurer les pointeurs de cadre ; cela permet
#aussi de libérer un registre supplémentaire qui pourra dès lors être disponible dans de nombreuses fonc-
#tions. Cela rend également l’utilisation d’un débogueur impossible sur certaines machines, notamment
#les machines de type x86, c’est la raison pour laquelle G95 n’active pas cette option lorsque l’on a
#uniquement recourt à -O2
#Sur certaines machines, telles que le VAX, ce fanion n’a aucun effet car la séquence standard d’appel
#traite automatiquement le pointeur de cadre et rien n’est sauvé en prétendant qu’il n’existe pas. La macro
#de description de machine « FRAME_POINTER_REQUIRED » contrôle si une machine cible supporte
#ce drapeau.

#-ffast-math réorganise les calculs lorsque cela ne change rien au résultat mathématique, c’est-à-dire
#qu’avec cette option le compilateur utilise les règles de l’arithmétique, typiquement les factorisations,
#pour minimiser le nombre de calcul. Cette réorganisation est purement mathématique, elle ne porte que
#sur les formules. Par exemple, dans le cadre spéciﬁque de ce stage, une telle option n’implique pas la
#réorganisation des basculements d’un mode d’arrondi à l’autre. Toutefois, il faut utiliser cette option
#avec parcimonie car les spéciﬁcités de l’arithmétique des nombres ﬂottants font que l’on peut préférer
#un ordre de calcul qui donne la meilleure précision.

#-march=pentium4 permet de réaliser un code spéciﬁque pour l’architecture spéciﬁée, ici le Pentium
#4 d’Intel et les processeurs compatibles.

#Les options de la forme -falign-* permettent d’aligner les données sur la largeur du bus du proces-
#seur. Intel (voir [17]) précise qu’un mauvais alignement peut, dans le pire des cas, créer une pénalité de
#plus de 100 cycles. Le core P6 (cœur des processeurs de type Pentium 4) dispose de 3 décodeurs. Ils re-
#çoivent ainsi le ﬂot d’instructions par paquets. Lors d’une boucle, si l’adresse de saut n’est pas alignée,
#au pire, un seul décodeur peut servir pour décoder l’instruction suivante. D’où l’intérêt de l’alignement
#pour garantir que les 3 décodeurs fonctionnent en même temps le plus souvent.
#Il existe le même problème pour les données. En effet, le bus faisant 64 bits de large (en prenant
#l’exemple des architectures de type x86), lors d’accès sur des données plus petites, le processeur doit
#faire un décalage après avoir chargé les données pour récupérer les bonnes valeurs. Pour éviter ce
#décalage de bits, qui fait perdre du temps, on force l’alignement. Pour faciliter la tâche du compilateur,
#il est utile de classer les données par ordre décroissant d’encombrement mémoire : les ﬂottants en
#double précision d’abord, puis ceux en simple précision, puis les entiers et ainsi de suite.
