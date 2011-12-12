#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# Of course, everything is not tested, but it is planed to test as many things as possible

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "21 Juillet 2011"
__version__ = "$Revision: 2.6.2 $"
__credits__ = """We run a test simulation and erase all the files created after the tests. The simulations files are thought to be 
in a "simu_test" subdirectory of the directory were are the sources (and binaries) of mercury (and this script)"""

from make import *
import os
import subprocess # To launch various process, get outputs et errors, returnCode and so on.
import difflib # To compare two strings
import pdb # To debug
import numpy as np # to calculate mean and standard deviation
from time import time
from mercury import * # In order to create a simulation via python

# We import all variables linked to original output files.
#from original_files import *

MERCURY_FILENAMES = ["info.out", "big.dmp", "small.dmp", "param.dmp", "restart.dmp", "big.tmp", "small.tmp", "param.tmp", "restart.tmp"]
ELEMENT_FILENAMES = ["APOLLO.aei", "JUPITER.aei", "MERCURY.aei", "ORPHEUS.aei", "TOUTATIS.aei", "EARTHMOO.aei", "KHUFU.aei", "MINOS.aei", 
  "PLUTO.aei", "URANUS.aei", "JASON.aei", "MARS.aei", "NEPTUNE.aei", "SATURN.aei", "VENUS.aei"]
CLOSE_FILENAMES = ["APOLLO.clo", "JUPITER.clo", "MERCURY.clo", "ORPHEUS.clo", "TOUTATIS.clo", "EARTHMOO.clo", "KHUFU.clo", "MINOS.clo", 
  "PLUTO.clo", "URANUS.clo", "JASON.clo", "MARS.clo", "NEPTUNE.clo", "SATURN.clo", "VENUS.clo"]

EXTENTION_ORIGINAL = ".ori"

MERCURY_FILENAMES_OLD = []
for file in MERCURY_FILENAMES:
  MERCURY_FILENAMES_OLD.append(file+EXTENTION_ORIGINAL)

ELEMENT_FILENAMES_OLD = []
for file in ELEMENT_FILENAMES:
  ELEMENT_FILENAMES_OLD.append(file+EXTENTION_ORIGINAL)
  
CLOSE_FILENAMES_OLD = []
for file in CLOSE_FILENAMES:
  CLOSE_FILENAMES_OLD.append(file+EXTENTION_ORIGINAL)

FOLDER = "simu_test"

##################
# Outputs of various binaries and tests to compare with the actual ones. 
# Theses outputs are those of the original version of mercury, that is, mercury6_2.for
##################

try:
  os.chdir(FOLDER)
except OSError:
  os.mkdir(FOLDER) # If the folder does not exist, we create it.
  
######################
# Génération de la simulation de base de mercury
######################
mercury = BodyCart("big",name="MERCURY", x=-3.83966017419175965E-01, y=-1.76865300855700736E-01, z=2.07959213998758705E-02,
vx=5.96286238644834141E-03, vy=-2.43281292146216750E-02, vz=-2.53463209848734695E-03, m=1.66013679527193009E-07,
r=20.0e0, d=5.43)
venus = BodyCart("big",name="VENUS", x=6.33469157915745540E-01, y=3.49855234102151691E-01, z=-3.17853172088953667E-02,
vx=-9.84258038001823571E-03, vy=1.76183746921837227E-02, vz=8.08822351013463794E-04, m=2.44783833966454430E-06,
r=20.0e0, d=5.24)
earthmoo = BodyCart("big", name="EARTHMOO", x=2.42093942183383037E-01, y=-9.87467766698604366E-01, z=-4.54276292555233496E-06,
vx=1.64294055023289365E-02, vy=4.03200725816140870E-03, vz=1.13609607260006795E-08, m=3.04043264264672381E-06,
r=20.0e0, d=5.52)
mars = BodyCart("big", name="MARS", x=2.51831018120174499E-01, y=1.52598983115984788E+00, z=2.57781137811807781E-02,
vx=-1.32744166042475433E-02, vy=3.46582959610421387E-03, vz=3.98930013246952611E-04, m=3.22715144505386530E-07,
r=20.0e0, d=3.94)
jupiter = BodyCart("big", name="JUPITER", x=4.84143144246472090E+00, y=-1.16032004402742839E+00, z=-1.03622044471123109E-01,
vx=1.66007664274403694E-03, vy=7.69901118419740425E-03, vz=-6.90460016972063023E-05, m=9.54791938424326609E-04,
r=3.0e0, d=1.33)
saturn = BodyCart("big", name="SATURN", x=8.34336671824457987E+00, y=4.12479856412430479E+00, z=-4.03523417114321381E-01,
vx=-2.76742510726862411E-03, vy=4.99852801234917238E-03, vz=2.30417297573763929E-05, m=2.85885980666130812E-04,
r=3.0e0, d=0.70)
uranus = BodyCart("big", name="URANUS", x=1.28943695621391310E+01, y=-1.51111514016986312E+01, z=-2.23307578892655734E-01,
vx=2.96460137564761618E-03, vy=2.37847173959480950E-03, vz=-2.96589568540237556E-05, m=4.36624404335156298E-05,
r=3.0e0, d=1.30)
neptune = BodyCart("big", name="NEPTUNE", x=1.53796971148509165E+01, y=-2.59193146099879641E+01, z=1.79258772950371181E-01,
vx=2.68067772490389322E-03, vy=1.62824170038242295E-03, vz=-9.51592254519715870E-05, m=5.15138902046611451E-05,
r=3.0e0, d=1.76)
pluto = BodyCart("big", name="PLUTO", x=-1.15095623952731607E+01, y=-2.70779438829451422E+01, z=6.22871533567077229E+00,
vx=2.97220056963797431E-03, vy=-1.69820233395912967E-03, vz=-6.76798264809371094E-04, m=7.39644970414201173E-09,
r=3.0e0, d=1.1)

apollo = BodyAst("small", name="APOLLO", a=1.4710345, e=.5600245, I=6.35621, 
g=285.63908, n=35.92313, M=15.77656, ep=2450400.5)
jason = BodyAst("small", name="JASON", a=2.2157309, e=.7644575, I=4.84834, 
g=336.49610, n=169.94137, M=293.37226, ep=2450400.5)
khufu = BodyAst("small", name="KHUFU", a=0.9894948, e=.4685310, I=9.91298, 
g=54.85927, n=152.64772, M=66.69818, ep=2450600.5)
minos = BodyAst("small", name="MINOS", a=1.1513383, e=.4127106, I=3.93863, 
g=239.50170, n=344.85893, M=8.93445, ep=2450400.5)
orpheus = BodyAst("small", name="ORPHEUS", a=1.2091305, e=.3226805, I=2.68180, 
g=301.55128, n=189.79654, M=28.31467, ep=2450400.5)
toutatis = BodyAst("small", name="TOUTATIS", a=2.5119660, e=.6335854, I=0.46976, 
g=274.82273, n=128.20968, M=50.00728, ep=2450600.5)

solarSystem = PlanetarySystem(bodies=[mercury, venus, earthmoo, mars, jupiter, saturn, uranus, 
neptune, pluto, apollo, jason, khufu, minos, orpheus, toutatis], m_star=1.0, epoch=2451000.5)

bigin = Big(solarSystem)
bigin.write()

smallin = Small(solarSystem)
smallin.write()

elementin = Element(format_sortie=" a8.5 e8.6 i8.4 g8.4 n8.4 l8.4 m13e ", coord="Cen", 
output_interval=365.2e1, time_format="years", relative_time="yes")
elementin.write()

closein = Close(time_format="years", relative_time="yes")
closein.write()

paramin = Param(algorithme="HYBRID", start_time=2451179.5, stop_time=2462502.5, output_interval=365.25e0, 
h=8, accuracy=1.e-12, stop_integration="no", collisions="no", fragmentation="no", 
time_format="years", relative_time="no", output_precision="medium", relativity="no", 
user_force="no", ejection_distance=100, radius_star=0.005, central_mass=1.0, 
J2=0, J4=0, J6=0, changeover=3., data_dump=500, periodic_effect=100)
paramin.write()

Files().write()
Message().write()

# We clean original files from mercury, element and close
clean(["ori"])

for algo in ["BS", "BS2", "MVS", "RADAU", "HYBRID"]:
	# We clean undesirable files. 
	clean(["out", "clo", "aei", "dmp", "tmp","ori"])
	
	paramin.set_algorithme(algo)
	paramin.write()
	print("##########################################")
	print("Running original binaries with "+algo+"...")
	(merc_or_stdout, merc_or_stderr) = run("../mercury_original/mercury")
	
	(clo_or_stdout, clo_or_stderr) = run("../mercury_original/close")

	(ele_or_stdout, ele_or_stderr) = run("../mercury_original/element")

	
	for file in MERCURY_FILENAMES:
		os.rename(file, file+EXTENTION_ORIGINAL)

	for file in CLOSE_FILENAMES:
	  os.rename(file, file+EXTENTION_ORIGINAL)
	  
	for file in ELEMENT_FILENAMES:
	  os.rename(file, file+EXTENTION_ORIGINAL)

	
	autiwa.printCR("Running new binaries with "+algo+"...")

	# We clean the .out files because we have not renamed them to .ori files
	clean(["out"])

	(merc_new_stdout, merc_new_stderr) = run("../mercury")

	(clo_new_stdout, clo_new_stderr) = run("../close")

	(ele_new_stdout, ele_new_stderr) = run("../element")
	print("Running new binaries with "+algo+"...ok")
	
	print("##########################################")


	diff = compare(merc_or_stdout, merc_new_stdout)
	if (diff != None):
	  print("\nTest of mercury with "+algo)
	  print("\tFor the Output of mercury")
	  print diff

	compare2file(MERCURY_FILENAMES_OLD, MERCURY_FILENAMES)

	diff = compare(clo_or_stdout, clo_new_stdout)
	if (diff != None):
	  print("\nTest of close with "+algo)
	  print("\tFor the Output of close")
	  print diff


	compare2file(CLOSE_FILENAMES_OLD, CLOSE_FILENAMES)

	diff = compare(ele_or_stdout, ele_new_stdout)
	if (diff != None):
	  print("\nTest of element with "+algo)
	  print("\tFor the Output of element")
	  print diff


	compare2file(ELEMENT_FILENAMES_OLD, ELEMENT_FILENAMES)


# CHANGELOG
# V2.6 : Je teste maintenant les 5 algorithmes usuels
# v2.6.1 : Rajout de lignes de séparation entre les tests des différents algorithmes pour plus de lisibilité.
# v2.6.2 : Suppression de la plupart des lignes si les fichiers n'ont pas de différences, ceci afin d'améliorer la lisibilité des 
#          différences justement.