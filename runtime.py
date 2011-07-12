#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# Of course, everything is not tested, but it is planed to test as many things as possible

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "12 Juillet 2011"
__version__ = "$Revision: 1.1.1 $"
__credits__ = """We run a test simulation several times using original binairies and new binaries in order to compare their running time."""

from make import *
import pdb # To debug
import numpy as np # to calculate mean and standard deviation
from time import time
import os
from mercury import *
import progressbar
nb_runs = 10
start_time = 2451179.5
delta_t = 2000

FOLDER = "simu_test"

def cputime(command, method='time'):
  """function that return the CPU time in second"""
  if (method == "time"):
    (stdout, stderr) = run("time "+command)
    temp = stderr.split("\n")[-2]
    temp = temp.split("\t")
    temp = temp[-1][0:-2].split("m")
    cpu_time = float(temp[1])
    if (temp[0] != "0"):
      cpu_time += float(temp[0]) * 60.
  elif (method == "python"):
    start = time()
    run(command)
    cpu_time = time() - start
  
  return cpu_time

os.chdir(FOLDER)

# Not representative, it depend really on the execution. We should do a mean value on several runs, or a longer run.
print("""
#############################
# Executing time comparison #
#############################""")
time_merc_new = []
time_ele_new = []
time_clo_new = []

time_merc_old = []
time_ele_old = []
time_clo_old = []

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

paramin = Param(algorithme="HYBRID", start_time=start_time, stop_time=start_time+delta_t, output_interval=365.25e0, 
h=8, accuracy=1.e-12, stop_integration="no", collisions="no", fragmentation="no", 
time_format="years", relative_time="no", output_precision="medium", relativity="no", 
user_force="no", ejection_distance=100, radius_star=0.005, central_mass=1.0, 
J2=0, J4=0, J6=0, changeover=3., data_dump=500, periodic_effect=100)
paramin.write()

Files().write()
Message().write()

print("Comparison of running time, for an average on "+str(nb_runs)+" runs")
widgets = [progressbar.Percentage(), ' ', progressbar.Bar(marker='=', progress_marker='>'),' ', progressbar.ETA()]
pbar = progressbar.ProgressBar(widgets=widgets, maxval=2*nb_runs)
pbar.start()
for i in range(nb_runs):
  #The two types of binaries are runed in the same loop in order to be as partial as possible for the calculation of runtime.
  #Old binaries
  clean(["out", "clo", "aei", "dmp", "tmp"])
  
  time_merc_old.append(cputime("../mercury_original/mercury"))
  
  time_clo_old.append(cputime("../mercury_original/close"))
  
  time_ele_old.append(cputime("../mercury_original/element"))
  pbar.update(2*i+1)
  # New binaries
  clean(["out", "clo", "aei", "dmp", "tmp"])
  
  time_merc_new.append(cputime("../mercury"))
  
  time_clo_new.append(cputime("../close"))
    
  time_ele_new.append(cputime("../element"))
  pbar.update(2*(i+1))
pbar.finish()
# We calculate mean values and standard deviation for running time of mercury, element and close
t_merc_old = np.mean(time_merc_old)
t_ele_old = np.mean(time_ele_old)
t_clo_old = np.mean(time_clo_old)

dt_merc_old = np.std(time_merc_old)
dt_ele_old = np.std(time_ele_old)
dt_clo_old = np.std(time_clo_old)
  
t_merc_new = np.mean(time_merc_new)
t_ele_new = np.mean(time_ele_new)
t_clo_new = np.mean(time_clo_new)

dt_merc_new = np.std(time_merc_new)
dt_ele_new = np.std(time_ele_new)
dt_clo_new = np.std(time_clo_new)
  
# We determine the pourcentage of difference from original and new binaries
pourcent = (t_merc_new - t_merc_old) / t_merc_old
if (pourcent>0):
  merc_pourcent = " (+"+str(round(pourcent*100.,2))+"%)"
else:
  merc_pourcent = " ("+str(round(pourcent*100.,2))+"%)"

pourcent = (t_ele_new - t_ele_old) / t_ele_old
if (pourcent>0):
  ele_pourcent = " (+"+str(round(pourcent*100.,2))+"%)"
else:
  ele_pourcent = " ("+str(round(pourcent*100.,2))+"%)"


pourcent = (t_clo_new - t_clo_old) / t_clo_old
if (pourcent>0):
  clo_pourcent = " (+"+str(round(pourcent*100.,2))+"%)"
else:
  clo_pourcent = " ("+str(round(pourcent*100.,2))+"%)"

print("Binary\tOld\tNew")
print("mercury\t("+str(round(t_merc_old,2))+" ± "+str(round(dt_merc_old,2))+") s\t("+str(round(t_merc_new,2))+" ± "+str(round(dt_merc_new,2))+") s "+merc_pourcent)
print("element\t("+str(round(t_ele_old,2))+"±"+str(round(dt_ele_old,2))+") s\t("+str(round(t_ele_new,2))+" ± "+str(round(dt_ele_new,2))+") s "+ele_pourcent)
print("close\t("+str(round(t_clo_old,2))+"±"+str(round(dt_clo_old,2))+") s\t("+str(round(t_clo_new,2))+" ± "+str(round(dt_clo_new,2))+") s "+clo_pourcent)


#~ #############################
#~ # Executing time comparison #
#~ #############################
#~ Comparison of running time, for an average on 10 runs
#~ Binary	Old	New
#~ mercury	1.08+-0.05	1.41+-0.02 (+30.82%)
#~ element	0.04+-0.01	0.04+-0.01 (+13.49%)
#~ close	0.02+-0.01	0.02+-0.01 (-3.15%)
#~ 
#~ #############################
#~ # Executing time comparison #
#~ #############################
#~ Comparison of running time, for an average on 100 runs
#~ Binary	Old	New
#~ mercury	2.9±1.88	3.46±2.09 (+19.54%)
#~ element	0.04±0.02	0.05±0.03 (+14.03%)
#~ close	0.02±0.01	0.02±0.04 (+15.75%)
#~ 
#~ #############################
#~ # Executing time comparison #
#~ #############################
#~ Comparison of running time, for an average on 100 runs
#~ Binary	Old	New
#~ mercury	7.35±2.8	8.36±2.83 (+13.84%)
#~ element	0.09±0.01	0.11±0.05 (+20.01%)
#~ close	0.03±0.01	0.03±0.01 (-2.08%)

#~ #############################
#~ # Executing time comparison #
#~ #############################
#~ Comparison of running time, for an average on 10 runs
#~ 100% |==========================================================| Time: 00:07:13
#~ Binary	Old	New
#~ mercury	(11.97 ± 0.2) s	(27.26 ± 0.27) s  (+127.76%)
#~ element	(0.56±0.05) s	(0.7 ± 0.04) s  (+24.91%)
#~ close	(0.17±0.0) s	(0.17 ± 0.01) s  (+0.57%)



#TODO utiliser progressbar pour voir combien de temps il reste avant la fin.