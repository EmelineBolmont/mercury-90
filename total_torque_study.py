#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that search for every gnuplot script available
from __future__ import print_function

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "9 août 2011"
__version__ = "$Revision: 1.0 $"
__credits__ = """Display each script available and allow us to choose which one we want to run"""

import os
import subprocess
import pdb
import mercury # for the object Diskin, that we create in the script
import shutil # to copy file from a destination to another one
import numpy as np

storeDir = "torque_parameter_space"

# Range for the parameters
SIGMA_0 = (300., 3000.) # g/cm^2
SIGMA_INDEX = (0.3, 1.5) # The slope of the surface density power law
ADIABATIC_INDEX = (1.2, 1.7) # the adiabatic index of the disk
VISCOSITY = (1.e13, 1.e17) # cm^2/s
B_H = (0.3, 0.7) # the smoothing length of the gravitationnal potential of the planet
MEAN_MOLECULAR_WEIGHT = (2., 5.) # The mean molecular weight in mass of the proton

NBVAL_SIGMA_0 = 10 # Number of values to test in the range
NBVAL_SIGMA_INDEX = 4 # Number of values to test in the range
NBVAL_ADIABATIC_INDEX = 2 # Number of values to test in the range
NBVAL_VISCOSITY = 3 # Number of values to test in the range
NBVAL_B_H = 2 # Number of values to test in the range
NBVAL_MEAN_MOLECULAR_WEIGHT = 2 # Number of values to test in the range

total_number = 1
values = [NBVAL_SIGMA_0, NBVAL_SIGMA_INDEX, NBVAL_ADIABATIC_INDEX, NBVAL_VISCOSITY, NBVAL_B_H, NBVAL_MEAN_MOLECULAR_WEIGHT]
for val in values:
  total_number *= val

max_length = len(str(total_number))

#----------------------------------
# Parameters for the disk
# None parameters will not be displayed in the parameter file and thus, default values of the code will be used.
sigma_0 = 450 # g/cm^2
sigma_index = 0.5 # The slope of the surface density power law
adiabatic_index = 1.4 # the adiabatic index of the disk
viscosity = 1.e15 # cm^2/s
b_h = 0.4 # the smoothing length of the gravitationnal potential of the planet
mean_molecular_weight = 2.35 # The mean molecular weight in mass of the proton
sample = 400
dissipation_type = 0
disk_exponential_decay = None # in years
inner_boundary_condition = 'open'
outer_boundary_condition = 'open'
torque_type = 'real'
torque_file = 'torque_profile.dat'
torque_profile_steepness = None
indep_cz = None
mass_dep_m_min = None
mass_dep_m_max = None
mass_dep_cz_m_min = None
mass_dep_cz_m_max = None

def run(command):
  """run a command that will be a string.
  The function return a tuple with the output, 
  the stderr and the return code. 
  It is fitted to run only command placed in the parent directory (of type "../command")"""
  
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  
  (process_stdout, process_stderr) = process.communicate()
  returnCode = process.poll()
  
  if (returnCode == 0):
    #print(command+" was runned successfully\n")
    
    # We delete the last \n introduced by subprocess and that corrupt the diff managed to see differences between outputs.
    return (process_stdout, process_stderr)
  else:
    print(command+" has terminated with an error code: "+str(returnCode))
    
    print("\nstderr:")
    if (type(process_stderr) == list):
      for line in process_stderr:
        print(line)
    else:
      print(process_stderr)
    
    return (process_stdout, process_stderr)


######################
# Début du programme #
######################

################## Preamble   ###############
# We create the directory where to store the outputs files
if not(os.path.isdir(os.path.join("unitary_tests", storeDir))):
  os.mkdir(os.path.join("unitary_tests", storeDir))

# We compile the fortran program that will create the datas for the gnuplot
print("Compiling the fortran code")
run("python maketest.py test_disk.f90")

print("Creating the gnuplot file")
gnuplotScript = open(os.path.join("unitary_tests", "torque_parameter_space.gnuplot"), 'w')
gnuplotScript.write('set terminal pngcairo crop enhanced size 1200, 1000\n')
gnuplotScript.write('set output "total_torque.png"\n')
gnuplotScript.write('set xlabel "semi major axis (AU)"\n')
gnuplotScript.write('set ylabel "Planet mass (m_{earth})" center\n')
gnuplotScript.write('set title "Evolution of the total torque {/Symbol G}_{tot}/{/Symbol G}_0 "\n')
gnuplotScript.write('set pm3d map\n')
gnuplotScript.write('set pm3d explicit\n')
gnuplotScript.write('set palette rgbformulae 22,13,-31\n')
gnuplotScript.write('set grid xtics ytics linetype 0\n')
gnuplotScript.write('set xrange [  9.99999977648258209E-003 :   50.000000000000000      ]\n')
gnuplotScript.write('set yrange [  0.10000000149011612      :   60.000000000000007      ]\n')
gnuplotScript.write('splot "test_total_torque.dat" with pm3d notitle, \\\n')
gnuplotScript.write('      "contour_total_torque.dat" with line linetype -1 linewidth 1 notitle, \\\n')
gnuplotScript.write('      "test_vector_total_torque.dat" with vector notitle head filled linestyle -1\n')
gnuplotScript.close()

dataBase = open(os.path.join("unitary_tests", storeDir, "correspondance.log"), 'w')
dataBase.write("%7s %7s %9s %21s %15s %3s\n" % ("test_ID", "sigma_0", "viscosity", "mean_molecular_weight", "adiabatic_index", "b_h"))
dataBase.close()

test_ID = 1

############## Beginning of the loop, if any ###################
for sigma_0 in np.linspace(*SIGMA_0, num=NBVAL_SIGMA_0):
  for sigma_index in np.linspace(*SIGMA_INDEX, num=NBVAL_SIGMA_INDEX):
    for viscosity in np.linspace(*VISCOSITY, num=NBVAL_VISCOSITY):
      for adiabatic_index in np.linspace(*ADIABATIC_INDEX, num=NBVAL_ADIABATIC_INDEX):
        for b_h in np.linspace(*B_H, num=NBVAL_B_H):
          for mean_molecular_weight in np.linspace(*MEAN_MOLECULAR_WEIGHT, num=NBVAL_MEAN_MOLECULAR_WEIGHT):

            print("Simulation %i/%i" % (test_ID, total_number))
            print("\tWriting 'disk.in'")
            diskin = mercury.Disk(b_over_h=b_h, adiabatic_index=adiabatic_index, mean_molecular_weight=mean_molecular_weight, 
                          surface_density=(sigma_0, sigma_index), 
                          disk_edges=(1., 100.), viscosity=viscosity, sample=sample, dissipation_type=dissipation_type, 
                          disk_exponential_decay=disk_exponential_decay, torque_type=torque_type,
                          inner_boundary_condition=inner_boundary_condition, outer_boundary_condition=outer_boundary_condition,
                          torque_profile_steepness=torque_profile_steepness, indep_cz=indep_cz, mass_dep_m_min=mass_dep_m_min, 
                          mass_dep_m_max=mass_dep_m_max, mass_dep_cz_m_min=mass_dep_cz_m_min, mass_dep_cz_m_max=mass_dep_cz_m_max)
            diskin.write()

            print("\tRunning test_disk")
            (stdout, stderr) = run("./test_disk")

            os.chdir("unitary_tests")


            print("\tRunning gnuplot")
            (stdout, stderr) = run("gnuplot torque_parameter_space.gnuplot")
            source = "total_torque.png"
            dest = "%0*i_total_torque.png" % (max_length, test_ID)
            print("\tCopying file")
            #~ pdb.set_trace()
            shutil.copy2(source,os.path.join(storeDir,dest))

            os.chdir(storeDir)

            dataBase = open("correspondance.log", 'a')
            dataBase.write("%-7i %-7.1f %-9.2e %-21.2f %-15.2f %-3.2f\n" % (test_ID, sigma_0, viscosity, mean_molecular_weight, adiabatic_index, b_h))
            dataBase.close()
            
            test_ID += 1
            
            os.chdir("../..")






