#!/usr/bin/env python
# script to show the evolution of the disk and especially the torque profile in function of some parameters
# Version 2.2

import mercury                  # module that contain object for each parameter file (and also for bodies and planetary system)
import pdb                      # to debug via pdb.set_trace()
import os                       # to create folder, change directory and so on
import subprocess               # to launch 'runjob'
from constants import *         # Several constants, including mass of planets and so on. All constants are in CAPSLOCK
import shutil                   # In particular to copy a file with .copy2()
import sys                      # To enable parameters in the script (in particular)
from math import *
import numpy as np

#-------------------------------------------------------------------------------
# MANUAL : 

# To define the range, you can use : 
#~ np.linspace(start, stop, nb_points)
#~ np.logspace(log(start), log(stop), nb_points)

# If FIXED_TOTAL_MASS = True, then it's the parameter TOTAL_MASS which is important, else, it's NB_PLANETS that is important.

density_constant = None
density_index = None

#-------------------------------------------------------------------------------
# Definition of functions

def set_default_parameters():
  """set the default values for various parameters"""
  global surface_density, adiabatic_index, viscosity_type, viscosity, alpha, b_h, torque_type, disk_edges, inner_smoothing_width
  global tau_viscous, tau_photoevap, dissipation_time_switch, is_irradiation, opacity_type
  global dissipation_type, disk_exponential_decay, sample, inner_boundary_condition, outer_boundary_condition
  global torque_profile_steepness, indep_cz, mass_dep_m_min, mass_dep_m_max, mass_dep_cz_m_min, mass_dep_cz_m_max
  global is_turbulence, turbulent_forcing, saturation_torque, density_index, density_constant
  global mean_molecular_weight
  
  #----------------------------------
  # Parameters for the disk
  # None parameters will not be displayed in the parameter file and thus, default values of the code will be used.
  density_constant = 1700.
  density_index = 0.5
  surface_density = (density_constant, density_index) # (g/cm^2, power law index)
  density_file = 'surface_density_profile.dat'
  adiabatic_index = 1.4 # the adiabatic index of the disk
  viscosity_type = 'constant' # constant, alpha
  viscosity = 1.e15 # cm^2/s
  alpha = None
  opacity_type = 'bell' # 'bell' or 'zhu' or 'hure' opacity table
  b_h = 0.4 # the smoothing length of the gravitationnal potential of the planet
  mean_molecular_weight = 2.35
  sample = 400
  disk_edges = (0.5, 100.) # (the inner and outer edge of the disk in AU)
  inner_smoothing_width = 0.05 # (in unit of the inner boundary radius) , the width of the region where the surface density decay to become 0 at the inner edge
  dissipation_type = 0
  disk_exponential_decay = None # in years
  tau_viscous = None # in years
  tau_photoevap = None # in years
  dissipation_time_switch = None # in years
  inner_boundary_condition = 'open'
  outer_boundary_condition = 'open'
  torque_type = 'real'
  torque_file = 'torque_profile.dat'
  torque_profile_steepness = None
  saturation_torque = None
  indep_cz = None
  mass_dep_m_min = None
  mass_dep_m_max = None
  mass_dep_cz_m_min = None
  mass_dep_cz_m_max = None
  is_turbulence = None
  turbulent_forcing = None
  is_irradiation = 1

def set_prefered_disk():
  """set the default values for various parameters"""
  global surface_density, adiabatic_index, viscosity_type, viscosity, alpha, b_h, torque_type, disk_edges, inner_smoothing_width
  global tau_viscous, tau_photoevap, dissipation_time_switch, is_irradiation, opacity_type
  global dissipation_type, disk_exponential_decay, sample, inner_boundary_condition, outer_boundary_condition
  global torque_profile_steepness, indep_cz, mass_dep_m_min, mass_dep_m_max, mass_dep_cz_m_min, mass_dep_cz_m_max
  global is_turbulence, turbulent_forcing, saturation_torque, density_index, density_constant
  global mean_molecular_weight
  
  #----------------------------------
  # Parameters for the disk
  # None parameters will not be displayed in the parameter file and thus, default values of the code will be used.
  density_constant = 1700.
  #~ density_index = 0.5
  density_index = 1.0
  surface_density = (density_constant, density_index) # (g/cm^2, power law index)
  density_file = 'surface_density_profile.dat'
  adiabatic_index = 1.4 # the adiabatic index of the disk
  viscosity_type = 'constant' # constant, alpha
  viscosity = 1.e15 # cm^2/s
  alpha = None
  #~ opacity_type = 'bell' # 'bell' or 'zhu' or 'hure' opacity table
  opacity_type = 'hure' # 'bell' or 'zhu' or 'hure' opacity table
  #~ b_h = 0.4 # the smoothing length of the gravitationnal potential of the planet
  b_h = 0.6 # the smoothing length of the gravitationnal potential of the planet
  mean_molecular_weight = 2.35
  sample = 400
  disk_edges = (0.5, 100.) # (the inner and outer edge of the disk in AU)
  inner_smoothing_width = 0.05 # (in unit of the inner boundary radius) , the width of the region where the surface density decay to become 0 at the inner edge
  dissipation_type = 0
  disk_exponential_decay = None # in years
  tau_viscous = None # in years
  tau_photoevap = None # in years
  dissipation_time_switch = None # in years
  inner_boundary_condition = 'open'
  outer_boundary_condition = 'open'
  torque_type = 'real'
  torque_file = 'torque_profile.dat'
  torque_profile_steepness = None
  saturation_torque = None
  indep_cz = None
  mass_dep_m_min = None
  mass_dep_m_max = None
  mass_dep_cz_m_min = None
  mass_dep_cz_m_max = None
  is_turbulence = None
  turbulent_forcing = None
  is_irradiation = 1

def generation_simulation_parameters():
  """the function generate simulations files in the current working directory, given the parameters on top of the script
  """
  
  diskin = mercury.Disk(b_over_h=b_h, adiabatic_index=adiabatic_index, mean_molecular_weight=mean_molecular_weight, surface_density=surface_density, 
                is_irradiation=is_irradiation,
                disk_edges=disk_edges, viscosity_type=viscosity_type, viscosity=viscosity, alpha=alpha, opacity_type=opacity_type, 
                sample=sample, dissipation_type=dissipation_type, 
                is_turbulence=is_turbulence, turbulent_forcing=turbulent_forcing, inner_smoothing_width=inner_smoothing_width,
                tau_viscous=tau_viscous, tau_photoevap=tau_photoevap, dissipation_time_switch=dissipation_time_switch, 
                disk_exponential_decay=disk_exponential_decay, torque_type=torque_type,
                inner_boundary_condition=inner_boundary_condition, outer_boundary_condition=outer_boundary_condition,
                torque_profile_steepness=torque_profile_steepness, indep_cz=indep_cz, mass_dep_m_min=mass_dep_m_min, 
                saturation_torque=saturation_torque,
                mass_dep_m_max=mass_dep_m_max, mass_dep_cz_m_min=mass_dep_cz_m_min, mass_dep_cz_m_max=mass_dep_cz_m_max)
  diskin.write()
  
def systemCommand(command):
  """will launch a command system, without paying any attention to the outputs and error, just waiting the job to finish"""
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  
  if (returncode != 0):
    pdb.set_trace()
  # there is .poll() or .wait() but I don't remember the difference. For some kind of things, one of the two was not working
  return (process_stdout, process_stderr, returncode)

def study_parameter_influence():
  global surface_density
  
  os.chdir("parameter_study")
  if os.path.exists(folder_name):
    raise NameError("The folder %s already exists" % folder_name)
  else:
    os.mkdir(folder_name)
    
    recap_file = open("%s/readme.txt" % folder_name, 'w')
    
    recap_file.write("parameter : %s\n" % parameter_name)
    recap_file.write("values : %s\n" % parameter_values)
    recap_file.close()
    
  os.chdir("..")

  index = 1
  for value in parameter_values:
    
    globals()[parameter_name] = value
    
    surface_density = (density_constant, density_index)
    
    generation_simulation_parameters()
    
    sys.stdout.write("%3d %s : %s                          \r" % (index, parameter_name, value))
    sys.stdout.flush()
    
    (process_stdout, process_stderr, returncode) = systemCommand("./test_disk")
    (process_stdout, process_stderr, returncode) = systemCommand("(cd unitary_tests && gnuplot total_torque.gnuplot)")
    shutil.copy2("unitary_tests/total_torque.png", "parameter_study/%s/total_torque%05d.png" % (folder_name, index))
    shutil.copy2("disk.out", "parameter_study/%s/disk.out" % (folder_name))
    
    index += 1
  
#    .-.     .-.     .-.     .-.     .-.     .-.     .-.     .-.     .-. 
#  .'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `.
# (    .     .-.     .-.     .-.     .-.     .-.     .-.     .-.     .    )
#  `.   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   .'
#    )    )                                                       (    (
#  ,'   ,'                                                         `.   `.
# (    (                     DEBUT DU PROGRAMME                     )    )
#  `.   `.                                                         .'   .' 
#    )    )                                                       (    (
#  ,'   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   `.
# (    '  _  `-'  _  `-'  _  `-'  _  `-'  _  `-'  _  `-'  _  `-'  _  `    )
#  `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .' `.   .'
#    `-'     `-'     `-'     `-'     `-'     `-'     `-'     `-'     `-'

scriptFolder = os.path.dirname(os.path.realpath(__file__)) # the folder in which the module is. 
binaryPath = os.path.join(scriptFolder, os.path.pardir)

set_default_parameters()

if not(os.path.exists("parameter_study")):
  os.mkdir("parameter_study")

#-------------------------------------------------------------------------------
# Definition of parameters of the script
#~ for opacity in ['bell', 'zhu', 'hure', 'chambers']:
  #~ set_default_parameters()
  #~ opacity_type = opacity
  #~ parameter_name = 'mean_molecular_weight'
  #~ folder_name = "%s_%s" % (opacity_type, parameter_name)
  #~ parameter_values = np.linspace(0.6, 2.35, 10)
  #~ study_parameter_influence()
  
  #~ set_default_parameters()
  #~ opacity_type = opacity
  #~ parameter_name = 'is_irradiation'
  #~ folder_name = "%s_%s" % (opacity_type, parameter_name)
  #~ parameter_values = [0,1]
  #~ study_parameter_influence()
  #~ 
  #~ set_default_parameters()
  #~ opacity_type = opacity
  #~ parameter_name = 'density_index'
  #~ folder_name = "%s_%s" % (opacity_type, parameter_name)
  #~ parameter_values = np.linspace(0.1, 1.9, 19)
  #~ study_parameter_influence()
#~ 
  #~ set_default_parameters()
  #~ opacity_type = opacity
  #~ parameter_name = 'density_constant'
  #~ folder_name = "%s_%s" % (opacity_type, parameter_name)
  #~ parameter_values = np.linspace(50., 2500., 21)
  #~ study_parameter_influence()
#~ 
  #~ set_default_parameters()
  #~ opacity_type = opacity
  #~ parameter_name = 'viscosity'
  #~ folder_name = "%s_%s" % (opacity_type, parameter_name)
  #~ parameter_values = np.logspace(12, 16, 15)
  #~ study_parameter_influence()
#~ 
  #~ set_default_parameters()
  #~ opacity_type = opacity
  #~ parameter_name = 'adiabatic_index'
  #~ folder_name = "%s_%s" % (opacity_type, parameter_name)
  #~ parameter_values = np.linspace(1, 2., 10)
  #~ study_parameter_influence()
#~ 
  #~ set_default_parameters()
  #~ opacity_type = opacity
  #~ viscosity_type = 'alpha'
  #~ parameter_name = 'alpha'
  #~ folder_name = "%s_%s" % (opacity_type, parameter_name)
  #~ parameter_values = np.logspace(-2, -4., 10)
  #~ study_parameter_influence()
  
  #~ set_default_parameters()
  #~ opacity_type = opacity
  #~ parameter_name = 'surface_density'
  #~ M_tot = 30. # dust total mass in earth mass
  #~ R_int = 0.5 # AU
  #~ R_out = 100. # AU
  #~ dtg = 0.01 # Dust to gas ratio
  #~ 
  #~ sigma_index = np.linspace(0.5, 1.6, 12)
  #~ 
  #~ sigma_0 = (M_tot / dtg) / ((2. * pi / (2. - sigma_index)) * (R_out**(2. - sigma_index) - R_int**(2. - sigma_index))) # in earth_mass / AU**2
  #~ sigma_0 = sigma_0 * (MT / AU**2)
  #~ folder_name = "%s_constant_total_mass" % (opacity_type)
  #~ parameter_values = zip(sigma_0, sigma_index)
  #~ study_parameter_influence()

# Study sigma_0 in hure case
set_prefered_disk()
parameter_name = 'density_constant'
folder_name = "%s_%s" % ("steepness_1", parameter_name)
parameter_values = np.linspace(500., 2000, 16)
study_parameter_influence()
  
# Study viscosity in hure case
set_prefered_disk()
parameter_name = 'viscosity'
folder_name = "%s_%s" % ("steepness_1", parameter_name)
parameter_values = np.logspace(12, 16, 15)
study_parameter_influence()

# Study alpha in hure case
set_prefered_disk()
viscosity_type = 'alpha'
parameter_name = 'alpha'
folder_name = "%s_%s" % ("steepness_1", parameter_name)
parameter_values = np.logspace(-2, -4., 10)
study_parameter_influence()
