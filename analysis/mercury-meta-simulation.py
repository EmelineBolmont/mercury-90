#!/usr/bin/env python
# script to test various pieces of python code
# Version 2.2

import mercury_utilities        # module that contain utilities for the mercury simulations
import simulations_utilities    # module that contain utilities to help launch simulations, regardless of the kind of simulations
import mercury                  # module that contain object for each parameter file (and also for bodies and planetary system)
import pdb                      # to debug via pdb.set_trace()
import os                       # to create folder, change directory and so on
import subprocess               # to launch 'runjob'
from constants import *         # Several constants, including mass of planets and so on. All constants are in CAPSLOCK
import random
import shutil                   # In particular to copy a file with .copy2()
import sys                      # To enable parameters in the script (in particular)
from math import *
import numpy as np

FOLDER_PREFIX = "simu" # the prefix for each sub simulation folder
SUB_FOLDER_LOG = "random_parameters.in" # the name of the log file where we will store meta simulation information to keep a trace.
toLaunch = True # Do we launch the simulation once the files are created?

#-------------------------------------------------------------------------------
# MANUAL : 

# To define parameters, you currently have 3 possibilities :
# _ parameters = 3 : all the planets will have the same value
# _ parameters = (1., 0.3, 'gaussian') : values will follow a gaussian law of mean value 1 and standard deviation 0.3
# _ parameters = (1., 3, 'uniform') : value will be randomly generated following a uniform law between 1 and 3

# If FIXED_TOTAL_MASS = True, then it's the parameter TOTAL_MASS which is important, else, it's NB_PLANETS that is important.

#-------------------------------------------------------------------------------
# Definition of parameters of the script

#----------------------------------
# Parameters for the meta simulation

# Number of simulations we want to launch with the same properties
NB_SIMULATIONS = 1

WALLTIME = 48 # estimated duration of the job (currently only used for avakas)

#----------------------------------
# Parameters for mercury

epoch = 0
time_format = "years"
relative_time = "yes"
nb_outputs = 1000
nb_dumps = 100 # Number of times we will create .dmp files and output dE/E
aei_outputs = 2000
user_force = "yes"
# integration_time = 1e6 * YEARS
timestep = 3 # days
m_star = 1.

EJECTION_DISTANCE = 100

#----------------------------------
# Parameters for the planetary systems

# If we want to set the planets for a fixed total mass or a fixed number of planets
#~ FIXED_TOTAL_MASS = True
#~ TOTAL_MASS = 50 # earth mass (only used if FIXED_TOTAL_MASS = True)
#~ NB_PLANETS = 10 # only used if FIXED_TOTAL_MASS = False
#~ mass_parameters = (1, 3, "uniform") # the mass (in earth mass)
#~ a_parameters = (1, 20, "uniform") # the semi major axis (in AU)
#~ e_parameters = (0.001, 0.5, "uniform") # the eccentricity
#~ I_parameters = (0.01, 3, "uniform") # The inclination (in degrees)
radius_star = 0.005 # The radius of the central star in AU
#----------------------------------
# Parameters for the disk
# None parameters will not be displayed in the parameter file and thus, default values of the code will be used.
surface_density = (1700, 0.5) # (g/cm^2, power law index)
density_file = 'surface_density_profile.dat'
adiabatic_index = 1.4 # the adiabatic index of the disk
viscosity_type = 'constant'
viscosity = 1.e15 # cm^2/s
alpha = None
opacity_type = 'bell' # 'bell' or 'zhu' opacity table
b_h = 0.4 # the smoothing length of the gravitationnal potential of the planet
sample = 400
disk_edges = (1., 100.) # (the inner and outer edge of the disk in AU)
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
is_irradiation = None

#-------------------------------------------------------------------------------
# Definition of functions

def planet_from_dust(surface_density, edges, dust_to_gas_ratio, iceline_prop, planet_mass, previous_position):
  """ The aim is to get the surface densitiy profile, and derive a dust density profile, integrated to have a linear density. 
  Based on the list of mass we want, we get the positions of the planets, integrating the dust surface profile. When the mass of 
  the first object is attained, the position wil be the position of the object, and so on, moving on the cumulated mass profile, 
  that is an integration of the linear density profile of dust

  Parameters : 
  surface_density : (Sigma_0, index) Sigma_0 in g/cm^2
  iceline_prop : a tuple of the shape : (radius_of_the_iceline [AU], factor). 
             After the given radius, the original dust to gas 
             ratio will have to be multiplied by 'factor'
  dust_to_gas_ratio : The ratio that allow us to get the dust density from the gas density
  planet_mass : the desired planet mass, in earth_mass
  previous_position : the position of the previous planet (or the inner edge of the disk, for the first planet)
  
  Distances are in AU, and masses in earth mass

  Surface density(R) = Sigma_0 * R**(-index)
  linear_density of dust = prefactor * DGR(R) * R**(1. - index)

  """
  
  (Sigma_0, index) = surface_density
  
  prefactor = 2 * pi * Sigma_0 * ((AU * 1e2)**2) / (MT * 1e3)
  
  # We generate a list of position and dust to gas ratio in the disk
  # The first element is the inner edge and the last the outer edge
  # Between each consecutive radius, the good dust to gas ratio is 
  # the one of the first inner radius of the couple
  (iceline, factor) = iceline_prop
  
  dust_radius = [edges[0]]
  dust = [dust_to_gas_ratio]
  
  dust_radius.append(iceline)
  dust.append(dust[0] * factor)
  
  dust_radius.append(edges[-1])
  dust.append(dust[-1])

  # We search for the position of the next planet, given its mass, 
  #  and the position of the previous planet
  # Since the dust to gas ratio varies inside the disk, we check in 
  #  which step of the step function the planet will lies
  
  # We search for the first step to start with
  step_index = np.searchsorted(dust_radius, previous_position)
  
  
  desired_mass = planet_mass
  current_position = previous_position
  current_dtg = dust[max(step_index - 1, 0)]
  for (dtg, iceline) in zip(dust[step_index:], dust_radius[step_index:]):
    mass_of_the_ring = prefactor * current_dtg / (2. - index) * \
                       (iceline**(2. - index) - current_position**(2. - index))
    
    # We substract the mass of the ring if the planet is going to be far away. 
    # But if the planet is inside this ring we stop the loop
    if (mass_of_the_ring < desired_mass):
      desired_mass -= mass_of_the_ring
      current_position = iceline
      current_dtg = dtg
    else:
      break
  new_position = (desired_mass * (2. - index) / (prefactor * current_dtg)
                + current_position**(2. - index))**(1./(2. - index))
  
  return new_position


def generate_meta_simulationin():
  """function that will generate a "meta_simulation.in" file to show how this works"""
  
  DEMO_FILE = \
"""#!/usr/bin/env python
## for a shortened parameter file, only "integration_time", 
## "mass_parameters", "a_parameters", "e_parameters", "I_parameters", 
## "FIXED_TOTAL_MASS" (and TOTAL_MASS or NB_PLANETS) are absolutely needed. 
## Other parameters contains default values inside the original script.
##-------------------------------------------------------------------------------
## Definition of parameters of the script

##----------------------------------
## Parameters for the meta simulation

## Number of simulations we want to launch with the same properties
NB_SIMULATIONS = 2

## The limited time above which the simulation will be stopped in avakas (useless for other servers)
WALLTIME = 119 # in hours

##----------------------------------
## Parameters for mercury

integration_time = 1e7 # in years
time_format = "years" # days or years
relative_time = "yes" # yes or no
nb_outputs = 100000 # number of outputs contained in xv.out (max possibility of accuracy in time)
nb_dumps = 1000 # number of outputs contained in generated .aei files by element
user_force = "yes" # yes or no, if we use user_module
timestep = 0.4 # days

##----------------------------------
## Parameters for the planetary systems

## If FIXED_TOTAL_MASS = True, then it's the parameter TOTAL_MASS which is important, else, it's NB_PLANETS that is important.
FIXED_TOTAL_MASS = True
TOTAL_MASS = 20 # earth mass, only used if FIXED_TOTAL_MASS = True
# NB_PLANETS = 5 # (number) only used if FIXED_TOTAL_MASS = False


## To define parameters, you currently have 3 possibilities :
## _ parameters = 3 : all the planets will have the same value
## _ parameters = (1., 0.3, 'gaussian') : values will follow a gaussian law of mean value 1 and standard deviation 0.3
## _ parameters = (1., 3, 'uniform') : value will be randomly generated following a uniform law between 1 and 3
## _ parameters = ([10.], (1., 3, 'uniform')) we copy the first list and complete it given the other parameters. 

##/!\ For mass, this only works if FIXED_TOTAL_MASS = False
mass_parameters = (0.1, 1, "uniform") # the mass (in earth mass)

##  For 'a', there is a special option that allow us to define a hill radii separation. (a_0, delta, 'rh')
## a_0 and delta can be calculated randomly for each planet if, instead of a value, a tuple of two values is given
## a_parameters = ((1, 1.5), (4., 6.), "rh")
##  Another one is 'from-dust' that will generate the positions of the planets using the masses
## (dust_to_gas_ratio, iceline, 'from-dust'), where 'icelines' is a tuple of 2 values (radius [AU], factor), 
## which means that from the original dtg (for the innermost region)
## , at this point you have to apply the given factor
#a_parameters = (0.01, (4., 2.), "from-dust")
a_parameters = ((1, 1.5), (4., 6.), "rh") # the semi major axis (in AU)

e_parameters = (1e-3, 1e-5, "uniform") # the eccentricity

I_parameters = (-1, 1, "uniform") # The inclination (in degrees)

radius_star = 0.05 # The radius of the central star in AU

##----------------------------------
## Parameters for the disk

## Surface density. We can define a power law through :
##    surface_density =  (sigma_0, alpha) : Sigma(R) = Sigma_0 * R**(-alpha)
##  But we can also define a manual surface density profile with 
##    surface_density =  "manual"
##  then, the file 'surface_density_profile.dat' must exist in the folder of meta_simulation.in so that it can be copied
surface_density = (500, 0.5) # g/cm^2 (sigma_0, alpha) help to define a power law : Sigma(R) = Sigma_0 * R**(-alpha)

adiabatic_index = 1.4 # The adiabatic index of the disk

viscosity_type = constant # constant, alpha
viscosity = 1.e15 # cm^2/s
# alpha = 1.e-2 # adim


opacity_type = bell # bell, chambers, zhu or hure

disk_edges = (0.1, 100.) # (the inner and outer edge of the disk in AU)
## The width of the region where the surface density decay to become 0 at the inner edge
inner_smoothing_width = 0.05 # (in unit of the inner boundary radius)

b/h = 0.4 # The smoothing width for gravitational effects

sample = 1000 # The number of points for the surface density profile

dissipation_type = 0 # The type of dissipation for the disk (0 for none)
#inner_boundary_condition = 'open' # 'open' or 'closed' (for dissipation_type=1)
#outer_boundary_condition = 'open' # 'open' or 'closed' (for dissipation_type=1)

#disk_exponential_decay = 1e6 # years (for dissipation_type=2)

#tau_viscous = 1e7 # years (for dissipation_type=3)
#tau_photoevap = 1e5 # years (for dissipation_type=3)
#dissipation_time_switch = 2e6 # years (for dissipation_type=3)

## We define the type of torque profile : 
## real : the profile defined by paardekooper formulas, nothing else to add in parameters here, everything is in the code
## manual : the file 'torque_profile.dat' must exist in the folder of meta_simulation.in and 
##          will contains two colums (first distance in AU, second the torque in units of Gamma_0)
## other options : linear_indep, tanh_indep, mass_dependant
torque_type = real # real, linear_indep, tanh_indep, mass_dependant, manual

#indep_cz = 3.0 # in AU, the position of the convergence zone (for linear or tanh_indep)
#saturation_torque = 1.0 # (in Gamma_0) the torque saturation value (for tanh_indep)
#torque_profile_steepness = 1.0 # the steepness for linear torque dependance (for linear_indep and mass_dependant)

## help to define a CZ(m) by defining two points (for mass_dependant)
#mass_dep_m_min = 5 # (earth mass)  help to define a CZ(m) by defining two points
#mass_dep_m_max = 20 # (earth mass) help to define a CZ(m) by defining two points
#mass_dep_cz_m_min = 50 # (AU)      help to define a CZ(m) by defining two points
#mass_dep_cz_m_max = 5 # (AU)       help to define a CZ(m) by defining two points

is_irradiation = 1 # is there stellar irradiation or not?
is_turbulence = 0 # is there turbulence or not?
#turbulent_forcing = 1e-4 # the turbulence forcing associated with turbulence, if any
"""

  demo_parameter_file = open("meta_simulation.in", 'w')
  demo_parameter_file.write(DEMO_FILE)
  demo_parameter_file.close()

def generate_diskin():
  """function that will generate a default 'disk.in' file"""
  
  diskin = mercury.Disk()
  diskin.demo()

def readParameterFile(parameter_file, COMMENT_CHARACTER="#", PARAMETER_SEPARATOR="="):
  """function that read the parameter file associated with meta 
  simulations to get values that define the meta simulation
  
  PARAMETER :
  parameter_file : the name of the file where are stored the parameters
  COMMENT_CHARACTER="#" : after this character in the parameter file, all the rest of the line will be ignored
  PARAMETER_SEPARATOR="=" : This character will separate the key and the value for a given parameter
  
  ACTION :
  This function doesn't return anything, but will store in global variable the values stored in the parameter file given in parameter.
  """
  
  global NB_SIMULATIONS, WALLTIME, isErase
  global integration_time, time_format, relative_time, nb_outputs, nb_dumps, user_force, timestep
  global FIXED_TOTAL_MASS, TOTAL_MASS, NB_PLANETS, mass_parameters, a_parameters, e_parameters, I_parameters, radius_star
  global surface_density, adiabatic_index, viscosity_type, viscosity, alpha, b_h, torque_type, disk_edges, inner_smoothing_width
  global tau_viscous, tau_photoevap, dissipation_time_switch, is_irradiation, opacity_type
  global dissipation_type, disk_exponential_decay, sample, inner_boundary_condition, outer_boundary_condition
  global torque_profile_steepness, indep_cz, mass_dep_m_min, mass_dep_m_max, mass_dep_cz_m_min, mass_dep_cz_m_max
  global is_turbulence, turbulent_forcing, saturation_torque
  global PARAMETERS
  
  
  f = open(parameter_file, 'r')
  lines = f.readlines()
  f.close()
  
  # Variable to store various information about the meta simulation that we will store in each sub simulation folder.
  PARAMETERS = ""

  for line in lines:
    # we get rid of any comments
    line = line.split(COMMENT_CHARACTER)[0]
    line = line.replace("\n", "")
    line = line.replace(" ", "")

    if (line.count(PARAMETER_SEPARATOR) == 1):
      (key, value) = line.split(PARAMETER_SEPARATOR)
      #----------------------------------
      # Meta simulation parameters
      if (key in ["nb_simulations", "NB_SIMULATIONS"]):
        NB_SIMULATIONS = int(value)
      elif (key in ["walltime", "WALLTIME"]):
        WALLTIME = int(value)
      #----------------------------------
      # Mercury parameters
      elif (key == "integration_time"):
        integration_time = float(value) * YEAR
      elif (key == "time_format"):
        time_format = simulations_utilities.str2str(value)
      elif (key == "relative_time"):
        relative_time = simulations_utilities.str2str(value)
      elif (key == "nb_outputs"):
        nb_outputs = int(value)
      elif (key == "nb_dumps"):
        nb_dumps = int(value)
      elif (key == "user_force"):
        user_force = simulations_utilities.str2str(value)
      elif (key == "timestep"):
        timestep = float(value)
      #----------------------------------
      # Planetary system parameters
      elif (key in ["fixed_total_mass", "FIXED_TOTAL_MASS"]):
        FIXED_TOTAL_MASS = simulations_utilities.str2bool(value)
      elif (key in ["total_mass", "TOTAL_MASS"]):
        
        TOTAL_MASS = float(value)
      elif (key in ["nb_planets", "NB_PLANETS"]):
        NB_PLANETS = int(value)
      elif (key == "mass_parameters"):
        mass_parameters = eval(value) # We must use 'eval' because else, the tuple will not be set correctly. I don't know better solution for now
      elif (key == "a_parameters"):
        a_parameters = eval(value)
      elif (key == "e_parameters"):
        e_parameters = eval(value)
      elif (key == "I_parameters"):
        I_parameters = eval(value)
      elif (key == "radius_star"):
        radius_star = float(value)
      #----------------------------------
      # Disk parameters
      elif (key == "surface_density"):
        if (value=='manual'):
          surface_density = value
        else:
          surface_density = eval(value)
      elif (key == "is_irradiation"):
        try:
          is_irradiation = int(value)
        except ValueError:
          raise ValueError("'is_irradiation' must be equal to 0 or 1")
      elif (key == "disk_edges"):
        disk_edges = eval(value)
      elif (key == "inner_smoothing_width"):
        inner_smoothing_width = float(value)
      elif (key == "adiabatic_index"):
        adiabatic_index = float(value)
      elif (key == "viscosity_type"):
        viscosity_type = simulations_utilities.str2str(value)
      elif (key == "viscosity"):
        viscosity = float(value)
      elif (key == "alpha"):
        alpha = float(value)
      elif (key == "opacity_type"):
        opacity_type = simulations_utilities.str2str(value)
      elif (key == "is_turbulence"):
        try:
          is_turbulence = int(value)
        except ValueError:
          raise ValueError("'is_turbulence' must be equal to 0 or 1")
      elif (key == "turbulent_forcing"):
        turbulent_forcing = float(value)
      elif (key == "b/h"):
        b_h = eval(value)
      elif (key == "dissipation_type"):
        dissipation_type = int(value)
      elif (key == "disk_exponential_decay"):
        disk_exponential_decay = float(value)
      elif (key == "tau_viscous"):
        tau_viscous = float(value)
      elif (key == "tau_photoevap"):
        tau_photoevap = float(value)
      elif (key == "dissipation_time_switch"):
        dissipation_time_switch = float(value)
      elif (key == "sample"):
        sample = int(value)
      elif (key == "inner_boundary_condition"):
        # We ensure that there is no extra quote in the string with the function str2str that remove extra ' and "
        # There extra quote will generate an error in the mercury.py class Disk
        inner_boundary_condition = simulations_utilities.str2str(value)
      elif (key == "outer_boundary_condition"):
        outer_boundary_condition = simulations_utilities.str2str(value)
      elif (key == "torque_type"):
        torque_type = simulations_utilities.str2str(value)
      elif (key == "torque_profile_steepness"):
        torque_profile_steepness = float(value)
      elif (key == "saturation_torque"):
        saturation_torque = float(value)
      elif (key == "indep_cz"):
        indep_cz = float(value)
      elif (key == "mass_dep_m_min"):
        mass_dep_m_min = float(value)
      elif (key == "mass_dep_m_max"):
        mass_dep_m_max = float(value)
      elif (key == "mass_dep_cz_m_min"):
        mass_dep_cz_m_min = float(value)
      elif (key == "mass_dep_cz_m_max"):
        mass_dep_cz_m_max = float(value)
      else:
        PARAMETERS += "the parameter %s is not known. His value was : %s\n" % (key, value)
        raise ValueError("no parameter is known for the key '%s'" % key)
      # Each parameter defined in the parameter file MUST BE a 'global()' variable, defined in the preamble of the function, in 
      # order to retrieve their values in the rest of the program
  
  if (FIXED_TOTAL_MASS==True):
    if ('TOTAL_MASS' not in globals()):
      print("You defined a Fixed total mass simulation but did not set the 'TOTAL_MASS' parameter")
      exit()
    elif ('NB_PLANETS' in globals()):
      print("Warning: NB_PLANETS has been defined but will not be used because FIXED_TOTAL_MASS is True")
  else:
    if ('NB_PLANETS' not in globals()):
      print("You defined a Fixed number of planets simulation but did not set the 'NB_PLANETS' parameter")
      exit()
    elif ('TOTAL_MASS' in globals()):
      print("Warning: TOTAL_MASS has been defined but will not be used because FIXED_TOTAL_MASS is False")
  
  # We prepare a variable we will use to store a log file of the parameters in each sub simulation folder
  PARAMETERS += "Here are the parameters used to generate the simulation.\n"
  PARAMETERS += "----------------------------------\nMercury Parameters\n\n"
  PARAMETERS += "integration time = %.2e years\n" % (integration_time/365.25)
  PARAMETERS += "number of outputs = "+str(nb_outputs)+"\n"
  PARAMETERS += "number of dumps = "+str(nb_dumps)+"\n"
  PARAMETERS += "user force = "+str(user_force)+"\n"
  PARAMETERS += "timestep = %f\n" % timestep
  PARAMETERS += "----------------------------------\nPlanetary System Parameters\n\n"
  PARAMETERS += "radius_star = %f\n" % radius_star
  PARAMETERS += "fixed total mass = "+str(FIXED_TOTAL_MASS)+"\n"
  if (FIXED_TOTAL_MASS==True):
    PARAMETERS += "total mass = "+str(TOTAL_MASS)+"\n"
  else:
    PARAMETERS += "number of planets = "+str(NB_PLANETS)+"\n"
  PARAMETERS += "mass parameters = "+str(mass_parameters)+"\n"
  PARAMETERS += "a parameters = "+str(a_parameters)+"\n"
  PARAMETERS += "e parameters = "+str(e_parameters)+"\n"
  PARAMETERS += "I parameters = "+str(I_parameters)+"\n"
  PARAMETERS += "----------------------------------\nDisk Parameters\n\n"
  if (surface_density == 'manual'):
    PARAMETERS += "surface_density = manual\n"
  else:
    PARAMETERS += "sigma_0 = %f g/cm^2 ; negative power law index = %f\n" % surface_density
  PARAMETERS += "is_irradiation = "+str(is_irradiation)+"\n"
  PARAMETERS += "sample = "+str(sample)+"\n"
  PARAMETERS += "disk_edges = "+str(disk_edges)+" AU\n"
  PARAMETERS += "inner_smoothing_width = "+str(inner_smoothing_width)+" AU\n"
  PARAMETERS += "adiabatic_index = "+str(adiabatic_index)+"\n"
  PARAMETERS += "viscosity_type = "+str(viscosity_type)+"\n"
  if (viscosity_type == 'constant'):
    PARAMETERS += "viscosity = "+str(viscosity)+" cm^2/s\n"
  if (viscosity_type == 'alpha'):
    PARAMETERS += "alpha = "+str(alpha)+"\n"
  PARAMETERS += "is_turbulence = "+str(is_turbulence)+"\n"
  PARAMETERS += "turbulent_forcing = "+str(turbulent_forcing)+"\n"
  PARAMETERS += "b/h = "+str(b_h)+"\n"
  PARAMETERS += "dissipation_type = "+str(dissipation_type)+"\n"
  if (dissipation_type == 2):
    PARAMETERS += "disk_exponential_decay = "+str(disk_exponential_decay)+"\n"
  if (dissipation_type == 3):
    PARAMETERS += "tau_viscous = "+str(tau_viscous)+" years\n"
    PARAMETERS += "tau photoevap = "+str(tau_photoevap)+" years\n"
    PARAMETERS += "switch time = "+str(dissipation_time_switch)+" years\n"
  PARAMETERS += "inner_boundary_condition = "+str(inner_boundary_condition)+"\n"
  PARAMETERS += "outer_boundary_condition = "+str(outer_boundary_condition)+"\n"
  PARAMETERS += "torque_type = "+str(torque_type)+"\n"
  PARAMETERS += "torque_profile_steepness = "+str(torque_profile_steepness)+"\n"
  PARAMETERS += "saturation_torque = "+str(saturation_torque)+"\n"
  PARAMETERS += "indep_cz = "+str(indep_cz)+"\n"
  PARAMETERS += "mass_dep_m_min = "+str(mass_dep_m_min )+"\n"
  PARAMETERS += "mass_dep_m_max = "+str(mass_dep_m_max )+"\n"
  PARAMETERS += "mass_dep_cz_m_min = "+str(mass_dep_cz_m_min )+"\n"
  PARAMETERS += "mass_dep_cz_m_max = "+str(mass_dep_cz_m_max )+"\n"




def generation_simulation_parameters():
  """the function generate simulations files in the current working directory, given the parameters on top of the script
  """
  
  output_interval = integration_time / float(nb_outputs)
  data_dump = abs(int(integration_time / (nb_dumps * timestep)))
  
  if (aei_outputs < nb_outputs):
    aei_time = integration_time / float(aei_outputs)
  else:
    aei_time = 0.0 # We want to see all the outputs of xv.out

  # We begin the random generation of parameters.
  if FIXED_TOTAL_MASS:
    total_mass = 0
    m = []
    if (type(mass_parameters[0]) == list):
      print("Warning: We can only specify a list of pre-planets in the case where the number of planet is fixed (and not the total mass)")
      exit()
    
    while (total_mass < TOTAL_MASS):
      m0 = simulations_utilities.setParameter(mass_parameters, 1, vmin=0.)[0] # The function return a list, and we want to have one element
      # the mass must be expressed relatively to the mass of the central body
      m.append(simulations_utilities.significativeRound(m0 * MT / MS, 4))
      total_mass += m0
    nb_planets = len(m)
  else:
    nb_planets = NB_PLANETS
    m = simulations_utilities.setParameter(mass_parameters, nb_planets, vmin=0.)
    m = [mi * MT / MS for mi in m]
  
  # By default, we do not specify density
  d = [None for mi in m]
  
  #~ # To add manually one massive planet in the set of planets.
  #~ random_index = random.randint(0, nb_planets)
  #~ random_mass = random.uniform(4, 5)
  #~ m.insert(random_index, random_mass * MT / MS)
  #~ nb_planets += 1

# This type is special and only for orbital distance. It define the various distances upon the mutual hill radii distance. 
# It is possible, by including two values in a tuple instead of a single value to define a random calculation of the distance 
# (whether the initial radius of the innermost planet of the hill separation.)
  if (type(a_parameters) in [list, tuple] and a_parameters[-1] == 'rh'):
    if (type(a_parameters[0]) == tuple):
      a = [random.uniform(a_parameters[0][0], a_parameters[0][1])]
    else:
      a = [a_parameters[0]]
      
    # If the value is a tuple, we generate a random value of delta between the two values, else, we take the only value.
    if (type(a_parameters[1]) == tuple):
      delta = [random.uniform(a_parameters[1][0], a_parameters[1][1]) for i in range(nb_planets)] # the first value will never be used
    else:
      delta = [a_parameters[1] for i in range(nb_planets)] # the first value will never be used
      
    for i in range(1,nb_planets):
      chi = delta[i] / 2. * ((m[i-1] + m[i]) / (3 * m_star))**(1/3.)
      ai = a[i-1] * (1 + chi) / (1 - chi)
      a.append(ai)
  elif (type(a_parameters) in [list, tuple] and a_parameters[-1] == 'from-dust'):
    a = []
    previous_position = disk_edges[0]
    dust_to_gas_ratio = a_parameters[0]
    iceline = a_parameters[1]
    HIGH_DENSITY = 3.0 # g/cm^3
    LOW_DENISTY  = 1.0 # g/cm^3
    
    # We will define densities
    d = []
    
    for mi in m:
      planet_mass = mi * (MS/MT) # m in solar mass, we want it in earth mass
      ai = planet_from_dust(surface_density, disk_edges, dust_to_gas_ratio, iceline, planet_mass, previous_position)
      a.append(ai)

      if (ai < iceline[0]):
        di = HIGH_DENSITY
      else:
        di = LOW_DENISTY
      d.append(di)
      previous_position = ai
    #~ prefactor = 2 * pi * 500. * ((AU * 1e2)**2) / (MT * 1e3)
    #~ total_mass = prefactor / 1.5 * 0.01 * (iceline[0]**1.5 - 0.1**1.5)
    #~ total_mass += prefactor / 1.5 * 0.01 * iceline[1] * (ai**1.5 - iceline[0]**1.5)
    #~ print(np.sum(m) * (MS/MT))
    #~ print("total_mass=%f" % total_mass)
    
  else:
    a = simulations_utilities.setParameter(a_parameters, nb_planets, vmin=0.)
    
  e = simulations_utilities.setParameter(e_parameters, nb_planets, vmin=0., vmax=1.)
  I = simulations_utilities.setParameter(I_parameters, nb_planets)
  
  # If there is planet beyong the ejection distance, we remove them before the simulation start :
  m_to_del = []
  a_to_del = []
  e_to_del = []
  I_to_del = []
  for (mi, ai, ei, Ii) in zip(m, a, e, I):
    if (ai > EJECTION_DISTANCE):
      m_to_del.append(mi)
      a_to_del.append(ai)
      e_to_del.append(ei)
      I_to_del.append(Ii)
  for (mi, ai, ei, Ii) in zip(m_to_del, a_to_del, e_to_del, I_to_del):
    m.remove(mi)
    a.remove(ai)
    e.remove(ei)
    I.remove(Ii)
  
  len_m = len(m)
  len_a = len(a)
  len_e = len(e)
  len_I = len(I)
  
  if ((len_m != len_a) or (len_a != len_e) or (len_e != len_I)):
    raise TypeError("The size of {m,a,e,I} is not the same")

  system = mercury_utilities.definePlanetarySystem(m=m, a=a, e=e, I=I, d=d)

  # We write the files

  bigin = mercury.Big(system)
  bigin.write()

  smallin = mercury.Small(system)
  smallin.write()
  if not(os.path.exists("small.in")):
    pdb.set_trace()

  # Setting the output interval to 0 ensure that we will have every output written in the xv.out in the element files.
  elementin = mercury.Element(format_sortie=" a8.5 e8.6 i8.4 g8.4 n8.4 l8.4 m13e ", coord="Cen", 
  output_interval=aei_time, time_format=time_format, relative_time=relative_time)
  elementin.write()

  closein = mercury.Close(time_format=time_format, relative_time=relative_time)
  closein.write()

  paramin = mercury.Param(algorithme="HYBRID", start_time=0, stop_time=integration_time, output_interval=output_interval, 
  h=timestep, accuracy=1.e-12, stop_integration="no", collisions="yes", fragmentation="no", 
  time_format=time_format, relative_time=relative_time, output_precision="medium", relativity="no", 
  user_force=user_force, ejection_distance=EJECTION_DISTANCE, radius_star=radius_star, central_mass=1.0, 
  J2=0, J4=0, J6=0, changeover=3., data_dump=data_dump, periodic_effect=100)
  paramin.write()

  filesin = mercury.Files()
  filesin.write()

  messagein = mercury.Message()
  messagein.write()
  
  

  if (user_force == "yes"):
    diskin = mercury.Disk(b_over_h=b_h, adiabatic_index=adiabatic_index, mean_molecular_weight=2.35, surface_density=surface_density, 
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
    
    # If we want to use a manual torque profile, we copy the torque profile from the current working directory
    if (torque_type == 'manual'):
      if (os.path.isfile("../"+torque_file)):
        shutil.copy2("../"+torque_file, "./"+torque_file)
      else:
        raise NameError("The file "+torque_file+" does not exist in the parent directory\n keep in mind to delete the created folder.")
        
        
    # If we want to use a manual surface density profile, we copy the density profile from the current working directory
    if (surface_density == 'manual'):
      if (os.path.isfile("../"+density_file)):
        shutil.copy2("../"+density_file, "./"+density_file)
      else:
        raise NameError("The file "+density_file+" does not exist in the parent directory\n keep in mind to delete the created folder.")
        
  
  # We store a log file of the simulation parameters
  f = open(SUB_FOLDER_LOG, 'w')
  f.write(PARAMETERS)
  f.close()
  
  # We reset the counters for planet names
  mercury.Body.resetCounter()

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

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * norun : will create the various folders and file, but will not run the simulation" + "\n" + \
" * disk : will create a 'disk.in' demo file (without anything else)" + "\n" + \
" * demo : will create a 'meta_simulation.in' file " + "\n" + \
"   (needed by the current script) to show what can be defined"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'norun'):
    toLaunch = False
  elif (key == 'demo'):
    print("A demo file 'meta_simulation.in' is being generated...")
    generate_meta_simulationin()
    exit()
  elif (key == 'disk'):
    print("A demo file 'disk.in' is being generated...")
    generate_diskin()
    exit()
  elif (key == 'help'):
    isProblem = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

scriptFolder = os.path.dirname(os.path.realpath(__file__)) # the folder in which the module is. 
binaryPath = os.path.join(scriptFolder, os.path.pardir)

# We try to read parameters from the file
readParameterFile("meta_simulation.in")

# We list all the folders. If the folder name contain the prefix for simulation name, we add the conversion in integer of all 
# that is after the FOLDER_PREFIX in the list. This list will then contain all the indexes of the simulations. By searching the 
# maximum, we will be able to know what index we need for the next simulation.
simus = [dir for dir in os.listdir(".") if (os.path.isdir(dir) and dir.count(FOLDER_PREFIX))]

index_simu = []
for dir in simus:
  try:
    index_simu.append(int(dir.strip(FOLDER_PREFIX)))
  except:
    pass

# Just in case there is no simulation at all, we add '0', then the first simulation will be 0+1 = 1
if (len(index_simu) > 0):
  starting_index = max(index_simu) + 1
else:
  starting_index = 1

for index_simu in range(starting_index, starting_index+NB_SIMULATIONS):
  folder_name = "%s%05i" % (FOLDER_PREFIX, index_simu)
  
  if not(os.path.exists(folder_name)):
    os.mkdir(folder_name)

  os.chdir(folder_name)
  
  generation_simulation_parameters()
  
  mercury_utilities.prepareSubmission(BinaryPath=binaryPath, walltime=WALLTIME)
  
  # We launch the job
  if toLaunch:
    print("We launch the job in "+folder_name)
    job = subprocess.Popen("./runjob", shell=True)
    returncode = job.wait()
  
  
  # We get back in the parent directory
  os.chdir("..")


