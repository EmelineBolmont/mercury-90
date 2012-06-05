#!/usr/bin/env python
# script to test various pieces of python code
# Version 1.6

from mercury_utilities import * # module that contain utilities to help create a planetary system. 
import simulations_utilities    # module that contain utilities to help launch simulations
#~ import mercury                  # module that contain object for each parameter file (and also for bodies and planetary system)
import pdb                      # to debug via pdb.set_trace()
import os                       # to create folder, change directory and so on
import subprocess               # to launch 'runjob'
from constants import *         # Several constants, including mass of planets and so on. All constants are in CAPSLOCK
import random
import shutil                   # In particular to copy a file with .copy2()

FOLDER_PREFIX = "simu" # the prefix for each sub simulation folder
SUB_FOLDER_LOG = "random_parameters.in" # the name of the log file where we will store meta simulation information to keep a trace.

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

# name of the queue where we want to lauch simulations (either syntax is correct : "arguin*.q", "arguin1.q" or "arguin1.q,arguin2.q")
QUEUE = "arguin1.q,arguin2.q"

# The absolute path to the folder were are the binary files of mercury
BINARY_FOLDER = "/home/cossou/bin/mercury"

#----------------------------------
# Parameters for mercury

epoch = 0
time_format = "years"
relative_time = "yes"
nb_outputs = 1000
user_force = "yes"
# integration_time = 1e6 * YEARS
timestep = 8 # days
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

#----------------------------------
# Parameters for the disk
# None parameters will not be displayed in the parameter file and thus, default values of the code will be used.
sigma_0 = 1700 # g/cm^2
adiabatic_index = 1.4 # the adiabatic index of the disk
viscosity = 1.e15 # cm^2/s
b_h = 0.4 # the smoothing length of the gravitationnal potential of the planet
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
turbulent_forcing = None

#-------------------------------------------------------------------------------
# Definition of functions

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
	
	global NB_SIMULATIONS, QUEUE, BINARY_FOLDER, isErase
	global epoch, integration_time, time_format, relative_time, nb_outputs, user_force
	global FIXED_TOTAL_MASS, TOTAL_MASS, NB_PLANETS, mass_parameters, a_parameters, e_parameters, I_parameters
	global sigma_0, adiabatic_index, viscosity, b_h, torque_type
	global dissipation_type, disk_exponential_decay, sample, inner_boundary_condition, outer_boundary_condition
	global torque_profile_steepness, indep_cz, mass_dep_m_min, mass_dep_m_max, mass_dep_cz_m_min, mass_dep_cz_m_max
	global turbulent_forcing
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
			elif (key in ["queue", "QUEUE"]):
				QUEUE = simulations_utilities.str2str(value)
			elif (key in ["binary_folder", "BINARY_FOLDER"]):
				BINARY_FOLDER = simulations_utilities.str2str(value)
			#----------------------------------
			# Mercury parameters
			elif (key == "epoch"):
				epoch = float(value) * YEAR
			elif (key == "integration_time"):
				integration_time = float(value) * YEAR
			elif (key == "time_format"):
				time_format = simulations_utilities.str2str(value)
			elif (key == "relative_time"):
				relative_time = simulations_utilities.str2str(value)
			elif (key == "nb_outputs"):
				nb_outputs = int(value)
			elif (key == "user_force"):
				user_force = simulations_utilities.str2str(value)
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
			#----------------------------------
			# Disk parameters
			elif (key == "sigma_0"):
				sigma_0 = float(value)
			elif (key == "adiabatic_index"):
				adiabatic_index = float(value)
			elif (key == "viscosity"):
				viscosity = float(value)
			elif (key == "turbulent_forcing"):
				viscosity = float(value)
			elif (key == "b/h"):
				b_h = eval(value)
			elif (key == "dissipation_type"):
				dissipation_type = int(value)
			elif (key == "disk_exponential_decay"):
				disk_exponential_decay = float(value)
			elif (key == "sample"):
				sample = int(value)
			elif (key == "inner_boundary_condition"):
				inner_boundary_condition = value
			elif (key == "outer_boundary_condition"):
				outer_boundary_condition = value
			elif (key == "torque_type"):
				torque_type = value
			elif (key == "torque_profile_steepness"):
				torque_profile_steepness = float(value)
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
				PARAMETERS += "the parameter "+key+" is not known. His value was : "+value+"\n"
				print("Warning: no parameter is known for the key '"+key+"'")
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
	PARAMETERS += "epoch = "+str(epoch/365.25)+" years\n"
	PARAMETERS += "integration time = "+str(integration_time/365.25)+" years\n"
	PARAMETERS += "number of outputs = "+str(nb_outputs)+"\n"
	PARAMETERS += "user force = "+str(user_force)+"\n"
	PARAMETERS += "----------------------------------\nPlanetary System Parameters\n\n"
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
	PARAMETERS += "sigma_0 = "+str(sigma_0)+" g/cm^2\n"
	PARAMETERS += "sample = "+str(sample)+"\n"
	PARAMETERS += "adiabatic_index = "+str(adiabatic_index)+"\n"
	PARAMETERS += "viscosity = "+str(viscosity)+" cm^2/s\n"
	PARAMETERS += "b/h = "+str(b_h)+"\n"
	PARAMETERS += "dissipation_type = "+str(dissipation_type)+"\n"
	if (dissipation_type == 2):
		PARAMETERS += "disk_exponential_decay = "+str(disk_exponential_decay)+"\n"
	PARAMETERS += "inner_boundary_condition = "+str(inner_boundary_condition)+"\n"
	PARAMETERS += "outer_boundary_condition = "+str(outer_boundary_condition)+"\n"
	PARAMETERS += "torque_type = "+str(torque_type)+"\n"
	PARAMETERS += "torque_profile_steepness = "+str(torque_profile_steepness)+"\n"
	PARAMETERS += "indep_cz = "+str(indep_cz )+"\n"
	PARAMETERS += "mass_dep_m_min = "+str(mass_dep_m_min )+"\n"
	PARAMETERS += "mass_dep_m_max = "+str(mass_dep_m_max )+"\n"
	PARAMETERS += "mass_dep_cz_m_min = "+str(mass_dep_cz_m_min )+"\n"
	PARAMETERS += "mass_dep_cz_m_max = "+str(mass_dep_cz_m_max )+"\n"




def generation_simulation_parameters():
	"""the function generate simulations files in the current working directory, given the parameters on top of the script
	"""
	
	output_interval = integration_time / float(nb_outputs)
	data_dump = abs(int(integration_time / (nb_outputs * timestep)))

	# We begin the random generation of parameters.
	if FIXED_TOTAL_MASS:
		total_mass = 0
		m = []
		while (total_mass < TOTAL_MASS):
			m0 = simulations_utilities.setParameter(mass_parameters, 1)[0] # The function return a list, and we want to have one element
			# the mass must be expressed relatively to the mass of the central body
			m.append(simulations_utilities.significativeRound(m0 * MT / MS, 4))
			total_mass += m0
		nb_planets = len(m)
	else:
		nb_planets = NB_PLANETS
		
		m = simulations_utilities.setParameter(mass_parameters, nb_planets)
		m = [mi * MT / MS for mi in m]

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
	else:
		a = simulations_utilities.setParameter(a_parameters, nb_planets)
		
	e = simulations_utilities.setParameter(e_parameters, nb_planets)
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

	system = definePlanetarySystem(m=m, a=a, e=e, I=I)

	# We write the files

	bigin = mercury.Big(system)
	bigin.write()

	smallin = mercury.Small(system)
	smallin.write()
	if not(os.path.exists("small.in")):
		pdb.set_trace()

	# Setting the output interval to 0 ensure that we will have every output written in the xv.out in the element files.
	elementin = mercury.Element(format_sortie=" a8.5 e8.6 i8.4 g8.4 n8.4 l8.4 m13e ", coord="Cen", 
	output_interval=0.0, time_format=time_format, relative_time=relative_time)
	elementin.write()

	closein = mercury.Close(time_format=time_format, relative_time=relative_time)
	closein.write()

	paramin = mercury.Param(algorithme="HYBRID", start_time=epoch, stop_time=epoch+integration_time, output_interval=output_interval, 
	h=timestep, accuracy=1.e-12, stop_integration="no", collisions="yes", fragmentation="no", 
	time_format=time_format, relative_time=relative_time, output_precision="medium", relativity="no", 
	user_force=user_force, ejection_distance=EJECTION_DISTANCE, radius_star=0.005, central_mass=1.0, 
	J2=0, J4=0, J6=0, changeover=3., data_dump=data_dump, periodic_effect=100)
	paramin.write()

	filesin = mercury.Files()
	filesin.write()

	messagein = mercury.Message()
	messagein.write()

	if (user_force == "yes"):
		diskin = mercury.Disk(b_over_h=b_h, adiabatic_index=adiabatic_index, mean_molecular_weight=2.35, surface_density=(sigma_0, 0.5), 
		              disk_edges=(1., 100.), viscosity=viscosity, sample=sample, dissipation_type=dissipation_type, 
		              turbulent_forcing=turbulent_forcing, 
		              disk_exponential_decay=disk_exponential_decay, torque_type=torque_type,
		              inner_boundary_condition=inner_boundary_condition, outer_boundary_condition=outer_boundary_condition,
	                torque_profile_steepness=torque_profile_steepness, indep_cz=indep_cz, mass_dep_m_min=mass_dep_m_min, 
	                mass_dep_m_max=mass_dep_m_max, mass_dep_cz_m_min=mass_dep_cz_m_min, mass_dep_cz_m_max=mass_dep_cz_m_max)
		diskin.write()
		
		# If we want to use a manual torque profile, we copy the torque profile from the current working directory
		if (torque_type == 'manual'):
			if (os.path.isfile("../"+torque_file)):
				shutil.copy2("../"+torque_file, "./"+torque_file)
			else:
				raise NameError("The file "+torque_file+" does not exist in the parent directory\n keep in mind to delete the created folder.")
	
	# We store a log file of the simulation parameters
	f = open(SUB_FOLDER_LOG, 'w')
	f.write(PARAMETERS)
	f.close()
	
	# We reset the counters for planet names
	mercury.Body.resetCounter()

#-------------------------------------------------------------------------------
# Beginning of the code

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
	folder_name = FOLDER_PREFIX+simulations_utilities.number_fill(index_simu, 5)
	
	if not(os.path.exists(folder_name)):
		os.mkdir(folder_name)
	#~ else:
		#~ if isErase:
			#~ process = subprocess.Popen("rm -f "+folder_name+"/*", shell=True)
			#~ (process_stdout, process_stderr) = process.communicate()
			#~ returnCode = process.wait()
			#~ if (process_stdout != None):
				#~ print(process_stdout)
			#~ if (process_stderr != None):
				#~ print(process_stderr)
			#~ if (returnCode != 0):
				#~ print("Warning: Unable to delete previous files")
				#~ pdb.set_trace()
		#~ else:
			#~ print("Warning: a simulation already exists")
	os.chdir(folder_name)
	
	generation_simulation_parameters()
	
	command = BINARY_FOLDER+"/mercury\n" + \
	          BINARY_FOLDER+"/element\n" + \
	          "echo `pwd` ': Done'>>~/qsub.log\n"
	script = simulations_utilities.SimpleJob(command) # For arguin
	#~ script = simulations_utilities.Job_PBS(command, walltime=48) # For avakas
	script.write()
	
	
	simulations_utilities.setExecutionRight("simulation.sh")
	
	# We define a bash script to launch the simulation in a queue
	simulations_utilities.writeRunjobSGE("simulation.sh", QUEUE) # For arguin
	#~ simulations_utilities.writeRunjobPBS("simulation.sh") # For avakas
	
	# We launch the job
	print("We launch the job in "+folder_name)
	job = subprocess.Popen("./runjob", shell=True)
	returncode = job.wait()
	
	
	# We get back in the parent directory
	os.chdir("..")


