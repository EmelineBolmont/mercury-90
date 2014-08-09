#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""module mercury qui permet de lancer des simulations via python"""
# Classe qui définit un environnement permettant de lancer et configurer une simulation mercury. Dans la pratique, il faut définir un script python pour une meta-simulation. Vous pouvez ensuite lancer ce script autant de fois que vous voulez, il va chercher les dossiers existants et en créer un en suivant pour la simulation suivante. Ainsi, on peut lancer plusieurs fois le script à la suite dans la queue du serveur (au hasard venus) pour en faire plusieurs en parallèle.
__author__ = "Autiwa <autiwa@gmail.com>"
__date__ = "2010-09-05"
__version__ = "1.0"

#########################################@
###### module et classe obsolète, je ne m'en sers plus. J'utilise mercury-meta-simulation.py à la place.
##########################################

# RMQ : J'ai modifié pour que print soit sous la forme print() en prévision de la v3 de python

import os   # used to change directories and create folders
import re   # used to determine the name of the next simulation
import subprocess   # Usefull to run the simulation
import time   # For the display of the running time of the simulation and the display of the current time in the logs
import pdb

#import sys
## We import my module. For that aim, we must add the path of its location in the python paths.
## in fact, no need to add because mercury.py and autiwa.py must be in the same directory. 
# sys.path.append(LOCATION_MODULES)
from autiwa import AutiwaObject, Temps  # We only import what interests us.

## CONSTANTS
# (chemin absolu) dossier où sont situées les programmes, notamment mercury6, element6 et close6
LOCATION_PRGM="/home/cossou/bin/mercury6"

# (chemin absolu) dossier où sont placés les meta-simu
LOCATION_METASIMU="/home/cossou/metasimu"

# dossier où seront stockées toutes les données des simulations (en particulier les .aei), sous la forme de sous-dossier pour chaque simulation appartenant à une meta-simulation donnée.
LOCATION_DATASIMU="/sse/cossou/Simulations"


# Pour des tests à la maison
LOCATION_PRGM="/home/autiwa/documents/travail/bin/mercury6"
LOCATION_METASIMU="/home/autiwa/documents/travail/Tests/meta_simu"
LOCATION_DATASIMU="/home/autiwa/documents/travail/Tests/data_meta_simu"


class Simulation(AutiwaObject):
  """Class that define a mercury simulation. It allows us to create the input files (param.in, element.in, big.in) and to run the 
  simulation when everything is defined. Each instance of that class is associated to a meta-simulation, 
  that is to say a group of simulation in wich we will run several simulation with the same specifications (to make statistics)

  Methods
  write_element : écrit element.in
  write_big : écrit big.in
  write_param : écrit param.in
   
  self.runing : (True, False) set to False if no instance is currently runing, but set to True while a process is runing into the directory. 
  system : an object of the class PlanetarySystem that contains all the orbital elements necessary for writing big.in
   
  Error : if two instance of the class use the same directory, there's nothing to test if one of the two is currently runing. 
  Each instance must be associated to a folder. Each time we run a simulation, we cut the output files in an output-directory which 
  can be a sub-directory of the current directory or a defined directory in a different location.
   
  à compléter
  """
  
  #FICHIERS_DEPLACER = ['big.dmp', 'small.dmp', 
            #'param.dmp', 'restart.dmp', 'big.tmp', 'small.tmp', 'param.tmp', 
            #'restart.tmp', 'param.in', 'big.in', 'element.in', 'xv.out', 
            #'ce.out', 'info.out', "mercury_stdout.txt", "mercury_stderr.txt"]
  #FICHIERS_EFFACER = ['big.dmp', 'small.dmp', 'param.dmp', 'restart.dmp', 
            #'big.tmp', 'small.tmp', 'param.tmp', 'restart.tmp', 'xv.out', 
            #'ce.out', 'info.out']
  
  # name of the sub folder of each meta simulation where will be stored every simulation runned under this meta-simulation
  FOLDER_SIMULATIONS = "simulations"
  
  #name for the file where will be stored the stdout output of mercury
  MERCURY_STDOUT = "mercury_stdout.txt"
  #name of the file where will be stored the stderr output of mercury
  MERCURY_STDERR = "mercury_stderr.txt"
  
  # name of the log file, stored in the meta-simu folder, where will be stored each simulation that is runned under this meta-simulation
  LOG_NAME = "simulations.log"
  
  # name of the log file, stored in the simulation folder, where is stored each important action regards to the simulation.
  LOG_RUNNING = "running.log"
  
  # Files that will be removed from the directory at the end of the simulation.
  REMOVE_FILES = ["files.in", "message.in", "small.in"]
  
  # list of files, without extension, that we exchance, from .tmp to .dmp in case of the .dmp are corrupted. (explained in the section (6) of the mercury manual
  LIST_EXCHANGE = ["big", "small", "element", "param", "restart"]
  
  
  # The number of digits we want to have in the numerotation of the folders for simulations in a meta_simulation folder (the 9th simulation will appear as '00009' simulation if NB_DIGITS is set to 5 for example
  NB_DIGITS = 5
   
  def __init__(self, element, param, system, meta_simu="test"):
    """
    On initialise la simulation. Pour cela on génère des 
    fichiers par défaut
    
    Parameters :
    element : an object of type 'Element'
    param : an object of type 'Param'
    system : an object whose type is 'PlanetarySystem'
    meta_simu="test" : the name of the meta-simu. It will be the name of the directory. In this directory will be stored all the simulation that have the parameters of this meta-simulation (i.e all the simulation that we will run with one particuliar instance of the class. BTW several class can have the same directory. Nothing in the class forbid to define different meta-simulation and run them in the same directory. The user MUST be aware of that and be cautious as a consequence.)
    """
    
    AutiwaObject.__init__(self)
    
    self.meta_simu = meta_simu
    
    # We change the directory to be in the meta simu folder
    prevDir = self.__changeDirectory(LOCATION_DATASIMU)
    
    self.param = param
    self.element = element
    self.system = system
    
    
    # we uniformize theses parameters that exists both in param.in and element.in. That's why it doesn't matter what values are set in the Element instance
    self.element.set_output_interval(self.param.get_output_interval())
    self.element.set_relative_time(self.param.get_relative_time())
    self.element.set_time_format(self.param.get_time_format())
    
    # We test and create the folders if necessary
    if not(os.path.exists(self.meta_simu)):
      os.mkdir(self.meta_simu)
    
    os.chdir(self.meta_simu)
    
    if not(os.path.exists(Simulation.FOLDER_SIMULATIONS)):
      os.mkdir(Simulation.FOLDER_SIMULATIONS)
    
    # We create the parameter files that will be the same for all the simulation in the parent directory, that is to say, the meta_simu folder.

    
    # We return in the previous directory
    self.__changeDirectory(prevDir)
    
  def __changeDirectory(self, directory):
    """method that change the directory to another one.
    
    Parameters :
    directory : the path, either relative or absolute where we want to be
    
    return : the path for the actual directory before the chdir
    """
    
    precedent_directory = os.getcwd()
    
    os.chdir(directory)
    
    return precedent_directory

  def __getFolderName(self):
    """method that return a string that represent the name of the simulation to run. 
    
    This function will search for all the names that match the regular expression we have defined for the names of the simulations and determine the next name to apply.
    
    Return : the name of the folder to create for the next simulation as a string
    """
    regular_expression = r"[0-9]+_"+self.param.algorithme
    
    dirs = [dir for dir in os.listdir(".") if os.path.isdir(dir) and re.search(regular_expression, dir)]

    numbers = []
    for dir in dirs:
      # We erase successively all the '0' to the left of the expression and take the left side of the expression, regardless of the character "_"
      numbers.append(eval(dir.lstrip("0").split("_")[0]))
    
    # We get the max number and increase its value by one in order to create a new folder
    if (numbers):
      number = max(numbers)+1
    else:
      number = 1

    folder_name = str(number).zfill(Simulation.NB_DIGITS)+"_"+self.param.algorithme

    return folder_name

  def __writeLog(self, *texts):
    """method that writes a log in a file whose name is determinated by LOG_NAME. This file will be stored in the meta_simu folder. When we write something, we must be able to identify to which simulation the line correspond. (give the name for example)
    
    Parameters
    texts : a sequence of elements that will be converted to strings. 
    
    Exemple :
    self.__writeLog("mercury has returned", integer, "and nothing went wrong")
    with, for instance, integer = 2
    
    return : nothing for the moment"""
    
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu))
    
    file_log = open(Simulation.LOG_NAME, 'a')
    temp_string = "["+str(time.strftime('%d/%m/%Y %H:%M:%S'))+"]"
    for object in texts:
      temp_string += " "+str(object)
    file_log.write(temp_string+"\n")
    file_log.close()
    
    self.__changeDirectory(prevDir)

  def __writeRunning(self, *texts):
    """method that writes a log in a file whose name is determinated by LOG_NAME. This file will be stored in the simulation folder.
    
    Parameters
    texts : a sequence of elements that will be converted to strings. 
    
    Exemple :
    self.__writeRunning("mercury has returned", integer, "and nothing went wrong")
    with, for instance, integer = 2
    
    return : nothing for the moment"""
    
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, self.folder_simulation))
    
    file_log = open(Simulation.LOG_RUNNING, 'a')
    temp_string = "["+str(time.strftime('%d/%m/%Y %H:%M:%S'))+"]"
    for object in texts:
      temp_string += " "+str(object)
    file_log.write(temp_string+"\n")
    file_log.close()
    
    print(temp_string)
    
    self.__changeDirectory(prevDir)
  
  def __extraSimu(self):
    """method that execute functions and things at the end of the simulation. Like sending en e-mail to advise the customer that the simulation has ended"""
    
    self.__writeRunning("Warning: The method __extraSimu() is currently not implemented. At first, you may add a function to send an e-mail at the end of the simulation")
    
    pass
  
  def __writeElement(self, element):
    """génère le fichier "element.in" pour que l'on ait les variables souhaitées ainsi que l'intervalle entre deux outputs voulus
    
    Parameters
    element : an object of the type 'Element'
    
    Return : doesn't return anything, but write the file element.in
    """
    # Il semble qu'en notation exponentielle, le chiffre soit le nombre 
    # total de caractères, incluant les signes '-', le "E" de l'exposant et cie.
    
    if (type(element) != Element):
      self.__writeRunning("'element' must be an objet of type 'Element'")
      raise TypeError("'element' must be an objet of type 'Element'")
    
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, self.folder_simulation))

    element.write()
    
    self.__changeDirectory(prevDir)
  
  def __writeParam(self, param):
    """Génère le fichier param.in avec les paramètres passés
    
    Parameters
    param : an object of the type 'Param'
    
    Return : doesn't return anything, but write the file param.in
    """
    # Version 1.1
    # Quand on définit un texte, toujours faire un retour chariot à la 
    # dernière ligne. mercury est très pointilleux sur le nombre de lignes du fichiers.
    
    if (type(param) != Param):
      self.__writeRunning("'param' must be an objet of type 'Param'")
      raise TypeError("'param' must be an objet of type 'Param'")
    
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, self.folder_simulation))

    param.write()
    
    self.__changeDirectory(prevDir)
  
  def __writeBig(self, big):
    """ Génère le fichier big.in à partir des paramètres d'entrée. 
    
    Parameters
    big : an object of the class Big that represent the file big.in
    
    Return : doesn't return anything, but write the file big.in
    """
    # Version 0.4
    
    
    if (type(big) != Big):
      self.__writeRunning("'big' must be an objet of type 'Big'")
      raise TypeError("'big' must be an objet of type 'Big'")
    
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, self.folder_simulation))
    
    # we define an 'Big' object that we write directly on a file
    big.write()
    
    
    self.__changeDirectory(prevDir)
    
  def __writeParameterFiles(self):
    """method that writes all the needed files in the current directory"""
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, self.folder_simulation))
    
    parameter_files.write_filesin()
    parameter_files.write_messagein()
    Small().write()
    
    self.__changeDirectory(prevDir)
  
  def __removeParameterFiles(self):
    """method that remove all the needed files to launch the simulation, once the simulation had terminated"""
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, self.folder_simulation))
    
    for file in Simulation.REMOVE_FILES:
      if os.path.exists(file):
        os.remove(file)
    
    self.__changeDirectory(prevDir)
  
  def startSimu(self):
    """method that start a new simulation of the current metasimulation
    
    This method will search all the simulation already finished or in progress for this meta-simulation and find the closest name of directory that match the norm of name we have fixed
    
    We are in the path of the meta-simulation in the partition where we run and store the datas. We MUST still change the directory to the sub-directory for the simulation that we will run. 
    """
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS))
    
    ## On crée un sous répertoire pour la simulation qu'on va lancer
    self.folder_simulation = self.__getFolderName()
    self.__writeLog(self.folder_simulation, ": We create the folder")
    os.mkdir(self.folder_simulation)  # faut créer le dossier maintenant. On est dans le bon dossier normalement
    
    self.__writeParam(self.param)
    self.__writeElement(self.element)
    
    # In case the system being partially random, we re-generate its values (in order to have different system if we launch several simulations with the same script. If the system is not random, then nothing will change.
    self.system.generateOrbitals()
    self.__writeBig(Big(self.system))
    self.__writeParameterFiles()
    
    os.chdir(self.folder_simulation)
    self.__writeLog(self.folder_simulation, ": We launch the simulation")
    self.__writeRunning("We launch the simulation")
    temps_debut = time.time()
    process = subprocess.Popen(LOCATION_PRGM+"/mercury6", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (instance_sortie, instance_erreur) = process.communicate()
    returnCode = process.poll()
    temps_fin = time.time()
    temps_exec = Temps(temps_fin - temps_debut)
    
    if (returnCode!=0):
      self.__writeRunning("mercury has returned an error "+str(returnCode))
    else:
      self.__writeRunning("mercury has been successfully executed (in "+str(temps_exec)+").")
    
    self.__writeRunning("we write stdout/stderr")
    # On écrit un fichier dans lequel on stocke les valeurs données 
    # par mercury. Le "f" noté devant est pour file, c'est l'objet correspondant au nom de fichier.
    fmercury_stdout = open(Simulation.MERCURY_STDOUT, 'w')
    for ligne in instance_sortie:
      fmercury_stdout.write(ligne)
    fmercury_stdout.close()
    
    fmercury_stderr = open(Simulation.MERCURY_STDERR, 'w')
    for ligne in instance_erreur:
      fmercury_stderr.write(ligne)
    fmercury_stderr.close()
    
    
    
    self.__writeRunning("we write the .aei")
    self.generateOutputFiles(self.folder_simulation)
    
    # we clean the directory
    self.__removeParameterFiles()
    
    
    self.__writeLog(self.folder_simulation, ": The simulation has terminated in", temps_exec)
    self.__extraSimu()
    self.__changeDirectory(prevDir)
  
  def restartSimu(self, folder):
    """method that allow to continue an existing simulation that crashed via .dmp
    
    folder : the name of the folder of the simulation we want to re-launch (to continue the integration where it crashed)
    """
    
    self.folder_simulation = folder
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, self.folder_simulation))
    
    self.__writeParameterFiles()
    
    self.__writeLog(self.folder_simulation, ": We re-launch the simulation")
    self.__writeRunning("We re-launch the simulation")
    temps_debut = time.time()
    process = subprocess.Popen(LOCATION_PRGM+"/mercury6", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (instance_sortie, instance_erreur) = process.communicate()
    returnCode = process.poll()
    temps_fin = time.time()
    temps_exec = Temps(temps_fin - temps_debut)
    
    if (returnCode!=0):
      self.__writeRunning("mercury has returned an error "+str(returnCode))
      bool_exchange_dmp_tmp = True
    else:
      self.__writeRunning("mercury has been successfully executed (in "+str(temps_exec)+").")
      bool_exchange_dmp_tmp = False
    
    if bool_exchange_dmp_tmp:
      self.__writeRunning("We copy *.tmp to *.dmp in case the .dmp's are corrupted")
      for file in Simulation.LIST_EXCHANGE:
        if os.path.exists(file):
          os.remove(file+".dmp")
      
      for file in Simulation.LIST_EXCHANGE:
        process = subprocess.Popen("cp "+file+".tmp "+file+".dmp", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        (instance_sortie, instance_erreur) = process.communicate()
        returnCode = process.poll()
        
        if (returnCode!=0):
          self.__writeRunning("cp has returned an error "+str(returnCode)+" while copying "+file+".tmp to "+file+".dmp")
          self.__writeRunning(instance_erreur)
      
      self.__writeRunning("We re-launch mercury")
      temps_debut = time.time()
      process = subprocess.Popen(LOCATION_PRGM+"/mercury6", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      (instance_sortie, instance_erreur) = process.communicate()
      returnCode = process.poll()
      temps_fin = time.time()
      temps_exec = Temps(temps_fin - temps_debut)
      
      if (returnCode!=0):
        self.__writeRunning("mercury has returned an error "+str(returnCode))
      else:
        self.__writeRunning("mercury has been successfully executed (in "+str(temps_exec)+").")
    
    self.__writeRunning("we write stdout/stderr")
    # On écrit un fichier dans lequel on stocke les valeurs données 
    # par mercury. Le "f" noté devant est pour file, c'est l'objet correspondant au nom de fichier.
    fmercury_stdout = open(Simulation.MERCURY_STDOUT, 'a')
    for ligne in instance_sortie:
      fmercury_stdout.write(ligne)
    fmercury_stdout.close()
    
    fmercury_stderr = open(Simulation.MERCURY_STDERR, 'a')
    for ligne in instance_erreur:
      fmercury_stderr.write(ligne)
    fmercury_stderr.close()
    
    
    self.__writeRunning("we write the .aei")
    self.generateOutputFiles(self.folder_simulation)
    
    # we clean the directory
    self.__removeParameterFiles()
    
    
    self.__writeLog(self.folder_simulation, ": The simulation has terminated in", temps_exec)
    self.__extraSimu()
    self.__changeDirectory(prevDir)
  
  def extendSimu(self, folder, suptime):
    """method that re-run the simulation for a defined time in years
    
    Parameters
    folder : the name of the folder of the simulation we want to extend
    suptime : the amount of time, in years, that we want to extend the simulation
    
    Do : will copy all the files of folder in a new one, suffixed by "_ext_by_suptime" then launch the simulation
    """
    self.folder_simulation = folder
    
    # We go in the folder of the simulation that interest us
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, self.folder_simulation))
    # We get the existing param.in
    paramin = readParam("param.in")
    
    # We modifie the stop_time
    paramin.set_stop_time(paramin.get_stop_time() + suptime)
    
    
    folder_extended_simulation = self.folder_simulation+"_ext_by_"+str(suptime)
    
    # We go back in the folder where all the simulations are contained, in order to create the new folder and copy all the data
    self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS))
    
    
    process = subprocess.Popen("cp -R "+self.folder_simulation+" "+folder_extended_simulation, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (instance_sortie, instance_erreur) = process.communicate()
    returnCode = process.poll()
    
    if (returnCode!=0):
      self.__writeRunning("cp has returned an error "+str(returnCode))
      self.__writeRunning(instance_erreur)
    
    # We go in the new directory to write the new param.in parameter file.
    self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, folder_extended_simulation))
    
    # We write the file in this new directory (also in 'param.tmp' because if dmp are corrupted, we can erase the .dmp and rename the .tmp in .dmp)
    paramin.write("param.dmp")
    paramin.write("param.tmp")
    
    self.__changeDirectory(prevDir)
    #prépare le dossier, modifie param.in (en le lisant?)
    
    self.restartSimu(folder_extended_simulation)
  
  def unitaryTests(self):
    """method that allow me to test all the private methods of the class"""
    
    print("unitaryTests is not implemented for the moment")
  
  def generateOutputFiles(self, folder):
    """method that generate outputs files in the current working directory. Before generating the .aei files, the method remove all the .aei files in the current directory, if they exists.
    
    Parameters :
    folder : the name of a simulation in the meta-simu directory defined for the instance
    """
    
    prevDir = self.__changeDirectory(os.path.join(LOCATION_DATASIMU, self.meta_simu, Simulation.FOLDER_SIMULATIONS, folder))
    
    # We list the .aei files that may exists in the directory
    process = subprocess.Popen("ls *.aei", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (instance_sortie, instance_erreur) = process.communicate()
    returnCode = process.poll()
    
    if (returnCode == 0):
      list_aei = instance_sortie.split("\n")
      for file in list_aei:
        if os.path.exists(file):
          os.remove(file)
    
    temps_debut = time.time()
    
    process = subprocess.Popen(LOCATION_PRGM+"/element6", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (instance_sortie, instance_erreur) = process.communicate()
    returnCode = process.poll()

    temps_fin = time.time()
    temps_exec = Temps(temps_fin - temps_debut)
    
    if (returnCode!=0):
      self.__writeRunning("element6 has returned an error "+str(returnCode))
      self.__writeRunning(instance_erreur)
    else:
      self.__writeRunning("element6 has been successfully executed (in "+str(temps_exec)+").")
