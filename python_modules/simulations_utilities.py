#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that generate parameter files for genesis simulation code.
# In particuliar, it helps generate random planetary systems.
from __future__ import print_function

__author__ = "Autiwa <autiwa@gmail.com>"
__date__ = "24 août 2011"
__version__ = "$Revision: 0.3 $"
__credits__ = """Module that contains some utilities to help generate simulations, no matter the code. Function to generate random parameters
function to round values with significative round and so on."""

from random import uniform, gauss
from constants import *
from math import *
import autiwa
import pdb
import subprocess
import os # at least to have access to the os.path.isfile() method.

def setExecutionRight(doc_name):
  """function that set the right for the file to be executed. This file must be in the current working directory

  Parameters
  doc_name = the name of the file

  Return : the return code of the chmod command. If 0, then everything went good.
  """

  import subprocess

  command = "chmod u+x "+doc_name

  process = subprocess.Popen(command, shell=True)
  returncode = process.wait()

  return returncode

def str2bool(value):
  """function that convert a string into a boolean given various ways to say 'True'. 
  The normal behaviour in python is to set to 
  True when the string is not empty, which is not what I want.
  
  Return : True (boolean) if value is 'True', 'yes', 1 or '.true.' 
  (no matter the case, i.e True, TRUE, true or tRuE will work)"""
  return value.lower() in ("yes", "true", "1", '.true.')
  
def str2str(value):
  """function that convert a string into a string. That might be weird to say, 
  but this function aim to suppress any '"' or "'" in string
  
  Return a string without external quotes that might exist
  """
  
  # We suppress extremal white spaces
  value = value.strip()
  
  value = value.strip('"').strip("'")
  
  return value


def writeRunjobSGE(command, queue="", nb_proc=1):
  """function that creates a script named 'runjob' that
  will run a job on a queue. If the number of processor exceed 1, then
   the function will try to launch the job on every queue. If not, it
   will launch the job on all the parallel queues.
   
  This script is done for the SGE (Sun Grid Engine) environment.

  Parameters
  nb_proc=1 : (integer) number of processor we want to use. By default, it will be 1
  queue : the queue you want to use to launch your job. You can use the various syntaxes allowed by the job scheduler. 
  command : The command you want the job to launch. 
  
  Example : 
  writeRunjob("./mercury", "arguin1.q,arguin2.q")
  will create a "runjob" script that will look like : 
  qsub -q arguin1.q,arguin2.q ./mercury

  Return : Nothing
  """

  NAME_SCRIPT = "runjob"

  
  
  if (queue == ""):
    queue_append = ""
  else:
    queue_append = " -q "+queue
  
  
  if (nb_proc > 1):
    qsub = "qsub -pe mpi "+str(nb_proc)+queue_append+" "+command
  else:
    qsub = "qsub"+queue_append+" "+command
  
  script = open(NAME_SCRIPT, 'w')
  script.write("stdout=$("+qsub+")\n")
  script.write("echo $stdout\n")
  script.write("echo `date '+%d-%m-%Y at %H:%M:%S'` `pwd` ':' $stdout>>~/qsub.log\n")
  script.close()

  setExecutionRight(NAME_SCRIPT)
  

def writeRunjobPBS(command):
  """function that creates a script named 'runjob' that
  will run a job on a queue. If the number of processor exceed 1, then
   the function will try to launch the job on every queue. If not, it
   will launch the job on all the parallel queues.
   
  This script is done for the Torque+Maui task scheduler, alias PBS. For the moment, the command is really simple because all the
  parameters

  Parameters
  command : The command you want the job to launch. 
  
  Example : 
  writeRunjob("./mercury")
  will create a "runjob" script that will look like : 
  qsub ./mercury

  Return : Nothing
  """

  NAME_SCRIPT = "runjob"
  
  
  qsub = "qsub "+command
  
  script = open(NAME_SCRIPT, 'w')
  script.write("stdout=$("+qsub+") # execute the command and store the output in '$stdout'\n")
  script.write("echo `date '+%d-%m-%Y at %H:%M:%S'` `pwd` ': launched'>>~/qsub.log\n")
  script.write("echo $stdout # display the output of the qsub\n")
  script.write("echo `date '+%d-%m-%Y at %H:%M:%S'` `pwd` ':' $stdout>>~/qsub.log\n")
  script.close()

  setExecutionRight(NAME_SCRIPT)
  
class Job_PBS(object):
  """class that define an object equivalent to a script needed to run a job on a PBS job engine
  
  For the moment, the class only allows us to define non parallel jobs that need only 1 proc.
  
  For jobs that only needs one proc, we define prologue and epilogue to transfert to the node all 
  the needed files and taking them back at the end of the job. This is done to prevent surcharge 
  of the hard disk. This is supposed to be completely transparent to the user on the express 
  condition that all the needed files for the simulation are in the 
  current working directory or with absolute paths.
  
  Parameters:
  command : the name of the command (or script) to run
  directory="." : the name of the directory where the job will be launched. By default, it will be the current working directory.
  walltime='00:10:00' : the maximum expected length of the job. If it's a number, then this will be the length in hour 
                        (converted in hour minutes seconds so you can put decimal point values). Else, you must specify the walltime 
                         with the following form 'hh:mm:ss'. 
  name="simulation.sh" : the name of the script that will contains the commands to launch the job
  """
  
  
  def __init__(self, command, directory='.', walltime='00:10:00', name="simulation.sh"):
    """initialisation of the class"""
    
    self.command = command
    
    self.directory = directory
    
    self.name = name
    
    self.proc_per_node = 1
    self.nodes = 1
    
    if (self.nodes == 1):
      self.isPrologEpilog = True
    else:
      self.isPrologEpilog = False
    
    if (type(walltime) == str):
      self.walltime = walltime
    elif (type(walltime) in [int, float]):
      rest = walltime
      hour = int(rest)
      
      rest = (rest - hour) * 60.
      minutes = int(rest)
      
      rest = (rest - minutes) * 60.
      seconds = int(rest)
      self.walltime = str(hour)+":"+number_fill(minutes,2)+":"+number_fill(seconds,2)
    else:
      raise TypeError("The argument walltime must be a number of hour or a string of the form 'hh:mm:ss'")
    
    
  def write(self):
    """write all the data in a file named self.name in the current working directory"""
    
    script = open(self.name, 'w')
    script.write("#!/bin/sh\n")
    script.write("\n")
    script.write("#############################\n")
    script.write("\n")
    script.write("# Your job name\n")
    script.write("#PBS -N "+str(self.name)+"\n")
    script.write("\n")
    script.write("# Specify the working directory\n")
    script.write("#PBS -d "+str(self.directory)+"\n")
    script.write("\n")
    # When all the calculations are done in one peculiar node, we create prolog and epilog to limit the use of the /home harddrive.
    if (self.isPrologEpilog):
      script.write("# prologue\n")
      script.write("#PBS -l prologue=./prolog.sh\n")
      script.write("\n")
      script.write("# epilogue\n")
      script.write("#PBS -l epilogue=./epilog.sh\n")
    script.write("# walltime (hh:mm::ss)\n")
    script.write("#PBS -l walltime="+self.walltime+"\n")
    script.write("\n")
    script.write("# Specify the number of nodes(nodes=) and the number of cores per nodes(ppn=) to be used\n")
    script.write("#PBS -l nodes="+str(self.nodes)+":ppn="+str(self.proc_per_node)+"\n")
    script.write("\n")
    script.write("#############################\n")
    script.write("\n")
    script.write("# modules cleaning\n")
    script.write("module purge\n")
    script.write("module add torque\n")
    script.write("module add gcc\n")
    script.write("\n")
    script.write("# useful informations to print\n")
    script.write("echo \"#############################\" \n")
    script.write("echo \"User:\" $USER\n")
    script.write("echo \"Date:\" `date`\n")
    script.write("echo \"Host:\" `hostname`\n")
    script.write("echo \"Directory:\" `pwd`\n")
    script.write("echo \"PBS_JOBID:\" $PBS_JOBID\n")
    script.write("echo \"PBS_O_WORKDIR:\" $PBS_O_WORKDIR\n")
    script.write("echo \"PBS_NODEFILE: \" `cat $PBS_NODEFILE | uniq`\n")
    script.write("echo \"#############################\" \n")
    script.write("\n")
    script.write("#############################\n")
    script.write("\n")
    if (self.isPrologEpilog):
      script.write("# With prologue and epilogue, we must now change the working derectory before running the simulation\n")
      script.write("myTmpRep=/tmp/$USER/job_$PBS_JOBID\n")
      script.write("cd $myTmpRep\n")
      script.write("echo \"pwd:\" `pwd`\n")
    script.write("# The job of the script is launched here\n")
    script.write(self.command)
    script.write("\n")
    script.write("# all done\n")
    script.write("echo `date` \"Job finished\" \n")
    script.close()
    
    if (self.isPrologEpilog):
      prolog = prolog_PBS()
      prolog.write()
      
      epilog = epilog_PBS()
      epilog.write()

class prolog_PBS(object):
  """
  This define an object and then a file 'prolog.sh' that will contains various commands for a
  prolog file in the PBS environment (like avakas).
  
  By default, this file will create a temporary directory on the targeted node, 
  and transfert all the files from the current directory to this tmp directory.
  
  (The same will then have to be done in the opposite way in an epilog script).
  
  Parameters :
  files="*" : By default, we copy all the files in the current working directory, but with this option, 
              we can specify to copy only peculiar files. All the specified files MUST NOT have 'spaces' in their names.
  
  Examples : 
  prolog_PBS(files="*.in")
  or
  prolog_PBS(files="param.in big.in")
  """
  
  def __init__(self, files="*"):
    """initialisation of the class"""
    
    self.files = files
    
    self.name = "prolog.sh"
  
  def __giveRights(self):
    """method to run after creating the file (after the 'write' method then) because to run the job, 
    the prolog and epilog must have some execution rights"""
    
    if os.path.isfile(self.name):
      setExecutionRight(self.name)
    else:
      raise NameError("the file '"+self.name+"' doesn't exist.")
    
  def write(self):
    """write all the data in a file named self.name in the current working directory"""
    
    script = open(self.name, 'w')
    script.write("#!/bin/sh\n")
    script.write("\n")
    script.write("#\n")
    script.write("echo \"-- begin prologue\"\n")
    script.write("\n")
    script.write("# prologue gets arguments:\n")
    script.write("# argv[1]   job id\n")
    script.write("# argv[2]   job execution user name\n")
    script.write("# argv[3]   job execution group name\n")
    script.write("# argv[4]   job name (TORQUE 1.2.0p4 and higher only)\n")
    script.write("# argv[5]   list of requested resource limits (TORQUE 1.2.0p4 and higher only)\n")
    script.write("# argv[6]   job execution queue (TORQUE 1.2.0p4 and higher only)\n")
    script.write("# argv[7]   job account (TORQUE 1.2.0p4 and higher only)jobid=$1\n")
    script.write("jobid=$1\n")
    script.write("user=$2\n")
    script.write("\n")
    script.write("#\n")
    script.write("echo \"User:\" $USER\n")
    script.write("echo \"Host:\" `hostname`\n")
    script.write("echo \"Directory:\" `pwd`\n")
    script.write("echo \"PBS_O_WORKDIR:\" $PBS_O_WORKDIR\n")
    script.write("\n")
    script.write("echo \"jobid:\" $jobid\n")
    script.write("echo \"user:\" $user\n")
    script.write("\n")
    script.write("# \n")
    script.write("myTmpRep=/tmp/$user/job_$jobid\n")
    script.write("mkdir -p $myTmpRep\n")
    for file in self.files.split():
      script.write("cp $PBS_O_WORKDIR/"+file+" $myTmpRep/.\n")
    script.write("\n")
    script.write("#\n")
    script.write("echo \"-- end prologue\"\n")
    script.write("\n")
    script.write("# Successful completion => Job will run\n")
    script.write("exit 0\n")
    script.close()
    
    # We need to give certain rights to the prologue script.
    self.__giveRights()
    
class epilog_PBS(object):
  """
  This define an object and then a file 'epilog.sh' that will contains various commands for an
  epilog file in the PBS environment (like avakas).
  
  By default, this file will copy back all the selected files form a tmp directory on the 
  node back to the working directory when we submitted the job
  
  (The same must have been done in the opposite way at the beginning of the job, in a prologue script).
  
  Parameters :
  files="*" : By default, we copy all the files in the current working directory, but with this option, 
              we can specify to copy only peculiar files. All the specified files MUST NOT have 'spaces' in their names.
  
  Examples : 
  epilog_PBS(files="*.out")
  or
  epilog_PBS(files="param.out big.out")
  """
  
  def __init__(self, files="*"):
    """initialisation of the class"""
    
    self.files = files
    
    self.name = "epilog.sh"
  
  def __giveRights(self):
    """method to run after creating the file (after the 'write' method then) because to run the job, 
    the prolog and epilog must have some execution rights"""
    
    if os.path.isfile(self.name):
      setExecutionRight(self.name)
    else:
      raise NameError("the file '"+self.name+"' doesn't exist.")
    
  def write(self):
    """write all the data in a file named self.name in the current working directory"""
    
    script = open(self.name, 'w')
    script.write("#!/bin/sh\n")
    script.write("\n")
    script.write("#\n")
    script.write("echo \"-- begin epilogue\"\n")
    script.write("\n")
    script.write("# epilogue gets arguments:\n")
    script.write("# argv[1]   job id\n")
    script.write("# argv[2]   job execution user name\n")
    script.write("# argv[3]   job execution group name\n")
    script.write("# argv[4]   job name\n")
    script.write("# argv[5]   session id\n")
    script.write("# argv[6]   list of requested resource limits\n")
    script.write("# argv[7]   list of resources used by job\n")
    script.write("# argv[8]   job execution queue\n")
    script.write("# argv[9]   job account\n")
    script.write("# argv[10]   job exit code\n")
    script.write("jobid=$1\n")
    script.write("user=$2\n")
    script.write("\n")
    script.write("#\n")
    script.write("echo \"User:\" $USER\n")
    script.write("echo \"Host:\" `hostname`\n")
    script.write("echo \"Directory:\" `pwd`\n")
    script.write("echo \"PBS_O_WORKDIR:\" $PBS_O_WORKDIR\n")
    script.write("\n")
    script.write("echo \"jobid:\" $jobid\n")
    script.write("echo \"user:\" $user\n")
    script.write("\n")
    script.write("#\n")
    script.write("myTmpRep=/tmp/$user/job_$jobid\n")
    for file in self.files.split():
      script.write("cp $myTmpRep/"+file+" $PBS_O_WORKDIR/\n")
    script.write("\n")
    script.write("#\n")
    script.write("echo \"-- end epilogue\"\n")
    script.write("\n")
    script.write("# Successful completion\n")
    script.write("exit 0\n")
    script.close()
    
    # We need to give certain rights to the epilogue script
    self.__giveRights()


class SimpleJob(object):
  """Class that allow us to define a script for simple jobs that only need a command line as argument
  
  Parameter:
  command : the commandes you want to launch
  """
  
  def __init__(self, command, name="simulation.sh"):
    """initialisation of the class"""
    
    self.name = name
    
    self.command = str(command)
    
  def write(self):
    """write all the data in a file named self.name in the current working directory"""
    script = open("simulation.sh", 'w')
    script.write("#!/bin/bash\n")
    script.write("\n")
    script.write(self.command)
    script.close()

def setParameter(parameter, nb_planets, vmin=None, vmax=None):
  """This function will generate the parameters list given a tuple a values.. 
  
  Parameter : 
  parameter : If parameter is a number, then it will be assumed that all the planets have the same value
        If parameter is a list or a tuple, then so far, it must have 3 elements. The third element must be either 'uniform' or 'gauss'.
        If it's 'uniform', then the first two elements will be assumed to be min and max for uniform random generation. 
        If it's 'gauss' then the first will be the mean value and the second the standard deviation.
        !!! IF IT'S A LIST OR A TUPLE, IT MUST HAVE 3 VALUES
  
  Exemple : 
  (if nb_planets=3)
  >>> a = setParameter(2.0, 3)
  >>> print(a)
  [2.0, 2.0, 2.0]
  
  >>> a = setParameter(2, 3)
  >>> print(a)
  [2, 2, 2]
  
  >>> a = setParameter((1, 2, 'uniform'), 3)
  >>> print(a)
  [1.0422813370661363, 1.5658739594149007, 1.3426454303851822]
  
  >>> a = setParameter((0.05, 0.01, 'gaussian'), 3)
  >>> print(a)
  [0.061185035654398909, 0.049134828792707953, 0.034602504519382703]
  
  
  return : the list of planet parameters.
  """
  # Significative numbers for the figures
  SIGNIFICATIVE_NUMBERS = 4
  
  if (type(parameter) in (int, float)):
    output_parameters = [significativeRound(parameter,SIGNIFICATIVE_NUMBERS) for i in range(nb_planets)]
    return output_parameters
  elif (type(parameter) == list):
    # We assume in this case that we specified explicitely the values for all the planet and do not change anything.
    return parameter[0:nb_planets]
  elif (type(parameter) == tuple):
    output_parameters = []
    if (type(parameter[0]) == list):
      output_parameters = list(parameter[0]) # Otherwise, the element itself is modified when output_parameters is modified.
      nb_planets -= len(output_parameters)
      parameter = parameter[1]
    if (parameter[2] == 'uniform'):
      output_parameters.extend([significativeRound(uniform(parameter[0], parameter[1]),SIGNIFICATIVE_NUMBERS) for i in range(nb_planets)])
      return output_parameters
    elif (parameter[2] == 'gaussian'):
      while (len(output_parameters) < nb_planets):
        tmp_value = significativeRound(gauss(parameter[0], parameter[1]),SIGNIFICATIVE_NUMBERS)
        
        above_min_value  = ((vmin == None) or ((vmin != None) and (vmin < tmp_value)))
        inside_max_value = ((vmax == None) or ((vmax != None) and (vmax > tmp_value)))
        
        if (above_min_value and inside_max_value):
          output_parameters.append(tmp_value)
        
      return output_parameters
  else:
    raise ValueError("The way you're trying to set the parameter doesn't seems to be implemented so far")

def significativeRound(real,roundNumber):
  """function that return a truncative floating number of the number 
  in parameter, given the roundNumber we want
  Version 2.0
  
  Parameters
  real : The floating point number we want to truncate
  roundNumber : the number of significative figures we want
  """
  #print "Warning:pour test, ne fait rien!!"
  #return real
  import math
  
  if (roundNumber == 0):
    print("'roundNumber' is set to 0")
    exit
  elif (type(roundNumber) != int):
    print("'roundNumber' is not an integer")
    exit
  
  # If an integer, we'll have problem for divisions.
  if (type(real) != float):
    real = float(real)
  
  if (real == 0.):
    return 0.
  
  # in case we have negative number
  if (real < 0.):
    negative = True
    real = - real
  else:
    negative = False
  
  significativeNumbers = []
  
  # First we initialise various parameters for the loop. The principle
  # is to extract the integer part of the logarithm. This way, we 
  # retrieve the first significative number. By substracting this 
  # number to the original number, we can successively retrieve all 
  # the significative figures.
  
  # We first get the order of magnitude of the number to be able to 
  # retrieve all the significative numbers in the loop after that.
  logNumber = int(math.log10(real))

  significativeNumbers = [0] # we force the 'while' loop to do at least one turn
  # For numbers under 1, we search for the first significative number
  while (significativeNumbers[-1] == 0):
    ordMagn = 10**logNumber
    significativeNumbers = [int(real/ordMagn)]
    decimalPart = real/ordMagn - significativeNumbers[-1]
    logNumber = logNumber - 1 # if this is the last turn, this line will do nothing. This avoid to substract one in the first turn (if we had done it at the beginning)

  if (significativeNumbers[-1] == 0):
    ordMagn = ordMagn / 10
    try:
      significativeNumbers = [int(real/ordMagn)]
    except:
      print(real)
      pdb.set_trace()
    decimalPart = real/ordMagn - significativeNumbers[-1]
    
  for i in range(roundNumber-2):
    significativeNumbers.append(int(decimalPart * 10))
    decimalPart = decimalPart * 10 - significativeNumbers[-1]
  
  # For the last number we must round.
  significativeNumbers.append(int(round(decimalPart * 10 + 0.5)))
    
  numberWithoutMagn = 0.
  for (i, si) in enumerate(significativeNumbers):
    numberWithoutMagn += si * 10**(-i)
  
  truncatedNumber = ordMagn * numberWithoutMagn
  
  if not(negative):
    return truncatedNumber
  else:
    return -truncatedNumber

def number_fill(number, fill):
  """function that create a string given a number in parameter and the length 
  of the final string we want. It is used, for instance, to generate name for 
  folders were a counter is used. We get the maximal length of the string with 
  the maximal number of folders, and with this, we will know how many "0" 
  we need to display the folders. 
  
  Parameters 
  number : the number we want to display
  fill : the total length of the string we want. Zero will be added on the left to get this length
  
  Return : A string that display the number. 
  
  RMQ : An error is returned if the length of the string is not sufficient to display the number.
  
  Example : 
  display = number_fill(18, 4) 
  print(display)
  0018
  """
  
  if (type(number) != int) or (type(fill) != int):
    raise TypeError("'number' and 'fill' must be integers")
  
  number = str(number)
  size = len(number)
  
  if (size > fill):
    raise ValueError("The desired length of the string is not sufficient to display the number")
  else:
    number = "0"*(fill - size) + number
  
  # While the total length is not sufficient, we add '0' to the left.
  while (len(number) < fill):
    number = "0"+number
  
  return number

def getHostname():
  """This function will return the machine name as a string. 
  For instance, if the machine name is arguin, and that I get arguin.obs.u-bordeaux1.fr, I can test that with :
  if ('arguin' in hostname):
  """
  
  (stdout, stderr, returnCode) = autiwa.lancer_commande("hostname")
  hostname = stdout.split("\n")[0]
  
  return hostname

if __name__=='__main__':
  autiwa.printCR("Test of str2bool...")
  tests_outputs = ""
  for input in ["1", 'true', '.true.', 'True', 'TRUE', 'TrUe']:
    output = str2bool(input)
    if (output != True):
      tests_outputs += "str2bool doesn't return 'True' when we give in parameter : "+input+"\n"
      
  
  for input in ["0", 'false', '.false.', 'False', 'FALSE', 'FaLse']:
    output = str2bool(input)
    if (output != False):
      tests_outputs += "str2bool doesn't return 'False' when we give in parameter : "+input+"\n"
  
  if (tests_outputs == ""):
    print("Test of str2bool : OK")
  else:
    print("Test of str2bool : FAIL")
    print(tests_outputs)
  
  autiwa.printCR("Test of Job_PBS...")
  job = Job_PBS("ls", directory='.', walltime=1, name="simulation.sh")
  job.write()
  print("Test of Job_PBS : OK")
  
  autiwa.printCR("Test of SimpleJob...")
  job = SimpleJob("ls", name="simulation.sh")
  job.write()
  print("Test of SimpleJob : OK")


