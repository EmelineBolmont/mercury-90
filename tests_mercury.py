#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# Of course, everything is not tested, but it is planed to test as many things as possible

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "07 juin 2011"
__version__ = "$Revision: 1.0 $"
__credits__ = """We run a test simulation and erase all the files created after the tests. The simulations files are thought to be 
in a "simu_test" subdirectory of the directory were are the sources (and binaries) of mercury (and this script)"""

from make import clean
import os
import subprocess

def run(command):
    """run a command that will be a string.
    The function return a tuple with the output, 
    the stderr and the return code. 
    It is fitted to run only command placed in the parent directory (of type "../command")"""

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
    (process_stdout, process_stderr) = process.communicate()
    returnCode = process.poll()
    
    if (returnCode == 0):
        print(command+" was runned successfully")
        return (process_stdout, process_stderr)
    else:
        print(command+" has terminated with an error code: "+str(returnCode))
        
        temp = command.lstrip("../")
        
        errlogname = temp+".err"
        print("Writing the error of the compilation in '"+errlogname+"'")
        f = open(errlogname,'w')
        f.write(process_stderr)
        f.close()
        
        outlogname = temp+".out"
        print("Writing the output of the compilation in '"+outlogname+"'")
        f = open(outlogname,'w')
        f.write(process_stdout)
        f.close()
        
        return 1

    

FOLDER = "simu_test"


os.chdir(FOLDER)

(process_stdout, process_stderr) = run("../mercury")

(process_stdout, process_stderr) = run("../close")

(process_stdout, process_stderr) = run("../element")





# We clean undesirable files. Indeed, we will compile everything everytime.
clean(["out", "clo", "aei", "dmp", "tmp"])