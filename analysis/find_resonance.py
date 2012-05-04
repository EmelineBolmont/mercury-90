#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script that searches for resonances close to the ratio given in paremeter."""
from __future__ import print_function

__author__ = "Autiwa <autiwa@gmail.com>"
__date__ = "14 fevrier 2012"
__version__ = "$Revision: 1.1"
__credits__ = """Thanks to Bastien for his help when I was learning Python and his usefull tutorial"""



from fractions import Fraction
import sys
import pdb

# TODO
# _ 

###################
## Configuration ##
###################
WRITE_TO_FILE = False

################
## Parameters ##
################
DENOMINATOR_LIMIT = 20 # Maximum value allowed of the denominator when we want to get a fraction from a decimal value
NUMBER_OF_VALUES = 100 # sampling for period ratio around the given value
UNCERTAINTY = 5 # In percentage
OUTPUT_FILE = "resonances.txt"

###############################################
## Beginning of the program
###############################################

uncertainty = 0.01 * float(UNCERTAINTY)

isProblem = False
isElement = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * ratio (the period ratio we want to investigate)" + "\n" + \
" * uncertainty (the tolerance over the period ratio in % ; the default is "+str(UNCERTAINTY)+"%)" + "\n" + \
" * element (to get possible resonances for every planet in an element.out file" + "\n" + \
" * help (to get this message)" + "\n" + \
"\nExample: find_resonances.py ratio=1.2"

# We get arguments from the script
for arg in sys.argv[1:]:
	try:
		(key, value) = arg.split("=")
	except:
		key = arg
	if (key == 'ratio'):
		periodRatio = float(value)
	elif (key == 'uncertainty'):
		uncertainty = 0.01 * float(value)
	elif (key == 'element'):
		isElement = True
	elif (key == 'help'):
		isProblem = True
	else:
		print("the key '"+key+"' does not match")
		isProblem = True

if ('periodRatio' not in vars().keys()):
	if not(isProblem):
		print("Error: The periodRatio is not set")
		isProblem = True

if isProblem:
	print(problem_message)
	exit()

periodMin = periodRatio * (1 - uncertainty)
periodMax = periodRatio * (1 + uncertainty)
periodWidth = periodMax - periodMin
deltaPeriod = periodWidth/NUMBER_OF_VALUES

periods = [periodMin + deltaPeriod * i for i in range(NUMBER_OF_VALUES)]

resonances = []
for period in periods:
  fraction = Fraction(period).limit_denominator(DENOMINATOR_LIMIT)
  #print(fraction)
  resonances.append(fraction)

# We exclude all values that appears several time to only keep one occurence of each value
resonances = list(set(resonances))
resonances = [(round(float(fraction),3), str(fraction)) for fraction in resonances]
resonances.sort()

print("Possible resonances for a period ratio close to "+str(periodRatio)+" from "+str(periodMin)+" to "+str(periodMax)+" : ")
print(*resonances,sep=" ; ")


# We write our results to a file if the option is set.
if WRITE_TO_FILE:
  outputFile = open(OUTPUT_FILE, 'w')
  outputFile.write("Possible resonances for a period ratio close to "+str(periodRatio)+" from "+str(periodMin)+" to "+str(periodMax)+"\n")
  for res in resonances:
    outputFile.write(str(res)+"\n")
  outputFile.close()

#pdb.set_trace()
