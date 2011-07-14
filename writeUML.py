#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# Of course, everything is not tested, but it is planed to test as many things as possible

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "14 Juillet 2011"
__version__ = "$Revision: 1.0 $"
__credits__ = """We run a test simulation and erase all the files created after the tests. The simulations files are thought to be 
in a "simu_test" subdirectory of the directory were are the sources (and binaries) of mercury (and this script)"""

from fortranSource import *
import pdb

excluded = []

sources_filename = lister("*.f90")

# We define an object for each source file in order to read its source
# code, once it is done, subroutines, functions and programs are already
# defined. All we have to do is to write the tree. If there are subroutine
# passed as argument, we must define manually the possible values of
# these externals in order to have a complete tree.
objects = []
for source in sources_filename:
	objects.append(FortranSource(source))

# Since I could not define externals directly in the classes, we must
# define them here manually
externals = {'mdt_hkce.force':["mfo_hkce"],
			'mal_hcon.coord':["mco_h2mvs", "mco_iden", "mco_h2dh"],
			'mal_hcon.bcoord':["mco_mvs2h", "mco_iden", "mco_dh2h"],
			'mal_hcon.onestep':["mdt_mvs", "mdt_hy"],
			'mal_hvar.onestep':["mdt_bs1", "mdt_bs2", "mdt_ra15"],
			'mdt_bs1.force':["mfo_all"],
			'mdt_bs2.force':["mfo_all"],
			'mdt_ra15.force':["mfo_all"]
			}

setExternals(externals)

for prgm in ProgramSource.findProgram.values():
	if (prgm.name == "mercury"):
		prgm.writeTree(prgm.name+"_UML.svg", excluded=excluded)
	

# Uncomment the following line to know the externals you will have to define in the dictionary
#~ print(SubroutineSource.findExternal.keys())



