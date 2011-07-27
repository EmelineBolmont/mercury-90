#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# Of course, everything is not tested, but it is planed to test as many things as possible

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "22 Juillet 2011"
__version__ = "$Revision: 2.6.1 $"
__credits__ = """We run a test simulation and erase all the files created after the tests. The simulations files are thought to be 
in a "simu_test" subdirectory of the directory were are the sources (and binaries) of mercury (and this script)"""

from make import sourceFile,lister

sources_filename = lister("*.f90")

objects = []
for source in sources_filename:
	objects.append(sourceFile(source))

mercury = sourceFile.findSource["mercury.f90"]
element = sourceFile.findSource["element.f90"]
close = sourceFile.findSource["close.f90"]

excluded = ["types_numeriques", "physical_constant", "mercury_constant"]

mercury.writeArchitecture("mercury_graph.svg", excluded=excluded, direction="topbottom")
element.writeArchitecture("element_graph.svg", excluded=excluded, direction="leftright")
close.writeArchitecture("close_graph.svg", excluded=excluded, direction="leftright")
