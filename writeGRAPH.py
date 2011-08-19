#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that generate a .svg file for each binaries of mercury that display the tree of modules used

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "19 août 2011"
__version__ = "$Revision: 2.6.2 $"

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
