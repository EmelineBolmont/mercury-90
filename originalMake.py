#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for the original mercury files

from make import *

FOLDER = "mercury_original"

os.chdir(FOLDER)


# We clean undesirable files. Indeed, we will compile everything everytime.
clean(["o"])

sources_filename = lister("*.for")

mains = {"mercury6_2.for":"mercury", "element6.for":"element", "close6.for":"close"}

# We create the binaries
make_binaries(sources_filename, mains)
