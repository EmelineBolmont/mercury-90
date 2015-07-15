Welcome to the git repository of Mercury-T

Here are a few things to know:
- In order to compile the model, you need to execute the python script Makefile.py
you can modify the compilation options to your liking. Right now, it's for ifort.
- There are some files with the time evolution of the radius and some other quantities
for brown-dwarfs (mass_xx.xxxx.dat, where xx.xxxx is the mass of the BD in Jupiter mass, 
i.e., 40.0000 corresponds to a 40 Mjup BD, or 0.04 Msun as you prefer;
and rg2BD.dat with the moment of inertia informations), for a 0.1 Msun Mdwarf
(01Msun.dat), for a 1 Sun-mass star (SRad_Spli_M-1_0000.dat).
- user_module.f90 is the place where tides are implemented
- tides_constant_GR.f90 is the file where the tidal parameters are and where the 
initialization is done. If you modify this file, you have to re-compile.
- All the rest can be used as the normal Mercury code
