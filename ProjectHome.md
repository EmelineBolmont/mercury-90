# Disclaimer #
The **master** contains only the code and should be used as his. Unfortunately, I worked during almost 2 years in master, adding my own modifications, that are now in the **migration** branch. As a consequence, the first useable version of **master**, for public, is the revision **97a6070cfc2c** launched in _2013 february the 5th_.

# Aim #
Introduction of fortran 90 modules, modification of arguments. The aim is to keep the original behaviour of mercury. The present code must respect :
  * Fully compatible with mercury configuration files (param.in, big.in, small.in and so on)
  * MUST return exactly the same results than the original 'mercury' binaries. A python script has been made to be sure of this. That means that even -fastmath compiling option is not allowed to retrieve the original outputs of mercury.

**Important Note : The file 'element.out' now display masses in earth mass instead of solar masses.**


# How To Use #
A folder containing an example simulation is available in example\_simulation/

To use the code, my advice would be to have one folder for the code, in a specific location. Then, in another location, you can
create as many folders as you want, with one specific simulation in each folder.

Be sure to run the installation script (to add useful shortcuts to your .bash\_profile):
```
configure_mercury.py
```

You can jump to your program directory via :
```
cd $mercury
```

To launch Mercury wherever you want in your server via :
```
mercury
```

Keep in mind that typing
```
> cd -
```
in your Terminal, you will go back to the previous directory (not the parent one, the previous one). Typing this command several
times allow you to switch easily between two folders.

# Basic commands #

Compile Mercury, Element and Close
```
Makefile.py
```

See all compilation options
```
Makefile.py help
```