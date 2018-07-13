; Constants 

; Units
AU        =  1.49597870700d11         ;m
yr        =  365.25*24*3600           ;s
day       =  24*3600.                 ;s
hr        =  3600.                    ;s

; Sun stuff
Msun      =  1.9818d30                ;kg
Rsun      =  6.957d8                  ;m
Lsun      =  3.939d26                 ;W

; Masses of some planets
Mjup      =  9.547922048d-4 * Msun    ;kg
Rjup      =  69.911d6                 ;m

m2earth   =  332946.050895d0          ;Msun/Mearth
Mearth    =  1.d0/m2earth * Msun      ;kg
Rearth    =  6.3781d6                 ;m

Mmars     =  6.4185d23                ;kg
Rmars     =  3.3895d6                 ;m

; Some physical constants
G         =  6.67428d-11              ;m^3.kg^-1.s^-2
cc        =  299792458.d0             ;m/s
hp        =  6.626d-34                ;J.s
sig       =  5.67d-8                  ;W.m-2.K-4
kb        =  1.3806503d-23            ;m2.kg.s-2.K-1


;! moment of inertia (I/mR^2)
rg2_sun    =  5.9d-2
rg2_jup    =  2.54d-1
rg2_earth  =  3.308d-1
