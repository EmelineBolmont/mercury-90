!******************************************************************************
! MODULE: physical_constant
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that contain all physical and mathematical constants 
!! needed for Mercury
!
!******************************************************************************

module physical_constant
  use types_numeriques

implicit none
!------------------------------------------------------------------------------

! Constants:

! DR = conversion factor from degrees to radians
! K2 = Gaussian gravitational constant 
! AU = astronomical unit in cm
! MSUN = mass of the Sun in g

real(double_precision), parameter :: PI = 3.141592653589793d0 
real(double_precision), parameter :: TWOPI = PI * 2.d0 
real(double_precision), parameter :: PIBY2 = PI * .5d0 
real(double_precision), parameter :: DEG2RAD = PI / 180.d0 !< DEG2RAD = conversion factor from degrees to radians
real(double_precision), parameter :: RAD2DEG = 180.d0 / PI !< RAD2DEG = conversion factor from radian to degrees

! Values in CGS
real(double_precision), parameter :: AU_cgs = 1.4959787e13 !< AU = astronomical unit in [cm]
real(double_precision), parameter :: MSUN_cgs = 1.9891e33 !< MSUN = mass of the Sun in [g]
real(double_precision), parameter :: DAY = 86400.d0 !< amount of second in one day. [s]

! Values in SI
! The IAU 2009 system of astronomical constants: the report of the IAU working group on numerical standards for Fundamental Astronomy
! https://en.wikipedia.org/wiki/Astronomical_constant
real(double_precision), parameter :: AU_SI = 1.49597870700d11
real(double_precision), parameter :: M_SUN = 1.9818e30 ! kg
real(double_precision), parameter :: G_SI = 6.67428d-11 ! m^3.kg^-1.s^-2 (S.I. units)

! values in numerical Units. Time is in Days, lengths are in AU and masses are in MSUN
! real(double_precision), parameter :: K2 = 2.959122082855911d-4 !< K2 = Gaussian gravitational constant [AU^3.MSUN-1.DAY-2]*/
real(double_precision), parameter :: K2 = G_SI/(AU_SI**3)*M_SUN*(DAY*DAY) !< K2 = Gaussian gravitational constant [AU^3.MSUN-1.DAY-2]

! The IAU 2010 system of astronomical constants: the report of the IAU working group on numerical standards for Fundamental Astronomy
! https://en.wikipedia.org/wiki/Astronomical_constant
! Ratio between mass of the Sun and mass of Earth
real(double_precision), parameter :: m2earth = 332946.050895d0
real(double_precision), parameter :: AU    = AU_SI
real(double_precision), parameter :: minau = 1.d0/AU

! Resolution B3 on recommended nominal conversion constants for selected solar and planetary properties
! http://adsabs.harvard.edu/abs/2015arXiv151007674M
! Radius of Earth in AU
real(double_precision), parameter :: rearth = 6.3781d6*minau
! Radius of the Sun in AU
real(double_precision), parameter :: rsun = 6.957d8*minau

! Various constants in (MSUN, AU, DAY) units
! real(double_precision), parameter :: EARTH_MASS = 3.00374072d-6 !< the mass of the earth in solar mass [MSUN]
real(double_precision), parameter :: EARTH_MASS = 1.d0/332946.050895d0 !< the mass of the earth in solar mass [MSUN]


! Numerical Constants
real(double_precision), parameter :: THIRD = .33333333333333333d0
real(double_precision), parameter :: TWOTHIRD = 0.66666666666666666d0

end module physical_constant
