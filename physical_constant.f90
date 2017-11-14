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

! values in numerical Units. Time is in Days, lengths are in AU and masses are in MSUN
real(double_precision), parameter :: K2 = 2.959122082855911d-4 !< K2 = Gaussian gravitational constant [AU^3.MSUN-1.DAY-2]

! Various constants in (MSUN, AU, DAY) units
real(double_precision), parameter :: EARTH_MASS = 3.00374072d-6 !< the mass of the earth in solar mass [MSUN]

! Numerical Constants
real(double_precision), parameter :: THIRD = .33333333333333333d0
real(double_precision), parameter :: TWOTHIRD = 0.66666666666666666d0

end module physical_constant
