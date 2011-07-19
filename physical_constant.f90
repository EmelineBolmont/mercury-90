module physical_constant
  use types_numeriques

implicit none
save
!------------------------------------------------------------------------------

! Constants:

! DR = conversion factor from degrees to radians
! K2 = Gaussian gravitational constant squared
! AU = astronomical unit in cm
! MSUN = mass of the Sun in g

real(double_precision), parameter :: PI = 3.141592653589793d0 
real(double_precision), parameter :: TWOPI = PI * 2.d0 
real(double_precision), parameter :: PIBY2 = PI * .5d0 
real(double_precision), parameter :: DEG2RAD = PI / 180.d0 ! DEG2RAD = conversion factor from degrees to radians
real(double_precision), parameter :: RAD2DEG = 180.d0 / PI ! RAD2DEG = conversion factor from radian to degrees

real(double_precision), parameter :: K2 = 2.959122082855911d-4 ! K2 = Gaussian gravitational constant squared
real(double_precision), parameter :: AU = 1.4959787e13 ! AU = astronomical unit in cm
real(double_precision), parameter :: MSUN = 1.9891e33 ! MSUN = mass of the Sun in g

! Numerical Constants
real(double_precision), parameter :: THIRD = .3333333333333333d0

end module physical_constant