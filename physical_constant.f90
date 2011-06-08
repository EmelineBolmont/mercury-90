module physical_constant
  use types_numeriques

implicit none
!
!------------------------------------------------------------------------------
!
! Constants:
!
! DR = conversion factor from degrees to radians
! K2 = Gaussian gravitational constant squared
! AU = astronomical unit in cm
! MSUN = mass of the Sun in g
!
real(double_precision), parameter :: PI = 3.141592653589793d0 
real(double_precision), parameter :: TWOPI = PI * 2.d0 
real(double_precision), parameter :: PIBY2 = PI * .5d0 
real(double_precision), parameter :: DR = PI / 180.d0 ! DR = conversion factor from degrees to radians

real(double_precision), parameter :: K2 = 2.959122082855911d-4 ! K2 = Gaussian gravitational constant squared
real(double_precision), parameter :: AU = 1.4959787e13 ! AU = astronomical unit in cm
real(double_precision), parameter :: MSUN = 1.9891e33 ! MSUN = mass of the Sun in g

end module physical_constant