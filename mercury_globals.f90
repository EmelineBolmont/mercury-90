module mercury_globals

!*************************************************************
!** Modules that contains all the globals variables of mercury
!**
!** Version 1.0 - june 2011
!*************************************************************
  use types_numeriques

  implicit none
  
  integer, dimension(8) :: opt = (/0,1,1,2,0,1,0,0/) ! Default options (can be overwritten later in the code) for mercury.



end module mercury_globals
