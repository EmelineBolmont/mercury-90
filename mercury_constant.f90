module mercury_constant
  use types_numeriques

implicit none
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MERCURY.INC    (ErikSoft   4 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Parameters that you may want to alter at some point:
!
! NMAX  = maximum number of bodies
! CMAX  = maximum number of close-encounter minima monitored simultaneously
! NMESS = maximum number of messages in message.in
! HUGE  = an implausibly large number
! NFILES = maximum number of files that can be open at the same time
!

integer, parameter :: NMAX = 2000 ! NMAX  = maximum number of bodies
integer, parameter :: CMAX = 50 ! CMAX  = maximum number of close-encounter minima monitored simultaneously
integer, parameter :: NMESS = 200 ! NMESS = maximum number of messages in message.in
integer, parameter :: NFILES = 50 ! NFILES = maximum number of files that can be open at the same time
real(double_precision), parameter :: HUGE = 9.9d29 ! HUGE  = an implausibly large number

end module mercury_constant