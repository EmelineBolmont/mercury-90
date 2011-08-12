module mercury_constant
  use types_numeriques

implicit none
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MERCURY.INC    (ErikSoft   4 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Parameters that you may want to alter at some point:

! NMAX  = maximum number of bodies
! CMAX  = maximum number of close-encounter minima monitored simultaneously
! NMESS = maximum number of messages in message.in
! HUGE  = an implausibly large number
! NFILES = maximum number of files that can be open at the same time


integer, parameter :: NMAX = 2000 ! NMAX  = maximum number of bodies
integer, parameter :: CMAX = 50 ! CMAX  = maximum number of close-encounter minima monitored simultaneously
integer, parameter :: NMESS = 200 ! NMESS = maximum number of messages in message.in
integer, parameter :: NFILES = 50 ! NFILES = maximum number of files that can be open at the same time
real(double_precision), parameter :: HUGE = 9.9d29 ! HUGE  = an implausibly large number

!Part were we define the parameters that were previously in swift.inc

!...   Maximum array size
integer, parameter :: NPLMAX = 202 ! max number of planets, including the Sun 
integer, parameter :: NTPMAX = 2000 ! max number of test particles

!...   Size of the test particle status flag
integer, parameter :: NSTAT = 3

!...   convergence criteria for danby
real(double_precision), parameter :: DANBYAC= 1.0d-14
real(double_precision), parameter :: DANBYB = 1.0d-13

!...    loop limits in the Laguerre attempts
integer, parameter :: NLAG1 = 50 
integer, parameter :: NLAG2 = 400

!...    A small number
real(double_precision), parameter :: TINY=4.D-15

end module mercury_constant