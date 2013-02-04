module user_module

!*************************************************************
!** Modules that contains user defined modules. 
!** Only mfo_user will be public.
!**
!** Version 1.0 - june 2011
!*************************************************************

! The user_module is divided in several parts. First of all, some parameters are calculated at the beginning of the run, in 
! init_globals. This subroutine must be executed before the rest. Since mfo_user is used by mercury, we can't execute init_globals 
! only once, we must include an 'if' test that  run init_globals only if it's the first time we pass in the routine. Another thing 
! to keep in mind is that the temperature profile of the disk is calculated manually. Since the density profile evolve in time, You
! must re-calculate the temperature profile each time you need it. A mistake can be made because the surface density will evolve 
! in time automatically. BUT, the temperature profile will not correspond to it. Another set of routine are prefixed with 'test_', 
! they are all coded to test something, either a routine or a physical value (and plot it). Each plot consist in a data file (*
! .dat) and a gnuplot file (*.gnuplot). You only have to run the gnuplot file with gnuplot with a command "gnuplot file.gnuplot". 

! Most of the plot are stored in a subdirectory "unitary_tests" of the current directory Theses tests are used ONLY in the 
! subroutine "unitary_tests" that can be used outside the module. For instance, I made a fortran source code that use this module 
! and call the routine 'unitary_tests'. That's how I test my module before running it in mercury. The problem is to prepare 
! variable in the same way, both in mfo_user and my tests. That's why init_globals exists, thus I can run it in my tests also.

  use types_numeriques
  use physical_constant

  implicit none
  
  private
  
  public :: mfo_user
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MFO_USER.FOR    (ErikSoft   2 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Applies an arbitrary force, defined by the user.

! If using with the symplectic algorithm MAL_MVS, the force should be
! small compared with the force from the central object.
! If using with the conservative Bulirsch-Stoer algorithm MAL_BS2, the
! force should not be a function of the velocities.

! Code Units are in AU, days and solar mass * K2 (mass are in solar mass, but multiplied by K2 earlier in the code).

! N.B. All coordinates and velocities must be with respect to central body
! ===
!------------------------------------------------------------------------------

subroutine mfo_user (time,jcen,n_bodies,n_big_bodies,mass,position,velocity,acceleration)
!  mass          = mass (in solar masses * K2)
!  position      = coordinates (x,y,z) with respect to the central body [AU]
!  velocity      = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  n_bodies      = current number of bodies (INCLUDING the central object)
!  n_big_bodies  =    "       "    " big bodies (ones that perturb everything else)
!  time          = current epoch [days]

  use mercury_constant
  
  implicit none

  
  ! Input
  integer, intent(in) :: n_bodies,n_big_bodies
  real(double_precision),intent(in) :: time,jcen(3),mass(n_bodies),position(3,n_bodies),velocity(3,n_bodies)
  
  ! Output
  real(double_precision),intent(out) :: acceleration(3,n_bodies)
  
  !------------------------------------------------------------------------------ 
  !------Local-------

  ! loop integers
  integer :: planet
  
  do planet = 1, n_bodies
    acceleration(1,planet) = 0.d0
    acceleration(2,planet) = 0.d0
    acceleration(3,planet) = 0.d0
  end do
  
  !------------------------------------------------------------------------------
  return
end subroutine mfo_user

end module user_module

