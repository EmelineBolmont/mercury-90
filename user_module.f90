!******************************************************************************
! MODULE: user_module
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contain user defined function. This module can call other module and subroutine. 
!! The only public routine is mfo_user, that return an acceleration that 
!! mimic a random force that depend on what the user want to model.
!
!******************************************************************************

module user_module

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
  integer, intent(in) :: n_bodies !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
  integer, intent(in) :: n_big_bodies !< [in] current number of big bodies (ones that perturb everything else)
  real(double_precision), intent(in) :: time !< [in] current epoch (days)
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: mass(n_bodies) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: position(3,n_bodies)
  real(double_precision), intent(in) :: velocity(3,n_bodies)
  
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

