module user_module

!*************************************************************
!** Modules that contains user defined modules. 
!** Only mfo_user will be public.
!**
!** Version 1.0 - june 2011
!*************************************************************
  use types_numeriques

  implicit none
  
  private
  
  public :: mfo_user
  
  real(double_precision), parameter :: b_over_h = 0.4 ! the smoothing length for the planet's potential
  real(double_precision), parameter :: adiabatic_index = 1.4 ! the adiabatic index for the gas equation of state
  
  ! Here we define the power law for surface density sigma(R) = sigma_0 * R^sigma_index
  real(double_precision), parameter :: sigma_0 = 1700 ! the surface density at (R=1) [g/cm^2]
  real(double_precision), parameter :: sigma_index = 0.5! the slope of the surface density power law
  
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

! N.B. All coordinates and velocities must be with respect to central body
! ===
!------------------------------------------------------------------------------

subroutine mfo_user (time,jcen,n_bodies,n_big_bodies,mass,position,velocity,acceleration)
!  mass      = mass (in solar masses)
!  position     = coordinates (x,y,z) with respect to the central body [AU]
!  speed     = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  n_bodies  = current number of bodies (INCLUDING the central object)
!  n_big_bodies  =    "       "    " big bodies (ones that perturb everything else)
!  time  = current epoch [days]

  use physical_constant
  use mercury_constant  

  implicit none

  
  ! Input
  integer, intent(in) :: n_bodies,n_big_bodies
  real(double_precision),intent(in) :: time,jcen(3),mass(n_bodies),position(3,n_bodies),velocity(3,n_bodies)
  
  ! Output
  real(double_precision),intent(out) :: acceleration(3,n_bodies)
  
  ! Local
  integer :: planet ! index for the planet loop
  real(double_precision) :: gamma_0 ! Normalization factor for all the torques.
  
  ! Properties of the planet
  real(double_precision), dimension(n_bodies) :: radius_p ! the radial position of the planet [AU]
  real(double_precision), dimension(n_bodies) :: velocity_p ! the norm of the speed [AU/day]
  real(double_precision), dimension(n_bodies) :: omega_p ! the angular rotation [day-1]
  !------------------------------------------------------------------------------
  ! Setup
  acceleration(:,:) = 0.d0
  !------------------------------------------------------------------------------
  
  ! Not for all the bodies because there might be problem to divide by radial position of the star.
  radius_p(2:) = sqrt(position(1,2:)**2 + position(2,2:)**2)
  velocity_p(2:) = sqrt(velocity(1,2:)**2 + velocity(2,2:)**2)
  
  ! \vect{v} = \vect{\omega} \wedge \vect{x}
  ! By assuming that the planet is in the plane of the disk (which is false) we get =
  omega_p(2:) = velocity_p(2:) / radius_p(2:)
  
  ! For every object except the central one
  do planet = 2,n_bodies
    
    
  !TODO loop for planets, but we might not apply the acceleration to all the bodies, we must also do a test regarding their masses for exemple.
    gamma_0 = (mass(planet) * get_scaleheight(radius_p(planet), time))**2 * get_sigma(radius_p(planet), time) &
            * radius_p(planet)**4 * omega_p(planet)**2
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_user

function get_sigma(radius, time)
  ! Function that return the surface density at the given radius.
  ! Maybe a next version will include dissipation of the disk with time.
  
  ! Parameter:
  !  radius : the radius at which we want the surface density [AU]
  !  time  = current epoch [days]
  
  ! Return:
  ! get_sigma: the surface density at R=radius [solar masses / AU^2]

  implicit none
  
  real(double_precision), intent(in) :: radius
  real(double_precision), intent(in) :: time
  real(double_precision) :: get_sigma
  
  get_sigma = sigma_0 * radius ** sigma_index
  
  return
end function get_sigma

function get_scaleheight(radius, time)
  ! Function that return the surface density at the given radius.
  ! Maybe a next version will include dissipation of the disk with time.
  
  ! Parameter:
  !  radius : the radius at which we want the surface density [AU]
  !  time  = current epoch [days]
  
  ! Return:
  ! scaleheight: the scaleheight at R=radius [AU]

  implicit none
  
  real(double_precision), intent(in) :: radius
  real(double_precision), intent(in) :: time
  
  real(double_precision) :: get_scaleheight
  
  get_scaleheight = 0.05 ! TODO for the moment, constant, but must be changed afterwards
  
  return
end function get_scaleheight

end module user_module

! TODO utiliser la masse des objets pour ne pas faire le calcul si trop massif, il faut respecter le domaine de validité des formules des couples
! TODO routine générale de conversion des couples en accélération afin de pouvoir réutiliser ailleurs
