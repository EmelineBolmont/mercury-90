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
  real(double_precision), parameter :: scaleheight = 0.05 ! the scaleheight of the disk. Is used by the function get_scaleheight
  
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

subroutine mfo_user (time,jcen,n_bodies,n_big_bodies,mass,position,velocity,acceleration)!,spin)
!  mass      = mass (in solar masses)
!  position     = coordinates (x,y,z) with respect to the central body [AU]
!  velocity     = velocities (vx,vy,vz) with respect to the central body [AU/day]
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
  
  !------------------------------------------------------------------------------ 
  !------Local-------
  
  ! Properties of the planet
  real(double_precision), dimension(n_bodies-1) :: radius_p ! the radial position of the planet [AU]
  real(double_precision), dimension(n_bodies-1) :: velocity_p ! the norm of the speed [AU/day]
  real(double_precision), dimension(n_bodies-1) :: omega_p ! the angular rotation [day-1]
  real(double_precision), dimension(n_bodies-1) :: gamma_0 ! Normalization factor for all the torques. [?]
  real(double_precision), dimension(n_bodies-1) :: torque ! the torque exerted by the disk on the planet [?]
  real(double_precision), dimension(n_bodies-1) :: sigma_p ! the surface density of the gas disk at the planet location [?]
  real(double_precision), dimension(n_bodies-1) :: scaleheight_p ! the scaleheight of the gas disk at the location of the planet [no dim]
  real(double_precision), dimension(n_bodies-1) :: time_mig ! The migration timescale for the planet [day]
  real(double_precision), dimension(n_bodies-1) :: angular_momentum ! the angular momentum of the planet [?]
  !------------------------------------------------------------------------------
  ! Setup
  acceleration(:,:) = 0.d0
  !------------------------------------------------------------------------------
  
  ! WE CALCULATE PROPERTIES OF THE PLANETS
  ! Not for all the bodies because there might be problem to divide by radial position of the star.
  radius_p(:) = sqrt(position(1,2:)**2 + position(2,2:)**2)
  velocity_p(:) = sqrt(velocity(1,2:)**2 + velocity(2,2:)**2)
  
  ! \vect{v} = \vect{\omega} \wedge \vect{x}
  ! By assuming that the planet is in the plane of the disk (which is false) we get =
  omega_p(:) = velocity_p(:) / radius_p(:) ! TODO UNITS!!!!!!!!!
  write(*,*) "Warning: Units for omega_p are currently not verified"
  
  sigma_p(:) = sigma_0 * radius_p(:) ** sigma_index ! need to be in numerical units, not the case for the moment
  scaleheight_p(:) = scaleheight
  
  !------------------------------------------------------------------------------
  ! WE CALCULATE TOTAL TORQUE EXERTED BY THE DISK ON THE PLANET
  gamma_0(:) = (mass(2:) / (mass(1) * scaleheight_p(:)))**2 * sigma_p(:) * radius_p(:)**4 * omega_p(:)**2
  
  !------------------------------------------------------------------------------

  
  time_mig(:) = angular_momentum(:) / torque(:)
  
  acceleration(1,:) = - (0.5d0 / time_mig(:)) * velocity(1,:)
  acceleration(2,:) = - (0.5d0 / time_mig(:)) * velocity(2,:)
  acceleration(3,:) = - (0.5d0 / time_mig(:)) * velocity(3,:)
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_user


end module user_module

! TODO utiliser la masse des objets pour ne pas faire le calcul si trop massif, il faut respecter le domaine de validité des formules des couples
! TODO routine générale de conversion des couples en accélération afin de pouvoir réutiliser ailleurs

!TODO : 
!_ amortissement de l'eccentricité
!_ amortissement de l'inclinaison
!_rajouter le spin comme variable d'entrée de la routine mfo_user
