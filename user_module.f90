module user_module

!*************************************************************
!** Modules that contains user defined modules. 
!** Only mfo_user will be public.
!**
!** Version 1.0 - june 2011
!*************************************************************
  use types_numeriques
  use physical_constant

  implicit none
  
!~   private
!~   
!~   public :: mfo_user
  
  real(double_precision), parameter :: b_over_h = 0.4 ! the smoothing length for the planet's potential
  real(double_precision), parameter :: adiabatic_index = 1.4 ! the adiabatic index for the gas equation of state
  real(double_precision), parameter :: aspect_ratio = 0.05 ! the aspect_ratio of the disk. Is used by the function get_aspect_ratio
  
  ! Here we define the power law for surface density sigma(R) = sigma_0 * R^sigma_index
  real(double_precision), parameter :: sigma_0 = 1700 ! the surface density at (R=1AU) [g/cm^2]
  real(double_precision), parameter :: sigma_index = 0.5! the slope of the surface density power law (alpha in the paper)
  real(double_precision), parameter :: sigma_0_num = sigma_0 * AU**2 / MSUN ! the surface density at (R=1AU) [Numerical Units]
  
  ! Here we define the power law for temperature T(R) = temperature_0 * R^temperature_index
  real(double_precision), parameter :: temperature_0 = 150 ! the temperature at (R=1AU) [K]
  real(double_precision), parameter :: temperature_index = 1.! the slope of the temperature power law (beta in the paper)
  
  !prefactors
  real(double_precision) :: x_s_prefactor
  real(double_precision) :: chi_p_prefactor
  real(double_precision) :: prefactor_lindblad ! prefactor for the lindblad torque
  
  real(double_precision) :: torque_hs_baro ! barotropic part of the horseshoe drag
  real(double_precision) :: torque_c_lin_baro ! barotropic part of the linear corotation torque
  
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

! Code Units are in AU, days and solar mass.

! N.B. All coordinates and velocities must be with respect to central body
! ===
!------------------------------------------------------------------------------

subroutine mfo_user (time,jcen,n_bodies,n_big_bodies,mass,position,velocity,acceleration)
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
  
  ! loop integers
  integer :: planet
  
  
  real(double_precision), dimension(n_bodies) :: torque ! the torque exerted by the disk on the planet [?]
  real(double_precision), dimension(n_bodies) :: time_mig ! The migration timescale for the planet [day]
  real(double_precision), dimension(n_bodies) :: angular_momentum ! the angular momentum of the planet [?]

  !------------------------------------------------------------------------------
  ! Setup
  acceleration(:,:) = 0.d0
  !------------------------------------------------------------------------------
  
  x_s_prefactor = 1.1 * (b_over_h / 0.4)**0.25 / mass(1) ! mass(1) is here for the ratio of mass q

  chi_p_prefactor = 16. * adiabatic_index * (adiabatic_index - 1) * SIGMA_STEFAN  / 3.

  prefactor_lindblad = -2.5 - 1.7 * temperature_index + 0.1 * sigma_index
  torque_hs_baro = 1.1 * (3.d0/2.d0 - sigma_index)
  torque_c_lin_baro = 0.7 * (3.d0/2.d0 - sigma_index)
  
  do planet=2,n_big_bodies
    torque(planet) = get_total_torque(mass(1), mass(planet), position(3,planet), velocity(3,planet))
    
    
  end do
  

  !------------------------------------------------------------------------------
  
  ! We calculate the angular momentum, that is, the z-composant. We assume that the x and y composant are negligible (to be tested)
  angular_momentum(:) = mass(:) * (position(1,:) * velocity(2,:) - position(2,:) * velocity(1,:))
!~   angular_momentum_x(:) = mass(2:) * (position(2,2:) * velocity(3,2:) - position(3,2:) * velocity(2,2:))
!~   angular_momentum_y(:) = mass(2:) * (position(3,2:) * velocity(1,2:) - position(1,2:) * velocity(3,2:))
  
  time_mig(:) = 0.5d0 * angular_momentum(:) / torque(:)
  
  acceleration(1,:) = - velocity(1,:) / time_mig(:)
  acceleration(2,:) = - velocity(2,:) / time_mig(:)
  acceleration(3,:) = - velocity(3,:) / time_mig(:)
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_user

function get_total_torque(stellar_mass, mass, position, velocity)
! function that return the total torque exerted by the disk on the planet 
!
! Parameter:
! mass : an array will all the masses (including the central one as the first element [solar mass]
! radius_p : the distance between the central body and each planet
  use orbital_elements, only : mco_x2a
  implicit none
  real(double_precision), intent(in) :: stellar_mass, mass, position(3), velocity(3)
  
  real(double_precision) :: get_total_torque

  !Local
  
  ! Meaningless parameters (intermediate constant and so on that doesn't mean anything physically)
  real(double_precision) :: Q_p ! parameter for gamma_eff (equation (45) of Paardekooper, Baruteau, 2010 II. Effects of diffusion)
  real(double_precision) :: k_p ! parameter for p_nu and p_chi for example  !!! This is not 'k' from equation (15)!
  real(double_precision) :: gm ! parameter to calculate the semi major axis. passed to mco_x2a
  
  ! Properties of the planet
  real(double_precision) :: radius_p ! the radial position of the planet [AU]
  real(double_precision) :: velocity_p, vel_squared ! the norm of the speed [AU/day] and the velocity squared.
  real(double_precision) :: omega_p ! the angular rotation [day-1]
  real(double_precision) :: semi_major_axis ! semi major axis of the planet [AU]
  
  !Properties of the disk at the location of the planet

  real(double_precision) :: Gamma_0 ! Normalization factor for all the torques. [?](equation (8) of Paardekooper, Baruteau, 2009)
  real(double_precision) :: sigma_p ! the surface density of the gas disk at the planet location [?]
  real(double_precision) :: aspect_ratio_p ! the aspect_ratio of the gas disk at the location of the planet [no dim]
  real(double_precision) :: x_s ! semi-width of the horseshoe region [radius_p (in unity of position of the planet)]
  real(double_precision) :: gamma_eff ! effective adiabatic index depending on several parameters
  real(double_precision) :: zeta_eff ! effective entropy index depending on gamma_eff
  real(double_precision) :: chi_p ! the thermal diffusion coefficient at the location of the planet [?]
  real(double_precision) :: nu_p ! the viscosity of the disk at the location of the planet [?]
  real(double_precision) :: opacity_p ! the opacity of the disk at the location of the planet [?]
  real(double_precision) :: temperature_p ! the temperature of the disk at the location of the planet [K]
  real(double_precision) :: p_nu ! parameter for saturation due to viscosity at the location of the planet
  real(double_precision) :: p_chi ! parameter for saturation due to thermal diffusion at the location of the planet
  
  !Torques (some depends of the planet)
  real(double_precision) :: torque_hs_ent ! entropy related part of the horseshoe drag
  real(double_precision) :: torque_c_lin_ent ! entropy related part of the linear corotation torque
  
  ! WE CALCULATE PROPERTIES OF THE PLANETS
  gm = K2 * (stellar_mass + mass)
  call mco_x2a (gm,position(1), position(2), position(3), velocity(1),velocity(2), velocity(3),semi_major_axis,radius_p,vel_squared)
  velocity_p = sqrt(vel_squared)
  
  !------------------------------------------------------------------------------
  omega_p = sqrt(gm / (semi_major_axis * semi_major_axis * semi_major_axis)) ! [day-1]
  write(*,*) "Warning: Units for omega_p are currently not verified"
  
  sigma_p = sigma_0_num * radius_p ** sigma_index ! [N.U.]
  temperature_p = temperature_0 * radius_p**temperature_index ! [K]
  aspect_ratio_p = aspect_ratio
  
  !------------------------------------------------------------------------------
  opacity_p = get_opacity(temperature_p, sigma_p)
  write(*,*) "Warning: opacity currently not set, just put to 1 to avoid compilation errors"
  
  !------------------------------------------------------------------------------
  ! WE CALCULATE TOTAL TORQUE EXERTED BY THE DISK ON THE PLANET
  Gamma_0 = (mass / (stellar_mass * aspect_ratio_p))**2 * sigma_p * radius_p**4 * omega_p**2
  
  chi_p = chi_p_prefactor * temperature_p**4 / (opacity_p * sigma_p**2 * omega_p**2)
  
  ! Q is needed by the lindblad torque. We set Q for m ~ 2 /3 h : 
  Q_p = 2 * chi_p / (3 * aspect_ratio_p**3 * radius_p**2 * omega_p)
  !------------------------------------------------------------------------------
  
  gamma_eff = 2 * Q_p * adiabatic_index / (adiabatic_index * Q_p + 0.5d0 * &
  sqrt(2.d0 * sqrt((adiabatic_index * adiabatic_index * Q_p * Q_p + 1)**2 + 16.d0 * Q_p * Q_p * (adiabatic_index - 1)) &
  + 2.d0 * adiabatic_index * adiabatic_index * Q_p * Q_p -2))
  
  !------------------------------------------------------------------------------
  zeta_eff = temperature_index - (gamma_eff - 1) * sigma_index
  
  x_s = x_s_prefactor / gamma_eff**0.25 * sqrt(mass / aspect_ratio_p)
  
  ! k_p is defined to limit the number of operation and to have a value independant from chi_p or nu_p
  k_p = radius_p * radius_p * omega_p * x_s * x_s * x_s / (2 * PI)
  
  !------------------------------------------------------------------------------
  
  p_nu = (2/3.d0) * sqrt(k_p / nu_p)
  
  p_chi = sqrt(k_p / chi_p)
  
  torque_hs_ent = 7.9 * zeta_eff / gamma_eff
  torque_c_lin_ent = (2.2 - 1.4 / gamma_eff) * zeta_eff
  
  !------------------------------------------------------------------------------

  get_total_torque = (Gamma_0 / gamma_eff) * (prefactor_lindblad &
    + torque_hs_baro * get_F(p_nu) * get_G(p_nu) + torque_c_lin_baro * (1 - get_K(p_nu)) &
    + torque_hs_ent * get_F(p_nu) * get_G(p_chi) * sqrt(get_G(p_nu) * get_G(p_chi)) &
    + torque_c_lin_ent * sqrt((1 - get_K(p_nu)) * (1 - get_K(p_chi))))

  
  return
end function get_total_torque

function get_F(p)
! F function (22) of the paper : "A torque formula for non-isothermal Type I planetary migration - II Effects of diffusion"
! By Paardekooper, baruteau and Kley, (2010)
! Equation (22) can be approximated within 5% by a much more simpler equation which will be the one we use. 
!
! Parameter:
! p : parameter
  implicit none
  real(double_precision), intent(in) :: p
  
  real(double_precision) :: get_F

  !Local
  !------------------------------------------------------------------------------

  
  get_F = 1 / (1 + (p / 1.3) * (p / 1.3))
  
  return
end function get_F

function get_G(p)
! G function (30) of the paper : "A torque formula for non-isothermal Type I planetary migration - II Effects of diffusion"
! By Paardekooper, baruteau and Kley, (2010)
! Equation (30) is a particuliar case of the (29) equation with l=3/4 and tau_0=2
!
! Parameter:
! p : parameter
  implicit none
  real(double_precision), intent(in) :: p
  
  real(double_precision) :: get_G

  !Local
  real(double_precision), parameter :: p_0 = sqrt(8. / (45. * PI)) ! With tau_0 = 2
  !------------------------------------------------------------------------------

    
  if (p.lt.p_0) then
    get_G = (16. / 25.) * (p / p_0)**(3./2.)
  else
    get_G = 1 - (9. / 25.) * (p_0 / p)**(8./3.)
  end if
  
  return
end function get_G

function get_K(p)
! G function (30) of the paper : "A torque formula for non-isothermal Type I planetary migration - II Effects of diffusion"
! By Paardekooper, baruteau and Kley, (2010)
! Equation (30) is a particuliar case of the (29) equation with l=3/4 and tau_0=2
!
! Parameter:
! p : parameter
  implicit none
  real(double_precision), intent(in) :: p
  
  real(double_precision) :: get_K

  !Local
  real(double_precision), parameter :: p_0 = sqrt(28. / (45. * PI)) ! With tau_0 = 7
  !------------------------------------------------------------------------------

  
  if (p.lt.p_0) then
    get_K = (16. / 25.) * (p / p_0)**(3./2.)
  else
    get_K = 1. - (9. / 25.) * (p_0 / p)**(8./3.)
  end if
  
  return
end function get_K

function get_opacity(temperature, sigma)
! subroutine that return the opacity of the disk at the location of the planet given various parameters
  implicit none
  
  ! Inputs 
  real(double_precision), intent(in) :: temperature & ! temperature of the disk [K]
                                        , sigma ! surface density of the disk [MSUN/AU^2]
  
  ! Output
  real(double_precision) :: get_opacity
  
  ! Local
  
  get_opacity = 1.
  write(*,*) "Warning: The opacity is currently not implemented"
  
  !------------------------------------------------------------------------------

  return
end function get_opacity

end module user_module

! TODO utiliser la masse des objets pour ne pas faire le calcul si trop massif, il faut respecter le domaine de validité des formules des couples
! TODO routine générale de conversion des couples en accélération afin de pouvoir réutiliser ailleurs

!TODO : 
!_ amortissement de l'eccentricité
!_ amortissement de l'inclinaison
!_rajouter le spin comme variable d'entrée de la routine mfo_user
