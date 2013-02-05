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
  
  private
  
  public :: mfo_user, unitary_tests ! unitary_tests must be used only to test function, you should not use it out of the debug phase
  
  !------------------------------------------------------------------------------
  ! Default values for parameters that are to be read in the parameter file 'disk.in'
  real(double_precision) :: b_over_h = 0.4 ! the smoothing length for the planet's potential
  real(double_precision) :: adiabatic_index = 1.4 ! the adiabatic index for the gas equation of state
  real(double_precision) :: mean_molecular_weight = 2.35 ! the mean molecular weight in mass of a proton
  
  ! Here we define the power law for surface density sigma(R) = sigma_0 * R^(-sigma_index)
  real(double_precision) :: sigma_0 = 1700 ! the surface density at (R=1AU) [g/cm^2]
  real(double_precision) :: sigma_index = 0.5! the negative slope of the surface density power law (alpha in the paper)
  real(double_precision) :: sigma_0_num ! the surface density at (R=1AU) [Msun/AU^2]
  
  ! Here we define the power law for viscosity viscosity(R) = alpha * c_s * H
  real(double_precision) :: alpha = 0.005 ! alpha parameter for an alpha prescription of the viscosity [No dim]
  
  ! Here we define the power law for temperature T(R) = temperature_0 * R^(-temperature_index)
  real(double_precision) :: temperature_0 = 400. ! the temperature at (R=1AU) [K]
  real(double_precision) :: temperature_index = 1.6! the negativeslope of the temperature power law (beta in the paper)
  
  !------------------------------------------------------------------------------
  ! prefactors
  real(double_precision) :: x_s_prefactor ! prefactor for the half width of the corotation region
  real(double_precision) :: chi_p_prefactor ! prefactor for the thermal diffusivity
  real(double_precision) :: scaleheight_prefactor ! prefactor for the scaleheight
  real(double_precision) :: lindblad_prefactor ! prefactor for the lindblad torque
  real(double_precision) :: migration_acc_prefactor ! prefactor for the migration acceleration
  real(double_precision) :: eccentricity_acc_prefactor ! prefactor for the eccentricity acceleration
  
  real(double_precision) :: torque_hs_baro ! barotropic part of the horseshoe drag
  real(double_precision) :: torque_c_lin_baro ! barotropic part of the linear corotation torque
  
  !------------------------------------------------------------------------------
  ! We define a new type for the properties of the planet
  type PlanetProperties
    ! Properties of the planet
    real(double_precision) :: angular_momentum ! the angular momentum of the planet [Ms.AU^2.day^-1]
    real(double_precision) :: radius ! the radial position of the planet [AU]
    real(double_precision) :: velocity ! the norm of the speed [AU/day]
    real(double_precision) :: omega ! the angular rotation [day-1]
    real(double_precision) :: semi_major_axis ! semi major axis of the planet [AU]
    real(double_precision) :: eccentricity ! the eccentricity of the planet
    real(double_precision) :: inclination ! the inclination of the planet [rad]
    
    ! Properties of the disk at the location of the planet
    real(double_precision) :: sigma ! the surface density of the gas disk at the planet location [MSUN.AU^-2]
    real(double_precision) :: scaleheight ! the scaleheight of the disk at the location of the planet [AU]
    real(double_precision) :: aspect_ratio ! the aspect_ratio of the gas disk at the location of the planet [no dim]
    real(double_precision) :: chi ! the thermal diffusion coefficient at the location of the planet [AU^2.day^-1]
    real(double_precision) :: nu ! the viscosity of the disk at the location of the planet [AU^2.day^-1]
    real(double_precision) :: temperature ! the temperature of the disk at the location of the planet [K] 
    real(double_precision) :: bulk_density ! the bulk_density of the disk at the location of the planet [MSUN.AU^-3]
    real(double_precision) :: opacity ! the opacity of the disk at the location of the planet [?]   
  end type PlanetProperties

  
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
  
  
  real(double_precision) :: torque ! the torque exerted by the disk on the planet [Ms.AU^2]
  real(double_precision) :: corotation_torque ! the corotation torque exerted by the disk on the planet in unit of torque_ref [No dim]
  real(double_precision) :: lindblad_torque ! the lindblad torque exerted by the disk on the planet in unit of torque_ref [No dim]
  real(double_precision) :: torque_ref ! a ref torque that depends on the properties of the planet [Ms.AU^2]
  real(double_precision) :: time_mig ! The migration timescale for the planet [day]
  real(double_precision) :: time_wave ! A timescale for the planet that I don't understand for the moment [day]
  real(double_precision) :: time_ecc ! The eccentricity damping timescale for the planet [day]
  real(double_precision) :: time_inc ! The inclination damping timescale for the planet [day]
  
  !local temporary parameters
  logical, save :: FirstCall = .True.
  type(PlanetProperties) :: p_prop ! various properties of a planet
  real(double_precision) :: e_h ! the ratio between the eccentricity and the aspect ratio for a given planet [no dim]
  real(double_precision) :: i_h ! the ratio between the inclination and the aspect ratio for a given planet [no dim]
  
  real(double_precision), dimension(3) :: migration_acceleration
  real(double_precision), dimension(3) :: eccentricity_acceleration
  real(double_precision) :: inclination_acceleration_z
  
!~   integer, dimension(n_bodies) :: calculation_flag=0 ! The calculation will be done only if this flag is equal to 0. Else, nothing will be done.

!~   
!~   ! modif pour amortissement de e et I
!~   real(double_precision) :: rad, ang, angold
  
  ! counter to check how many time the routine have been called
  integer :: counter = 0
  
  !------------------------------------------------------------------------------
  ! Setup
  
  counter = counter + 1
  
  do planet=1,n_bodies
    acceleration(1,planet) = 0.d0
    acceleration(2,planet) = 0.d0
    acceleration(3,planet) = 0.d0
  end do
  !------------------------------------------------------------------------------
  call init_globals(stellar_mass=mass(1))
  
!~   if (FirstCall) then
!~     FirstCall = .False.
!~     write(*,*) 'Warning: the torque is specified manually'
!~   end if
  
  
  do planet=2,n_big_bodies
    ! because ongoing deletion of planets put their mass to 0 for a few steps, we must check. Else, we will have an error "NaN".
    if (mass(planet).gt.TINY) then
      !------------------------------------------------------------------------------
      ! we store in a structure, various properties of the planet usefull in all the 
      ! subroutine to avoid multiple calculation of the same parameters
      call get_planet_properties(stellar_mass=mass(1), mass=mass(planet), & ! input
      position=position(1:3, planet), velocity=velocity(1:3,planet),& ! input
      p_prop=p_prop) ! Output
      
      !------------------------------------------------------------------------------
      ! prefactor calculation for eccentricity and inclination damping
      e_h = p_prop%eccentricity / p_prop%aspect_ratio
      i_h = p_prop%inclination / p_prop%aspect_ratio
      time_wave = mass(1)**2 * p_prop%aspect_ratio**4 / (mass(planet) * K2 * p_prop%sigma * &
                  p_prop%semi_major_axis**2 * p_prop%omega)
!~       time_wave = abs(time_mig) * p_prop%aspect_ratio**2

      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to migration
      call get_corotation_torque(mass(1), mass(planet), p_prop, & ! input
      corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output

      torque = torque_ref * (lindblad_torque + corotation_torque)
!~       torque = torque_ref * (torque_at_zero / z_c * (p_prop%semi_major_axis - z_c))
      
      
      time_mig = 0.5d0 * p_prop%angular_momentum / torque
      
      migration_acc_prefactor = 1.d0 / time_mig
      
      migration_acceleration(1) = migration_acc_prefactor * velocity(1,planet)
      migration_acceleration(2) = migration_acc_prefactor * velocity(2,planet)
      migration_acceleration(3) = migration_acc_prefactor * velocity(3,planet)
      

      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to eccentricity damping
      
      time_ecc = time_wave / 0.780d0 * (1.d0 - 0.14d0 * e_h**2 + 0.06 * e_h**3 + 0.18 * e_h * i_h**2)
      
      eccentricity_acc_prefactor = -2.d0 * (position(1,planet) * velocity(1,planet) + position(2,planet) * velocity(2,planet) + &
      position(3,planet) * velocity(3,planet)) / (p_prop%radius**2 * time_ecc)
      
      eccentricity_acceleration(1) = eccentricity_acc_prefactor * position(1,planet)
      eccentricity_acceleration(2) = eccentricity_acc_prefactor * position(2,planet)
      eccentricity_acceleration(3) = eccentricity_acc_prefactor * position(3,planet)
      
      !------------------------------------------------------------------------------
      
      ! Calculation of the acceleration due to the inclination damping
      
      time_inc = time_wave / 0.544d0 * (1.d0 - 0.30d0 * i_h**2 + 0.24 * i_h**3 + 0.14 * e_h**2 * i_h)
      
      inclination_acceleration_z = - velocity(3,planet) / time_inc
      
      !------------------------------------------------------------------------------
      ! Calculation of the total acceleration on the planet
      
      acceleration(1,planet) = migration_acceleration(1) + eccentricity_acceleration(1)
      acceleration(2,planet) = migration_acceleration(2) + eccentricity_acceleration(2)
      acceleration(3,planet) = migration_acceleration(3) + eccentricity_acceleration(3) + inclination_acceleration_z
  !~ 
  !~     ! Les commandes en dessous sont la deuxième façon de décrire l'amortissement de e et I (pour limiter les effets sur a)
  !~     !------------------------------------------------------------------------------
  !~     ! prefactor calculation for eccentricity and inclination damping
  !~     e_h = p_prop%eccentricity / p_prop%aspect_ratio
  !~     i_h = p_prop%inclination / p_prop%aspect_ratio
  !~     time_wave = mass(1)**2 * p_prop%aspect_ratio**4 / (mass(planet) * K2 * p_prop%sigma * p_prop%semi_major_axis**2 * p_prop%omega)
  !~     
  !~     !------------------------------------------------------------------------------
  !~     ! Calculation of the acceleration due to eccentricity damping
  !~     
  !~     time_ecc = time_wave / 0.780d0 * (1.d0 - 0.14d0 * e_h**2 + 0.06 * e_h**3 + 0.18 * e_h * i_h**2)
  !~     
  !~     rad = (position(1,planet) * velocity(1,planet) + position(2,planet) * velocity(2,planet)) / p_prop%radius**2
  !~     
  !~     ang = (position(1,planet) * velocity(2,planet) - position(2,planet) * velocity(1,planet)) / p_prop%radius**2
  !~     
  !~     angold = ang - sqrt((mass(1) + mass(planet)) * p_prop%radius) / p_prop%radius**2
  !~     
  !~     eccentricity_acceleration(1) = - (rad * position(1,planet) - angold * position(2,planet)) / time_ecc
  !~     eccentricity_acceleration(2) = - (rad * position(2,planet) + angold * position(1,planet)) / time_ecc
  !~     
  !~     !------------------------------------------------------------------------------
  !~     
  !~     ! Calculation of the acceleration due to the inclination damping
  !~     
  !~     time_inc = time_wave / 0.544d0 * (1.d0 - 0.30d0 * i_h**2 + 0.24 * i_h**3 + 0.14 * e_h**2 * i_h)
  !~     
  !~     inclination_acceleration_z = - velocity(3,planet) / time_inc
  !~     
  !~     !------------------------------------------------------------------------------
  !~     ! Calculation of the total acceleration on the planet
  !~     
  !~     acceleration(1,planet) = migration_acceleration(1) + eccentricity_acceleration(1)
  !~     acceleration(2,planet) = migration_acceleration(2) + eccentricity_acceleration(2)
  !~     acceleration(3,planet) = migration_acceleration(3) + inclination_acceleration_z

!~ !      if (counter.gt.807934) then
!~      if ((counter.gt.78400).and.(counter.lt.78500)) then
!~ !       if (isnan(acceleration(1,planet))) then
!~         open(15, file="leak.out", status="old", access="append")
!~         write(15,*) counter, time, planet, mass(planet)
!~         write(15,*) 'position:', position(1,planet), position(2,planet), position(3,planet)
!~         write(15,*) 'velocity:', velocity(1,planet), velocity(2,planet), velocity(3,planet)
!~         write(15,*) 'acceleration:', acceleration(1,planet), acceleration(2,planet), acceleration(3,planet)
!~         write(15,*) 'times:', time_mig, time_ecc, time_inc
!~         write(15,*) 'torques:', torque, lindblad_torque, corotation_torque
!~         write(15,*) 'angular momentum:', p_prop%angular_momentum, 'radius:', p_prop%radius, 'velocity:', p_prop%velocity
!~         write(15,*) 'omega:', p_prop%omega, 'semi_major_axis:', p_prop%semi_major_axis, 'eccentricity:', p_prop%eccentricity
!~         write(15,*) 'inclination:', p_prop%inclination, 'sigma:', p_prop%sigma, 'scaleheight:', p_prop%scaleheight
!~         write(15,*) 'aspect_ratio:', p_prop%aspect_ratio, 'chi:', p_prop%chi, 'nu:', p_prop%nu
!~         write(15,*) 'temperature:', p_prop%temperature, 'bulk_density:', p_prop%bulk_density, 'opacity:', p_prop%opacity
!~         close(15)
!~       end if
    end if
  end do
  
  
!~   open(15, file="migration_time.out", status="old", access="append")
!~   write(15,*) time, time_mig, time_ecc, time_inc
!~   close(15)
  
  !------------------------------------------------------------------------------
  
!~   if (FirstCall) then
!~     FirstCall = .False.
!~     write (*,*) 'time=', time
!~     write (*,*) 'jcen=', jcen
!~     write (*,*) 'mass=', mass
!~     write (*,*) 'position=', position
!~     write (*,*) 'velocity=', velocity
!~     write (*,*) '--------------------'
!~     write (*,*) 'time_mig=', time_mig
!~     write (*,*) 'angular_momentum=', angular_momentum
!~     write (*,*) 'torque=', torque
!~     write (*,*) 'acceleration=', acceleration
!~   end if
  !------------------------------------------------------------------------------
  return
end subroutine mfo_user

subroutine get_parameter_value(line, isParameter, id, value)
! subroutine that try to split the line in two part, given a separator value (set in parameter of the subroutine)
! The routine return 3 values : 
!
! Return
! isParameter : is a boolean to say whether or not there is a parameter on this line. 
!  i.e if there is an occurence of the separator in the input line
! id : a string that contain the name of the parameter
! value : a string that contains the value(s) associated with the parameter name

  implicit none
  
  ! Input
  character(len=80), intent(in) :: line
  
  ! Output
  logical, intent(out) :: isParameter
  character(len=80), intent(out) :: id, value
  
  ! Local
  character(len=1), parameter :: sep = '=' ! the separator of a parameter line

  integer :: sep_position ! an integer to get the position of the separator

  !------------------------------------------------------------------------------

  sep_position = index(line, sep)
  
  if (sep_position.ne.0) then
    isParameter = .true.
    id = line(1:sep_position-1)
    value = line(sep_position+1:)
  else
    isParameter = .false.
  end if

end subroutine get_parameter_value

subroutine read_disk_properties()
! subroutine that read the 'disk.in' file to retrieve disk properties. Default value exist, if a parameter is not defined

  implicit none
  
  character(len=80) :: line
  character(len=1) :: comment_character = '!' ! character that will indicate that the reste of the line is a comment
  integer :: comment_position ! the index of the comment character on the line. If zero, there is none on the current string
  integer :: error ! to store the state of a read instruction
  
  logical :: isParameter, isDefined
  character(len=80) :: identificator, value
  !------------------------------------------------------------------------------
  
  inquire(file='disk.in', exist=isDefined)
  if (isDefined) then
  
    open(10, file='disk.in', status='old')
    
    do
      read(10, '(a80)', iostat=error), line
      if (error /= 0) exit
        
      ! We get only what is on the left of an eventual comment parameter
        comment_position = index(line, comment_character)
      
      ! If there are comments on the current line, we get rid of them
      if (comment_position.ne.0) then
        line = line(1:comment_position - 1)
      end if
      
      call get_parameter_value(line, isParameter, identificator, value)
        
      if (isParameter) then
        select case(identificator)
        case('b/h')
          read(value, *) b_over_h
        
        case('adiabatic_index')
          read(value, *) adiabatic_index
          
        case('mean_molecular_weight')
          read(value, *) mean_molecular_weight
          
        case('surface_density')
          read(value, *) sigma_0, sigma_index

        case('temperature')
          read(value, *) temperature_0, temperature_index
          
        case('alpha')
          read(value, *) alpha
          
        case default
          write(*,*) 'Warning: An unknown parameter has been found'
          write(*,*) "identificator='", trim(identificator), "' ; value(s)='", trim(value),"'"
        end select
      end if
    end do
    
    close(10)
  else
    write (*,*) 'Warning: The file "disk.in" does not exist. Default values have been used'
  end if
  
  sigma_0_num = sigma_0 * AU**2 / MSUN ! the surface density at (R=1AU) [Msun/AU^2]
  
end subroutine read_disk_properties

subroutine get_planet_properties(stellar_mass, mass, position, velocity, p_prop)

! subroutine that return numerous properties of the planet and its environment given its mass, position and velocity
! Note that some parameters are global and accessed directly by the subroutine

! Parameter:
! stellar_mass : the mass of the central star in [solar mass * K2]
! mass : the planet mass in [solar mass * K2]
! position(3) : position of the planet relatively to the central star [AU]
! velocity(3) : velocity of the planet relatively to the central star [AU/day]

! return : 
! An object of type 'PlanetProperties'. 

! BEWARE : the angular velocity is calculated from the radius of the planet and NOT from his semi major axis.


  use orbital_elements, only : mco_x2ae
  
  implicit none
  real(double_precision), intent(in) :: stellar_mass ! the mass of the central body [Msun * K2]
  real(double_precision), intent(in) :: mass ! the mass of the planet [Msun * K2]
  real(double_precision), dimension(3), intent(in) :: position ! Cartesian position of the planet [AU]
  real(double_precision), dimension(3), intent(in) :: velocity ! Cartesian velocity of the planet [AU/day]
  type(PlanetProperties), intent(out) :: p_prop
  
  ! Local
  real(double_precision) :: gm ! sum of mass (since mass are multiplied implicitely by K2)
  real(double_precision) :: h_p ! the angular momentum given by the calculation of orbital elements. i.e without the mass in it.
  real(double_precision) :: velocity2_p ! the norm of the velocity squared [AU^2 day^-2]

  gm = stellar_mass + mass
  call mco_x2ae(gm,position(1),position(2),position(3),velocity(1),velocity(2),velocity(3),&
                p_prop%semi_major_axis,p_prop%eccentricity,p_prop%inclination,p_prop%radius,velocity2_p,h_p)
  
  
  !------------------------------------------------------------------------------
  p_prop%sigma = sigma_0_num / p_prop%radius ** sigma_index ! [Msun/AU^3]
  p_prop%temperature = temperature_0 / p_prop%radius**temperature_index ! [K]
  
  ! We calculate the angular momentum
  p_prop%angular_momentum = (mass / K2) * h_p  
  p_prop%velocity = sqrt(velocity2_p) ! [AU/day]
  p_prop%omega = sqrt(gm / (p_prop%radius**3)) ! [day-1]
  
  !------------------------------------------------------------------------------
  ! H = sqrt(k_B * T / (omega^2 * mu * m_H))
  p_prop%scaleheight = scaleheight_prefactor * sqrt(p_prop%temperature) / p_prop%omega
!~   p_prop%scaleheight = 0.05 * p_prop%radius
  !------------------------------------------------------------------------------
!~   p_prop%nu = alpha * p_prop%omega * p_prop%scaleheight**2 ! [AU^2.day-1]
  p_prop%nu = 1.d15 * DAY / AU**2
  p_prop%aspect_ratio = p_prop%scaleheight / p_prop%radius
!~     write(*,'(e12.4)') p_prop%nu * AU**2 / DAY 
  !------------------------------------------------------------------------------
  p_prop%bulk_density = 0.5d0 * p_prop%sigma / p_prop%scaleheight
  !------------------------------------------------------------------------------
  p_prop%opacity = get_opacity(p_prop%temperature, p_prop%bulk_density)
  
  p_prop%chi = 1.d-5 * p_prop%radius**2 * p_prop%omega 
!~   p_prop%chi = chi_p_prefactor * p_prop%temperature**4 / (p_prop%opacity * p_prop%sigma**2 * p_prop%omega**2)


  
end subroutine get_planet_properties

subroutine get_corotation_torque(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, Gamma_0)
! function that return the total torque exerted by the disk on the planet 
!

  implicit none
  real(double_precision), intent(in) :: stellar_mass ! the mass of the central body [Msun * K2]
  ! Properties of the planet
  real(double_precision), intent(in) :: mass ! the mass of the planet [Msun * K2]
  type(PlanetProperties), intent(in) :: p_prop ! various properties of the planet
  
  
  real(double_precision), intent(out) :: corotation_torque
  real(double_precision), intent(out) :: lindblad_torque !  lindblad torque exerted by the disk on the planet [\Gamma_0]
  real(double_precision), intent(out) :: Gamma_0 ! canonical torque value [Ms.AU^2](equation (8) of Paardekooper, Baruteau, 2009)


  !Local
  
  ! Meaningless parameters (intermediate constant and so on that doesn't mean anything physically)
  real(double_precision) :: Q_p ! parameter for gamma_eff (equation (45) of Paardekooper, Baruteau, 2010 II. Effects of diffusion)
  real(double_precision) :: k_p ! parameter for p_nu and p_chi for example  !!! This is not 'k' from equation (15)!
  real(double_precision) :: lindblad_parameter ! intermediate of calculation for the lindblad torque. 
  

  
  !Properties of the disk at the location of the planet
  real(double_precision) :: x_s ! semi-width of the horseshoe region [radius_p (in unity of position of the planet)]
  real(double_precision) :: zeta_eff ! effective entropy index depending on gamma_eff [no dim]
  real(double_precision) :: p_nu ! parameter for saturation due to viscosity at the location of the planet [no dim]
  real(double_precision) :: p_chi ! parameter for saturation due to thermal diffusion at the location of the planet [no dim]
  real(double_precision) :: gamma_eff ! effective adiabatic index depending on several parameters [no dim]
  
  !Torques (some depends of the planet)
  real(double_precision) :: torque_hs_ent ! entropy related part of the horseshoe drag
  real(double_precision) :: torque_c_lin_ent ! entropy related part of the linear corotation torque
  
  !------------------------------------------------------------------------------
  ! WE CALCULATE TOTAL TORQUE EXERTED BY THE DISK ON THE PLANET
  Gamma_0 = (mass / (stellar_mass * p_prop%aspect_ratio))**2 * p_prop%sigma * p_prop%radius**4 * p_prop%omega**2
  
  !------------------------------------------------------------------------------
  ! Q is needed by the lindblad torque. We set Q for m ~ 2 /3 h (45): 
  Q_p = TWOTHIRD * p_prop%chi / (p_prop%aspect_ratio * p_prop%scaleheight**2 * p_prop%omega) ! p_prop%aspect_ratio**3 * p_prop%radius**2 = aspect_ratio * scaleheight**2
  !------------------------------------------------------------------------------
  
  gamma_eff = 2.d0 * Q_p * adiabatic_index / (adiabatic_index * Q_p + 0.5d0 * &
  sqrt(2.d0 * sqrt((adiabatic_index * adiabatic_index * Q_p * Q_p + 1.d0)**2 - 16.d0 * Q_p * Q_p * (adiabatic_index - 1.d0)) &
  + 2.d0 * adiabatic_index * adiabatic_index * Q_p * Q_p - 2.d0))
  
  !------------------------------------------------------------------------------
  zeta_eff = temperature_index - (gamma_eff - 1.d0) * sigma_index
  
  x_s = x_s_prefactor / gamma_eff**0.25d0 * sqrt(mass / p_prop%aspect_ratio)
  
  !------------------------------------------------------------------------------
  ! k_p is defined to limit the number of operation and to have a value independant from chi_p or nu_p
  k_p = p_prop%radius * p_prop%radius * p_prop%omega * x_s * x_s * x_s / (2.d0 * PI)
  
  !------------------------------------------------------------------------------
  p_nu = TWOTHIRD * sqrt(k_p / p_prop%nu)
  
  p_chi = sqrt(k_p / p_prop%chi)

  lindblad_parameter = sqrt(p_prop%chi / ( 2.d0 * p_prop%omega * p_prop%scaleheight**2))
!~   lindblad_torque = lindblad_prefactor * (lindblad_parameter + 1.d0 / gamma_eff) / (lindblad_parameter + 1.d0) ! lindblad torque from masset & casoli 2010
  lindblad_torque = lindblad_prefactor / gamma_eff ! lindblad torque formulae from pardekooper, 2010
  torque_hs_ent = 7.9d0 * zeta_eff / gamma_eff
  torque_c_lin_ent = (2.2d0 - 1.4d0 / gamma_eff) * zeta_eff

  
  !------------------------------------------------------------------------------

  corotation_torque = (1.d0 / gamma_eff) * torque_hs_baro * get_F(p_nu) * get_G(p_nu) + torque_c_lin_baro * (1 - get_K(p_nu)) &
    + torque_hs_ent * get_F(p_nu) * get_F(p_chi) * sqrt(get_G(p_nu) * get_G(p_chi)) &
    + torque_c_lin_ent * sqrt((1 - get_K(p_nu)) * (1 - get_K(p_chi)))
  
!~   write (*,*) lindblad_torque, corotation_torque, Gamma_0
!~   write (*,*) 'stellar_mass=',stellar_mass
!~   write (*,*) 'mass=',mass
!~   write (*,*) 'position=',position
!~   write (*,*) 'velocity=',velocity

  return
end subroutine get_corotation_torque

subroutine init_globals(stellar_mass)
! subroutine that initialize global values that define prefactors or values for torque that does not depend on the planet properties
! Parameters
! stellar_mass : the mass of the central object in solar mass
  
  implicit none
  real(double_precision), intent(in) :: stellar_mass
  logical, save :: FirstCall = .True.
  
  if (FirstCall) then
    FirstCall = .False.
    
    call read_disk_properties()
    
    open(15, file="leak.out", status="replace")
    write(15,*) "various data for debuguing that are output when the acceleration get NaN in mfo_user"
    close(15)
!~     open(15, file="migration_time.out", status="replace")
!~     write(15,*) "# The migration time in days, to be compared with real migration time mesured on the plot."
!~     close(15)
    
    x_s_prefactor = 1.1d0 * (b_over_h / 0.4d0)**0.25d0 / sqrt(stellar_mass) ! mass(1) is here for the ratio of mass q

    chi_p_prefactor = (16.d0  / 3.d0) * adiabatic_index * (adiabatic_index - 1.d0) * SIGMA_STEFAN
    
    ! AU is in cm, so we must turn into meter before doing the conversion
    ! division of k_B by m_H is done separately for exponant and value to have more precision
    ! sqrt(k_B/m_H) in numerical units, knowing that [k_B]=[m^2.kg.s^-2K^-1] and [m_H]=[kg]. 
    scaleheight_prefactor = sqrt(1.3806503d0/(1.67262158d0 * mean_molecular_weight) * 1.d4) * DAY / (AU * 1.d-2) 
    
  !~   lindblad_prefactor = -(2.3d0 + 0.4d0 * temperature_index - 0.1d0 * sigma_index) ! masset & casoli 2010
    lindblad_prefactor = -(2.5d0 + 1.7d0 * temperature_index - 0.1d0 * sigma_index) ! paardekooper, baruteau & kley 2010
    
    write(*,*) 'Warning: lindblad torque prefactor need to be calculated at each step if the temperature index change'
    
    torque_hs_baro = 1.1d0 * (1.5d0 - sigma_index)
    torque_c_lin_baro = 0.7d0 * (1.5d0 - sigma_index)
    
    call calculate_temperature_profile()
    
  !~   write(*,*) 'Warning: il y a un offset au couple de corotation'
!~     write(*,*) 'Warning: h est fixé à 0.05'
    write(*,*) 'Warning: nu est fixé à la main à 10^15'
    write(*,*) 'Warning: chi est fixé à la main à 10^-6'
  endif
  
end subroutine init_globals

  function get_F(p)
  ! F function (22) of the paper : "A torque formula for non-isothermal Type I planetary migration - II Effects of diffusion"
  ! By Paardekooper, baruteau and Kley, (2010)
  ! Equation (22) can be approximated within 5% by a much more simpler equation which will be the one we use. 
  !
  ! Parameter:
  ! p : parameter
    implicit none
    real(double_precision), intent(in) :: p
    real(double_precision), parameter :: tmp = 1.d0 / (1.3d0**2)
    
    real(double_precision) :: get_F

    !Local
    !------------------------------------------------------------------------------

    
!~     get_F = 1.d0 / (1.d0 + (p / 1.3d0) * (p / 1.3d0))
    get_F = 1.d0 / (1.d0 + p * p * tmp)
    
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
    real(double_precision), parameter :: p_0 = sqrt(8.d0 / (45.d0 * PI)) ! With tau_0 = 2
    real(double_precision), parameter :: idx_lt = 1.5d0
    real(double_precision), parameter :: idx_gt = 8.d0/3.d0 ! index for the case p greater than p_0
    
    real(double_precision), parameter :: prefactor_lt = 16.d0 / (25.d0 * p_0**idx_lt)
    real(double_precision), parameter :: prefactor_gt = 9.d0 / 25.d0 * p_0**idx_gt
    !------------------------------------------------------------------------------
    ! 16/25 = 0.64
    ! 9/25 = 0.36
      
    if (p.lt.p_0) then
      get_G = prefactor_lt * p**idx_lt
    else
      get_G = 1.d0 - prefactor_gt / p**(idx_gt)
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
    real(double_precision), parameter :: p_0 = sqrt(28.d0 / (45.d0 * PI)) ! With tau_0 = 7
    real(double_precision), parameter :: idx_lt = 1.5d0
    real(double_precision), parameter :: idx_gt = 8.d0/3.d0 ! index for the case p greater than p_0
    
    real(double_precision), parameter :: prefactor_lt = 16.d0 / (25.d0 * p_0**idx_lt)
    real(double_precision), parameter :: prefactor_gt = 9.d0 / 25.d0 * p_0**idx_gt
    !------------------------------------------------------------------------------
    ! 16/25 = 0.64
    ! 9/25 = 0.36
    
    if (p.lt.p_0) then
      get_K = prefactor_lt * p**idx_lt
    else
      get_K = 1.d0 - prefactor_gt / p**(idx_gt)
    end if
    
    return
  end function get_K

  function get_opacity(temperature, num_bulk_density)
  ! subroutine that return the opacity of the disk at the location of the planet given various parameters
    implicit none
    
    ! Inputs 
    real(double_precision), intent(in) :: temperature & ! temperature of the disk [K]
                                          , num_bulk_density ! bulk density of the gas disk [MSUN/AU^3] (in numerical units)
    
    ! Output
    real(double_precision) :: get_opacity
    
    ! Local
    real(double_precision) :: temp34, temp45, temp56, temp67
    real(double_precision) :: bulk_density ! [g/cm^3] (in physical units needed for the expression of the opacity)
    real(double_precision), parameter :: num_to_phys_bulk_density = MSUN / AU**3
    real(double_precision), parameter :: phys_to_num_opacity = MSUN / AU**2

    ! we convert the bulk_density from numerical units(AU, MS, DAY) to physical units (CGS)
    bulk_density = num_to_phys_bulk_density * num_bulk_density
    
    
    ! We get the transition point between the various regimes
    temp34 = 2286.787d0 * bulk_density**(2.d0/49.d0)
    temp45 = 2029.764d0 * bulk_density**(1.d0/81.d0)
    temp56 = 10000.d0 * bulk_density**(1.d0/21.d0)
    temp67 = 30000.d0 * bulk_density**(4.d0/75.d0)
    
    if (temperature.le.166.81d0) then ! regime 1
      get_opacity = 2d-4 * temperature * temperature
    elseif ((temperature.gt.166.81d0).and.(temperature.le.202.677d0)) then ! regime 2
      get_opacity = 2.d16 / temperature**7.d0
    elseif ((temperature.gt.202.677d0).and.(temperature.le.temp34)) then ! regime 3
      get_opacity = 0.1d0 * sqrt(temperature)
    elseif ((temperature.gt.temp34).and.(temperature.le.temp45)) then ! regime 4
      get_opacity = 2.d81 * bulk_density / temperature**24.d0
    elseif ((temperature.gt.temp45).and.(temperature.le.temp56)) then ! regime 5
      get_opacity = 1.d-8 * bulk_density**TWOTHIRD * temperature**3
    elseif ((temperature.gt.temp56).and.(temperature.le.temp67)) then ! regime 6
      get_opacity = 1.d-36 * bulk_density**THIRD * temperature**10
    else ! regime 7
      get_opacity = 1.5d20 * bulk_density / temperature**2.5d0
    endif
    
    ! we change the opacity from physical units to numerical units
    get_opacity = phys_to_num_opacity * get_opacity
    
    !------------------------------------------------------------------------------

    return
  end function get_opacity

  subroutine unitary_tests()
    ! public subroutine to test various function and subroutine of the module user_module. 
    ! For each procedure, the tests should return a data file for values and an associated gnuplot file. 
    ! 
    ! Example : 
    ! To get the plot, you just have to enter in a terminal :
    !> gnuplot "name.gnuplot"
    
    ! For the moment, the function that are tested :
    ! _ get_F
    
    implicit none
    
!~     call test_functions_FGK()
!~     call test_get_opacity()
!~     call test_torques()
!~     call test_torques_fixed_a()
!~     call test_torques_fixed_m()
!~     call test_paardekooper_corotation()
    call test_function_zero_temperature()
    call test_temperature_profile()
    
  end subroutine unitary_tests

  subroutine test_functions_FGK
  ! subroutine that test the functions 'get_F', 'get_G' and 'get_K' and 
  
  ! Return:
  !  a data file 'test_functions_FGK.dat' 
  ! and an associated gnuplot file 'functions_FGK.gnuplot' that display values for get_F, get_G and get_K for a range of p values.
    implicit none
    
    real(double_precision) :: p, f_p, g_p, k_p
    
    real(double_precision), parameter :: p_min = 0.1
    real(double_precision), parameter :: p_max = 100.
    integer, parameter :: nb_points = 100
    real(double_precision), parameter :: p_step = (p_max/p_min) ** (1/(nb_points-1.d0))
    
    integer :: i ! for loops
    
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/test_functions_FGK.dat')
    write(10,*) "Correspond to the figure 2 of II. Effects of diffusion"
    write(10,*) '# p ; F(p) ; G(p) ; K(p)'
    
    
    do i=1,nb_points
      p = p_min * p_step ** (i-1)
      f_p = get_F(p)
      g_p = get_G(p)
      k_p = get_K(p)
      
      write(10,*) p, f_p, g_p, k_p
    end do
    
    close(10)
    
    open(10, file="unitary_tests/functions_FGK.gnuplot")

    write(10,*) 'set xlabel "p"'
    write(10,*) 'set ylabel "function"'
    write(10,*) 'set logscale x'
    write(10,*) 'set mxtics 10'
    write(10,*) 'set grid xtics ytics mxtics'
    write(10,*) 'set xrange [', p_min, ':', p_max, ']'
    write(10,*) 'plot "test_functions_FGK.dat" using 1:2 with lines title "F(p)",\'
    write(10,*) "     '' using 1:3 with lines title 'G(p)',\"
    write(10,*) "     '' using 1:4 with lines title 'K(p)'"
    write(10,*) '#pause -1 # wait until a carriage return is hit'
    write(10,*) 'set terminal pdfcairo enhanced'
    write(10,*) 'set output "functions_FGK.pdf"'
    write(10,*) 'replot' 
    
    close(10)
      
  end subroutine test_functions_FGK

  subroutine test_get_opacity
  ! subroutine that test the function 'get_opacity'
  
  ! Return:
  !  a data file 'test_opacity.dat' 
  ! and an associated gnuplot file 'opacity.gnuplot' that display values for get_opacity for a range of p values.
    implicit none
    
    real(double_precision), dimension(5) :: bulk_density = (/ 1.d-5, 1.d-6, 1.d-7, 1.d-8, 1.d-9/)
    real(double_precision), dimension(5) :: opacity
    real(double_precision) :: temperature
    
    real(double_precision), parameter :: T_min = 80.
    real(double_precision), parameter :: T_max = 4.5d5
    integer, parameter :: nb_points = 400
    real(double_precision), parameter :: T_step = (T_max/T_min) ** (1/(nb_points-1.d0))
    
    integer :: i,j ! for loops
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/test_opacity.dat')
    write(10,*) "Correspond to the figure 9a of Bell & Lin 1994"
    write(10,*) '# Temperature (K) ; Opacity for bulk density from 1e-5 to 1e-9 by power of ten'
    
    do i=1, nb_points
      temperature = T_min * T_step ** (i-1)
      
      do j=1,5
        opacity(j) = get_opacity(temperature, bulk_density(j)) * MSUN / AU**3
      end do
      
      write(10,*) temperature, opacity(1), opacity(2), opacity(3), opacity(4), opacity(5)
    end do
    close(10)
    
    open(10, file="unitary_tests/opacity.gnuplot")
    write(10,*) 'set terminal wxt enhanced'
    write(10,*) 'set xlabel "Temperature T"'
    write(10,*) 'set ylabel "Opacity {/Symbol k}"'
    write(10,*) 'set logscale x'
    write(10,*) 'set logscale y'
    write(10,*) 'set grid'
    write(10,*) 'set xrange [', T_min, ':', T_max, ']'
    write(10,*) "plot 'test_opacity.dat' using 1:2 with lines title '{/Symbol r}=10^{-5}',\"
    write(10,*) "     '' using 1:3 with lines title '{/Symbol r}=10^{-6}',\"
    write(10,*) "     '' using 1:4 with lines title '{/Symbol r}=10^{-7}',\"
    write(10,*) "     '' using 1:5 with lines title '{/Symbol r}=10^{-8}',\"
    write(10,*) "     '' using 1:6 with lines title '{/Symbol r}=10^{-9}'"
    write(10,*) "#pause -1 # wait until a carriage return is hit"
    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) "set output 'opacity.pdf'"
    write(10,*) "replot # pour générer le fichier d'output"  
    
    close(10)
    
  end subroutine test_get_opacity
  
  subroutine test_torques
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    integer, parameter :: nb_mass = 100
    real(double_precision), parameter :: mass_min = 1. * EARTH_MASS
    real(double_precision), parameter :: mass_max = 60. * EARTH_MASS
    real(double_precision), parameter :: mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    
    integer, parameter :: nb_points = 200
    real(double_precision), parameter :: a_min = 0.01
    real(double_precision), parameter :: a_max = 50.
    ! step for log sampling
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_points-1.d0)
    
    real(double_precision) :: a, mass, total_torque, total_torque_units, corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: stellar_mass, position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! stellar mass
    stellar_mass = 1.d0 * K2
    
    call init_globals(stellar_mass)
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/test_corotation_torque.dat')
    open(11, file='unitary_tests/test_total_torque.dat')
    open(12, file='unitary_tests/test_total_torque_units.dat')
    open(13, file='unitary_tests/test_lindblad_torque.dat')
    open(14, file='unitary_tests/test_ref_torque.dat')
    
    
    write(10,*) '# semi major axis (AU) ; mass in earth mass ; corotation torque (no dim)'
    write(11,*) '# semi major axis (AU) ; mass in earth mass ; total torque (no dim)'
    write(12,*) '# semi major axis (AU) ; mass in earth mass ; total torque in M_s.AU^2.day^{-2}'
    write(13,*) '# semi major axis (AU) ; mass in earth mass ; lindblad torque (no dim)'
    write(14,*) '# semi major axis (AU) ; mass in earth mass ; reference torque in M_s.AU^2.day^{-2}'
    
    
    do i=1, nb_points ! loop on the position
      a = a_min + a_step * (i-1)
      
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a
      
      
      do j=1,nb_mass
        mass = (mass_min + mass_step * (j - 1.d0)) * K2
        
        ! We generate cartesian coordinate for the given mass and semi major axis
        velocity(2) = sqrt((stellar_mass + mass) / position(1))
        
        ! we store in global parameters various properties of the planet
        call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! input
         p_prop=p_prop) ! Output
        call get_corotation_torque(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        total_torque = lindblad_torque + corotation_torque
        total_torque_units = torque_ref * total_torque
        
                
        write(10,*) a, mass / (EARTH_MASS*K2), corotation_torque
        write(11,*) a, mass / (EARTH_MASS*K2), total_torque
        write(12,*) a, mass / (EARTH_MASS*K2), total_torque_units
        write(13,*) a, mass / (EARTH_MASS*K2), lindblad_torque
        write(14,*) a, mass / (EARTH_MASS*K2), torque_ref
        
      end do
      
      do j=10,14
        write(j,*) ""! we write a blank line to separate them in the data file, else, gnuplot doesn't want to make the surface plot
      end do
    end do
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    
    
    open(10, file="unitary_tests/corotation_torque.gnuplot")
    open(11, file="unitary_tests/total_torque.gnuplot")
    open(12, file="unitary_tests/total_torque_units.gnuplot")
    open(13, file="unitary_tests/lindblad_torque.gnuplot")
    open(14, file="unitary_tests/ref_torque.gnuplot")
    
    
    ! --------------------------------------------
    ! Part for contour for only a part of the plots. We want an isoline for zero torque
    do j=11,12
      write(j,*) 'set contour base; set cntrparam levels discret 0.'
      write(j,*) 'unset surface'
      write(j,*) 'set table "contour.dat"'
      write(j,*) 'set dgrid3d 30,30,10'
    end do
    
    write(11,*) "splot 'test_total_torque.dat'"
    write(12,*) "splot 'test_total_torque_units.dat'"
    
    do j=11,12
      write(j,*) 'unset table'

      write(j,*) '# Draw the plot'
      write(j,*) 'reset'
    end do
    ! --------------------------------------------
    
    do j=10,14
      write(j,*) 'set terminal wxt enhanced'
      write(j,*) 'set xlabel "semi major axis (AU)"'
      write(j,*) 'set ylabel "Planet mass (m_{earth})" center'
    end do
    
    write(10,*) 'set title "Evolution of the corotation torque {/Symbol G}_c/{/Symbol G}_0"'
    write(11,*) 'set title "Evolution of the total torque {/Symbol G}_{tot}/{/Symbol G}_0 "'
    write(13,*) 'set title "Evolution of the lindblad torque {/Symbol G}_L/{/Symbol G}_0 "'
    write(14,*) 'set title "Evolution of the reference torque {/Symbol G}_0 [M_s.AU^2.day^{-2}]"'
    write(12,*) 'set title "Evolution of the total torque {/Symbol G}_{tot} [M_s.AU^2.day^{-2}]"'
    
    do j=10,14
      write(j,*) 'set pm3d map'
      write(j,*) 'set pm3d explicit'
      write(j,*) 'set palette rgbformulae 22,13,-31'
      write(j,*) 'set mxtics 5'
      write(j,*) 'set mytics 5'
      write(j,*) 'set grid xtics ytics mxtics mytics linetype -1, linetype 0'
      write(j,*) 'set xrange [', a_min, ':', a_max, ']'
      write(j,*) 'set yrange [', mass_min / EARTH_MASS, ':', mass_max / EARTH_MASS, ']'
    end do

    write(10,*) "splot 'test_corotation_torque.dat' with pm3d notitle"
    write(11,*) "splot 'test_total_torque.dat' with pm3d notitle, 'contour.dat' with line linetype -1 title '{/Symbol G}=0'"
    write(12,*) "splot 'test_total_torque_units.dat' with pm3d notitle, 'contour.dat' with line linetype -1 title '{/Symbol G}=0'"
    write(13,*) "splot 'test_lindblad_torque.dat' with pm3d notitle"
    write(14,*) "splot 'test_ref_torque.dat' with pm3d notitle"

    
    do j=10,14
      write(j,*) "#pause -1 # wait until a carriage return is hit"
      write(j,*) "set terminal pngcairo enhanced"
    end do
    
    write(10,*) "set output 'corotation_torque.png'"
    write(11,*) "set output 'total_torque.png'"
    write(12,*) "set output 'total_torque_units.png'"
    write(13,*) "set output 'lindblad_torque.png'"
    write(14,*) "set output 'ref_torque.png'"
        
    do j=10,14
      write(j,*) "replot # pour générer le fichier d'output"
    end do
    
    do j=11,12
      write(j,*) '!rm contour.dat' ! to delete the temporary file
    end do
    
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    
  end subroutine test_torques

  subroutine test_torques_fixed_a
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    integer, parameter :: nb_mass = 100
    real(double_precision), parameter :: mass_min = 1. ! in earth mass
    real(double_precision), parameter :: mass_max = 75. ! in earth mass
    real(double_precision), parameter :: mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    
    real(double_precision), parameter :: a = 5.2
    
    real(double_precision) :: mass, total_torque, corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: lindblad_torque_units, corotation_torque_units, total_torque_units
    real(double_precision) :: stellar_mass, position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! stellar mass
    stellar_mass = 1.d0 * K2
    
    call init_globals(stellar_mass)
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/test_torques_fixed_a.dat')
    open(11, file='unitary_tests/test_ref_torque_fixed_a.dat')
    open(12, file='unitary_tests/test_torques_fixed_a_units.dat')
    open(13, file='unitary_tests/test_specific_torque_fixed_a.dat')
    
    write(10,*) '# mass in earth mass ; corotation torque (no dim), lindblad torque (no dim), total torque (no dim)'
    write(11,*) '# mass in earth mass ; reference torque in M_s.AU^2.day^{-2}'
    write(12,*) '# mass in earth mass ; corotation torque (M_s.AU^2.day^{-2}), lindblad torque (M_s.AU^2.day^{-2})&
                 &, total torque (M_s.AU^2.day^{-2})'
    write(13,*) '# mass in earth mass ; specific torque (a_jup^2.Omega_jup^2)'             

    ! We generate cartesian coordinate for the given semi major axis
    position(1) = a
    
    
    do j=1,nb_mass
      mass = (mass_min + mass_step * (j - 1.d0)) * EARTH_MASS  * K2
      
      ! We generate cartesian coordinate for the given mass and semi major axis
      velocity(2) = sqrt((stellar_mass + mass) / position(1))
      
      ! we store in global parameters various properties of the planet
      call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! input
       p_prop=p_prop) ! Output
       
      call get_corotation_torque(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
      
      total_torque = lindblad_torque + corotation_torque

      lindblad_torque_units = torque_ref * lindblad_torque
      corotation_torque_units = torque_ref * corotation_torque
      total_torque_units = torque_ref * total_torque
                    
      write(10,*) mass / (EARTH_MASS * K2), corotation_torque, lindblad_torque, total_torque
      write(11,*) mass / (EARTH_MASS * K2), torque_ref
      write(12,*) mass / (EARTH_MASS * K2), corotation_torque_units, lindblad_torque_units, total_torque_units
      write(13,*) mass / (EARTH_MASS * K2), total_torque_units * (K2 / (5.2d0 * mass * p_prop%omega))**2
    end do
    
    close(10)
    close(11)
    close(12)
    close(13)
    
    
    open(10, file="unitary_tests/torques_fixed_a.gnuplot")
    open(11, file="unitary_tests/ref_torque_fixed_a.gnuplot")
    open(12, file="unitary_tests/torques_fixed_a_units.gnuplot")
    open(13, file="unitary_tests/specific_torque_fixed_a.gnuplot")

    
    do j=10,13
      write(j,*) 'set terminal wxt enhanced'
      write(j,*) 'set xlabel "Planet mass (m_{earth})"'
      write(j,*) 'set title "for a semi major axis a=',a,' AU"'
    end do
    
    write(10,*) 'set ylabel "torque [{/Symbol G}_0]"'
    
    write(11,*) 'set ylabel "reference torque {/Symbol G}_0 [M_s.AU^2.day^{-2}]"'
    write(11,*) 'set nokey'
    
    write(12,*) 'set ylabel "torque [M_s.AU^2.day^{-2}]"'
    
    write(13,*) 'set ylabel "specific torque [a_{jup}^2.{/Symbol W}_{jup}^2]"'
    write(13,*) 'set nokey'
        
    do j=10,13
!~     write(j,*) 'set xrange [', a_min, ':', a_max, '] noreverse'
      write(j,*) 'set grid'
      write(j,*) 'set xrange [', mass_min, ':', mass_max, ']'
    end do

    write(10,*) "plot 'test_torques_fixed_a.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(10,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(10,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"
    
    write(11,*) "plot 'test_ref_torque_fixed_a.dat' using 1:2 with lines"
    
    write(12,*) "plot 'test_torques_fixed_a_units.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(12,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(12,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"
    
    write(13,*) "plot 'test_specific_torque_fixed_a.dat' using 1:2 with lines"
    
    do j=10,13
      write(j,*) "#pause -1 # wait until a carriage return is hit"
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) "set output 'torques_fixed_a.pdf'"
    write(11,*) "set output 'ref_torque_fixed_a.pdf'"
    write(12,*) "set output 'torques_fixed_a_units.pdf'"
    write(13,*) "set output 'specific_torque_fixed_a.pdf'"
    
    do j=10,13
      write(j,*) "replot # pour générer le fichier d'output"
    end do
    
    close(10)
    close(11)
    close(12)
    close(13)

    
  end subroutine test_torques_fixed_a

  subroutine test_torques_fixed_m
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    integer, parameter :: nb_a = 400
    real(double_precision), parameter :: a_min = 0.1 ! in AU
    real(double_precision), parameter :: a_max = 60. ! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision), parameter :: mass = 20. * EARTH_MASS * K2
    
    real(double_precision) :: a, total_torque, corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: lindblad_torque_units, corotation_torque_units, total_torque_units
    real(double_precision) :: stellar_mass, position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! stellar mass
    stellar_mass = 1.d0 * K2
    
    call init_globals(stellar_mass)
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/test_torques_fixed_m.dat')
    open(11, file='unitary_tests/test_ref_torque_fixed_m.dat')
    open(12, file='unitary_tests/test_torques_fixed_m_units.dat')
    
    write(10,*) '# a in AU ; corotation torque (no dim), lindblad torque (no dim), total torque (no dim)'
    write(11,*) '# a in AU ; reference torque in M_s.AU^2.day^{-2}'
    write(12,*) '# a in AU ; corotation torque (M_s.AU^2.day^{-2}), lindblad torque (M_s.AU^2.day^{-2})&
                 &, total torque (M_s.AU^2.day^{-2})'

    ! We generate cartesian coordinate for the given semi major axis
    position(1) = a
    
    
    do j=1,nb_a
      a = (a_min + a_step * (j - 1.d0))
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a
      
      ! We generate cartesian coordinate for the given mass and semi major axis
      velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
      
      ! we store in global parameters various properties of the planet
      call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
       p_prop=p_prop) ! Output
      
      call get_corotation_torque(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
      
      total_torque = lindblad_torque + corotation_torque
      
      lindblad_torque_units = torque_ref * lindblad_torque
      corotation_torque_units = torque_ref * corotation_torque
      total_torque_units = torque_ref * total_torque
      
              
      write(10,*) a, corotation_torque, lindblad_torque, total_torque
      write(11,*) a, torque_ref
      write(12,*) a, corotation_torque_units, lindblad_torque_units, total_torque_units
    end do
    
    close(10)
    close(11)
    close(12)
    
    
    open(10, file="unitary_tests/torques_fixed_m.gnuplot")
    open(11, file="unitary_tests/ref_torque_fixed_m.gnuplot")
    open(12, file="unitary_tests/torques_fixed_m_units.gnuplot")
    
    do j=10,12
      write(j,*) 'set terminal wxt enhanced'
      write(j,*) 'set xlabel "semi major axis a (in AU)"'
      write(j,*) 'set title "for planet mass =',mass / (EARTH_MASS * K2),' m_{earth}"'
    end do
    
    write(10,*) 'set ylabel "torque [{/Symbol G}_0]"'
    
    write(11,*) 'set ylabel "reference torque {/Symbol G}_0 [M_s.AU^2.day^{-2}]"'
    write(11,*) 'set nokey'
    
    write(12,*) 'set ylabel "torque [M_s.AU^2.day^{-2}]"'
    
    do j=10,12
!~     write(j,*) 'set xrange [', a_min, ':', a_max, '] noreverse'
      write(j,*) 'set grid'
      write(j,*) 'set xrange [', a_min, ':', a_max, ']'
    end do

    write(10,*) "plot 'test_torques_fixed_m.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(10,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(10,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"
    
    write(11,*) "plot 'test_ref_torque_fixed_m.dat' using 1:2 with lines"
        
    write(12,*) "plot 'test_torques_fixed_m_units.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(12,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(12,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"

    
    do j=10,12
      write(j,*) "#pause -1 # wait until a carriage return is hit"
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) "set output 'torques_fixed_m.pdf'"
    write(11,*) "set output 'ref_torque_fixed_m.pdf'"
    write(12,*) "set output 'gamma_eff_fixed_m.pdf'"
    
    do j=10,12
      write(j,*) "replot # pour générer le fichier d'output"
    end do
    
    close(10)
    close(11)
    close(12)

    
  end subroutine test_torques_fixed_m
  
  subroutine calculate_temperature_profile()
! Subroutine that test the finding of the temperature profile and store a plot of the temperature profile of the disk
! A gnuplot file and a data file are created to display the temperature profile.
! TODO calculer aussi l'exposant local du profil de température

    implicit none
    
    integer, parameter :: nb_a = 400
    real(double_precision), parameter :: a_min = 1 ! in AU
    real(double_precision), parameter :: a_max = 60. ! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision), parameter :: mass = 20. * EARTH_MASS * K2
    
    real(double_precision) :: a
    real(double_precision) :: stellar_mass, position(3), velocity(3), temperature, exponant
    type(PlanetProperties) :: p_prop
    ! value for the precedent step of the loop. In order to calculate the index of the local temperature power law.
    real(double_precision) :: a_old, temperature_old 
    
    integer :: i,j ! for loops
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! stellar mass
    stellar_mass = 1.d0 * K2  
  
    ! We open the file where we want to write the outputs
    open(10, file='temperature_profile.dat', status='replace')
    
    write(10,*) '# exponant ; a in AU ; temperature (in K)'    
    ! We generate cartesian coordinate for the given semi major axis
    position(1) = a
    
    a = 0.9d0 * a_min
    ! We generate cartesian coordinate for the given semi major axis
    position(1) = a
    
    ! We generate cartesian coordinate for the given mass and semi major axis
    velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
    
    ! we store in global parameters various properties of the planet
    call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
     p_prop=p_prop) ! Output
    
    temperature = zero_finding_zbrent(x_min=1.d-5, x_max=1.d4, tolerance=1d-4, p_prop=p_prop)

    do j=1,nb_a
      a_old = a
      temperature_old = temperature
      
      a = (a_min + a_step * (j - 1.d0))
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a
      
      ! We generate cartesian coordinate for the given mass and semi major axis
      velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
      
      ! we store in global parameters various properties of the planet
      call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
       p_prop=p_prop) ! Output
      
      temperature = zero_finding_zbrent(x_min=1.d-5, x_max=1.d4, tolerance=1d-4, p_prop=p_prop)
      
      exponant = - (log(temperature) - log(temperature_old)) / (log(a) - log(a_old))
      
      write(10,*) log(a), log(temperature), exponant, a, temperature
    end do
    
    close(10)
  end subroutine calculate_temperature_profile
  
  subroutine test_temperature_profile()
! Subroutine that test the finding of the temperature profile and store a plot of the temperature profile of the disk
! A gnuplot file and a data file are created to display the temperature profile.
! TODO calculer aussi l'exposant local du profil de température

    implicit none
    
    integer, parameter :: nb_a = 400
    real(double_precision), parameter :: a_min = 1 ! in AU
    real(double_precision), parameter :: a_max = 60. ! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision), parameter :: mass = 20. * EARTH_MASS * K2
    
    real(double_precision) :: a
    real(double_precision) :: stellar_mass, position(3), velocity(3), temperature, exponant
    type(PlanetProperties) :: p_prop
    ! value for the precedent step of the loop. In order to calculate the index of the local temperature power law.
    real(double_precision) :: a_old, temperature_old 
    
    integer :: i,j ! for loops
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! stellar mass
    stellar_mass = 1.d0 * K2
    
    call init_globals(stellar_mass)
    
    
    open(10, file="unitary_tests/temperature_profile.gnuplot")
    open(11, file="unitary_tests/temperature_index.gnuplot")
    
    do j=10,11
      write(j,*) 'set terminal wxt enhanced'
      write(j,*) 'set xlabel "semi major axis a (in AU)"'
      write(j,*) 'set nokey'
    end do
    
    write(10,*) 'set ylabel "Temperature [K]"'
    
    write(11,*) 'set ylabel "Temperature law index"'
    
    do j=10,11
      write(j,*) 'set grid'
      write(j,*) 'cd ".."'
!~       write(j,*) 'set xrange [', a_min, ':', a_max, ']'
    end do

    write(10,*) "plot 'temperature_profile.dat' using 4:5 with lines notitle"
    
    write(11,*) "plot 'temperature_profile.dat' using 4:3 with lines notitle"
    

    
    do j=10,11
      write(j,*) "#pause -1 # wait until a carriage return is hit"
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) '!rm "unitary_tests/temperature_profile.pdf"'
    write(10,*) "set output 'unitary_tests/temperature_profile.pdf'"
    write(11,*) '!rm "unitary_tests/temperature_index.pdf"'
    write(11,*) "set output 'unitary_tests/temperature_index.pdf'"
    
    
    do j=10,11
      write(j,*) "replot # pour générer le fichier d'output"
    end do
    
    close(10)
    close(11)
  
  end subroutine test_temperature_profile

  subroutine test_function_zero_temperature
  ! subroutine that test the function 'zero_finding_temperature'
  
  ! Return:
  !  a data file and an associated gnuplot file.
    implicit none
    
    real(double_precision) :: temperature
    
    real(double_precision), parameter :: T_min = 0.
    real(double_precision), parameter :: T_max = 1000.
    integer, parameter :: nb_points = 2000
!~     real(double_precision), parameter :: T_step = (T_max/T_min) ** (1/(nb_points-1.d0))
    real(double_precision), parameter :: T_step = (T_max - T_min) / (nb_points - 1.d0)
    
    real(double_precision), parameter :: mass = 20. * EARTH_MASS * K2
    
    real(double_precision) :: zero_function ! value that we want to output
    
    integer :: i,j ! for loops
    real(double_precision) :: stellar_mass, position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    real(double_precision) :: prefactor ! prefactor for the calculation of the function of the temperature whose zeros are searched

  !------------------------------------------------------------------------------


    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! stellar mass
    stellar_mass = 1.d0 * K2
    
    call init_globals(stellar_mass)
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/test_function_zero_temperature.dat')
    
    ! We generate cartesian coordinate for the given semi major axis
    position(1) = 10.d0
    
    ! We generate cartesian coordinate for the given mass and semi major axis
    velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
    
    ! we store in global parameters various properties of the planet
    call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
     p_prop=p_prop) ! Output
    
    ! We calculate this value outside the function because we only have to do this once per step (per radial position)
    prefactor = - (9.d0 * p_prop%nu * p_prop%sigma * p_prop%omega**2 / 16.d0)
    
    write(10,*) '# properties of the disk at the location of the planet hat influence the value of the temperature'
    write(10,*) '# radial position of the planet (in AU) :', p_prop%radius
    write(10,*) '# viscosity :', p_prop%nu
    write(10,*) '# surface density :', p_prop%sigma
    write(10,*) '# bulk density :', p_prop%bulk_density
    write(10,*) '# angular velocity :', P_prop%omega
    write(10,*) '# Temperature (K) ; value of the function. The right temperature is when the function is 0'
    
    do i=1, nb_points
!~       temperature = T_min * T_step ** (i-1)
      temperature = (T_min + T_step * (i - 1.d0))
      
      zero_function = zero_finding_temperature(temperature=temperature, rho=p_prop%bulk_density, &
                                               omega=p_prop%omega, prefactor=prefactor)
      
      write(10,*) temperature, zero_function
    end do
    close(10)
    
    open(10, file="unitary_tests/function_zero_temperature.gnuplot")
    write(10,*) 'set terminal wxt enhanced'
    write(10,*) 'set xlabel "Temperature T"'
    write(10,*) 'set ylabel "zero function"'
    write(10,*) 'set grid'
    write(10,*) 'set xrange [', T_min, ':', T_max, ']'
    write(10,*) "plot 'test_function_zero_temperature.dat' using 1:2 with lines notitle"
    write(10,*) "#pause -1 # wait until a carriage return is hit"
    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) "set output 'function_zero_temperature.pdf'"
    write(10,*) "replot # pour gÃ©nÃ©rer le fichier d'output"  
    
    close(10)
    
  end subroutine test_function_zero_temperature
  
subroutine paardekooper_corotation(stellar_mass, mass, position, velocity)
! function that return the total torque exerted by the disk on the planet 
!
! Parameter:
! mass : an array will all the masses (including the central one as the first element [solar mass]
! p_prop%radius : the distance between the central body and each planet
  use orbital_elements, only : mco_x2a
  implicit none
  real(double_precision), intent(in) :: stellar_mass, mass, position(3), velocity(3)

  !Local
  
  ! Meaningless parameters (intermediate constant and so on that doesn't mean anything physically)
  real(double_precision) :: Q_p ! parameter for gamma_eff (equation (45) of Paardekooper, Baruteau, 2010 II. Effects of diffusion)
  real(double_precision) :: k_p ! parameter for p_nu and p_chi for example  !!! This is not 'k' from equation (15)!
  
  !Properties of the disk at the location of the planet
  type(PlanetProperties) :: p_prop

  real(double_precision) :: x_s ! semi-width of the horseshoe region [radius_p (in unity of position of the planet)]
  real(double_precision) :: zeta_eff ! effective entropy index depending on gamma_eff [no dim]
  real(double_precision) :: p_nu ! parameter for saturation due to viscosity at the location of the planet [no dim]
  real(double_precision) :: p_chi ! parameter for saturation due to thermal diffusion at the location of the planet [no dim]
  
  !Torques (some depends of the planet)
  real(double_precision) :: torque_hs_ent ! entropy related part of the horseshoe drag
  real(double_precision) :: torque_c_lin_ent ! entropy related part of the linear corotation torque

  real(double_precision), parameter :: p_nu_min = 0.1d0
  real(double_precision), parameter :: p_nu_max = 20.d0
  integer, parameter :: nb_points = 400
  real(double_precision), parameter :: p_nu_step = (p_nu_max/p_nu_min) ** (1/(nb_points-1.d0))
  
  real(double_precision), parameter, dimension(4) :: list_chi = (/1.d-5, 2.d-6, 1.d-6, 1.d-7/)
  real(double_precision) :: corotation_torque(4)
  real(double_precision) :: Gamma_0 ! canonical torque value [Ms.AU^2](equation (8) of Paardekooper, Baruteau, 2009)
  real(double_precision) :: gamma_eff ! effective adiabatic index depending on several parameters [no dim]
  
  integer :: i,j ! for loops
  
  
  ! we store in global parameters various properties of the planet
  call get_planet_properties(stellar_mass=stellar_mass, mass=mass, & ! input
  position=position(1:3), velocity=velocity(1:3),& ! input
  p_prop=p_prop) ! Output
  
  !------------------------------------------------------------------------------
  ! WE CALCULATE TOTAL TORQUE EXERTED BY THE DISK ON THE PLANET
  Gamma_0 = (mass / (stellar_mass * p_prop%aspect_ratio))**2 * p_prop%sigma * p_prop%radius**4 * p_prop%omega**2
  
  open(10, file='unitary_tests/test_paardekooper_figure6.dat')
  write(10,*) "Correspond to the figure 6 of Paardekooper 2010"
  write(10,*) '# p_nu ; corotation torque for chi=(1e-5, 2e-6, 1e-6, 1e-7)'
  
  do i=1,nb_points
    do j=1,4 ! loop on chi
      p_prop%chi = list_chi(j) * p_prop%radius**2 * p_prop%omega
      
      !------------------------------------------------------------------------------
      ! Q is needed by the lindblad torque. We set Q for m ~ 2 /3 h (45): 
      Q_p = TWOTHIRD * p_prop%chi / (p_prop%aspect_ratio * p_prop%scaleheight**2 * p_prop%omega) ! p_prop%aspect_ratio**3 * p_prop%radius**2 = aspect_ratio * scaleheight**2
      !------------------------------------------------------------------------------
      
      gamma_eff = 2.d0 * Q_p * adiabatic_index / (adiabatic_index * Q_p + 0.5d0 * &
      sqrt(2.d0 * sqrt((adiabatic_index * adiabatic_index * Q_p * Q_p + 1.d0)**2 - 16.d0 * Q_p * Q_p * (adiabatic_index - 1.d0)) &
      + 2.d0 * adiabatic_index * adiabatic_index * Q_p * Q_p - 2.d0))
      
      !------------------------------------------------------------------------------
      zeta_eff = temperature_index - (gamma_eff - 1.d0) * sigma_index
      
      x_s = x_s_prefactor / gamma_eff**0.25d0 * sqrt(mass / p_prop%aspect_ratio)
      
      !------------------------------------------------------------------------------
      ! k_p is defined to limit the number of operation and to have a value independant from chi_p or nu_p
      k_p = p_prop%radius * p_prop%radius * p_prop%omega * x_s * x_s * x_s / (2.d0 * PI)
      
      !------------------------------------------------------------------------------
    !~   p_nu = TWOTHIRD * sqrt(k_p / p_prop%nu)
      p_nu = p_nu_min * p_nu_step ** (i-1)
      
      p_chi = sqrt(k_p / p_prop%chi)
      
      torque_hs_ent = 7.9d0 * zeta_eff / gamma_eff
      torque_c_lin_ent = (2.2d0 - 1.4d0 / gamma_eff) * zeta_eff
      
      !------------------------------------------------------------------------------

      corotation_torque(j) = (1.d0 / gamma_eff) * &
        (torque_hs_baro * get_F(p_nu) * get_G(p_nu) + torque_c_lin_baro * (1 - get_K(p_nu)) &
        + torque_hs_ent * get_F(p_nu) * get_F(p_chi) * sqrt(get_G(p_nu) * get_G(p_chi)) &
        + torque_c_lin_ent * sqrt((1 - get_K(p_nu)) * (1 - get_K(p_chi))))
    end do
    write (10,*) p_nu, corotation_torque
  end do
  close(10)
  
  open(10, file="unitary_tests/paardekooper_figure6.gnuplot")
  write(10,*) 'set terminal wxt enhanced'
  write(10,*) 'set xlabel "p_{/Symbol n}"'
  write(10,*) 'set ylabel "{/Symbol G}_c/{/Symbol G}_0"'
  write(10,*) 'set logscale x'
  write(10,*) 'set mxtics 10'
  write(10,*) 'set grid xtics ytics mxtics'
  write(10,*) 'set xrange [', p_nu_min, ':', p_nu_max, ']'
  write(10,*) "plot 'test_paardekooper_figure6.dat' using 1:2 with lines title '{/Symbol c}_p=10^{-5}',\"
  write(10,*) "     '' using 1:3 with lines title '{/Symbol c}_p=2*10^{-6}',\"
  write(10,*) "     '' using 1:4 with lines title '{/Symbol c}_p=10^{-6}',\"
  write(10,*) "     '' using 1:5 with lines title '{/Symbol c}_p=10^{-7}'"
  write(10,*) "#pause -1 # wait until a carriage return is hit"
  write(10,*) "set terminal pdfcairo enhanced"
  write(10,*) "set output 'paardekooper_figure6.pdf'"
  write(10,*) "replot # pour générer le fichier d'output"  

  close(10)
  
  return
end subroutine paardekooper_corotation

  subroutine test_paardekooper_corotation
  ! subroutine that try to reproduce the figure 6 of paardekooper, 2010
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    
    real(double_precision) :: a, mass
    real(double_precision) :: stellar_mass, position(3), velocity(3)
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! stellar mass
    stellar_mass = 1.d0 * K2
    
    mass = 1.26d-5  * K2
    
    call init_globals(stellar_mass)
    
    a = 5.d0
    position(1) = a
    
    velocity(2) = sqrt((stellar_mass + mass) / position(1))
    
    
    call paardekooper_corotation(stellar_mass, mass, position, velocity)
    
  end subroutine test_paardekooper_corotation

function zero_finding_zbrent(x_min, x_max, tolerance, p_prop)
! Using Brent's method, find the root of a function 'func' known to lie between 'x_min' and 'x_max'. 
! The root, returned as 'zero_finding_zbrent', will be refined until its accuray is 'tolerance'. 

! Parameters :
! ITMAX : maximum allowed number of iterations
! EPS : machine floating-point precision.

! REMARK : This function is based on the zbrent function in fortran 90 of numerical recipes

! Output
real(double_precision) :: zero_finding_zbrent

! Input 
real(double_precision), intent(in) :: tolerance, x_min, x_max
type(PlanetProperties), intent(in) :: p_prop ! various properties of a planet

! Parameters
! the routine zbrent works best when PES is exactly the machine precision. 
! The fortran 90 intrinsic function epsilon allows us to code this in a portable fashion.
real(double_precision), parameter :: EPS=epsilon(x_min) 

integer, parameter :: ITMAX=100

! Locals
integer :: iter
real(double_precision) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
real(double_precision) :: prefactor ! prefactor for the calculation of the function of the temperature whose zeros are searched

!------------------------------------------------------------------------------

! We calculate this value outside the function because we only have to do this once per step (per radial position)
prefactor = - (9.d0 * p_prop%nu * p_prop%sigma * p_prop%omega**2 / 16.d0)

a = x_min
b = x_max
fa = zero_finding_temperature(temperature=a, rho=p_prop%bulk_density, omega=p_prop%omega, prefactor=prefactor)
fb = zero_finding_temperature(temperature=b, rho=p_prop%bulk_density, omega=p_prop%omega, prefactor=prefactor)

if (((fa.gt.0.).and.(fb.gt.0.)).or.((fa.lt.0.).and.(fb.lt.0.))) then
  write(*,*) 'root must be bracketed for zbrent'
  write(*,*) '  T_min =', x_min, 'f(T_min) =', fa
  write(*,*) '  T_max =', x_max, 'f(T_max) =', fb
  write(*,*) 'properties of the disk at the location of the planet hat influence the value of the temperature'
  write(*,*) '  radial position of the planet (in AU) :', p_prop%radius
  write(*,*) '  viscosity :', p_prop%nu
  write(*,*) '  surface density :', p_prop%sigma
  write(*,*) '  bulk density :', p_prop%bulk_density
  write(*,*) '  angular velocity :', P_prop%omega
  stop
endif

c = b
fc = fb
do iter=1,ITMAX
  if (((fb.gt.0.).and.(fc.gt.0.)).or.((fb.lt.0.).and.(fc.lt.0.))) then
    ! rename a, b, c and adjust bouding interval d.
    c = a
    fc = fa
    d = b - a
    e = d
  endif
  
  if (abs(fc).lt.abs(fb)) then
    a = b
    b = c
    c = a
    fa = fb
    fb = fc
    fc = fa
  endif
  
  ! convergence check
  tol1 = 2. * EPS * abs(b) + 0.5 * tolerance
  xm = .5 * (c - b)
  
  if (abs(xm).le.tol1 .or. fb.eq.0.) then
    zero_finding_zbrent = b
    return
  endif
  
  if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
    ! attempt inverse quadratic interpolation
    s = fb / fa
    if(a.eq.c) then
      p = 2. * xm * s
      q = 1. - s
    else
      q = fa / fc
      r = fb / fc
      p = s * (2. * xm * q * (q - r) - (b - a) * (r - 1.))
      q = (q - 1.) * (r - 1.) * (s - 1.)
    endif
    
    if(p.gt.0.) q=-q ! check whether in bounds
    p=abs(p)
    if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
      ! accept interpolation
      e = d
      d = p / q
    else
      ! interpolation failed, use bisection
      d = xm
      e = d
    endif
  else
    ! bound decreasing too slowly, use bisection
    d = xm
    e = d
  endif
  
  ! move last best guest to a
  a = b
  fa = fb
  
  ! evaluate new trial root
  b = b + merge(d, sign(tol1,xm), abs(d) .gt. tol1)
  
  fb = zero_finding_temperature(temperature=b, rho=p_prop%bulk_density, omega=p_prop%omega, prefactor=prefactor)
end do
write(*,*) 'Warning: zbrent exceeding maximum iterations'
zero_finding_zbrent=b
return
end function zero_finding_zbrent

function zero_finding_temperature(temperature, rho, omega, prefactor)
! function that is thought to be equal to zero when the good temperature is retrieved. For that purpose, various parameters are needed. 
! This f(x) = 0 function is obtained by using (37) in (36) (paardekooper, baruteau & kley 2010). 
! We also use the opacity given in Bell & lin 1994. 

! REMARKS : The scaleheight of the disk is determined directly in the function, because it depends on the temperature


! Output
real(double_precision) :: zero_finding_temperature ! the value of the function

! Input
real(double_precision), intent(in) :: temperature ! the temperature at a given position (in K)
real(double_precision), intent(in) :: rho ! the bulk density at a given position (in MS/AU**3)
real(double_precision), intent(in) :: omega ! the angular velocity of the disk at a given position
real(double_precision), intent(in) :: prefactor ! = - (9.d0 * nu * sigma * omega**2 / 16.d0)

! Local
real(double_precision) :: optical_depth ! the optical depth at a given position

!------------------------------------------------------------------------------

optical_depth = get_opacity(temperature, rho) * rho * scaleheight_prefactor * sqrt(temperature) / omega

! 1.7320508075688772d0 = sqrt(3)
zero_finding_temperature = SIGMA_STEFAN * temperature**4 + prefactor * &
                           (1.5d0 * optical_depth  + 1.7320508075688772d0 + 1 / (optical_depth))
return
end function zero_finding_temperature

subroutine get_temperature(ln_x, ln_y, radius, temperature)
! subroutine that interpolate a value of the temperature at a given radius with input arrays of radius (x) and temperature (y)

! Warning : 
! the 'x' array must contains equally spaced 'r' values in linear basis (but not thei logarithm values of course). 
! BUT 'x' and 'y' MUST be respectively log(r) and log(T).
! 'x' and 'y' must have the same size !

! If the given radius is out of the radius boundaries of the temperature profile, 
! then the temperature of the closest bound of the temperature profile will be given.

! Return : 
! temperature : the temperature (in K) at the radius 'radius'

real(double_precision), dimension(:), intent(in) :: ln_x, ln_y
real(double_precision), intent(in) :: radius

real(double_precision), intent(out) :: temperature

! Local
real(double_precision) :: radius_step ! the step between each radial values in the 'x' array
integer :: id_max ! the last index of arrays
integer :: closest_low_id ! the index of the first closest lower value of radius regarding the radius value given in parameter of the subroutine. 
real(double_precision) :: x_min, x_max
real(double_precision) :: ln_x1, ln_x2, ln_y1, ln_y2

id_max = ubound(ln_x,1)
x_min = exp(ln_x(1))
x_max = exp(ln_x(id_max))


if ((radius .gt. x_min) .and. (radius .gt. x_max)) then
  ! in the range
  radius_step = exp(ln_x(2)) - x_min ! Only valid if the separation between radial values is linear. (equal space between the values, but not their logarithm)
  closest_low_id = 1 + int((radius - x_min) / radius_step)
  
  ln_x1 = ln_x(closest_low_id)
  ln_x2 = ln_x(closest_low_id + 1)
  ln_y1 = ln_y(closest_low_id)
  ln_y2 = ln_y(closest_low_id + 1)
  
  temperature = exp(ln_y2 + (ln_y1 - ln_y2) * (log(radius) - ln_x2) / (ln_x1 - ln_x2))
else if (radius .lt. x_min) then
  temperature = exp(ln_y(1))
else if (radius .gt. x_max) then
  temperature = exp(ln_y(id_max))
end if

end subroutine get_temperature
end module user_module

! TODO utiliser la masse des objets pour ne pas faire le calcul si trop massif, il faut respecter le domaine de validité des formules des couples
! TODO routine générale de conversion des couples en accélération afin de pouvoir réutiliser ailleurs

!TODO : 
!_ amortissement de l'eccentricité
!_ amortissement de l'inclinaison
!_rajouter le spin comme variable d'entrée de la routine mfo_user
