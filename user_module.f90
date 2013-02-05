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
  
  public :: mfo_user, unitary_tests ! unitary_tests must be used only to test function, you should not use it out of the debug phase
  
  ! If we don't put a cutoff, then the simulation will crash if the inclination become exactly equal to zero, 
  ! which in addition is not really physically possible.
  ! The other idea behind this cutoff is to allow planets to come close, and pass one in front of the other without collision. 
  ! Hence, the idea is to have a cutoff sufficiently high to allow, at a given position, 
  ! 2 planet radius in the apperture of this cutoff (a triangle with angle, orbital distance and 2 planet radius)
  real(double_precision), parameter :: INCLINATION_CUTOFF = 5.d-4 ! (in rad) the value below whom there will be no inclination damping anymore.
  
  !------------------------------------------------------------------------------
  ! Default values for parameters that are to be read in the parameter file 'disk.in'
  real(double_precision) :: B_OVER_H = 0.4 ! the smoothing length for the planet's potential
  real(double_precision) :: ADIABATIC_INDEX = 1.4 ! the adiabatic index for the gas equation of state
  real(double_precision) :: MEAN_MOLECULAR_WEIGHT = 2.35 ! the mean molecular weight in mass of a proton
  
  ! Here we define the power law for surface density sigma(R) = INITIAL_SIGMA_0 * R^(-INITIAL_SIGMA_INDEX)
  real(double_precision) :: INITIAL_SIGMA_0 = 450 ! the surface density at (R=1AU) [g/cm^2]
  real(double_precision) :: INITIAL_SIGMA_INDEX = 0.5! the negative slope of the surface density power law (alpha in the paper)
  real(double_precision) :: INITIAL_SIGMA_0_NUM ! the surface density at (R=1AU) [Msun/AU^2]
  logical :: IS_DISSIPATION = .True. ! boolean to tell if there is dissipation of the disk or not.
  real(double_precision) :: dissipation_timestep ! the timestep between two computation of the disk [in days]
  character(len=80) :: INNER_BOUNDARY_CONDITION = 'open' ! 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid
  character(len=80) :: OUTER_BOUNDARY_CONDITION = 'open' ! 'open' or 'closed'. If open, gas can cross the outer edge. If closed, nothing can escape the grid
  
  ! Here we define the constant value of the viscosity of the disk
  real(double_precision) :: viscosity = 1.d15 ! viscosity of the disk [cm^2/s]
  
  !------------------------------------------------------------------------------
  ! prefactors
  real(double_precision) :: X_S_PREFACTOR ! prefactor for the half width of the corotation region
  real(double_precision) :: SCALEHEIGHT_PREFACTOR ! prefactor for the scaleheight
  
  !------------------------------------------------------------------------------
  ! Here we define properties common to the profiles
  real(double_precision) :: INNER_BOUNDARY_RADIUS = 1d0
  real(double_precision) :: OUTER_BOUNDARY_RADIUS = 100.d0
  integer :: NB_SAMPLE_PROFILES = 200 ! number of points for the sample of radius of the temperature profile
  real(double_precision) :: X_SAMPLE_STEP ! the constant step for the x_sample. Indeed, due to diffusion equation, the sample must be constant in X, and not in r. 

  real(double_precision), dimension(:), allocatable :: distance_log_sample ! values of 'a' in log()
  real(double_precision), dimension(:), allocatable :: x_sample ! values of 'x' with x = 2*sqrt(r) (used for the diffusion equation)
  real(double_precision), dimension(:), allocatable :: surface_density_profile ! values of the density in log() for each value of the 'a' sample
  real(double_precision), dimension(:), allocatable :: surface_density_index ! values of the local negative slope of the surface density profile
  real(double_precision), dimension(:), allocatable :: temperature_profile ! values of the temperature in log() for each value of the 'a' sample
  real(double_precision), dimension(:), allocatable :: temp_profile_index ! values of the local negative slope of the temperature profile
  real(double_precision), dimension(:), allocatable :: chi_profile ! thermal diffusivity
  real(double_precision), dimension(:), allocatable :: tau_profile ! optical depth 
  
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
    real(double_precision) :: sigma_index ! the negative slope of the surface density profile at the location of the planet.
    real(double_precision) :: scaleheight ! the scaleheight of the disk at the location of the planet [AU]
    real(double_precision) :: aspect_ratio ! the aspect_ratio of the gas disk at the location of the planet [no dim]
    real(double_precision) :: chi ! the thermal diffusion coefficient at the location of the planet [AU^2.day^-1]
    real(double_precision) :: nu ! the viscosity of the disk at the location of the planet [AU^2.day^-1]
    real(double_precision) :: temperature ! the temperature of the disk at the location of the planet [K] 
    real(double_precision) :: temperature_index ! the negative temperature index of the disk at the location of the planet [no dim] 
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

! Global parameters
! IS_DISSIPATION : boolean to tell if there is dissipation of the disk or not.
! dissipation_timestep : the timestep between two computation of the disk [in days]
! X_SAMPLE_STEP : the constant step for the x_sample. Indeed, due to diffusion equation, the sample must be constant in X, and not in r. 
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

  real(double_precision) :: migration_acc_prefactor ! prefactor for the migration acceleration
  real(double_precision) :: eccentricity_acc_prefactor ! prefactor for the eccentricity acceleration

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
  type(PlanetProperties) :: p_prop ! various properties of a planet
  real(double_precision) :: e_h ! the ratio between the eccentricity and the aspect ratio for a given planet [no dim]
  real(double_precision) :: i_h ! the ratio between the inclination and the aspect ratio for a given planet [no dim]
  
  real(double_precision), dimension(3) :: migration_acceleration
  real(double_precision), dimension(3) :: eccentricity_acceleration
  real(double_precision) :: inclination_acceleration_z
  real(double_precision), save :: next_step = -1.d0 ! next time at which we will compute the thermal properties of the disk?
  
  !------------------------------------------------------------------------------
  ! Setup
  
  do planet=1,n_bodies
    acceleration(1,planet) = 0.d0
    acceleration(2,planet) = 0.d0
    acceleration(3,planet) = 0.d0
  end do
  
  call init_globals(stellar_mass=mass(1))
    
  !------------------------------------------------------------------------------
  ! If it's time (depending on the timestep we want between each calculation of the disk properties)
  ! The first 'next_step' is set to '-1' to force the calculation for the first timestep. In fact, the first timestep will be done fornothing, but we need this in order to have a clean code.
  if (IS_DISSIPATION) then
    if (time.gt.next_step) then
      dissipation_timestep = 0.5d0 * X_SAMPLE_STEP**2 / (4 * get_viscosity(1.d0)) ! a correction factor of 0.5 has been applied. No physical reason to that, just intuition and safety
      ! TODO if the viscosity is not constant anymore, the formulae for the dissipation timestep must be changed
      next_step = time + dissipation_timestep
      
      ! we get the density profile.
      call calculate_density_profile()
      
      ! we get the temperature profile.
      call calculate_temperature_profile()
      
      ! we store in a .dat file the temperature profile
      call store_temperature_profile(filename='temperature_profile.dat')
      call store_density_profile(filename='density_profile.dat')
      call store_scaleheight_profile()
    end if
  end if
  !------------------------------------------------------------------------------

  
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

      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to migration
      call get_corotation_torque(mass(1), mass(planet), p_prop, & ! input
      corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output

      torque = torque_ref * (lindblad_torque + corotation_torque)      
      
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
      
      if (p_prop%inclination.gt.INCLINATION_CUTOFF) then
        time_inc = time_wave / 0.544d0 * (1.d0 - 0.30d0 * i_h**2 + 0.24 * i_h**3 + 0.14 * e_h**2 * i_h)
        
        inclination_acceleration_z = - velocity(3,planet) / time_inc
      else
        time_inc = 0.d0
        inclination_acceleration_z = 0.d0
      end if
      
      !------------------------------------------------------------------------------
      ! Calculation of the total acceleration on the planet
      
      acceleration(1,planet) = migration_acceleration(1) + eccentricity_acceleration(1)
      acceleration(2,planet) = migration_acceleration(2) + eccentricity_acceleration(2)
      acceleration(3,planet) = migration_acceleration(3) + eccentricity_acceleration(3) + inclination_acceleration_z
      
    end if
  end do
  
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

! Global Parameters
! B_OVER_H : the smoothing length for the planet's potential
! ADIABATIC_INDEX : the adiabatic index for the gas equation of state
! MEAN_MOLECULAR_WEIGHT : the mean molecular weight in mass of a proton
! INITIAL_SIGMA_0 : the surface density at (R=1AU) [g/cm^2]
! INITIAL_SIGMA_INDEX : the negative slope of the surface density power law (alpha in the paper)
! INITIAL_SIGMA_0_NUM : the surface density at (R=1AU) [Msun/AU^2]
! INNER_BOUNDARY_RADIUS : the inner radius of the various profiles (all based on the radius profile)
! OUTER_BOUNDARY_RADIUS : the outer radius of the various profiles (all based on the radius profile)
! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
! viscosity : the viscosity of the disk in [cm^2/s]
! IS_DISSIPATION : boolean to tell if there is dissipation of the disk or not.
! INNER_BOUNDARY_CONDITION : 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid
! OUTER_BOUNDARY_CONDITION : 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid

  implicit none
  
  character(len=80) :: line
  character(len=1), parameter :: comment_character = '!' ! character that will indicate that the rest of the line is a comment
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
          read(value, *) B_OVER_H
        
        case('ADIABATIC_INDEX')
          read(value, *) ADIABATIC_INDEX
          
        case('MEAN_MOLECULAR_WEIGHT')
          read(value, *) MEAN_MOLECULAR_WEIGHT
          
        case('surface_density')
          read(value, *) INITIAL_SIGMA_0, INITIAL_SIGMA_INDEX

        case('temperature')
          read(value, *) INNER_BOUNDARY_RADIUS, OUTER_BOUNDARY_RADIUS
        
        case('sample')
          read(value, *) NB_SAMPLE_PROFILES
          
        case('viscosity')
          read(value, *) viscosity
          
        case('is_dissipation')
          read(value, *) IS_DISSIPATION
          
        case('INNER_BOUNDARY_CONDITION')
          read(value, *) INNER_BOUNDARY_CONDITION
        
        case('OUTER_BOUNDARY_CONDITION')
          read(value, *) OUTER_BOUNDARY_CONDITION
          
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
  
  INITIAL_SIGMA_0_NUM = INITIAL_SIGMA_0 * AU**2 / MSUN ! the surface density at (R=1AU) [Msun/AU^2]
  
end subroutine read_disk_properties

subroutine get_planet_properties(stellar_mass, mass, position, velocity, p_prop)

! subroutine that return numerous properties of the planet and its environment given its mass, position and velocity
! Note that some parameters are global and accessed directly by the subroutine

! Global parameters
! SCALEHEIGHT_PREFACTOR : prefactor for the scaleheight

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
!~   p_prop%sigma = get_surface_density(radius=p_prop%radius) ! [Msun/AU^3]
  call get_surface_density(radius=p_prop%radius, sigma=p_prop%sigma, sigma_index=p_prop%sigma_index)
!~   call print_planet_properties(p_prop)
  call get_temperature(radius=p_prop%radius, & ! Input
                       temperature=p_prop%temperature, temperature_index=p_prop%temperature_index, chi=p_prop%chi) ! Output
  
  ! We calculate the angular momentum
  p_prop%angular_momentum = (mass / K2) * h_p  
  p_prop%velocity = sqrt(velocity2_p) ! [AU/day]
  p_prop%omega = sqrt(gm / (p_prop%radius**3)) ! [day-1]
  
  !------------------------------------------------------------------------------
  ! H = sqrt(k_B * T / (omega^2 * mu * m_H))
  p_prop%scaleheight = SCALEHEIGHT_PREFACTOR * sqrt(p_prop%temperature) / p_prop%omega
!~   p_prop%scaleheight = 0.05 * p_prop%radius

  !------------------------------------------------------------------------------
!~   p_prop%nu = alpha * p_prop%omega * p_prop%scaleheight**2 ! [AU^2.day-1]
  p_prop%nu = get_viscosity(p_prop%radius)

  p_prop%aspect_ratio = p_prop%scaleheight / p_prop%radius
!~     write(*,'(e12.4)') p_prop%nu * AU**2 / DAY 
  
  !------------------------------------------------------------------------------
!~   p_prop%chi = 1.d-5 * p_prop%radius**2 * p_prop%omega ! comment if you want to use the thermal diffusivity calculated from the temperature profile
  
end subroutine get_planet_properties

subroutine get_corotation_torque(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, Gamma_0)
! function that return the total torque exerted by the disk on the planet 
!
! Global parameters
! ADIABATIC_INDEX : the adiabatic index for the gas equation of state
! X_S_PREFACTOR : prefactor for the half width of the corotation region

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
  real(double_precision) :: lindblad_prefactor ! prefactor for the lindblad torque
  
  !Properties of the disk at the location of the planet
  real(double_precision) :: x_s ! semi-width of the horseshoe region [radius_p (in unity of position of the planet)]
  real(double_precision) :: zeta_eff ! effective entropy index depending on gamma_eff [no dim]
  real(double_precision) :: p_nu ! parameter for saturation due to viscosity at the location of the planet [no dim]
  real(double_precision) :: p_chi ! parameter for saturation due to thermal diffusion at the location of the planet [no dim]
  real(double_precision) :: gamma_eff ! effective adiabatic index depending on several parameters [no dim]
  
  !Torques (some depends of the planet)
  real(double_precision) :: torque_hs_ent ! entropy related part of the horseshoe drag
  real(double_precision) :: torque_c_lin_ent ! entropy related part of the linear corotation torque
  real(double_precision) :: torque_hs_baro ! barotropic part of the horseshoe drag
  real(double_precision) :: torque_c_lin_baro ! barotropic part of the linear corotation torque
  
  !------------------------------------------------------------------------------
  ! WE CALCULATE TOTAL TORQUE EXERTED BY THE DISK ON THE PLANET
  Gamma_0 = (mass / (stellar_mass * p_prop%aspect_ratio))**2 * p_prop%sigma * p_prop%radius**4 * p_prop%omega**2
  
  !------------------------------------------------------------------------------
  ! Q is needed by the lindblad torque. We set Q for m ~ 2 /3 h (45): 
  Q_p = TWOTHIRD * p_prop%chi / (p_prop%aspect_ratio * p_prop%scaleheight**2 * p_prop%omega) ! p_prop%aspect_ratio**3 * p_prop%radius**2 = aspect_ratio * scaleheight**2
  !------------------------------------------------------------------------------
  
  gamma_eff = 2.d0 * Q_p * ADIABATIC_INDEX / (ADIABATIC_INDEX * Q_p + 0.5d0 * &
  sqrt(2.d0 * sqrt((ADIABATIC_INDEX * ADIABATIC_INDEX * Q_p * Q_p + 1.d0)**2 - 16.d0 * Q_p * Q_p * (ADIABATIC_INDEX - 1.d0)) &
  + 2.d0 * ADIABATIC_INDEX * ADIABATIC_INDEX * Q_p * Q_p - 2.d0))
  
  !------------------------------------------------------------------------------
  zeta_eff = p_prop%temperature_index - (gamma_eff - 1.d0) * p_prop%sigma_index
  
  x_s = X_S_PREFACTOR / gamma_eff**0.25d0 * sqrt(mass / p_prop%aspect_ratio)
  
  !------------------------------------------------------------------------------
  ! k_p is defined to limit the number of operation and to have a value independant from chi_p or nu_p
  k_p = p_prop%radius * p_prop%radius * p_prop%omega * x_s * x_s * x_s / (2.d0 * PI)
  
  !------------------------------------------------------------------------------
  p_nu = TWOTHIRD * sqrt(k_p / p_prop%nu)
  
  p_chi = sqrt(k_p / p_prop%chi)
  
  lindblad_prefactor = -(2.5d0 + 1.7d0 * p_prop%temperature_index - 0.1d0 * p_prop%sigma_index) ! paardekooper, baruteau & kley 2010
  lindblad_torque = lindblad_prefactor / gamma_eff ! lindblad torque formulae from pardekooper, 2010
  
  torque_hs_ent = 7.9d0 * zeta_eff / gamma_eff
  torque_c_lin_ent = (2.2d0 - 1.4d0 / gamma_eff) * zeta_eff
  
  torque_hs_baro = 1.1d0 * (1.5d0 - p_prop%sigma_index)
  torque_c_lin_baro = 0.7d0 * (1.5d0 - p_prop%sigma_index)
  
  !------------------------------------------------------------------------------

  corotation_torque = (1.d0 / gamma_eff) * torque_hs_baro * get_F(p_nu) * get_G(p_nu) + torque_c_lin_baro * (1 - get_K(p_nu)) &
    + torque_hs_ent * get_F(p_nu) * get_F(p_chi) * sqrt(get_G(p_nu) * get_G(p_chi)) &
    + torque_c_lin_ent * sqrt((1 - get_K(p_nu)) * (1 - get_K(p_chi)))
  

  return
end subroutine get_corotation_torque

subroutine init_globals(stellar_mass)
! subroutine that initialize global values that define prefactors or values for torque that does not depend on the planet properties
!
! Global Parameters
! B_OVER_H : the smoothing length for the planet's potential
! MEAN_MOLECULAR_WEIGHT : the mean molecular weight in mass of a proton
! INNER_BOUNDARY_RADIUS : the inner radius of the various profiles (all based on the radius profile)
! OUTER_BOUNDARY_RADIUS : the outer radius of the various profiles (all based on the radius profile)
! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
! X_SAMPLE_STEP : the constant step for the x_sample. Indeed, due to diffusion equation, the sample must be constant in X, and not in r. 
! X_S_PREFACTOR : prefactor for the half width of the corotation region
! SCALEHEIGHT_PREFACTOR : prefactor for the scaleheight
! distance_log_sample : values of 'a' in log()
! x_sample : values of 'x' with x = 2*sqrt(r) (used for the diffusion equation)
! surface_density_profile : values of the density in log() for each value of the 'a' sample
! surface_density_index : values of the local negative slope of the surface density profile
! temperature_profile : values of the temperature in log() for each value of the 'a' sample
! temp_profile_index : values of the local negative slope of the temperature profile
! chi_profile : thermal diffusivity
! tau_profile : optical depth 
!
! Parameters
! stellar_mass : the mass of the central object in solar mass (times K2)


  
  implicit none
  real(double_precision), intent(in) :: stellar_mass
  logical, save :: FirstCall = .True.
  integer :: i  
  
  if (FirstCall) then
    FirstCall = .False.
    
    call read_disk_properties()
    
    allocate(distance_log_sample(NB_SAMPLE_PROFILES))
    distance_log_sample(1:NB_SAMPLE_PROFILES) = 0.d0
    
    allocate(x_sample(NB_SAMPLE_PROFILES))
    x_sample(1:NB_SAMPLE_PROFILES) = 0.d0
    
    allocate(surface_density_profile(NB_SAMPLE_PROFILES))
    allocate(surface_density_index(NB_SAMPLE_PROFILES))
    surface_density_profile(1:NB_SAMPLE_PROFILES) = 0.d0
    surface_density_index(1:NB_SAMPLE_PROFILES) = 0.d0
    
    allocate(temperature_profile(NB_SAMPLE_PROFILES))
    allocate(temp_profile_index(NB_SAMPLE_PROFILES))
    temperature_profile(1:NB_SAMPLE_PROFILES) = 0.d0
    temp_profile_index(1:NB_SAMPLE_PROFILES) = 0.d0
    
    allocate(chi_profile(NB_SAMPLE_PROFILES))
    allocate(tau_profile(NB_SAMPLE_PROFILES))
    chi_profile(1:NB_SAMPLE_PROFILES) = 0.d0
    tau_profile(1:NB_SAMPLE_PROFILES) = 0.d0
    
    ! We calculate the initial surface density profile.
    ! First, we want a constant spaced x_sample (which is propto sqrt(r)). Because it is important for diffusion equation which is solved depending on X and not R
    x_sample(1) = 2.d0 * sqrt(INNER_BOUNDARY_RADIUS)
    distance_log_sample(1) = log(INNER_BOUNDARY_RADIUS)
    
    x_sample(NB_SAMPLE_PROFILES) = 2.d0 * sqrt(OUTER_BOUNDARY_RADIUS)
    distance_log_sample(NB_SAMPLE_PROFILES) = log(OUTER_BOUNDARY_RADIUS)
    
    ! We initialize the global variable (in the module) for the constant step of x_sample
    X_SAMPLE_STEP = (x_sample(NB_SAMPLE_PROFILES) - x_sample(1)) / (NB_SAMPLE_PROFILES - 1.d0)
    
    do i=2, NB_SAMPLE_PROFILES - 1
      x_sample(i) = x_sample(1) + X_SAMPLE_STEP * (i - 1.d0)
      distance_log_sample(i) = log(0.25d0 * x_sample(i)**2)
    end do
    
    ! The x_s value is corrected from (paardekooper, 2010). The expression used is the one from (paardekooper, 2009a)
    X_S_PREFACTOR = 1.1d0 * (0.4d0 / B_OVER_H)**0.25d0 / sqrt(stellar_mass) ! mass(1) is here for the ratio of mass q
    
    ! AU is in cm, so we must turn into meter before doing the conversion
    ! division of k_B by m_H is done separately for exponant and value to have more precision
    ! sqrt(k_B/m_H) in numerical units, knowing that [k_B]=[m^2.kg.s^-2K^-1] and [m_H]=[kg]. 
    SCALEHEIGHT_PREFACTOR = sqrt(1.3806503d0/(1.67262158d0 * MEAN_MOLECULAR_WEIGHT) * 1.d4) * DAY / (AU * 1.d-2) 
    
    call initial_density_profile()
    
    ! we get the temperature profile, but we need the surface density profile before.
    call calculate_temperature_profile() ! WARNING : SCALEHEIGHT_PREFACTOR must exist before the temperature profile is computed !
    
    ! we store in a .dat file the temperature profile
    call store_temperature_profile(filename='temperature_profile.dat')
    call store_density_profile(filename='density_profile.dat')
    call store_scaleheight_profile()
    
    ! Here we display various warning for specific modification of the code that must be kept in mind (because this is not the normal behaviour of the code)

  endif
end subroutine init_globals

subroutine initial_density_profile()
  ! subroutine that store in the global parameters the value of the initial surface density profile. 
  !
  ! Global parameters
  ! INITIAL_SIGMA_0 : the surface density at (R=1AU) [g/cm^2]
  ! INITIAL_SIGMA_INDEX : the negative slope of the surface density power law (alpha in the paper)
  ! distance_log_sample : values of 'a' in log()
  ! surface_density_profile : values of the density in log() for each value of the 'a' sample
  ! surface_density_index : values of the local negative slope of the surface density profile
  
  implicit none
  
  integer :: i ! for loops
  
  do i=1,NB_SAMPLE_PROFILES
    surface_density_profile(i) = log(INITIAL_SIGMA_0_NUM * exp(distance_log_sample(i))**(-INITIAL_SIGMA_INDEX))
    surface_density_index(i) = INITIAL_SIGMA_INDEX
  end do
end subroutine initial_density_profile

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
  
  subroutine get_surface_density(radius, sigma, sigma_index)
    ! function that interpolate the value of the surface density at the given radius
    
    ! Global parameter
    ! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
    ! X_SAMPLE_STEP : the constant step for the x_sample. Indeed, due to diffusion equation, the sample must be constant in X, and not in r. 
    ! INNER_BOUNDARY_RADIUS : the inner edge of the different profiles
    ! OUTER_BOUNDARY_RADIUS : the outer edge of the different profiles
    ! distance_log_sample : values of 'a' in log()
    ! x_sample : values of 'x' with x = 2*sqrt(r) (used for the diffusion equation)
    ! surface_density_profile : values of the density in log() for each value of the 'a' sample
    ! surface_density_index : values of the local negative slope of the surface density profile
    
    ! Parameters : 
    ! radius : the orbital distance [in AU]

    ! Warning : ! the surface density profile is a global parameter of the module. So nothing is given in parameter because it's 
    ! this global array that change whenever needed. 

    ! Return : 
    ! temperature : the temperature (in K) at the radius 'radius'
    ! If the given radius is out of the radius boundaries of the temperature profile, 
    ! then the temperature of the closest bound of the temperature profile will be given.

    real(double_precision), intent(in) :: radius

    real(double_precision), intent(out) :: sigma ! the surface density at 'radius' in [MSUN/AU^2]
    real(double_precision), intent(out) :: sigma_index ! the negative slope of the surface density profile at the location of the planet.

    ! Local
    integer :: closest_low_id ! the index of the first closest lower value of radius regarding the radius value given in parameter of the subroutine. 
    real(double_precision) :: ln_x1, ln_x2, ln_y1, ln_y2
    real(double_precision) :: x_radius ! the corresponding 'x' value for the radius given in parameter of the routine. we can retrieve the index of the closest values in this array in only one calculation.

    if ((radius .ge. INNER_BOUNDARY_RADIUS) .and. (radius .lt. OUTER_BOUNDARY_RADIUS)) then
      
      x_radius = 2.d0 * sqrt(radius)
      ! in the range
      closest_low_id = 1 + int((x_radius - x_sample(1)) / X_SAMPLE_STEP) ! X_SAMPLE_STEP being a global value, x_sample also
      
      ln_x1 = distance_log_sample(closest_low_id)
      ln_x2 = distance_log_sample(closest_low_id + 1)
      ln_y1 = surface_density_profile(closest_low_id)
      ln_y2 = surface_density_profile(closest_low_id + 1)

      sigma = exp(ln_y2 + (ln_y1 - ln_y2) * (log(radius) - ln_x2) / (ln_x1 - ln_x2))
      sigma_index = surface_density_index(closest_low_id) ! for the temperature index, no interpolation.
    else if (radius .lt. INNER_BOUNDARY_RADIUS) then
      sigma = exp(surface_density_profile(1))
      sigma_index = surface_density_index(1)
    else if (radius .gt. OUTER_BOUNDARY_RADIUS) then
      sigma = exp(surface_density_profile(NB_SAMPLE_PROFILES))
      sigma_index = surface_density_index(NB_SAMPLE_PROFILES)
    end if
  end subroutine get_surface_density
  
  function get_viscosity(radius)
  ! function that return the viscosity of the disk in [AU^2.day^-1]
  
  ! Global parameters
  ! density : the viscosity of the disk in [cm^2/s]
  
  implicit none
  
  real(double_precision) :: get_viscosity ! the viscosity of the disk in [AU^2.day^-1]
  
  real(double_precision), intent(in) :: radius ! the distance from the central object in AU
  
  get_viscosity = viscosity * DAY / AU**2
  ! TODO if the viscosity is not constant anymore, the formulae for the dissipation timestep must be changed
  
  end function get_viscosity
  
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

  subroutine calculate_temperature_profile()
! subroutine that calculate the temperature profile of the disk given various parameters including the surface density profile.
! 
! Global parameters
! ADIABATIC_INDEX : the adiabatic index for the gas equation of state
! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
! distance_log_sample : values of 'a' in log()
! temperature_profile : values of the temperature in log() for each value of the 'a' sample
! temp_profile_index : values of the local negative slope of the temperature profile
! chi_profile : thermal diffusivity
! tau_profile : optical depth 
!
! Return
! Nothing, but store in the associated global variable the temperature profile

    implicit none

    ! Local
    real(double_precision), parameter :: mass = 20. * EARTH_MASS * K2
    
    real(double_precision) :: a
    real(double_precision) :: stellar_mass, position(3), velocity(3), temperature, exponant
    type(PlanetProperties) :: p_prop

    ! value for the precedent step of the loop. In order to calculate the index of the local temperature power law.
    real(double_precision) :: a_old, temperature_old, tmp
    
    integer :: i,j ! for loops
    
    
    position(1:3) = 0.d0
    velocity(1:3) = 0.d0
    
    ! stellar mass
    stellar_mass = 1.d0 * K2  
  
    
    a = exp(distance_log_sample(1))
    ! We generate cartesian coordinate for the given semi major axis
    position(1) = a
    
    ! We generate cartesian coordinate for the given mass and semi major axis
    velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
    
    ! we store in global parameters various properties of the planet
    call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
     p_prop=p_prop) ! Output
     
    call zbrent(x_min=1.d-5, x_max=1.d4, tolerance=1d-4, p_prop=p_prop, & ! Input
                            temperature=temperature, optical_depth=tau_profile(1)) ! Output    

    temperature_profile(1) = log(temperature)
    chi_profile(1) = 1.5d0 * p_prop%nu * ADIABATIC_INDEX * (ADIABATIC_INDEX - 1.d0) * &
                      (1.5d0 + sqrt(3.d0) / tau_profile(1) + 1 / tau_profile(1)**2)
    do j=2,NB_SAMPLE_PROFILES
      a_old = a
      temperature_old = temperature
      
      a = exp(distance_log_sample(j)) ! Be carefull, the step between 'a' values is not constant !
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a
      
      ! We generate cartesian coordinate for the given mass and semi major axis
      velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
      
      call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
       p_prop=p_prop) ! Output
      
      call zbrent(x_min=1.d-5, x_max=1.d4, tolerance=1d-4, p_prop=p_prop, & ! Input
                              temperature=temperature, optical_depth=tau_profile(j)) ! Output
      
      temperature_profile(j) = log(temperature)
      temp_profile_index(j) = - (temperature_profile(j) - temperature_profile(j-1)) / &
                                (distance_log_sample(j) - distance_log_sample(j-1))
      chi_profile(j) = 1.5d0 * p_prop%nu * ADIABATIC_INDEX * (ADIABATIC_INDEX - 1.d0) * &
                      (1.5d0 + sqrt(3.d0) / tau_profile(j) + 1 / tau_profile(j)**2)
      
    end do
    
    temp_profile_index(2) = temp_profile_index(3)! We do this because the first two values are polluted by boundary conditions
    temp_profile_index(1) = temp_profile_index(2) ! the first element of the index array is forced to be equal to the second one because else, we can't calculate it.
    
  end subroutine calculate_temperature_profile
  
  subroutine calculate_density_profile()
! subroutine that calculate the temperature profile of the disk given various parameters including the surface density profile.
! 
! Global Parameters 
! INNER_BOUNDARY_CONDITION : 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid
! OUTER_BOUNDARY_CONDITION : 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid
! dissipation_timestep : the timestep between two computation of the disk [in days]
! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
! x_sample : values of 'x' with x = 2*sqrt(r) (used for the diffusion equation)
! surface_density_profile : values of the density in log() for each value of the 'a' sample
! surface_density_index : values of the local negative slope of the surface density profile

    use mercury_constant

    implicit none

    ! value for the precedent step of the loop. In order to calculate the index of the local temperature power law.
    real(double_precision) :: a_im1, a_i, a_ip1, a_ip2 ! ip1 == (i+1) ; im1 == (i-1)
    real(double_precision) :: f_i, f_ip1, f_ip2
    real(double_precision) :: flux_im1, flux_i, flux_ip1, flux_ip2 ! ip1 == (i+1) ; im1 == (i-1)
    integer :: i ! for loops
    logical, save :: FirstCall = .True.
    real(double_precision) :: tmp
    
    ! The dissipation_timestep is a global variable that is the timestep between two calculation of the surface density profile 
    ! This value is calculated in mfo_user, in the 'if' loop where the re-calculation is done. This value is needed here in the 
    ! part that dissipate the disk, but is not needed for the first run of this routine.
    
    ! On the first run, we compute the initial profile, but on the other ones, we run the diffusion equation to make it dissipate.
    
    if (.not.FirstCall) then
      ! Warning : X_SAMPLE_STEP is not constant because equally spaced a values does not mean equally spaced X values (due to the square root)
            
      ! for i=1
      ! surface_density_index(1) is to be defined later, using surface_density_index(2)
      
      ! We prepare values for the first step of the loop 
      ! (i.e i=1 because we simulate the ghost step of the loop to be able to shift the temp values)
      a_i = 0.25d0   * x_sample(1) * x_sample(1)
      a_ip1 = 0.25d0 * x_sample(2) * x_sample(2)
      
      f_i   = 1.5d0 * x_sample(1) * exp(surface_density_profile(1)) 
      f_ip1 = 1.5d0 * x_sample(2) * exp(surface_density_profile(2))
      f_ip2 = 1.5d0 * x_sample(3) * exp(surface_density_profile(3))
      
      flux_i = 3.d0 * get_viscosity(a_i) / a_i * (f_ip1 - f_i) / X_SAMPLE_STEP ! not a centered scheme, simple step, so 3 instead of 1.5 (because we divide by (X_SAMPLE_STEP) instead of (2*X_SAMPLE_STEP)
      flux_ip1 = 1.5d0 * get_viscosity(a_ip1) / a_ip1 * (f_ip2 - f_i) / X_SAMPLE_STEP
      
      ! CONDITION AT THE INNER EDGE
      select case(INNER_BOUNDARY_CONDITION)
      case('closed') 
      ! correspond to the case where the velocity at the inner edge is forced to be zero. This is equivalent to a 0-flux condition
        flux_i = 0.d0
        surface_density_profile(1) = log((f_i + dissipation_timestep * flux_ip1 / X_SAMPLE_STEP) / (1.5d0 * x_sample(i)))
      case('open') 
      ! correspond to the case where the surface density is forced to be zero. This is equivalent to an accretion 
      ! condition. Density is free to dissipate outside the grid.
        surface_density_profile(1) = log(TINY) ! in log, '0' is not defined, so we put a huge negative value to get close to 0
        f_i = 0.d0 ! We impose the condition sigma=0 so in practice it should be equal to 0
      case default
        write(*,*) 'Warning: An unknown inner boundary condition has been found'
        write(*,*) "inner_boundary_condition=", trim(INNER_BOUNDARY_CONDITION)
      end select
      
      ! We store the previous values in order to avoid the use of an array. 
      ! This way, it should more efficient, because we don't have to create several arrays with thousands of elements
      do i=2,NB_SAMPLE_PROFILES-3
        
        ! We shift the indexes by 1
        a_im1 = a_i
        a_i = a_ip1
        a_ip1 = 0.25d0 * x_sample(i+1) * x_sample(i+1)
        
        f_i = f_ip1
        f_ip1 = f_ip2
        f_ip2 = 1.5d0 * x_sample(i+2) * exp(surface_density_profile(i+2))
        
        flux_im1 = flux_i
        flux_i = flux_ip1
        flux_ip1 = 1.5d0 * get_viscosity(a_ip1) / a_ip1 * (f_ip2 - f_i) / X_SAMPLE_STEP
        
        tmp = (f_i + dissipation_timestep * (flux_ip1 - flux_im1) / (2 * X_SAMPLE_STEP)) / (1.5d0 * x_sample(i))
        if (tmp.lt.0.) then
          write(*,*) 'ERROR: tmp is negative!!!!'
          write(*,*) a_i, a_ip1, f_i, f_ip1, f_ip2, flux_i, flux_ip1
        end if
        surface_density_profile(i) = log(tmp) ! the (1.5d0 * x_i) is here to convert from 'f' to Sigma
        surface_density_index(i) = - (surface_density_profile(i) - surface_density_profile(i-1)) &
                                      / (distance_log_sample(i) - distance_log_sample(i-1))
              
      end do
      !------------------------------------------------------------------------------
      i=NB_SAMPLE_PROFILES-2
      ! We shift the indexes by 1
      a_im1 = a_i
      a_i = a_ip1
      a_ip1 = 0.25d0 * x_sample(i+1) * x_sample(i+1)
      
      f_i = f_ip1
      f_ip1 = f_ip2
      f_ip2 = 1.5d0 * x_sample(i+2) * exp(surface_density_profile(i+2))
      
      flux_im1 = flux_i
      flux_i = flux_ip1
      flux_ip1 = 1.5d0 * get_viscosity(a_ip1) / a_ip1 * (f_ip2 - f_i) / X_SAMPLE_STEP
      
      
      tmp = (f_i + dissipation_timestep * (flux_ip1 - flux_im1) / (2 * X_SAMPLE_STEP)) / (1.5d0 * x_sample(i))
      surface_density_profile(i) = log(tmp) ! the (1.5d0 * x_i) is here to convert from 'f' to Sigma
      surface_density_index(i) = - (surface_density_profile(i) - surface_density_profile(i-1)) &
                                    / (distance_log_sample(i) - distance_log_sample(i-1))
      
      ! The boundary condition is computed here because some values are needed by "i=NB_SAMPLE_PROFILES-1 surface density value".
      select case(OUTER_BOUNDARY_CONDITION)
      case('closed') 
      ! here i=NB_SAMPLE_PROFILES-2
      ! correspond to the case where the velocity at the inner edge is forced to be zero. This is equivalent to a 0-flux condition
        flux_ip2 = 0.d0
        
        surface_density_profile(NB_SAMPLE_PROFILES) = log((f_ip2 - dissipation_timestep * flux_ip1 / X_SAMPLE_STEP) &
                                                          / (1.5d0 * x_sample(i+2)))
      case('open') 
      ! correspond to the case where the surface density is forced to be zero. This is equivalent to an accretion 
      ! condition. Density is free to dissipate outside the grid.
        surface_density_profile(NB_SAMPLE_PROFILES) = -HUGE ! in log, '0' is not defined, so we put a huge negative value to get close to 0
        
        a_ip2 = 0.25d0 * x_sample(i+2) * x_sample(i+2)
        flux_ip2 = - 3.d0 * get_viscosity(a_ip2) / a_ip2 * f_ip1 / X_SAMPLE_STEP ! simple step, so 3 instead of 1.5 (because we divide by (X_SAMPLE_STEP) instead of (2*X_SAMPLE_STEP)
      case default
        write(*,*) 'Warning: An unknown outer boundary condition has been found'
        write(*,*) "outer_boundary_condition=", trim(OUTER_BOUNDARY_CONDITION)
      end select    
      
      !------------------------------------------------------------------------------
      i=NB_SAMPLE_PROFILES-1
      ! We shift the indexes by 1
      f_i = f_ip1
      f_ip1 = f_ip2
      
      flux_im1 = flux_i
      flux_i = flux_ip1
      flux_ip1 = flux_ip2 ! The outer boundary condition determiner flux_ip2
      
      surface_density_profile(i) = log((f_i + dissipation_timestep * (flux_ip1 - flux_im1) / (2 * X_SAMPLE_STEP)) &
                                            / (1.5d0 * x_sample(i))) ! the (1.5d0 * x_i) is here to convert from 'f' to Sigma
      surface_density_index(i) = - (surface_density_profile(i) - surface_density_profile(i-1)) &
                                        / (distance_log_sample(i) - distance_log_sample(i-1))
      
      !------------------------------------------------------------------------------
      ! And the last index NB_SAMPLE_PROFILES 
      ! We DO NOT compute the surface density since it's ruled by the boundary condition. But we compute the index
      
      surface_density_index(i) = - (surface_density_profile(i) - surface_density_profile(i-1)) &
                                      / (distance_log_sample(i) - distance_log_sample(i-1))
  
      !------------------------------------------------------------------------------
      ! We copy paste the index for the first value because however it is not defined.
      surface_density_index(2) = surface_density_index(3) ! We do this because the first two values are polluted by boundary conditions
      surface_density_index(1) = surface_density_index(2)
    else ! If it's the first call, we compute the initial profile
      FirstCall = .False.
      
      call initial_density_profile()
    end if
  
  end subroutine calculate_density_profile
  
  subroutine store_temperature_profile(filename)
  ! subroutine that store in a '.dat' file the temperature profile and negative index of the local power law
  
  ! Global parameters
  ! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
  ! temperature_profile : values of the temperature in log() for each value of the 'a' sample
  ! temp_profile_index : values of the local negative slope of the temperature profile
  ! chi_profile : thermal diffusivity
  ! tau_profile : optical depth 
  
  implicit none
  
  character(len=*), intent(in) :: filename
  
  integer :: j ! for loops
  
  ! We open the file where we want to write the outputs
  open(10, file=filename, status='replace')
  write(10,*) '# a in AU            ;    temperature (in K)    ;       exponant   &
              &; chi (thermal diffusivity) ;    tau (optical depth)'

  do j=1,NB_SAMPLE_PROFILES
    write(10,*) exp(distance_log_sample(j)), exp(temperature_profile(j)), temp_profile_index(j), tau_profile(j), chi_profile(j)!, distance_log_sample(j), temperature_profile(j)
  end do
  
  close(10)
  
  end subroutine store_temperature_profile
  
  subroutine store_density_profile(filename)
  ! subroutine that store in a '.dat' file the temperature profile and negative index of the local power law
  
  ! Global parameters
  ! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
  ! surface_density_index : values of the local negative slope of the surface density profile
  
  
  implicit none
  
  character(len=*), intent(in) :: filename
  
  real(double_precision), parameter :: NUM2PHYS = MSUN / AU**2
  
  integer :: j ! for loops
  
  ! We open the file where we want to write the outputs
  open(10, file=filename, status='replace')
  write(10,*) '#       a in AU       ; surface density (in g/cm^2) ;    exponant'

  do j=1,NB_SAMPLE_PROFILES
    write(10,*) exp(distance_log_sample(j)), exp(surface_density_profile(j)) * NUM2PHYS, surface_density_index(j)
  end do
  
  close(10)
  
  end subroutine store_density_profile
  
  subroutine store_scaleheight_profile()
  ! subroutine that store in a '.dat' file the scaleheight profile
  
  ! Global parameters
  ! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
  
  implicit none
  
  integer :: j ! for loops
  real(double_precision) :: a, scaleheight
  real(double_precision) :: position(3), velocity(3), stellar_mass, mass
  type(PlanetProperties) :: p_prop

  position(1:3) = 0.d0
  velocity(1:3) = 0.d0

  ! stellar mass
  stellar_mass = 1.d0 * K2
  ! planet mass
  mass = 20. * EARTH_MASS * K2
  
  ! We open the file where we want to write the outputs
  open(10, file='scaleheight_profile.dat', status='replace')
  write(10,*) '# a in AU ; scaleheight (AU) ; aspect ratio'
  do j=1,NB_SAMPLE_PROFILES
    a = exp(distance_log_sample(j))
    ! We generate cartesian coordinate for the given semi major axis
    position(1) = a
    
    ! We generate cartesian coordinate for the given mass and semi major axis
    velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
    
    ! we store in global parameters various properties of the planet
    call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
     p_prop=p_prop) ! Output
    write(10,*) a, p_prop%scaleheight, p_prop%scaleheight/a
  end do
  
  close(10)

  
  end subroutine store_scaleheight_profile

subroutine zbrent(x_min, x_max, tolerance, p_prop, temperature, optical_depth)
! Using Brent's method, find the root of a function 'func' known to lie between 'x_min' and 'x_max'. 
! The root, returned as 'zero_finding_zbrent', will be refined until its accuray is 'tolerance'. 

! Parameters :
! ITMAX : maximum allowed number of iterations
! EPS : machine floating-point precision.

! REMARK : This function is based on the zbrent function in fortran 90 of numerical recipes

! Output
real(double_precision), intent(out) :: temperature
real(double_precision), intent(out) :: optical_depth

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
real(double_precision) :: tau_a, tau_b

if (isnan(p_prop%sigma)) then
  write(*,*) 'Error: the surface density is equal to NaN when we want to calculate the temperature profile'
end if 

if (isnan(p_prop%nu)) then
  write(*,*) 'Error: the viscosity is equal to NaN when we want to calculate the temperature profile'
end if 

if (isnan(p_prop%omega)) then
  write(*,*) 'Error: the angular velocity is equal to NaN when we want to calculate the temperature profile'
end if 

if (isnan(p_prop%radius)) then
  write(*,*) 'Error: the distance is equal to NaN when we want to calculate the temperature profile'
end if

!------------------------------------------------------------------------------

! We calculate this value outside the function because we only have to do this once per step (per radial position)
prefactor = - (9.d0 * p_prop%nu * p_prop%sigma * p_prop%omega**2 / 16.d0)

a = x_min
b = x_max
call zero_finding_temperature(temperature=a, sigma=p_prop%sigma, omega=p_prop%omega, prefactor=prefactor,& ! Input
                              funcv=fa, optical_depth=tau_a) ! Output
call zero_finding_temperature(temperature=b, sigma=p_prop%sigma, omega=p_prop%omega, prefactor=prefactor,& ! Input
                              funcv=fb, optical_depth=tau_b) ! Output

!~ fa = zero_finding_temperature(temperature=a, sigma=p_prop%sigma, omega=p_prop%omega, prefactor=prefactor)
!~ fb = zero_finding_temperature(temperature=b, sigma=p_prop%sigma, omega=p_prop%omega, prefactor=prefactor)

if (((fa.gt.0.).and.(fb.gt.0.)).or.((fa.lt.0.).and.(fb.lt.0.))) then
  write(*,*) 'root must be bracketed for zbrent'
  write(*,*) '  T_min =', x_min, 'f(T_min) =', fa
  write(*,*) '  T_max =', x_max, 'f(T_max) =', fb
  write(*,*) 'properties of the disk at the location of the planet hat influence the value of the temperature'
  write(*,*) '  radial position of the planet (in AU) :', p_prop%radius
  write(*,*) '  viscosity :', p_prop%nu
  write(*,*) '  surface density :', p_prop%sigma
  write(*,*) '  angular velocity :', P_prop%omega
  stop
endif

! these values force the code to go into the first 'if' statement. 
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
    temperature = b
    optical_depth = tau_b
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
  
  call zero_finding_temperature(temperature=b, sigma=p_prop%sigma, omega=p_prop%omega, prefactor=prefactor,& ! Input
                                funcv=fb, optical_depth=tau_b) ! Output
end do
write(*,*) 'Warning: zbrent exceeding maximum iterations'
temperature = b
optical_depth = tau_b
return
end subroutine zbrent

subroutine zero_finding_temperature(temperature, sigma, omega, prefactor, funcv, optical_depth)
! function that is thought to be equal to zero when the good temperature is retrieved. For that purpose, various parameters are needed. 
! This f(x) = 0 function is obtained by using (37) in (36) (paardekooper, baruteau & kley 2010). 
! We also use the opacity given in Bell & lin 1994. 

! REMARKS : The scaleheight of the disk is determined directly in the function, because it depends on the temperature

! Global parameters
! SCALEHEIGHT_PREFACTOR : prefactor for the scaleheight


! Output
real(double_precision), intent(out) :: funcv ! the value of the function
real(double_precision), intent(out) :: optical_depth ! the optical depth at a given position

! Input
real(double_precision), intent(in) :: temperature ! the temperature at a given position (in K)
real(double_precision), intent(in) :: sigma ! the surface density at a given position (in MS/AU**2)
real(double_precision), intent(in) :: omega ! the angular velocity of the disk at a given position
real(double_precision), intent(in) :: prefactor ! = - (9.d0 * nu * sigma * omega**2 / 16.d0)

! Local
real(double_precision) :: scaleheight ! the scaleheight of the disk at a given position
real(double_precision) :: rho ! the bulk density of the disk at a given position
!------------------------------------------------------------------------------
scaleheight = SCALEHEIGHT_PREFACTOR * sqrt(temperature) / omega
rho = 0.5d0 * sigma / scaleheight
optical_depth = get_opacity(temperature, rho) * rho * scaleheight ! even if there is scaleheight in rho, the real formulae is this one. The formulae for rho is an approximation.

! 1.7320508075688772d0 = sqrt(3)
funcv = 2.d0 * SIGMA_STEFAN * temperature**4 + prefactor * &
                           (1.5d0 * optical_depth  + 1.7320508075688772d0 + 1 / (optical_depth))

return
end subroutine zero_finding_temperature

subroutine get_temperature(radius, temperature, temperature_index, chi)
! subroutine that interpolate a value of the temperature at a given radius with input arrays of radius (x) and temperature (y)

! Global parameters
! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
! X_SAMPLE_STEP : the constant step for the x_sample. Indeed, due to diffusion equation, the sample must be constant in X, and not in r.
! INNER_BOUNDARY_RADIUS : the inner edge of the different profiles
! OUTER_BOUNDARY_RADIUS : the outer edge of the different profiles
! distance_log_sample : values of 'a' in log()
! temperature_profile : values of the temperature in log() for each value of the 'a' sample
! temp_profile_index : values of the local negative slope of the temperature profile
! chi_profile : thermal diffusivity

! Warning : 
! the 'x' array must contains equally spaced 'r' values in linear basis (but not thei logarithm values of course). 
! BUT 'x' and 'y' MUST be respectively log(r) and log(T).
! 'x' and 'y' must have the same size !

! If the given radius is out of the radius boundaries of the temperature profile, 
! then the temperature of the closest bound of the temperature profile will be given.

! Return : 
! temperature : the temperature (in K) at the radius 'radius'

real(double_precision), intent(in) :: radius

real(double_precision), intent(out) :: temperature
real(double_precision), intent(out) :: temperature_index
real(double_precision), intent(out) :: chi

! Local
integer :: closest_low_id ! the index of the first closest lower value of radius regarding the radius value given in parameter of the subroutine. 
real(double_precision) :: ln_x1, ln_x2, ln_y1, ln_y2
real(double_precision) :: x_radius ! the corresponding 'x' value for the radius given in parameter of the routine. we can retrieve the index of the closest values in this array in only one calculation.


if ((radius .ge. INNER_BOUNDARY_RADIUS) .and. (radius .lt. OUTER_BOUNDARY_RADIUS)) then
  
  x_radius = 2.d0 * sqrt(radius)
  ! in the range
  closest_low_id = 1 + int((x_radius - x_sample(1)) / X_SAMPLE_STEP) ! X_SAMPLE_STEP being a global value, x_sample also
  
  ln_x1 = distance_log_sample(closest_low_id)
  ln_x2 = distance_log_sample(closest_low_id + 1)
  ln_y1 = temperature_profile(closest_low_id)
  ln_y2 = temperature_profile(closest_low_id + 1)

  temperature = exp(ln_y2 + (ln_y1 - ln_y2) * (log(radius) - ln_x2) / (ln_x1 - ln_x2))
  temperature_index = temp_profile_index(closest_low_id) ! for the temperature index, no interpolation.
  chi = chi_profile(closest_low_id)
else if (radius .lt. INNER_BOUNDARY_RADIUS) then
  temperature = exp(temperature_profile(1))
  temperature_index = temp_profile_index(1)
  chi = chi_profile(1)
else if (radius .gt. OUTER_BOUNDARY_RADIUS) then
  temperature = exp(temperature_profile(NB_SAMPLE_PROFILES))
  temperature_index = temp_profile_index(NB_SAMPLE_PROFILES)
  chi = chi_profile(NB_SAMPLE_PROFILES)
end if

end subroutine get_temperature

subroutine print_planet_properties(p_prop)
! subroutine that display in the terminal all the values 
! contained in the instance of planet properties given in parameters
!
! Parameters
! p_prop : an object of type 'PlanetProperties'
  implicit none
  type(PlanetProperties), intent(in) :: p_prop
  
  write (*,*) 'angular_momentum :', p_prop%angular_momentum 
  write (*,*) 'radius :', p_prop%radius 
  write (*,*) 'velocity :', p_prop%velocity 
  write (*,*) 'omega :', p_prop%omega 
  write (*,*) 'semi_major_axis :', p_prop%semi_major_axis 
  write (*,*) 'eccentricity :', p_prop%eccentricity 
  write (*,*) 'inclination :', p_prop%inclination 
  write (*,*) 'sigma :', p_prop%sigma 
  write (*,*) 'sigma_index :', p_prop%sigma_index
  write (*,*) 'scaleheight :', p_prop%scaleheight 
  write (*,*) 'aspect_ratio :', p_prop%aspect_ratio 
  write (*,*) 'chi :', p_prop%chi 
  write (*,*) 'nu :', p_prop%nu
  write (*,*) 'temperature :', p_prop%temperature 
  write (*,*) 'temperature_index:', p_prop%temperature_index
end subroutine print_planet_properties

! ##################################
! TESTS OF THE ROUTINES
! ##################################
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
    
    real(double_precision) :: stellar_mass

    stellar_mass = 1.d0 * K2
    
    write(*,*) 'Initialisation'
    call init_globals(stellar_mass)
    ! Note that the initial density profile and temperature profile are calculated inside the 'init_globals' routine.
    
    ! we store in a .dat file the temperature profile
    call store_temperature_profile(filename='temperature_profile.dat')
    call store_density_profile(filename='density_profile.dat')
    call store_scaleheight_profile()
    
    call test_functions_FGK()
    call test_opacity_profile()
    call test_torques(stellar_mass)
    call test_torques_fixed_a(stellar_mass)
    call test_torques_fixed_m(stellar_mass)
    call test_function_zero_temperature(stellar_mass)
    call test_temperature_profile(stellar_mass)
    call test_temperature_interpolation()
    call test_density_interpolation()
    call test_optical_depth_profile(stellar_mass)
    call test_thermal_diffusivity_profile(stellar_mass)
    call test_scaleheight_profile()
!~     
!~     call test_dissipation_of_the_disk(stellar_mass)
    
  end subroutine unitary_tests

! %%% unitary tests of some routines %%%
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
    
    write(*,*) 'test of functions F, G, K from (paardekooper, 2010)'
    
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

  subroutine test_function_zero_temperature(stellar_mass)
  ! subroutine that test the function 'zero_finding_temperature'
  
  ! Return:
  !  a data file and an associated gnuplot file.
    implicit none
    real(double_precision), intent(in) :: stellar_mass
    
    real(double_precision) :: temperature
    
    real(double_precision), parameter :: T_min = 0.
    real(double_precision), parameter :: T_max = 1000.
    integer, parameter :: nb_points = 2000
!~     real(double_precision), parameter :: T_step = (T_max/T_min) ** (1/(nb_points-1.d0))
    real(double_precision), parameter :: T_step = (T_max - T_min) / (nb_points - 1.d0)
    
    real(double_precision), parameter :: mass = 20. * EARTH_MASS * K2
    
    real(double_precision) :: zero_function, tmp ! value that we want to output and a dummy argument 'tmp'
    
    integer :: i,j ! for loops
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    real(double_precision) :: prefactor ! prefactor for the calculation of the function of the temperature whose zeros are searched

  !------------------------------------------------------------------------------
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    write(*,*) 'Test of the zero function used to calculate the temperature at a given radius'
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/test_function_zero_temperature.dat')
    
    ! We generate cartesian coordinate for the given semi major axis
    position(1) = 1.d0
    
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
    write(10,*) '# angular velocity :', P_prop%omega
    write(10,*) '# Temperature (K) ; value of the function. The right temperature is when the function is 0'
    
    do i=1, nb_points
!~       temperature = T_min * T_step ** (i-1)
      temperature = (T_min + T_step * (i - 1.d0))
      
      call zero_finding_temperature(temperature=temperature, sigma=p_prop%sigma, omega=p_prop%omega, prefactor=prefactor,& ! Input
                              funcv=zero_function, optical_depth=tmp) ! Output
!~       zero_function = zero_finding_temperature(temperature=temperature, sigma=p_prop%sigma, &
!~                                                omega=p_prop%omega, prefactor=prefactor)
      
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
    write(10,*) "replot # to generate the output file"  
    
    close(10)
    
  end subroutine test_function_zero_temperature

  subroutine test_temperature_interpolation()
  
    implicit none
    
    integer, parameter :: nb_a = 1000
    real(double_precision), parameter :: a_min = 0.d0 ! in AU
    real(double_precision), parameter :: a_max = 100.d0! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision) :: a, temperature, temperature_index, chi
    
    integer :: j ! for loops
    
    write(*,*) 'Test of the temperature interpolation'
    
    open(10, file='unitary_tests/test_temperature_interpolation.dat')
    do j=1,nb_a
      a = (a_min + a_step * (j - 1.d0))
      ! We generate cartesian coordinate for the given semi major axis
      
      call get_temperature(radius=a, & ! Input
                           temperature=temperature, temperature_index=temperature_index, chi=chi) ! Output
      
      
      write(10,*) a, temperature, temperature_index, chi
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="unitary_tests/temperature_interpolation.gnuplot")
    

    write(10,*) 'set terminal wxt enhanced'
    write(10,*) 'set xlabel "semi major axis a (in AU)"'

    
    write(10,*) 'set ylabel "temperature [K]"'
      
    write(10,*) 'set grid'


    write(10,*) 'plot "test_temperature_interpolation.dat" using 1:2 with lines title "Interpolation",\'
    write(10,*) '     "../temperature_profile.dat" using 1:2 with lines title "Profile"'
        

    
    write(10,*) "#pause -1 # wait until a carriage return is hit"
    write(10,*) "set terminal pdfcairo enhanced"

    
    write(10,*) '!rm "temperature_interpolation.pdf"'
    write(10,*) "set output 'temperature_interpolation.pdf'"

    
    
    write(10,*) "replot # to generate the output file"
    
    close(10)
  
  end subroutine test_temperature_interpolation
  
  subroutine test_density_interpolation()
  
    implicit none
    
    integer, parameter :: nb_a = 1000
    real(double_precision), parameter :: a_min = 0.d0 ! in AU
    real(double_precision), parameter :: a_max = 100.d0! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision) :: a, sigma, sigma_index
    
    integer :: j ! for loops
    
    write(*,*) 'Test of the surface density interpolation'
    
    open(10, file='unitary_tests/test_density_interpolation.dat')
    do j=1,nb_a
      a = (a_min + a_step * (j - 1.d0))
      ! We generate cartesian coordinate for the given semi major axis
      
      call get_surface_density(radius=a, sigma=sigma, sigma_index=sigma_index)
      
      
      write(10,*) a, sigma, sigma_index
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="unitary_tests/density_interpolation.gnuplot")
    

    write(10,*) 'set terminal wxt enhanced'
    write(10,*) 'set xlabel "semi major axis a (in AU)"'

    
    write(10,*) 'set ylabel "density [g/cm^2]"'
      
    write(10,*) 'set grid'


    write(10,*) 'plot "test_density_interpolation.dat" using 1:2 with lines title "Interpolation",\'
    write(10,*) '     "../density_profile.dat" using 1:2 with lines title "Profile"'
        

    
    write(10,*) "#pause -1 # wait until a carriage return is hit"
    write(10,*) "set terminal pdfcairo enhanced"

    
    write(10,*) '!rm "density_interpolation.pdf"'
    write(10,*) "set output 'density_interpolation.pdf'"

    
    
    write(10,*) "replot # to generate the output file"
    
    close(10)
  
  end subroutine test_density_interpolation

! %%% Physical behaviour %%%
  subroutine test_temperature_profile(stellar_mass)
! Subroutine that test the finding of the temperature profile and store a plot of the temperature profile of the disk
! A gnuplot file and a data file are created to display the temperature profile.

    implicit none
    
    real(double_precision), intent(in) :: stellar_mass
    
    integer :: j ! for loops
    
    write(*,*) 'Test of the temperature profile'
    
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
    end do

    write(10,*) "plot 'temperature_profile.dat' using 1:2 with lines notitle"
    
    write(11,*) "plot 'temperature_profile.dat' using 1:3 with lines notitle"
    

    
    do j=10,11
      write(j,*) "#pause -1 # wait until a carriage return is hit"
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) '!rm "unitary_tests/temperature_profile.pdf"'
    write(10,*) "set output 'unitary_tests/temperature_profile.pdf'"
    write(11,*) '!rm "unitary_tests/temperature_index.pdf"'
    write(11,*) "set output 'unitary_tests/temperature_index.pdf'"
    
    
    do j=10,11
      write(j,*) "replot # to generate the output file"
    end do
    
    close(10)
    close(11)
  
  end subroutine test_temperature_profile
  
  subroutine test_scaleheight_profile
    implicit none
    
    integer :: j
    
    write(*,*) 'Test of the scaleheight'
    
  open(10, file="unitary_tests/scaleheight_profile.gnuplot")
  open(11, file="unitary_tests/aspect_ratio.gnuplot")
  
  do j=10,11
    write(j,*) 'set terminal wxt enhanced'
    write(j,*) 'set xlabel "semi major axis a (in AU)"'
    write(j,*) 'set nokey'
  end do
  
  write(10,*) 'set ylabel "Scaleheight H [AU]"'
  
  write(11,*) 'set ylabel "Aspect ratio h=H/R"'
  
  do j=10,11
    write(j,*) 'set grid'
    write(j,*) 'cd ".."'
  end do

  write(10,*) "plot 'scaleheight_profile.dat' using 1:2 with line linetype -1 notitle, \"
  write(10,*) "     '' using 1:(-$2) with line linetype -1 notitle"
  
  write(11,*) "plot 'scaleheight_profile.dat' using 1:3 with line linetype -1 notitle, \"
  write(11,*) "     '' using 1:(-$3) with line linetype -1 notitle"
  

  
  do j=10,11
    write(j,*) "#pause -1 # wait until a carriage return is hit"
    write(j,*) "set terminal pdfcairo enhanced"
  end do
  
  write(10,*) '!rm "unitary_tests/scaleheight_profile.pdf"'
  write(10,*) "set output 'unitary_tests/scaleheight_profile.pdf'"
  write(11,*) '!rm "unitary_tests/aspect_ratio_profile.pdf"'
  write(11,*) "set output 'unitary_tests/aspect_ratio_profile.pdf'"
  
  
  do j=10,11
    write(j,*) "replot # to generate the output file"
  end do
  
  close(10)
  close(11)
  
  
  end subroutine test_scaleheight_profile
  
  subroutine test_optical_depth_profile(stellar_mass)
! Subroutine that test the finding of the optical depth profile and store a plot of the temperature profile of the disk
! A gnuplot file and a data file are created to display the temperature profile.

    implicit none
    
    real(double_precision), intent(in) :: stellar_mass
    
    integer :: j ! for loops
    
    
    open(10, file="unitary_tests/optical_depth_profile.gnuplot")
    
    write(*,*) 'Test of the optical depth'

    write(10,*) 'set terminal wxt enhanced'
    write(10,*) 'set xlabel "semi major axis a (in AU)"'
    write(10,*) 'set nokey'

    
    write(10,*) 'set ylabel "Optical depth {/Symbol t}"'
        
    write(10,*) 'set grid'
    write(10,*) 'cd ".."'


    write(10,*) "plot 'temperature_profile.dat' using 1:4 with lines notitle"

    
    write(10,*) "#pause -1 # wait until a carriage return is hit"
    write(10,*) "set terminal pdfcairo enhanced"
    
    write(10,*) '!rm "unitary_tests/optical_depth_profile.pdf"'
    write(10,*) "set output 'unitary_tests/optical_depth_profile.pdf'"
    
    
    write(10,*) "replot # to generate the output file"
    
    close(10)
  
  end subroutine test_optical_depth_profile
  
  subroutine test_thermal_diffusivity_profile(stellar_mass)
! Subroutine that test the finding of the optical depth profile and store a plot of the temperature profile of the disk
! A gnuplot file and a data file are created to display the temperature profile.

    implicit none
    
    real(double_precision), intent(in) :: stellar_mass
    
    integer :: j ! for loops
    
    write(*,*) 'Test of the thermal diffusivity'
    
    open(10, file="unitary_tests/thermal_diffusivity_profile.gnuplot")
    

    write(10,*) 'set terminal wxt enhanced'
    write(10,*) 'set xlabel "semi major axis a (in AU)"'
    write(10,*) 'set nokey'

    
    write(10,*) 'set ylabel "Thermal diffusivity {/Symbol c} [AU^2/day]"'
        
    write(10,*) 'set grid'
    write(10,*) 'cd ".."'


    write(10,*) "plot 'temperature_profile.dat' using 1:5 with lines notitle"

    
    write(10,*) "#pause -1 # wait until a carriage return is hit"
    write(10,*) "set terminal pdfcairo enhanced"
    
    write(10,*) '!rm "unitary_tests/thermal_diffusivity_profile.pdf"'
    write(10,*) "set output 'unitary_tests/thermal_diffusivity_profile.pdf'"
    
    
    write(10,*) "replot # to generate the output file"
    
    close(10)
  
  end subroutine test_thermal_diffusivity_profile

  subroutine test_opacity_profile
  ! subroutine that test the function 'get_opacity'
  
  ! Return:
  ! a data file 'test_opacity.dat' 
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
    
    write(*,*) 'Test of the opacity'
    
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
    write(10,*) "replot # to generate the output file"  
    
    close(10)
    
  end subroutine test_opacity_profile

  subroutine test_torques_fixed_a(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass
    
    integer, parameter :: nb_mass = 100
    real(double_precision), parameter :: mass_min = 1. ! in earth mass
    real(double_precision), parameter :: mass_max = 75. ! in earth mass
    real(double_precision), parameter :: mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    
    real(double_precision), parameter :: a = 5.2
    
    real(double_precision) :: mass, total_torque, corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: lindblad_torque_units, corotation_torque_units, total_torque_units
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    
    write(*,*) 'Evolution of the torque for a fixed distance "a"'
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
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
      write(j,*) "replot # to generate the output file"
    end do
    
    close(10)
    close(11)
    close(12)
    close(13)

    
  end subroutine test_torques_fixed_a

  subroutine test_torques_fixed_m(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass
    
    integer, parameter :: nb_a = 400
    real(double_precision), parameter :: a_min = 0.1 ! in AU
    real(double_precision), parameter :: a_max = 60. ! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision), parameter :: mass = 5. * EARTH_MASS * K2
    
    real(double_precision) :: a, total_torque, corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: lindblad_torque_units, corotation_torque_units, total_torque_units
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    
    write(*,*) 'Evolution of the torque for a fixed planet mass "m"'

    position(:) = 0.d0
    velocity(:) = 0.d0
    
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
      write(j,*) "replot # to generate the output file"
    end do
    
    close(10)
    close(11)
    close(12)

    
  end subroutine test_torques_fixed_m
  

  subroutine test_torques(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
  
    use contour
    
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass
    
    integer, parameter :: nb_mass = 150
    real(double_precision), parameter :: mass_min = 0.1 * EARTH_MASS
    real(double_precision), parameter :: mass_max = 60. * EARTH_MASS
    real(double_precision), parameter :: mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    real(double_precision), dimension(nb_mass) :: mass
    
    integer, parameter :: nb_points = 100
    real(double_precision), parameter :: a_min = 0.01
    real(double_precision), parameter :: a_max = 50.
    ! step for log sampling
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_points-1.d0)
    real(double_precision), dimension(nb_points) :: a
    
    real(double_precision), dimension(nb_points, nb_mass) :: total_torque, total_torque_units
    
    real(double_precision) :: corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    
    write(*,*) 'Evolution of the total, lindblad and corotation torques depending on the planet mass and distance'
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
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
      a(i) = a_min + a_step * (i-1)
      
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a(i)
      
      
      do j=1,nb_mass
        mass(j) = (mass_min + mass_step * (j - 1.d0)) * K2
        
        ! We generate cartesian coordinate for the given mass and semi major axis
        velocity(2) = sqrt((stellar_mass + mass(j)) / position(1))
        
        ! we store in global parameters various properties of the planet
        call get_planet_properties(stellar_mass=stellar_mass, & ! Input
         mass=mass(j), position=position(1:3), velocity=velocity(1:3),& ! Input
         p_prop=p_prop) ! Output
        call get_corotation_torque(stellar_mass, mass(j), p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        total_torque(i,j) = lindblad_torque + corotation_torque
        total_torque_units(i,j) = torque_ref * total_torque(i,j)
        
                
        write(10,*) a(i), mass(j) / (EARTH_MASS*K2), corotation_torque
        write(11,*) a(i), mass(j) / (EARTH_MASS*K2), total_torque(i,j)
        write(12,*) a(i), mass(j) / (EARTH_MASS*K2), total_torque_units(i,j)
        write(13,*) a(i), mass(j) / (EARTH_MASS*K2), lindblad_torque
        write(14,*) a(i), mass(j) / (EARTH_MASS*K2), torque_ref
        
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
    
    ! We want to get the contour for the torque equal to 0 for total torque both in physical dimension of units of gamma_0
    mass(1:nb_mass) = mass(1:nb_mass) / (EARTH_MASS*K2)
    call get_contour(total_torque, a, mass,'unitary_tests/contour_total_torque.dat', 0.d0)
!~     call get_contour(total_torque_units, a, mass,'unitary_tests/contour_total_torque_units.dat', 0.d0)
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
    write(11,*) "splot 'test_total_torque.dat' with pm3d notitle, \"
    write(11,*) "      'contour_total_torque.dat' with line linetype -1 title '{/Symbol G}=0'"
    write(12,*) "splot 'test_total_torque_units.dat' with pm3d notitle, \"
    write(12,*) "      'contour_total_torque_units.dat' with line linetype -1 title '{/Symbol G}=0'"
    write(13,*) "splot 'test_lindblad_torque.dat' with pm3d notitle"
    write(14,*) "splot 'test_ref_torque.dat' with pm3d notitle"

    
    do j=10,14
      write(j,*) "#pause -1 # wait until a carriage return is hit"
      write(j,*) "set terminal pngcairo enhanced size 1024, 768"
    end do
    
    write(10,*) "set output 'corotation_torque.png'"
    write(11,*) "set output 'total_torque.png'"
    write(12,*) "set output 'total_torque_units.png'"
    write(13,*) "set output 'lindblad_torque.png'"
    write(14,*) "set output 'ref_torque.png'"
        
    do j=10,14
      write(j,*) "replot # to generate the output file"
    end do
    
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    
  end subroutine test_torques
  
  subroutine test_dissipation_of_the_disk(stellar_mass)
  ! subroutine that plot several values depending on the properties of the disk during its dissipation
  
  ! Global parameters
  ! dissipation_timestep : the timestep between two computation of the disk [in days]
  ! INNER_BOUNDARY_RADIUS : the inner radius of the various profiles (all based on the radius profile)
  ! OUTER_BOUNDARY_RADIUS : the outer radius of the various profiles (all based on the radius profile)
  ! surface_density_profile : values of the density in log() for each value of the 'a' sample
  ! temperature_profile : values of the temperature in log() for each value of the 'a' sample
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
  
    use contour
    
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass
    
    ! mass sample
    integer, parameter :: nb_mass = 150
    real(double_precision), parameter :: mass_min = 0.1 * EARTH_MASS
    real(double_precision), parameter :: mass_max = 60. * EARTH_MASS
    real(double_precision), parameter :: mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    real(double_precision), dimension(nb_mass) :: mass, mass_earth
    
    ! orbital distance sample
    integer, parameter :: nb_points = 100
    real(double_precision), parameter :: a_min = 0.01
    real(double_precision), parameter :: a_max = 50.
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_points - 1.d0)
    real(double_precision), dimension(nb_points) :: a
    
    ! time sample
    real(double_precision), parameter :: t_min = 0. ! time in years
    real(double_precision), parameter :: t_max = 1.d7 ! time in years
    real(double_precision), dimension(:), allocatable :: time, time_temp ! time in days
    integer :: time_size ! the size of the array 'time'. 
    
    ! Array to store values of the torque. Used to find the zero torque line
    real(double_precision), dimension(nb_points, nb_mass) :: total_torque
    
    real(double_precision) :: corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    real(double_precision) :: temp_min, temp_max, density_min, density_max, torque_min, torque_max
    character(len=80) :: filename_torque, filename_density, filename_temperature, filename_contour
    character(len=80) :: output_torque, output_density, output_temperature, output_time, time_format, purcent_format
    integer :: time_length ! the length of the displayed time, usefull for a nice display
    
    integer :: i,j,k ! for loops
    integer :: nb_time ! The total number of 't' values. 
    integer :: error ! to retrieve error, especially during allocations
    !------------------------------------------------------------------------------
    write(*,*) 'Evolution of the total torque during the dissipation of the disk'
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
!~     call system("rm dissipation/*")
    
    do i=1, nb_points ! loop on the position
      a(i) = a_min + a_step * (i-1)
    end do
    
    do j=1,nb_mass
      mass(j) = (mass_min + mass_step * (j - 1.d0)) * K2
      mass_earth(j) = (mass_min + mass_step * (j - 1.d0)) / EARTH_MASS
    end do
    
    ! We want to know the max size of the time display in order to have a nice display, with filled spaces in the final plots
    write(output_time, '(f0.0)') t_max
    time_length = len(trim(output_time))
    write(time_format, *) '(f',time_length,'.0)'
    write(purcent_format, *) '(f',time_length,'.0,"/",f',time_length,'.0," years")'
    
    !------------------------------------------------------------------------------
    k = 1
    time_size = 512 ! the size of the array. 
    allocate(time(time_size), stat=error)
    time(1) = t_min * 365.25d0 ! days
    
!~     write(*,*) time(k), t_max*365.25d0
    do while (time(k).lt.(t_max*365.d0))
      ! We calculate the temperature profile for the current time (because the surface density change in function of time)
      
      ! We open the file where we want to write the outputs
      write(filename_torque, '(a,i0.5,a)') 'dissipation/total_torque',k,'.dat'
      write(filename_density, '(a,i0.5,a)') 'dissipation/surface_density',k,'.dat'
      write(filename_temperature, '(a,i0.5,a)') 'dissipation/temperature',k,'.dat'
      write(filename_contour, '(a,i0.5,a)') 'dissipation/contour',k,'.dat'
      
      open(11, file=filename_torque)
            
      ! we get the density profile.
      call calculate_density_profile()
      
      if (k.eq.time_size) then
        ! If the limit of the array is reach, we copy the values in a temporary array, allocate with a double size, et paste the old values in the new bigger array
        allocate(time_temp(time_size), stat=error)
        time_temp(1:time_size) = time(1:time_size)
        deallocate(time, stat=error)
        time_size = time_size * 2
        allocate(time(time_size), stat=error)
        time(1:time_size/2) = time_temp(1:time_size/2)
        deallocate(time_temp, stat=error)
      end if
      
      write(*,purcent_format) time(k)/365.25d0, t_max
      
      dissipation_timestep = 0.5d0 * X_SAMPLE_STEP**2 / (4 * get_viscosity(1.d0)) ! a correction factor of 0.5 has been applied. No physical reason to that, just intuition and safety
      ! TODO if the viscosity is not constant anymore, the formulae for the dissipation timestep must be changed
      time(k+1) = time(k) + dissipation_timestep * 365.25d0 ! days
      ! we get the temperature profile.
      call calculate_temperature_profile()
      
      ! we store in a .dat file the temperature profile
      call store_temperature_profile(filename=filename_temperature)
      call store_density_profile(filename=filename_density)
      
      if (k.eq.1) then
        temp_max = exp(temperature_profile(1))
        temp_min = 0.
        
        ! We want the extremum of the surface density during the dissipation of the disk in order to have nice plots
        density_min = 0.
        density_max = exp(surface_density_profile(1)) * MSUN / AU**2
      end if
      
      write(11,*) '# semi major axis (AU) ; mass in earth mass ; total torque (no dim)'
      
      do i=1, nb_points ! loop on the position
        
        ! We generate cartesian coordinate for the given semi major axis
        position(1) = a(i)
        

        
        do j=1,nb_mass
          mass(j) = (mass_min + mass_step * (j - 1.d0)) * K2
          
          ! We generate cartesian coordinate for the given mass and semi major axis
          velocity(2) = sqrt((stellar_mass + mass(j)) / position(1))
          
          ! we store in global parameters various properties of the planet
          call get_planet_properties(stellar_mass=stellar_mass, & ! Input
           mass=mass(j), position=position(1:3), velocity=velocity(1:3),& ! Input
           p_prop=p_prop) ! Output
          call get_corotation_torque(stellar_mass, mass(j), p_prop, corotation_torque, lindblad_torque, torque_ref)
          
          total_torque(i,j) = lindblad_torque + corotation_torque        
          
          write(11,*) a(i), mass_earth(j), total_torque(i,j)
          
          
          
        end do

        write(11,*) ""! we write a blank line to separate them in the data file, else, gnuplot doesn't want to make the surface plot

      end do

      close(11)
      ! We want to get the contour for the torque equal to 0 for total torque both in physical dimension of units of gamma_0
      call get_contour(total_torque, a, mass_earth,filename_contour, 0.d0)
      
      if (k.eq.1) then
        torque_min = minval(total_torque)
        torque_max = maxval(total_torque)
      end if
      
      k = k + 1 ! We increment the integer that point the time in the array (since it's a 'while' and not a 'do' loop)
    end do
    
    nb_time = k - 1 ! since for the last step with incremented 'k' by one step that is beyond the limit.
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the total torque
    open(11, file="dissipation/total_torque.gnuplot")
    write(11,*) "set terminal pngcairo enhanced size 1024, 768"
    write(11,*) 'set xlabel "semi major axis (AU)"'
    write(11,*) 'set ylabel "Planet mass (m_{earth})" center'
    write(11,*) 'set pm3d map'
    write(11,*) 'set pm3d explicit'
    write(11,*) 'set palette rgbformulae 22,13,-31'
    write(11,*) 'set mxtics 5'
    write(11,*) 'set mytics 5'
    write(11,*) 'set grid xtics ytics mxtics mytics linetype -1, linetype 0'
    write(11,*) 'set xrange [', a_min, ':', a_max, ']'
    write(11,*) 'set yrange [', mass_min / EARTH_MASS, ':', mass_max / EARTH_MASS, ']'
    write(11,*) 'set cbrange [', torque_min, ':', torque_max, ']'
    
    do k=1, nb_time
      write(filename_torque, '(a,i0.5,a)') 'total_torque',k,'.dat'
      write(filename_contour, '(a,i0.5,a)') 'contour',k,'.dat'
      write(output_torque, '(a,i0.5,a)') 'total_torque',k,'.png'
      write(output_time, time_format) time(k)/365.25
      
      write(11,*) "set output '",trim(output_torque),"'"
      write(11,*) 'set title "total torque {/Symbol G}_{tot}/{/Symbol G}_0 T=', trim(output_time),' years"'
      write(11,*) "splot '",trim(filename_torque),"' with pm3d notitle, \"
      write(11,*) "      '",trim(filename_contour),"' with line linetype -1 title '{/Symbol G}=0'"
      write(11,*) ""
    end do
    close(11)
    
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the temperature profile
    open(12, file="dissipation/temperature.gnuplot")
    write(12,*) "set terminal pngcairo enhanced size 1024, 768"
    write(12,*) 'set xlabel "semi major axis (AU)"'
    write(12,*) 'set ylabel "Temperature (K)"'
    write(12,*) 'set grid'
    write(12,*) 'set xrange [', a_min, ':', a_max, ']'
    write(12,*) 'set yrange [', temp_min, ':', temp_max, ']'
    
    do k=1, nb_time
      write(filename_temperature, '(a,i0.5,a)') 'temperature',k,'.dat'
      write(output_temperature, '(a,i0.5,a)') 'temperature',k,'.png'
      write(output_time, time_format) time(k)/365.25
      
      write(12,*) "set output '",trim(output_temperature),"'"
      write(12,*) 'set title "T=', trim(output_time),' years"'
      write(12,*) "plot '",trim(filename_temperature),"' using 1:2 with lines notitle"
      write(12,*) ""
    end do
    close(12)
    
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the density
    open(13, file="dissipation/density.gnuplot")
    write(13,*) "set terminal pngcairo enhanced size 1024, 768"
    write(13,*) 'set xlabel "semi major axis (AU)"'
    write(13,*) 'set ylabel "Surface density (g/cm^2)"'
    write(13,*) 'set grid'
    write(13,*) 'set xrange [', INNER_BOUNDARY_RADIUS, ':', OUTER_BOUNDARY_RADIUS, ']'
    write(13,*) 'set yrange [', density_min, ':', density_max, ']'
    
    do k=1, nb_time
      write(filename_density, '(a,i0.5,a)') 'surface_density',k,'.dat'
      write(output_density, '(a,i0.5,a)') 'surface_density',k,'.png'
      write(output_time, time_format) time(k)/365.25
      
      write(13,*) "set output '",trim(output_density),"'"
      write(13,*) 'set title "T=', trim(output_time),' years"'
      write(13,*) "plot '",trim(filename_density),"' using 1:2 with lines notitle"
      write(13,*) ""
    end do
    close(13)
  
  end subroutine test_dissipation_of_the_disk
end module user_module
  
! TODO utiliser la masse des objets pour ne pas faire le calcul si trop massif, il faut respecter le domaine de validité des formules des couples
! TODO routine générale de conversion des couples en accélération afin de pouvoir réutiliser ailleurs

!TODO : 
!_rajouter le spin comme variable d'entrée de la routine mfo_user
