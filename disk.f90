module disk

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
  integer :: DISSIPATION_TYPE = 1 ! integer to tell if there is dissipation of the disk or not. 0 for no dissipation, 1 for viscous dissipation and 2 for exponential decay of the initial profile
  real(double_precision) :: TAU_DISSIPATION = -1.d0 ! the characteristic time for the exponential decay of the surface density (in years)
  real(double_precision) :: dissipation_timestep ! the timestep between two computation of the disk [in days]
  character(len=80) :: INNER_BOUNDARY_CONDITION = 'closed' ! 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid
  character(len=80) :: OUTER_BOUNDARY_CONDITION = 'closed' ! 'open' or 'closed'. If open, gas can cross the outer edge. If closed, nothing can escape the grid
  character(len=80) :: TORQUE_TYPE = 'real' ! values possible to change the properties of the torque. 'real', 'mass_independant', 'mass_dependant'
  ! Here we define the constant value of the viscosity of the disk
  real(double_precision) :: viscosity = 1.d15 ! viscosity of the disk [cm^2/s]
  
  !------------------------------------------------------------------------------
  ! Values for manual torque profiles (mass dependant or independant for instance)
  real(double_precision) :: TORQUE_PROFILE_STEEPNESS = 1. ! increase, in units of Gamma_0 of the torque per 10AU
  
  ! For mass independant convergence zone
  real(double_precision) :: INDEP_CZ = 8. ! Position of the mass independant convergence zone in AU
  
  ! For mass dependant convergence zone
  ! boundary mass values that have a zero torque zone
  real(double_precision) :: MASS_DEP_M_MIN = 1.  ! Minimum mass for the mass dependant convergence zone (in earth mass)
  real(double_precision) :: MASS_DEP_M_MAX = 60. ! Maximum mass for the mass dependant convergence zone (in earth mass)
  
  ! position of the zero torque zone for the minimum and maximum mass
  real(double_precision) :: MASS_DEP_CZ_M_MIN = 4.  ! position of the mass dependant convergence zone for the minimum mass (in AU)
  real(double_precision) :: MASS_DEP_CZ_M_MAX = 30. ! position of the mass dependant convergence zone for the maximum mass (in AU)
  
  !------------------------------------------------------------------------------
  ! prefactors
  real(double_precision) :: X_S_PREFACTOR ! prefactor for the half width of the corotation region
  real(double_precision) :: SCALEHEIGHT_PREFACTOR ! prefactor for the scaleheight
  
  !------------------------------------------------------------------------------
  ! Here we define properties common to the profiles
  real(double_precision) :: INNER_BOUNDARY_RADIUS = 1.d0
  real(double_precision) :: OUTER_BOUNDARY_RADIUS = 100.d0
  integer :: NB_SAMPLE_PROFILES = 200 ! number of points for the sample of radius of the temperature profile
  real(double_precision) :: X_SAMPLE_STEP ! the constant step for the x_sample. Indeed, due to diffusion equation, the sample must be constant in X, and not in r. 

  real(double_precision), dimension(:), allocatable :: distance_sample ! values of 'a' in AU
  real(double_precision), dimension(:), allocatable :: x_sample ! values of 'x' with x = 2*sqrt(r) (used for the diffusion equation)
  real(double_precision), dimension(:), allocatable :: surface_density_profile ! values of the density in MSUN/AU^2 for each value of the 'a' sample
  real(double_precision), dimension(:), allocatable :: surface_density_index ! values of the local negative slope of the surface density profile
  real(double_precision), dimension(:), allocatable :: temperature_profile ! values of the temperature in K for each value of the 'a' sample
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
  character(len=1), parameter :: SEP = '=' ! the separator of a parameter line

  integer :: sep_position ! an integer to get the position of the separator

  !------------------------------------------------------------------------------

  sep_position = index(line, SEP)
  
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
! DISSIPATION_TYPE : integer to tell if there is dissipation of the disk or not. 0 for no dissipation, 1 for viscous dissipation and 2 for exponential decay of the initial profile
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
        
        case('adiabatic_index')
          read(value, *) ADIABATIC_INDEX
          
        case('mean_molecular_weight')
          read(value, *) MEAN_MOLECULAR_WEIGHT
          
        case('surface_density')
          read(value, *) INITIAL_SIGMA_0, INITIAL_SIGMA_INDEX

        case('disk_edges')
          read(value, *) INNER_BOUNDARY_RADIUS, OUTER_BOUNDARY_RADIUS
        
        case('sample')
          read(value, *) NB_SAMPLE_PROFILES
          
        case('viscosity')
          read(value, *) viscosity
          
        case('dissipation_type')
          read(value, *) DISSIPATION_TYPE
        
        case('disk_exponential_decay')
          read(value, *) TAU_DISSIPATION
          DISSIPATION_TYPE = 2
          
        case('inner_boundary_condition')
          read(value, *) INNER_BOUNDARY_CONDITION
        
        case('outer_boundary_condition')
          read(value, *) OUTER_BOUNDARY_CONDITION
          
        case('torque_type')
          read(value, *) TORQUE_TYPE
          
        case default
          write(*,*) 'Warning: An unknown parameter has been found'
          write(*,*) "identificator='", trim(identificator), "' ; value(s)='", trim(value),"'"
        end select
      end if
    end do
    close(10)
    
    ! problem if exponential decay timescale is set, but the dissipation_type is not exponentiel decay
    if ((TAU_DISSIPATION.gt.0.d0).and.(DISSIPATION_TYPE.ne.2)) then
      write(*,*) 'Warning: Exponential decay of surface density must exist only if Dissipation method is "exponential decay"'
      DISSIPATION_TYPE = 2
    end if
    
    ! problem is we want exponential decay of the surface density but we did not set the exponential decay timescale
    if ((TAU_DISSIPATION.lt.0.d0).and.(DISSIPATION_TYPE.eq.2)) then
      write(*,*) 'Error: since dissipation_type=2, "disk_exponential_decay" is expected to be set.'
      stop
    end if
  else
    write (*,*) 'Warning: The file "disk.in" does not exist. Default values have been used'
  end if
  
  ! Whether we modified or not 'INITIAL_SIGMA_0' we calculate the value in numerical units of the initial sigma_0 value.
  INITIAL_SIGMA_0_NUM = INITIAL_SIGMA_0 * AU**2 / MSUN ! the surface density at (R=1AU) [Msun/AU^2]
  
end subroutine read_disk_properties

subroutine write_disk_properties()
! subroutine that write the parameters of the user_module into the file 'disk.out'

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
! DISSIPATION_TYPE : integer to tell if there is dissipation of the disk or not. 0 for no dissipation, 1 for viscous dissipation and 2 for exponential decay of the initial profile
! INNER_BOUNDARY_CONDITION : 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid
! OUTER_BOUNDARY_CONDITION : 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid

  implicit none
  
  open(10, file='disk.out')
  
  write(10,'(a)') '------------------------------------'
  write(10,'(a)') '|   Properties of the disk         |'
  write(10,'(a)') '------------------------------------'
  write(10,'(a,f4.2)') 'b/h = ',B_OVER_H
  write(10,'(a,f4.2)') 'adiabatic index = ', ADIABATIC_INDEX
  write(10,'(a,f4.2)') 'mean molecular weight = ', MEAN_MOLECULAR_WEIGHT
  write(10,'(a,es10.1e2,a)') 'viscosity = ', viscosity, ' (cm^2/s)'
  write(10,'(a,f6.1,a,f4.2 ,a)') 'initial surface density = ',INITIAL_SIGMA_0, ' * R^(-',INITIAL_SIGMA_INDEX,') (g/cm^2)'
  write(10,*) ''
  write(10,'(a,f6.1,a)') 'inner edge of the disk = ',INNER_BOUNDARY_RADIUS, ' (AU)'
  write(10,'(a,f6.1,a)') 'outer edge of the disk = ',OUTER_BOUNDARY_RADIUS, ' (AU)'
  write(10,*) ''
  write(10,'(a)') 'Possible values : &
  0 for no dissipation, 1 for viscous dissipation and 2 for exponential decay of the initial profile'
  write(10,'(a,i1)') 'dissipation of the disk = ',DISSIPATION_TYPE
  write(10,'(a)') '  Case(1) : viscous dissipation'
  write(10,'(a)') "    Possible values : 'open', 'closed'"
  write(10,'(a, a)') '    inner boundary condition = ', trim(INNER_BOUNDARY_CONDITION)
  write(10,'(a, a)') '    outer boundary condition = ', trim(OUTER_BOUNDARY_CONDITION)
  write(10,'(a)') '  Case(2) : exponential decay'
  write(10,'(a,es10.1e2)') '    characteristic time for decay = ',TAU_DISSIPATION
  write(10,*) ''
  write(10,'(a)') '------------------------------------'
  write(10,'(a)') '|     Interactions disk/planets     |'
  write(10,'(a)') '------------------------------------'
  write(10,'(a)') 'When the inclination damping stops'
  write(10,'(a,es10.1e2,a)') 'inclination cutoff = ',INCLINATION_CUTOFF, ' (rad)'
  write(10,*) ''
  write(10,'(a)') "Possible values : 'real', 'mass_independant', 'mass_dependant'"
  write(10,'(a, a)') 'torque type = ', trim(TORQUE_TYPE)
  write(10,'(a)') '  Case(mass_independant) : manual torque profile with a mass independant convergence zone'
  write(10,'(a,f6.1,a)') '    steepness of the torque profile = ', TORQUE_PROFILE_STEEPNESS, ' (unit/10 AU)'
  write(10,'(a,f6.1,a)') '    position of the convergence zone = ', INDEP_CZ, ' (AU)'
  write(10,'(a)') '  Case(mass_dependant) : manual torque profile with a mass independant convergence zone'
  write(10,'(a,f6.1,a)') '    minimum mass = ', MASS_DEP_M_MIN, ' (earth mass)'
  write(10,'(a,f6.1,a)') '    maximum mass = ', MASS_DEP_M_MAX, ' (earth mass)'
  write(10,'(a,f6.1,a)') '    CZ for minimum mass = ', MASS_DEP_CZ_M_MIN, ' (AU)'
  write(10,'(a,f6.1,a)') '    CZ for maximum mass = ', MASS_DEP_CZ_M_MAX, ' (AU)'
  write(10,*) ''
  close(10)
  
end subroutine write_disk_properties

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
  
  ! Since the sigma_index is dependant on the location of the planet and the time (thanks to the dissipation), we need to calculate these values at each timestep
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
! distance_sample : values of 'a' in AU
! x_sample : values of 'x' with x = 2*sqrt(r) (used for the diffusion equation)
! surface_density_profile : values of the density in MSUN/AU^2 for each value of the 'a' sample
! surface_density_index : values of the local negative slope of the surface density profile
! temperature_profile : values of the temperature in K for each value of the 'a' sample
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
    ! we write all the values used by user_module, those given by the user, and the default ones, in 'disk.out' file
    call write_disk_properties() 
    
    allocate(distance_sample(NB_SAMPLE_PROFILES))
    distance_sample(1:NB_SAMPLE_PROFILES) = 0.d0
    
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
    distance_sample(1) = INNER_BOUNDARY_RADIUS
    
    x_sample(NB_SAMPLE_PROFILES) = 2.d0 * sqrt(OUTER_BOUNDARY_RADIUS)
    distance_sample(NB_SAMPLE_PROFILES) = OUTER_BOUNDARY_RADIUS
    
    ! We initialize the global variable (in the module) for the constant step of x_sample
    X_SAMPLE_STEP = (x_sample(NB_SAMPLE_PROFILES) - x_sample(1)) / (NB_SAMPLE_PROFILES - 1.d0)
    
    do i=2, NB_SAMPLE_PROFILES - 1
      x_sample(i) = x_sample(1) + X_SAMPLE_STEP * (i - 1.d0)
      distance_sample(i) = 0.25d0 * x_sample(i)**2
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
  ! distance_sample : values of 'a' in AU
  ! surface_density_profile : values of the density in MSUN/AU^2 for each value of the 'a' sample
  ! surface_density_index : values of the local negative slope of the surface density profile
  
  implicit none
  
  integer :: i ! for loops
!~   write(*,*) 'Warning: the initial profil is linear for tests of viscous dissipation!'
!~   do i=1,NB_SAMPLE_PROFILES
!~     surface_density_profile(i)=INITIAL_SIGMA_0_NUM*(1-distance_sample(i)/distance_sample(NB_SAMPLE_PROFILES)) ! linear decay of the surface density, for tests
!~   end do
  
  do i=1,NB_SAMPLE_PROFILES
    surface_density_profile(i) = INITIAL_SIGMA_0_NUM * distance_sample(i)**(-INITIAL_SIGMA_INDEX)
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
    ! function that interpolate the value of the surface density at the given radius. 
    ! A linear interpolation is used, in particuliar to avoid problems with log(sigma) when the surface density is 0.
    
    ! Global parameter
    ! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
    ! X_SAMPLE_STEP : the constant step for the x_sample. Indeed, due to diffusion equation, the sample must be constant in X, and not in r. 
    ! INNER_BOUNDARY_RADIUS : the inner edge of the different profiles
    ! OUTER_BOUNDARY_RADIUS : the outer edge of the different profiles
    ! distance_sample : values of 'a' in AU
    ! x_sample : values of 'x' with x = 2*sqrt(r) (used for the diffusion equation)
    ! surface_density_profile : values of the density in MSUN/AU^2 for each value of the 'a' sample
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
    real(double_precision) :: x1, x2, y1, y2
    real(double_precision) :: x_radius ! the corresponding 'x' value for the radius given in parameter of the routine. we can retrieve the index of the closest values in this array in only one calculation.

    if ((radius .ge. INNER_BOUNDARY_RADIUS) .and. (radius .lt. OUTER_BOUNDARY_RADIUS)) then
      
      x_radius = 2.d0 * sqrt(radius)
      ! in the range
      closest_low_id = 1 + int((x_radius - x_sample(1)) / X_SAMPLE_STEP) ! X_SAMPLE_STEP being a global value, x_sample also
      
      x1 = distance_sample(closest_low_id)
      x2 = distance_sample(closest_low_id + 1)
      y1 = surface_density_profile(closest_low_id)
      y2 = surface_density_profile(closest_low_id + 1)

      sigma = y2 + (y1 - y2) * (radius - x2) / (x1 - x2)
      sigma_index = surface_density_index(closest_low_id) ! for the temperature index, no interpolation.
    else if (radius .lt. INNER_BOUNDARY_RADIUS) then
      sigma = surface_density_profile(1)
      sigma_index = surface_density_index(1)
    else if (radius .gt. OUTER_BOUNDARY_RADIUS) then
      sigma = surface_density_profile(NB_SAMPLE_PROFILES)
      sigma_index = surface_density_index(NB_SAMPLE_PROFILES)
    end if
  end subroutine get_surface_density
  
  function get_viscosity(radius)
  ! function that return the viscosity of the disk in [AU^2.day^-1]
  
  ! Parameters
  ! radius : The orbital distance in AU
  
  ! Global parameters
  ! viscosity : the viscosity of the disk in [cm^2/s]
  
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
! distance_sample : values of 'a' in AU
! temperature_profile : values of the temperature in K for each value of the 'a' sample
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
  
    
    a = distance_sample(1)
    ! We generate cartesian coordinate for the given semi major axis
    position(1) = a
    
    ! We generate cartesian coordinate for the given mass and semi major axis
    velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
    
    ! we store in global parameters various properties of the planet
    call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
     p_prop=p_prop) ! Output
     
    call zbrent(x_min=1.d-5, x_max=1.d4, tolerance=1d-4, p_prop=p_prop, & ! Input
                            temperature=temperature, optical_depth=tau_profile(1)) ! Output    

    temperature_profile(1) = temperature
    chi_profile(1) = 1.5d0 * p_prop%nu * ADIABATIC_INDEX * (ADIABATIC_INDEX - 1.d0) * &
                      (1.5d0 + sqrt(3.d0) / tau_profile(1) + 1 / tau_profile(1)**2)
    do j=2,NB_SAMPLE_PROFILES
      a_old = a
      temperature_old = temperature
      
      a = distance_sample(j) ! Be carefull, the step between 'a' values is not constant !
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a
      
      ! We generate cartesian coordinate for the given mass and semi major axis
      velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
      
      call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
       p_prop=p_prop) ! Output
      
      call zbrent(x_min=1.d-5, x_max=1.d4, tolerance=1d-4, p_prop=p_prop, & ! Input
                              temperature=temperature, optical_depth=tau_profile(j)) ! Output
      
      temperature_profile(j) = temperature
      temp_profile_index(j) = - (log(temperature_profile(j)) - log(temperature_profile(j-1))) / &
                                (log(distance_sample(j)) - log(distance_sample(j-1)))
      chi_profile(j) = 1.5d0 * p_prop%nu * ADIABATIC_INDEX * (ADIABATIC_INDEX - 1.d0) * &
                      (1.5d0 + sqrt(3.d0) / tau_profile(j) + 1 / tau_profile(j)**2)
      
    end do
    
    temp_profile_index(2) = temp_profile_index(3)! We do this because the first two values are polluted by boundary conditions
    temp_profile_index(1) = temp_profile_index(2) ! the first element of the index array is forced to be equal to the second one because else, we can't calculate it.
    
  end subroutine calculate_temperature_profile
  
  subroutine exponential_decay_density_profile()
! subroutine that calculate new surface density profile given the old one and a timestep for the exponential decay
! 
! Global Parameters 
! dissipation_timestep : the timestep between two computation of the disk [in days]
! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
! surface_density_profile : values of the density in MSUN/AU^2 for each value of the 'a' sample
! surface_density_index : values of the local negative slope of the surface density profile
!
! Goal
! This routine will update the surface_density_index and surface_density_profile for the next timestep dissipation_timestep of the disk
! With exponential decay, the index of the slope is not supposed to change, so the values are not modified for surface_density_index
    use mercury_constant

    implicit none
    
    ! Local values
    integer :: i
    
    !------------------------------------------------------------------------------
    
    do i=1,NB_SAMPLE_PROFILES
      surface_density_profile(i) = surface_density_profile(i) * exp(- dissipation_timestep / (TAU_DISSIPATION * 365.25d0))
!~       surface_density_index(i) = - (log(surface_density_profile(i)) - log(surface_density_profile(i-1))) &
!~                                   / (log(distance_sample(i)) - log(distance_sample(i-1)))
    end do
  
  end subroutine exponential_decay_density_profile
  
  subroutine dissipate_density_profile()
! subroutine that calculate the temperature profile of the disk given various parameters including the surface density profile.
! 
! Global Parameters 
! INNER_BOUNDARY_CONDITION : 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid
! OUTER_BOUNDARY_CONDITION : 'open' or 'closed'. If open, gas can fall on the star. If closed, nothing can escape the grid
! dissipation_timestep : the timestep between two computation of the disk [in days]
! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
! x_sample : values of 'x' with x = 2*sqrt(r) (used for the diffusion equation)
! surface_density_profile : values of the density in MSUN/AU^2 for each value of the 'a' sample
! surface_density_index : values of the local negative slope of the surface density profile
!
! Goal
! This routine will update the surface_density_index and surface_density_profile for the next timestep dissipation_timestep of the disk

    use mercury_constant

    implicit none

    ! value for the precedent step of the loop. In order to calculate the index of the local temperature power law.
    real(double_precision) :: a_im12, a_ip12 ! ip12 == (i+1/2) ; im12 == (i-1/2)
    real(double_precision) :: f_i, f_ip1 ! ip1 == (i+1) ; im1 == (i-1)
    real(double_precision) :: flux_im12, flux_ip12 ! ip12 == (i+1/2) ; im12 == (i-1/2)
    integer :: i ! for loops
    real(double_precision) :: tmp ! for temporary values
    
    ! The dissipation_timestep is a global variable that is the timestep between two calculation of the surface density profile 
    ! This value is calculated in mfo_user, in the 'if' loop where the re-calculation is done. This value is needed here in the 
    ! part that dissipate the disk, but is not needed for the first run of this routine.
    
    ! for i=1
    ! surface_density_index(1) is to be defined later, using surface_density_index(2)
    
    ! We prepare values for the first step of the loop 
    ! (i.e i=1 because we simulate the ghost step of the loop to be able to shift the temp values)
    a_ip12 = 0.125d0 * (x_sample(1) * x_sample(1) + x_sample(2) * x_sample(2))
    
    f_i   = 1.5d0 * x_sample(1) * surface_density_profile(1)
    f_ip1 = 1.5d0 * x_sample(2) * surface_density_profile(2)
    
    flux_ip12 = 3.d0 * get_viscosity(a_ip12) / a_ip12 * (f_ip1 - f_i) / X_SAMPLE_STEP
    
    ! CONDITION AT THE INNER EDGE
    select case(INNER_BOUNDARY_CONDITION)
      case('closed') 
!~       ! correspond to the case where the velocity at the inner edge is forced to be zero. This is equivalent to a 0-flux condition
!~         flux_ip12 = 0.d0
        ! If closed condition, the surface density at the boundary 'surface_density_profile(1)' doesn't change. So we do not compute anything.
        tmp = (f_i + dissipation_timestep * flux_ip12 / X_SAMPLE_STEP) 
      
        surface_density_profile(1) = tmp / (1.5d0 * x_sample(1)) ! the (1.5d0 * x_i) is here to convert from 'f' to Sigma
      case('open') 
      ! correspond to the case where the surface density is forced to be zero. This is equivalent to an accretion 
      ! condition. Density is free to dissipate outside the grid.
        surface_density_profile(1) = 0.d0 ! in log, '0' is not defined, so we put a huge negative value to get close to 0
        f_i = 0.d0 ! We impose the condition sigma=0 so in practice it should be equal to 0
      case default
        write(*,*) 'Warning: An unknown inner boundary condition has been found'
        write(*,*) "inner_boundary_condition=", trim(INNER_BOUNDARY_CONDITION)
    end select
    
    ! We store the previous values in order to avoid the use of an array. 
    ! This way, it should more efficient, because we don't have to create several arrays with thousands of elements
    do i=2,NB_SAMPLE_PROFILES-2
      
      ! We shift the indexes by 1
      a_ip12 = 0.125d0 * (x_sample(i) * x_sample(i) + x_sample(i+1) * x_sample(i+1))
      
      f_i = f_ip1
      f_ip1 = 1.5d0 * x_sample(i+1) * surface_density_profile(i+1)
      
      flux_im12 = flux_ip12
      flux_ip12 = 3.d0 * get_viscosity(a_ip12) / a_ip12 * (f_ip1 - f_i) / X_SAMPLE_STEP
      
      tmp = (f_i + dissipation_timestep * (flux_ip12 - flux_im12) / X_SAMPLE_STEP)
      
      if (tmp.lt.0.) then
        write(*,*) 'ERROR: tmp and thus the surface density is negative!!!!'
        write(*,*) a_ip12, f_i, f_ip1, flux_im12, flux_ip12
        stop
      end if
      
      surface_density_profile(i) = tmp / (1.5d0 * x_sample(i)) ! the (1.5d0 * x_i) is here to convert from 'f' to Sigma
      surface_density_index(i) = - (log(surface_density_profile(i)) - log(surface_density_profile(i-1))) &
                                    / (log(distance_sample(i)) - log(distance_sample(i-1)))
            
    end do
    !------------------------------------------------------------------------------
    i=NB_SAMPLE_PROFILES-1
    ! We shift the indexes by 1
    f_i = f_ip1
    
    flux_im12 = flux_ip12
    
    a_ip12 = 0.125d0 * (x_sample(i) * x_sample(i) + x_sample(i+1) * x_sample(i+1))

    f_ip1 = 1.5d0 * x_sample(i+1) * surface_density_profile(i+1)

    flux_ip12 = 3.d0 * get_viscosity(a_ip12) / a_ip12 * (f_ip1 - f_i) / X_SAMPLE_STEP
        
    ! The boundary condition is computed here because some values are needed by "i=NB_SAMPLE_PROFILES-1 surface density value".
    select case(OUTER_BOUNDARY_CONDITION)
      case('closed') 
      ! here i=NB_SAMPLE_PROFILES-2
!~       ! correspond to the case where the velocity at the inner edge is forced to be zero. This is equivalent to a 0-flux condition
!~         flux_ip12 = 0.d0
        ! If closed condition, the surface density at the boundary 'surface_density_profile(1)' doesn't change. So we do not compute anything.
        tmp = (f_i - dissipation_timestep * flux_im12 / X_SAMPLE_STEP)
      
        surface_density_profile(NB_SAMPLE_PROFILES) = tmp / (1.5d0 * x_sample(NB_SAMPLE_PROFILES)) ! the (1.5d0 * x_i) is here to convert from 'f' to Sigma
      case('open') 
      ! correspond to the case where the surface density is forced to be zero. This is equivalent to an accretion 
      ! condition. Density is free to dissipate outside the grid.
        surface_density_profile(NB_SAMPLE_PROFILES) = 0.d0 ! in log, '0' is not defined, so we put a huge negative value to get close to 0
        
      case default
        write(*,*) 'Warning: An unknown outer boundary condition has been found'
        write(*,*) "outer_boundary_condition=", trim(OUTER_BOUNDARY_CONDITION)
    end select    
    
    tmp = (f_i + dissipation_timestep * (flux_ip12 - flux_im12) / X_SAMPLE_STEP) / (1.5d0 * x_sample(i))
      
    if (tmp.lt.0.) then
      write(*,*) 'ERROR: tmp is negative!!!!'
      write(*,*) a_ip12, f_i, f_ip1, flux_im12, flux_ip12
    end if
    
    surface_density_profile(i) = tmp ! the (1.5d0 * x_i) is here to convert from 'f' to Sigma
    surface_density_index(i) = - (log(surface_density_profile(i)) - log(surface_density_profile(i-1))) &
                                  / (log(distance_sample(i)) - log(distance_sample(i-1)))
    
    !------------------------------------------------------------------------------
    ! And the last index NB_SAMPLE_PROFILES 
    ! We DO NOT compute the surface density since it's ruled by the boundary condition. But we compute the index
    i=NB_SAMPLE_PROFILES
    surface_density_index(i) = - (log(surface_density_profile(i)) - log(surface_density_profile(i-1))) &
                                    / (log(distance_sample(i)) - log(distance_sample(i-1)))

    !------------------------------------------------------------------------------
    ! We copy paste the index for the first value because however it is not defined.
    surface_density_index(2) = surface_density_index(3) ! We do this because the first two values are polluted by boundary conditions
    surface_density_index(1) = surface_density_index(2)
  
  end subroutine dissipate_density_profile

  subroutine store_temperature_profile(filename)
  ! subroutine that store in a '.dat' file the temperature profile and negative index of the local power law
  
  ! Global parameters
  ! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
  ! temperature_profile : values of the temperature in K for each value of the 'a' sample
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
    write(10,*) distance_sample(j), temperature_profile(j), temp_profile_index(j), tau_profile(j), chi_profile(j)!, distance_sample(j), temperature_profile(j)
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
    write(10,*) distance_sample(j), surface_density_profile(j) * NUM2PHYS, surface_density_index(j)
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
    a = distance_sample(j)
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
  call print_planet_properties(p_prop)
  stop
end if 

if (isnan(p_prop%nu)) then
  write(*,*) 'Error: the viscosity is equal to NaN when we want to calculate the temperature profile'
  call print_planet_properties(p_prop)
  stop
end if 

if (isnan(p_prop%omega)) then
  write(*,*) 'Error: the angular velocity is equal to NaN when we want to calculate the temperature profile'
  call print_planet_properties(p_prop)
  stop
end if 

if (isnan(p_prop%radius)) then
  write(*,*) 'Error: the distance is equal to NaN when we want to calculate the temperature profile'
  call print_planet_properties(p_prop)
  stop
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
  write(*,*) 'properties of the disk at the location of the planet that influence the value of the temperature'
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
! distance_sample : values of 'a' in AU
! temperature_profile : values of the temperature in K for each value of the 'a' sample
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
  
  ln_x1 = log(distance_sample(closest_low_id))
  ln_x2 = log(distance_sample(closest_low_id + 1))
  ln_y1 = log(temperature_profile(closest_low_id))
  ln_y2 = log(temperature_profile(closest_low_id + 1))

  temperature = exp(ln_y2 + (ln_y1 - ln_y2) * (log(radius) - ln_x2) / (ln_x1 - ln_x2))
  temperature_index = temp_profile_index(closest_low_id) ! for the temperature index, no interpolation.
  chi = chi_profile(closest_low_id)
else if (radius .lt. INNER_BOUNDARY_RADIUS) then
  temperature = temperature_profile(1)
  temperature_index = temp_profile_index(1)
  chi = chi_profile(1)
else if (radius .gt. OUTER_BOUNDARY_RADIUS) then
  temperature = temperature_profile(NB_SAMPLE_PROFILES)
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
  
  ! %%% Local modifications of the code %%%
subroutine get_corotation_torque_mass_indep_CZ(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, Gamma_0)
! function that return the total torque exerted by the disk on the planet 
!
! Global parameters
! ADIABATIC_INDEX : the adiabatic index for the gas equation of state
! X_S_PREFACTOR : prefactor for the half width of the corotation region
! TORQUE_PROFILE_STEEPNESS : increase, in units of Gamma_0 of the torque per 10AU
! INDEP_CZ : Position of the mass independant convergence zone in AU


  implicit none
  real(double_precision), intent(in) :: stellar_mass ! the mass of the central body [Msun * K2]
  ! Properties of the planet
  real(double_precision), intent(in) :: mass ! the mass of the planet [Msun * K2]
  type(PlanetProperties), intent(in) :: p_prop ! various properties of the planet
  
  
  real(double_precision), intent(out) :: corotation_torque
  real(double_precision), intent(out) :: lindblad_torque !  lindblad torque exerted by the disk on the planet [\Gamma_0]
  real(double_precision), intent(out) :: Gamma_0 ! canonical torque value [Ms.AU^2](equation (8) of Paardekooper, Baruteau, 2009)
  
  !------------------------------------------------------------------------------
  ! WE CALCULATE TOTAL TORQUE EXERTED BY THE DISK ON THE PLANET
  Gamma_0 = (mass / (stellar_mass * p_prop%aspect_ratio))**2 * p_prop%sigma * p_prop%radius**4 * p_prop%omega**2
  
  
  lindblad_torque = 0.
  
  !------------------------------------------------------------------------------

  corotation_torque = TORQUE_PROFILE_STEEPNESS * 0.1d0 * (INDEP_CZ - p_prop%radius)
  

  return
end subroutine get_corotation_torque_mass_indep_CZ

subroutine get_corotation_torque_mass_dep_CZ(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, Gamma_0)
! function that return the total torque exerted by the disk on the planet 
!
! Global parameters
! ADIABATIC_INDEX : the adiabatic index for the gas equation of state
! X_S_PREFACTOR : prefactor for the half width of the corotation region
! TORQUE_PROFILE_STEEPNESS : increase, in units of Gamma_0 of the torque per 10AU
! MASS_DEP_M_MIN = 1.  ! Minimum mass for the mass dependant convergence zone (in earth mass)
! MASS_DEP_M_MAX = 60. ! Maximum mass for the mass dependant convergence zone (in earth mass)
! MASS_DEP_CZ_M_MIN = 4.  ! position of the mass dependant convergence zone for the minimum mass (in AU)
! MASS_DEP_CZ_M_MAX = 30. ! position of the mass dependant convergence zone for the maximum mass (in AU)

  implicit none
  real(double_precision), intent(in) :: stellar_mass ! the mass of the central body [Msun * K2]
  ! Properties of the planet
  real(double_precision), intent(in) :: mass ! the mass of the planet [Msun * K2]
  type(PlanetProperties), intent(in) :: p_prop ! various properties of the planet
  
  
  real(double_precision), intent(out) :: corotation_torque
  real(double_precision), intent(out) :: lindblad_torque !  lindblad torque exerted by the disk on the planet [\Gamma_0]
  real(double_precision), intent(out) :: Gamma_0 ! canonical torque value [Ms.AU^2](equation (8) of Paardekooper, Baruteau, 2009)
   
  ! coeff for the function that give the position of the convergence zone in function of mass
  real(double_precision) :: a ! The mass will be given in solar mass. So we add a corrective factor to get planet mass in earth mass
  real(double_precision) :: b
  
  ! Local 
  real(double_precision) :: planet_mass ! the mass of the current planet in earth mass
  !------------------------------------------------------------------------------
  ! position of zero torque zone in function of the mass : 
  ! r(m_min) = a_min
  ! r(m_max) = a_max
  !
  ! r(m) = a * m + b
  ! a = (a_max - a_min) / (m_max - m_min)
  ! b = a_min - m_min * a
  
  a = (MASS_DEP_CZ_M_MAX - MASS_DEP_CZ_M_MIN) / ((MASS_DEP_M_MAX - MASS_DEP_M_MIN))
  b = MASS_DEP_CZ_M_MIN - MASS_DEP_M_MIN * a
  
  planet_mass = mass / (3.00374072d-6 * K2)
  
  ! WE CALCULATE TOTAL TORQUE EXERTED BY THE DISK ON THE PLANET
  Gamma_0 = (mass / (stellar_mass * p_prop%aspect_ratio))**2 * p_prop%sigma * p_prop%radius**4 * p_prop%omega**2
  
  lindblad_torque = 0.
  
  !------------------------------------------------------------------------------
  corotation_torque = TORQUE_PROFILE_STEEPNESS * 0.1d0 * ((a * planet_mass + b) - p_prop%radius)


  return
end subroutine get_corotation_torque_mass_dep_CZ

end module disk