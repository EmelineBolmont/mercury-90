module disk_properties

!*************************************************************
!** Modules that contains disk parameters that can be used by several other modules
!**
!** Version 1.0 - june 2012
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
  ! values possible to change the properties of the torque. 'real', 'mass_independant', 'mass_dependant', 'manual'. 
  ! If 'manual' is chosen, the code will read the file 'torque_profile.dat' that must exist and the first column must be semi major axis in AU, and the second one is the torque (in units of \Gamma_0 for the moment)
  character(len=80) :: TORQUE_TYPE = 'real' 
  ! Here we define the constant value of the viscosity of the disk
  real(double_precision) :: viscosity = 1.d15 ! viscosity of the disk [cm^2/s]
  logical :: isTurbulence = .False. ! if there is turbulence or not inside the disk
  real(double_precision) :: TURBULENT_FORCING = 1.d-5 ! the turbulent forcing parameter, which controls the amplitude of the stochastic density perturbations.
  ! the value of the turbulent forcing gamma is related to the alpha parameter of the viscosity prescription by : alpha = 120 (gamma / h)^2 where h is the aspect ratio
  
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
  ! Here we define properties common to the profiles
  real(double_precision) :: INNER_BOUNDARY_RADIUS = 1.d0
  real(double_precision) :: OUTER_BOUNDARY_RADIUS = 100.d0
  integer :: NB_SAMPLE_PROFILES = 200 ! number of points for the sample of radius of the temperature profile
  
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
    real(double_precision) :: mass ! the mass of the planet [Msun * K2]
    
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

subroutine print_planet_properties(p_prop)
! subroutine that display in the terminal all the values 
! contained in the instance of planet properties given in parameters
!
! Parameters
! p_prop : an object of type 'PlanetProperties'
  implicit none
  type(PlanetProperties), intent(in) :: p_prop
  
  write (*,*) 'angular_momentum :', p_prop%angular_momentum, ' [Ms.AU^2.day^-1]'
  write (*,*) 'radius :', p_prop%radius, ' [AU]'
  write (*,*) 'velocity :', p_prop%velocity , ' [AU/day]'
  write (*,*) 'omega :', p_prop%omega , ' [day-1]'
  write (*,*) 'mass :', p_prop%mass, ' [Msun * K2]' 
  write (*,*) 'semi_major_axis :', p_prop%semi_major_axis, ' [AU]'
  write (*,*) 'eccentricity :', p_prop%eccentricity 
  write (*,*) 'inclination :', p_prop%inclination, ' [rad]'
  write (*,*) 'sigma :', p_prop%sigma , ' [Msun.AU^-2]'
  write (*,*) 'sigma_index :', p_prop%sigma_index
  write (*,*) 'scaleheight :', p_prop%scaleheight , ' [AU]'
  write (*,*) 'aspect_ratio :', p_prop%aspect_ratio 
  write (*,*) 'chi :', p_prop%chi , ' [AU^2.day^-1]'
  write (*,*) 'nu :', p_prop%nu, ' [AU^2.day^-1]'
  write (*,*) 'temperature :', p_prop%temperature , ' [K]'
  write (*,*) 'temperature_index:', p_prop%temperature_index
end subroutine print_planet_properties

end module disk_properties