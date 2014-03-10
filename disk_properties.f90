!******************************************************************************
! MODULE: disk_properties
!******************************************************************************
!
! DESCRIPTION: 
!> @brief contains disk parameters that can be used by several other 
!! modules, that is, global parameters. Defines also interfaces for subroutine pointers
!
!******************************************************************************

module disk_properties

  use types_numeriques
  use physical_constant

  implicit none
  
  !------------------------------------------------------------------------------
  ! Star parameters :
  real(double_precision) :: R_STAR = 4.6491d-3 ! Stellar radius in AU
  real(double_precision) :: T_STAR = 5700.d0 ! in K
  
  !------------------------------------------------------------------------------
  ! Default values for parameters that are to be read in the parameter file 'disk.in'
  real(double_precision) :: B_OVER_H_LINDBLAD = 0.4d0 ! the smoothing length for the planet's potential, related to the lindblad torque
  real(double_precision) :: B_OVER_H_COROTATION = 0.4d0 ! the smoothing length for the planet's potential, related to the corotation torque
  real(double_precision) :: ADIABATIC_INDEX = 1.4d0 ! the adiabatic index for the gas equation of state
  real(double_precision) :: MEAN_MOLECULAR_WEIGHT = 2.35d0 ! the mean molecular weight in mass of a proton
  
  ! Here we define the power law for surface density sigma(R) = INITIAL_SIGMA_0 * R^(-INITIAL_SIGMA_INDEX)
  logical :: IS_MANUAL_SURFACE_DENSITY = .false. ! If true, the surface density profile is read from the file 'surface_density_profile.dat' instead of being a power law following the 2 next parameters.
  real(double_precision) :: GROUND_SURFACE_DENSITY = 5.d0 ! (g/cm^2) the minimum surface density (for dissipation or for inner edge smoothing)
  real(double_precision) :: INITIAL_SIGMA_0 = 450.d0 ! the surface density at (R=1AU) [g/cm^2]
  real(double_precision) :: INITIAL_SIGMA_INDEX = 0.5d0! the negative slope of the surface density power law (alpha in the paper)
  real(double_precision), parameter :: SIGMA_CGS2NUM = AU**2 / MSUN ! The factor to convert surface density from g/cm^2 to MSUN/AU^2 (the numerical units)
  real(double_precision), parameter :: SIGMA_NUM2CGS = MSUN / AU**2 ! The factor to convert surface density from MSUN/AU^2 (numerical units) to g/cm^2 (CGS)
  
  !------------------------------------------------------------------------------
  ! Irradiation Parameters
  logical :: IS_IRRADIATION = .false. ! if there is irradiation or not for the calculation of the temperature profile
  real(double_precision) :: DISK_ALBEDO = 0.5d0
  
  integer :: DISSIPATION_TYPE = 0 ! integer to tell if there is dissipation of the disk or not. 0 for no dissipation, 1 for viscous dissipation and 2 for exponential decay of the initial profile. 3 for mixed dissipation, both viscously and with photoevaporation, with two timescales
  real(double_precision) :: next_dissipation_step = -1.d0 !< next time at which we will compute the thermal properties of the disk?
  real(double_precision) :: TAU_DISSIPATION = 1.d6 ! the characteristic time for the exponential decay of the surface density (in years) (dissipation_type=2)
  
  
  real(double_precision) :: DISSIPATION_TIME_SWITCH = 2d6 ! (years) the time at which we switch from viscous to photoevaporation exponential decrease (dissipatio_type=3)
  real(double_precision) :: TAU_VISCOUS = 1.d7 ! (years) the characteristic time for the viscous exponential decay
  real(double_precision) :: TAU_PHOTOEVAP = 3.d4 ! (years) the characteristic time for the photoevaporation exponential decay

  ! values possible to change the properties of the torque. 'real', 'mass_independant', 'mass_dependant', 'manual'. 
  ! If 'manual' is chosen, the code will read the file 'torque_profile.dat' that must exist and the first column must be semi major axis in AU, and the second one is the torque (in units of \Gamma_0 for the moment)
  character(len=80) :: TORQUE_TYPE = 'real' 
  character(len=80) :: OPACITY_TYPE = 'bell' 
  character(len=80) :: DAMPING_TYPE = 'cossou' 
  
  ! Here we define the constant value of the viscosity of the disk
  character(len=80) :: VISCOSITY_TYPE = 'constant' ! 'constant' or 'alpha', or 'alpha_dz'
  real(double_precision) :: VISCOSITY = 1.d15 ! viscosity of the disk [cm^2/s] (must not be used directly into the code, use get_viscosity(radius) instead)
  real(double_precision) :: ALPHA = 1.d-3 ! in case of alpha prescription, the alpha of the disk.
  real(double_precision), dimension(3) :: alpha_dz = (/ 1d-2, 1d-4, 1d-2/) ! in case of alpha_dz prescription, the alpha of the disk.
  real(double_precision), dimension(2) :: radius_dz = (/ 1.d0, 10.d0/)! [AU] in case of alpha_dz prescription, the two radius that separate the 3 different alpha regions.
  !------------------------------------------------------------------------------
  logical :: IS_TURBULENCE = .False. ! if there is turbulence or not inside the disk
  real(double_precision) :: TURBULENT_FORCING_MANUAL = 0.d0 ! We can specify manually the value of the turbulent forcing in disk.in
  real(double_precision) :: TURBULENT_FORCING ! the turbulent forcing parameter, which controls the amplitude of the stochastic density perturbations.
  ! the value of the turbulent forcing gamma is related to the alpha parameter of the viscosity prescription by : alpha = 120 (gamma / h)^2 where h is the aspect ratio
  !------------------------------------------------------------------------------
  ! Values for manual torque profiles (mass dependant or independant for instance)
  real(double_precision) :: TORQUE_PROFILE_STEEPNESS = 1.d0 ! increase, in units of Gamma_0 of the torque per 10AU
  real(double_precision) :: SATURATION_TORQUE = 1.d0
  
  ! For mass independant convergence zone
  real(double_precision) :: INDEP_CZ = 8.d0 ! Position of the mass independant convergence zone in AU
  
  ! For mass dependant convergence zone
  ! boundary mass values that have a zero torque zone
  real(double_precision) :: MASS_DEP_M_MIN = 1.d0  ! Minimum mass for the mass dependant convergence zone (in earth mass)
  real(double_precision) :: MASS_DEP_M_MAX = 60.d0 ! Maximum mass for the mass dependant convergence zone (in earth mass)
  
  ! position of the zero torque zone for the minimum and maximum mass
  real(double_precision) :: MASS_DEP_CZ_M_MIN = 4.d0  ! position of the mass dependant convergence zone for the minimum mass (in AU)
  real(double_precision) :: MASS_DEP_CZ_M_MAX = 30.d0 ! position of the mass dependant convergence zone for the maximum mass (in AU)
  
  !------------------------------------------------------------------------------
  ! Here we define properties common to the profiles
  real(double_precision) :: INNER_BOUNDARY_RADIUS = 0.2d0
  real(double_precision) :: OUTER_BOUNDARY_RADIUS = 100.d0
  real(double_precision) :: INNER_SMOOTHING_WIDTH = 0.05d0 ! the width (in unit of the inner boundary radius) of the inner region where the surface density is smoothed to 0
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
    ! the scaleheight and/or aspect ratio is not used in the calculation of the turbulence, where the value 0.05 is used directly into the code
    real(double_precision) :: aspect_ratio ! the aspect_ratio of the gas disk at the location of the planet [no dim]
    real(double_precision) :: chi ! the thermal diffusion coefficient at the location of the planet [AU^2.day^-1]
    real(double_precision) :: nu ! the viscosity of the disk at the location of the planet [AU^2.day^-1]
    real(double_precision) :: temperature ! the temperature of the disk at the location of the planet [K] 
    real(double_precision) :: temperature_index ! the negative temperature index of the disk at the location of the planet [no dim] 
  end type PlanetProperties
  
  procedure(get_torques_interface), pointer :: get_torques
  procedure(function_temperature_interface), pointer :: zero_finding_temperature
  procedure(get_opacity_interface), pointer :: get_opacity
  procedure(get_viscosity_interface), pointer :: get_temp_viscosity ! Only used to retrieve the temperature profile. 
  ! After that, we use a tabulated viscosity profile, fixed once and for all.
  procedure(get_corotation_damping_interface), pointer :: get_corotation_damping
  
  abstract interface 
  subroutine get_torques_interface(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, Gamma_0, ecc_corot)
    import 
    
    implicit none
    real(double_precision), intent(in) :: stellar_mass ! the mass of the central body [Msun * K2]
    ! Properties of the planet
    real(double_precision), intent(in) :: mass ! the mass of the planet [Msun * K2]
    type(PlanetProperties), intent(in) :: p_prop ! various properties of the planet
    
    
    real(double_precision), intent(out) :: corotation_torque
    real(double_precision), intent(out) :: lindblad_torque !  lindblad torque exerted by the disk on the planet [\Gamma_0]
    real(double_precision), intent(out) :: Gamma_0 ! canonical torque value [Ms.AU^2](equation (8) of Paardekooper, Baruteau, 2009)
    real(double_precision), intent(out) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
  end subroutine get_torques_interface
  end interface
  
  abstract interface 
  subroutine function_temperature_interface(temperature, sigma, omega, distance_new, scaleheight_old, distance_old, &
                                            funcv, optical_depth, nu)
    import
    
  ! Output
  real(double_precision), intent(out) :: funcv ! the value of the function
  real(double_precision), intent(out) :: optical_depth ! the optical depth at a given position
  real(double_precision), intent(out) :: nu ! the viscosity in numerical units

  ! Input
  real(double_precision), intent(in) :: temperature ! the temperature at a given position (in K)
  real(double_precision), intent(in) :: sigma ! the surface density at a given position (in MS/AU**2)
  real(double_precision), intent(in) :: distance_new ! current orbital distance [AU]
  real(double_precision), intent(in) :: omega ! the angular velocity of the disk at a given position
  real(double_precision), intent(in) :: scaleheight_old ! aspect ratio of the previous point
  real(double_precision), intent(in) :: distance_old ! orbital distance of the previous point [AU]
  end subroutine function_temperature_interface
  end interface
  
  abstract interface 
  function get_opacity_interface(temperature, num_bulk_density)
  ! subroutine that return the opacity of the disk at the location of the planet given various parameters
    import
    
    implicit none
    
    ! Inputs 
    real(double_precision), intent(in) :: temperature & ! temperature of the disk [K]
                                          , num_bulk_density ! bulk density of the gas disk [MSUN/AU^3] (in numerical units)
  
    ! Outputs
    real(double_precision) :: get_opacity_interface
  
  end function get_opacity_interface
  end interface
  
  abstract interface 
  function get_viscosity_interface(omega, scaleheight, radius)
  ! subroutine that return the opacity of the disk at the location of the planet given various parameters
    import
    
    implicit none
    
    ! Inputs 
    real(double_precision), intent(in) :: omega ! The angular speed in DAY-1
    real(double_precision), intent(in) :: scaleheight ! the scaleheight of the disk in AU
    real(double_precision), intent(in) :: radius ! the distance from the host star in AU
  
    ! Outputs
    real(double_precision) :: get_viscosity_interface ! in [AU^2.day-1]
  
  end function get_viscosity_interface
  end interface
  
  abstract interface 
  function get_corotation_damping_interface(e, x_s, h)
  ! subroutine that return the corotation damping for a planet given its eccentricity
    import
    
    implicit none
    
    ! Inputs 
    real(double_precision), intent(in) :: e ! eccentricity
    real(double_precision), intent(in) :: x_s  ! half-width of the horseshoe region
    real(double_precision), intent(in) :: h ! aspect ratio of the disk
  
    ! Outputs
    real(double_precision) :: get_corotation_damping_interface
  
  end function get_corotation_damping_interface
  end interface
contains

subroutine print_planet_properties(p_prop, output)
! subroutine that display in the terminal all the values 
! contained in the instance of planet properties given in parameters
!
! Parameters
! p_prop : an object of type 'PlanetProperties'
! output : the unit where to write the informations. By default, if nothing 
!        specified, the information are displayed on the screen
  implicit none
  type(PlanetProperties), intent(in) :: p_prop
  integer, optional :: output
  
! Local
integer :: unit


  if (present(output)) then
    unit = output
  else
    unit = 6
  end if

  write(unit,'(a)')            '________________________________________________'
  write(unit,'(a)')            '|####################################'
  write(unit,'(a)')            '|#     Properties of the planet     #'
  write(unit,'(a)')            '|####################################'
  write(unit,'(a,2(es10.2e2,a))') '| Mass : ', p_prop%mass, ' [Msun * K2] (', p_prop%mass/K2 / EARTH_MASS, ' earth mass)'
  write(unit,'(a,es14.4e2,a)')     '| Semi-major axis : ', p_prop%semi_major_axis, ' [AU]'
  write(unit,'(a,f10.7)')       '| Eccentricity : ', p_prop%eccentricity 
  write(unit,'(a,f7.3,a)')     '| Inclination : ', p_prop%inclination*180./PI, ' [degrees]'
  write(unit,'(a,es12.2e2,a)')     '| Radius : ', p_prop%radius, ' [AU]'
  write(unit,'(a,es10.2e2,a)') '| Velocity : ', p_prop%velocity , ' [AU/day]'
  write(unit,'(a,es10.2e2,a)') '| Omega : ', p_prop%omega , ' [day-1]'
  write(unit,'(a,es10.2e2,a)') '| Angular_momentum : ', p_prop%angular_momentum, ' [Ms.AU^2.day^-1]'
  write(unit,'(a)')            '|####################################'
  write(unit,'(a)')            '|#  Properties of the disk at the   #'
  write(unit,'(a)')            '|#     location of the planet       #'
  write(unit,'(a)')            '|####################################'
  write(unit,'(a,es10.2e2,a)') '| Sigma : ', p_prop%sigma , ' [Msun.AU^-2]'
  write(unit,'(a,f7.2)')       '| Sigma_index : ', p_prop%sigma_index
  write(unit,'(a,es10.2e2,a)')     '| Scaleheight : ', p_prop%scaleheight , ' [AU]'
  write(unit,'(a,f6.4)')       '| Aspect_ratio : ', p_prop%aspect_ratio 
  write(unit,'(a,es10.2e2,a)')     '| Chi : ', p_prop%chi , ' [AU^2.day^-1]'
  write(unit,'(a,es10.2e2,a)') '| Nu : ', p_prop%nu, ' [AU^2.day^-1]'
  write(unit,'(a,es12.2e2,a)')     '| Temperature : ', p_prop%temperature , ' [K]'
  write(unit,'(a,f8.3)')       '| Temperature_index : ', p_prop%temperature_index
  write(unit,'(a)')            '------------------------------------------------'

end subroutine print_planet_properties

end module disk_properties
