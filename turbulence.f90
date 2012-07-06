module turbulence

!*************************************************************
!** Modules that contains routines useful for turbulences
!**
!** Version 1.0 - mai 2012
!*************************************************************

! WARNING : You must call once the subroutine "init_turbulence", in the global code or whatever

  use types_numeriques
  use physical_constant
  use disk_properties

  implicit none
    
  integer, parameter :: nb_modes = 50 ! The total number of mode existing at the same time and that represent the turbulence in the disk at any time.
  integer :: wavenumber_min = 1
  integer :: wavenumber_max = 96
  integer :: wavenumber_cutoff = 6 ! If the wavenumber if greater than this number, then we set the gravitational potential of this mode to 0, to gain computational time.
  real(double_precision) :: lifetime_prefactor = 1.0 ! a security factor to limit the lifetime of a mode and fit more accurately to hydro simulations
  
  ! We define a new type for the properties of the planet
  type TurbulenceMode
    ! Properties of the mode
    real(double_precision) :: r ! the radial position of the mode [AU]
    real(double_precision) :: phi ! the azimuthal position of the mode [radians]
    integer :: wavenumber ! The mode number associated to the mode
    real(double_precision) :: t_init ! The time at which the mode was initialised [days]
    real(double_precision) :: lifetime ! the duration of the mode since its initialisation [days]
    real(double_precision) :: chi ! I don't know exactly what it is. Seems like a prefactor [dimensionless]
    real(double_precision) :: radial_extent ! The radial extent of the mode [AU]

  end type TurbulenceMode
  
  type(TurbulenceMode), dimension(:), allocatable :: turbulence_mode ! an array containing the values for the 50 modes of turbulence

  contains

subroutine print_turbulencemode_properties(mode)
! subroutine that display in the terminal all the values 
! contained in the instance of turbulenceMode given in parameters
!
! Parameters
! mode : an object of type 'TurbulenceMode'
  implicit none
  type(TurbulenceMode), intent(in) :: mode
  
  write(*,'(a)')               '________________________________________________'
  write(*,'(a,es10.3e2,a)')    '| Time of creation of the mode : ', mode%t_init, ' days'
  write(*,'(a,f9.2,a)')        '| Lifetime of the mode : ', mode%lifetime, ' days'
  write(*,'(a,I2)')            '| Wavenumber : ', mode%wavenumber
  write(*,'(a,f4.1,a)')        '| Radial extent : ', mode%radial_extent, ' AU'
  write(*,'(a,f4.1,a,f5.1,a)') '| Position in the disk (r,theta) : (', mode%r , ' AU, ', mode%phi * 180./PI, ' degrees)'
  write(*,'(a,f5.2)')          '| Chi : ', mode%chi 
  write(*,'(a)')               '------------------------------------------------'

end subroutine print_turbulencemode_properties

subroutine init_mode(time, mode)
  ! Routine that create a mode for turbulence, and return an object of type "TurbulenceMode"
  
  ! In this routine, the aspect ratio is supposed to be 0.05, and is not related to the 'real' aspect ratio of the disk. 
  ! This is only used for the lifetime of a mode.
  
  ! Parameter
  ! time : the time of creation of the mode
  !
  ! Return : mode
  
  implicit none
  ! Input
  real(double_precision), intent(in) :: time
  
  ! Outputs
  type(TurbulenceMode), intent(out) :: mode
  
  ! Locals
  real(double_precision) :: random_number_01, random_number_02, random_number_03, random_number_normal
  type(PlanetProperties) :: d_prop ! to get some properties of the disk at a given location.calculate
  real(double_precision) :: aspect_ratio = 0.05d0 ! the aspect ratio of the disk. This has nothing to do with the value computed elsewhere. 
  
  !------------------------------------------------------------------------------
  call random_number(random_number_01) ! By default between 0 and 1
  call random_number(random_number_02)
  call random_number(random_number_03)
  
  random_number_normal = dble(wavenumber_min) * exp(log(dble(wavenumber_max) / dble(wavenumber_min)) * random_number_01) ! pour le mode m, on passe d'une distrib uniforme a une  distrib normale, modemax vaut 96 et modemin vaut 1
  mode%wavenumber = dble(int(random_number_normal))! on force le fait qeu ce soit un entier

  mode%r = INNER_BOUNDARY_RADIUS + (OUTER_BOUNDARY_RADIUS - INNER_BOUNDARY_RADIUS) * random_number_02 ! c'est le rk du papier, c'est la position radiale de l'origine de la turbulence

  mode%phi = TWOPI * random_number_03 ! pareil pour phik
  call normal(mode%chi) ! ca genere xi_k, distribution normale de moyenne nulle et d'ecart type 1
  
  mode%radial_extent = 0.25d0 * PI * mode%r / mode%wavenumber ! c'est l'etendu du mode de la turbulence
  
  mode%t_init = time ! origine de la generation du mode 
  ! The factor lifetime_prefactor is here to have a better agreement with mhd simulations. lifemode value has 0.1
  mode%lifetime = lifetime_prefactor * TWOPI * mode%r**(1.5d0) / (mode%wavenumber * 0.05) ! lifetime of the mode
  
!~   call print_turbulencemode_properties(mode)
  
end subroutine init_mode

subroutine init_turbulence(time)
! routine that initialize all the parameters of the turbulence
  implicit none
  
  ! Inputs
  real(double_precision), intent(in) :: time
  
  ! Locals
  integer :: i
  
  !------------------------------------------------------------------------------
  
  allocate(turbulence_mode(nb_modes))
  
  ! We initialize the random seed
  call init_random_seed()
  
  do i=1, nb_modes
	call init_mode(time, turbulence_mode(i))
  end do
  
end subroutine init_turbulence

subroutine get_turbulence_acceleration(time, p_prop, position, turbulence_acceleration)
! Routine that return the acceleration feeled by the planet described by p_prop and position due to turbulence.
! 
! Global Parameter
! TURBULENT_FORCING : the turbulent forcing parameter, which controls the amplitude of the stochastic density perturbations.
!                     the value of the turbulent forcing gamma is related to the alpha parameter of 
!                     the viscosity prescription by : alpha = 120 (gamma / h)^2 where h is the aspect ratio



  use utilities, only : get_polar_coordinates

  implicit none
  
  ! Inputs
  real(double_precision), intent(in) :: time ! The current time [days]
  real(double_precision), dimension(3), intent(in) :: position ! the cartesian position [AU]
  type(PlanetProperties) :: p_prop ! An object to store various properties of a planet
  
  ! Outputs
  real(double_precision), dimension(3), intent(out) :: turbulence_acceleration ! The gravitational potential induced by the turbulence
  
  ! Locals
  real(double_precision), parameter :: acc_prefactor = (64.d0 / (K2 * PI * PI)) ! the prefactor used to calculate the acceleration
  real(double_precision) :: r, phi ! polar coordinates of the planet
  real(double_precision), dimension(3) :: shifted_position
  real(double_precision) :: single_prefactor ! common prefactor for a given mode, for each position
  real(double_precision) :: planet_prefactor ! prefactor including properties of the planet that must be calculated each time.
  real(double_precision) :: argument ! argument for cos and sin of the single potential (created by one mode)
  real(double_precision) :: lambda_cm, lambda_sm ! single potential created by one mode (and its derivative) used for the expression of the force.
  real(double_precision) :: sum_element_r, sum_element_theta ! result of the sum on all the modes for the r or theta part of the force
  real(double_precision) :: force_radius, force_theta ! force exerted by the turbulence, on r or theta direction
  real(double_precision) :: relative_time ! the current time, with respect to the creation time of a given mode
  integer :: k ! for loops  
  
  !------------------------------------------------------------------------------  
  
  ! initialisation
  ! You must call once the subroutine "init_turbulence", in the global code or whatever
  sum_element_r = 0.d0
  sum_element_theta  = 0.d0
  
  ! We get the polar coordinates of the planet
  call get_polar_coordinates(position(1), position(2), position(3), r, phi)

  ! la on calcule le potentiel turbulent, on fait la smme de tous les
  ! modes en tenant compte de leur evolution temporelle
  do k=1,nb_modes
	relative_time = time - turbulence_mode(k)%t_init ! c'est le temps relatif au temps d 'origine du mode
	
	! If needed, we replace an old mode by a new one
	if (relative_time.ge.turbulence_mode(k)%lifetime) then
	  call init_mode(time, turbulence_mode(k))
	  relative_time = 0.d0
	end if

	if (turbulence_mode(k)%wavenumber.le.wavenumber_cutoff) then
	  single_prefactor = turbulence_mode(k)%chi * sin(PI * relative_time / turbulence_mode(k)%lifetime) * &
						 exp(-((r - turbulence_mode(k)%r) / turbulence_mode(k)%radial_extent)**2)
	  argument = turbulence_mode(k)%wavenumber * phi - turbulence_mode(k)%phi - p_prop%omega * relative_time
	
	  lambda_cm = single_prefactor * cos(argument) ! equation (7) (Ogihara, 2007)
	  lambda_sm = single_prefactor * sin(argument) ! equation (10) (Ogihara, 2007)
	
	  sum_element_r = sum_element_r + (1.d0 + 2.d0 * r * (r - turbulence_mode(k)%r) / turbulence_mode(k)%radial_extent**2) * &
	                  lambda_cm
	  sum_element_theta = sum_element_theta + turbulence_mode(k)%wavenumber * lambda_sm
	  

	endif
  enddo

  ! We apply at the end the prefactor of the gravitational potential
  planet_prefactor = acc_prefactor * TURBULENT_FORCING * p_prop%radius**3 * p_prop%omega**2 * p_prop%sigma
  
  force_radius = planet_prefactor * sum_element_r
  force_theta  = planet_prefactor * sum_element_theta
  
  ! To get the turbulence acceleration, we use (28), (29) and (30) of (ogihara, 2007). Some constant calculation is putted in a 
  ! prefactor, some other calculation that depend upon the planet are calculed in another prefactor
  
  turbulence_acceleration(1) = cos(phi) * force_radius - sin(phi) * force_theta
  turbulence_acceleration(2) = sin(phi) * force_radius + cos(phi) * force_theta
  turbulence_acceleration(3) = 0.d0

end subroutine get_turbulence_acceleration

subroutine get_turbulence_potential(time, p_prop, position, full_gravitational_potential)
  use utilities, only : get_polar_coordinates

  implicit none
  
  ! Inputs
  real(double_precision), intent(in) :: time ! The current time [days]
  real(double_precision), dimension(3), intent(in) :: position ! the cartesian position [AU]
  type(PlanetProperties) :: p_prop ! An object to store various properties of a planet
  
  ! Outputs
  real(double_precision), intent(out) :: full_gravitational_potential ! The gravitational potential induced by the turbulence
  
  ! Locals
  integer k
  real(double_precision) :: single_potential ! the gravitational potential of a single turbulence mode
  real(double_precision) :: relative_time ! the current time, with respect to the creation time of a given mode
  real(double_precision) :: theta_planet ! the angle of the planet in polar coordinates
  real(double_precision) :: radius_planet ! the radius of the planet [AU] (currently not needed here since we have the value in the p_prop object)

  !------------------------------------------------------------------------------
  
  call get_polar_coordinates(position(1), position(2), position(3), radius_planet, theta_planet) 
  ! This radius must be used instead of p_prop%radius because p_prop is the same during the calculation of the derivative
  
  ! initialisation
  ! You must call once the subroutine "init_turbulence", in the global code or whatever

full_gravitational_potential = 0.0d0

! here we calculate the turbulent potential by calculating the potential of each mode, taking into account their dissipation through time.
do k=1,nb_modes
  relative_time = time - turbulence_mode(k)%t_init ! time expressed with respect to the creation time of the mode.
  
  ! If needed, we replace an old mode by a new one
  if (relative_time.ge.turbulence_mode(k)%lifetime) then
	call init_mode(time, turbulence_mode(k))
	relative_time = 0.d0
  end if

  if (turbulence_mode(k)%wavenumber.gt.wavenumber_cutoff) then
	single_potential = 0.0d0
  else

	single_potential =  turbulence_mode(k)%chi * & 
				   exp(-(radius_planet - turbulence_mode(k)%r)**2 / turbulence_mode(k)%radial_extent) * &
				   cos(turbulence_mode(k)%wavenumber * theta_planet - &
				   turbulence_mode(k)%phi - p_prop%omega * relative_time) * &
				   sin(pi * relative_time / turbulence_mode(k)%lifetime) ! expression 14 of pierens et al.

	full_gravitational_potential = full_gravitational_potential + single_potential
  endif
enddo

! We apply at the end the prefactor of the gravitational potential
full_gravitational_potential = turbulent_forcing * p_prop%radius**2 * p_prop%omega**2 * full_gravitational_potential

end subroutine get_turbulence_potential

SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE

SUBROUTINE normal(x)
  ! generate a random number through a gaussian distribution of mean=0 and standard deviation of 1
  implicit none
  
  ! Output
  real(double_precision), intent(out) :: x
  
  ! Locals
  real(double_precision) :: v1,v2,v11,v22
  real(double_precision) :: r

  r = 1.5d0
  do while ((r.gt.1.d0).OR.(r.eq.0.d0))
	call random_number(v1)
	call random_number(v2)
	v11 = 2.d0 * v1 - 1.d0
	v22 = 2.d0 * v2 - 1.d0
	r = v11**2 + v22**2
  enddo
  
  x = v11 * dsqrt(-2.d0 * dlog(r) / r)
END SUBROUTINE normal

end module turbulence