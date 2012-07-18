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

subroutine print_turbulencemode_properties(mode, output)
! subroutine that display in the terminal all the values 
! contained in the instance of turbulenceMode given in parameters
!
! Parameters
! mode : an object of type 'TurbulenceMode'
! output : the unit where to write the informations. By default, if nothing 
!        specified, the information are displayed on the screen
  implicit none
  type(TurbulenceMode), intent(in) :: mode
  integer, optional :: output
  
  ! Local
  integer :: unit
  !------------------------------------------------------------------------------
  
  if (present(output)) then
    unit = output
  else
    unit = 6
  end if
  
  write(unit,'(a)')               '________________________________________________'
  write(unit,'(a,es10.3e2,a)')    '| Time of creation of the mode : ', mode%t_init, ' days'
  write(unit,'(a,f9.2,a)')        '| Lifetime of the mode : ', mode%lifetime, ' days'
  write(unit,'(a,I2)')            '| Wavenumber : ', mode%wavenumber
  write(unit,'(a,f4.1,a)')        '| Radial extent : ', mode%radial_extent, ' AU'
  write(unit,'(a,f4.1,a,f5.1,a)') '| Position in the disk (r,theta) : (', mode%r , ' AU, ', mode%phi * 180./PI, ' degrees)'
  write(unit,'(a,f5.2)')          '| Chi : ', mode%chi 
  write(unit,'(a)')               '------------------------------------------------'

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
  real(double_precision) :: random_number_01, random_number_02, random_number_03, random_number_log
  type(PlanetProperties) :: d_prop ! to get some properties of the disk at a given location.calculate
  real(double_precision) :: aspect_ratio = 0.05d0 ! the aspect ratio of the disk. This has nothing to do with the value computed elsewhere. 
  
  !------------------------------------------------------------------------------
  call random_number(random_number_01) ! By default between 0 and 1
  call random_number(random_number_02)
  call random_number(random_number_03)
  
  random_number_log = dble(wavenumber_min) * exp(log(dble(wavenumber_max) / dble(wavenumber_min)) * random_number_01) ! pour le mode m, on passe d'une distrib uniforme a une  distrib normale, modemax vaut 96 et modemin vaut 1
  mode%wavenumber = dble(int(random_number_log))! on force le fait qeu ce soit un entier

  mode%r = INNER_BOUNDARY_RADIUS + (OUTER_BOUNDARY_RADIUS - INNER_BOUNDARY_RADIUS) * random_number_02 ! c'est le rk du papier, c'est la position radiale de l'origine de la turbulence

  mode%phi = TWOPI * random_number_03 ! pareil pour phik
  call normal(mode%chi) ! ca genere xi_k, distribution normale de moyenne nulle et d'ecart type 1
  
  mode%radial_extent = 0.25d0 * PI * mode%r / mode%wavenumber ! c'est l'etendu du mode de la turbulence
  
  mode%t_init = time ! origine de la generation du mode 

  mode%lifetime = TWOPI * mode%r**(1.5d0) / (mode%wavenumber * aspect_ratio) ! lifetime of the mode
  
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
  ! We initialize only if turbulence_mode is not allocated yet (which would mean that the initialisation has already been done).
  if (.not.allocated(turbulence_mode)) then
	allocate(turbulence_mode(nb_modes))
	
	! We initialize the random seed
	call init_random_seed()
	
	do i=1, nb_modes
	  call init_mode(time, turbulence_mode(i))
	end do
  end if
  
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
  real(double_precision), parameter :: acc_prefactor = (64.d0 / (PI * PI)) ! the prefactor used to calculate the acceleration
  ! The surface density of the disk is given in MSUN, with no K2, so the solar mass is given in MSUN too, without K2
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
	  
!~ 	  write(*,'(a,es12.3e2,a,i2)') 'time = ', time/365.25d0, ' years. Initialisation of mode k=',k 
	  
!~ 	  open(10, file='turbulence_modes.out', access='append')
!~ 	  write(10,*) 'time=', time, 'creation of mode ',k
!~ 	  call print_turbulencemode_properties(turbulence_mode(k), unit=10)
!~ 	  close(10)
	end if

	! if the mode is too faint, we neglect it, instead of calculating a very small number
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