! Program that calculate the total torque on one planet, that is the resultant torque from his torque, and the torque of all the 
! planet trapped in a resonant chain with it. 
program resultant_torque
  use types_numeriques
  use disk
    
  implicit none
  
  real(double_precision), dimension(:), allocatable :: semi_major_axis, eccentricity, inclination, mass
  real(double_precision), dimension(:), allocatable :: single_torque_adim, ref_torque, single_torque, merged_torque
  real(double_precision) :: torque_sum, total_torque
  integer :: nb_planets, i_start, i_stop
  logical :: isChangeOfSign
  
  ! For loops
  integer :: planet, i
  !-------------------------------------------------------------------
  
  call initialization()

  call read_element_out(semi_major_axis=semi_major_axis, eccentricity=eccentricity, inclination=inclination, mass=mass)
  nb_planets = size(mass)
  
  call calculate_torques(mass=mass, semi_major_axis=semi_major_axis, single_torque_adim=single_torque_adim, ref_torque=ref_torque)
    
  ! We calculate the single torque in numerical units for each planet (in AU, MSUN and DAYS)
  allocate(single_torque(nb_planets))
  single_torque(1:nb_planets) = single_torque_adim(1:nb_planets) * ref_torque(1:nb_planets)
  
  allocate(merged_torque(nb_planets))
  merged_torque(1:nb_planets) = 0.d0
  
  i = 1
  do while (i.le.nb_planets)
    isChangeOfSign = .False.
    total_torque = 0.
    i_start = i
    do while (.not.(isChangeOfSign.and.(single_torque(i_start)*single_torque(i).gt.0.)))
      total_torque = total_torque + single_torque(i)
      
      if (single_torque(i_start)*single_torque(i).lt.0.) then
        isChangeOfSign = .True.
      end if
      
      i = i + 1
      if (i.gt.nb_planets) then
        exit
      end if
    end do
  
  i_stop = i - 1
  
  merged_torque(i_start:i_stop) = total_torque
  end do
    
  ! We write the file 'torque.out' that store informations of the torques of the planets
  open(10, file='torque.out', status='replace')
  write(10,'(a)') '   a     mass  Gamma_p    Gamma_0  Gamma_res'
  do planet=1,nb_planets
  write(10, '(f6.2, f6.1, 3es11.2e2)') semi_major_axis(planet), mass(planet), single_torque(planet), &
                                       ref_torque(planet), merged_torque(planet)
  end do
  close(10)
  
  
  contains
  
  subroutine initialization()
  ! subroutine to initialize some values
  
    implicit none
    
    real(double_precision) :: stellar_mass

    stellar_mass = 1.d0 * K2
    
    ! We force the value to be interesting for our tests
    TORQUE_TYPE = 'mass_dependant' ! 'real', 'mass_independant', 'mass_dependant'
    
    call init_globals(stellar_mass)
    ! Note that the initial density profile and temperature profile are calculated inside the 'init_globals' routine.
  
  end subroutine initialization
  
  subroutine read_element_out(semi_major_axis, eccentricity, inclination, mass)
  ! subroutine that get properties of the planets from the file element.out
  
  implicit none
  
  real(double_precision), dimension(:), allocatable, intent(out) :: semi_major_axis, eccentricity, inclination, mass

  
  character(len=80) :: line
  character(len=8) :: name
  integer :: error ! to store the state of a read instruction
  real(double_precision) :: value1, value2, value3, value4, value5, value6
  integer :: nb_planets
  
  logical :: isDefined
  
  ! For loops
  integer :: i
  !------------------------------------------------------------------------------
  
  inquire(file='element.out', exist=isDefined)
  if (isDefined) then
  
    open(10, file='element.out', status='old')
    
    ! We read the header. There the time in there, but we do not need it for the moment
    do i=1,5
      read(10, '(a80)', iostat=error), line
    end do
    
    nb_planets = 0
    do 
      read(10, *, iostat=error), line
      if (error.eq.0) then
        nb_planets = nb_planets + 1
      else
        exit
      end if
    end do
    
    ! One we have the number of planets, we rewind the file and read again
    rewind(unit=10, iostat=error)
    
    allocate(semi_major_axis(nb_planets))
    allocate(eccentricity(nb_planets))
    allocate(inclination(nb_planets))
    allocate(mass(nb_planets))
    
    ! We skip the header
    do i=1,5
      read(10, '(a80)', iostat=error), line
    end do
    
    i=1
    do
      read(10, *, iostat=error), name, value1, value2, value3, value4, value5, value6
      if (error /= 0) exit
      semi_major_axis(i) = value1
      eccentricity(i) = value2
      inclination(i) = value3
      mass(i) = value4
      
      i = i + 1
    
    end do
  else
    write(*,*) 'Error: the file "element.out" does not exist'
  end if
  
  end subroutine read_element_out

  subroutine calculate_torques(mass, semi_major_axis, single_torque_adim, ref_torque)
    ! subroutine that return the torques in two lists, each element for one planet. 
    ! The first list for torque/gamma_0 and the second for Gamma_0
    ! 
    use physical_constant, only : EARTH_MASS, K2
    
    implicit none
    real(double_precision), dimension(:), intent(in) :: mass, semi_major_axis
    real(double_precision), dimension(:), allocatable, intent(out) :: single_torque_adim, ref_torque
    
    real(double_precision) :: corotation_torque, lindblad_torque ! in units of Gamma_0
    real(double_precision) :: stellar_mass, planet_mass ! mass in (Msun * K2)
    integer :: nb_planets
    
    type(PlanetProperties) :: p_prop ! various properties of a planet
    real(double_precision) :: position(3), velocity(3) ! respectively in AU and AU/day
    
    ! For loops
    integer :: planet
    !--------------------------------------------------------------
    stellar_mass = 1.d0 * K2
    position(1:3) = 0.d0
    velocity(1:3) = 0.d0
    
    nb_planets = size(mass)
    allocate(single_torque_adim(nb_planets))
    allocate(ref_torque(nb_planets))
    
    do planet=1,nb_planets
      planet_mass = mass(planet) * EARTH_MASS * K2
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = semi_major_axis(planet)
      
      ! We generate cartesian coordinate for the given mass and semi major axis
      velocity(2) = sqrt(K2 * (stellar_mass + planet_mass) / position(1))
      
      ! we store in global parameters various properties of the planet
      call get_planet_properties(stellar_mass=stellar_mass, mass=planet_mass, position=position(1:3), velocity=velocity(1:3),& ! Input
       p_prop=p_prop) ! Output
       
!~       call print_planet_properties(p_prop)
      
      call get_corotation_torque_mass_dep_CZ(stellar_mass, planet_mass, p_prop, & ! input
          corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=ref_torque(planet)) ! Output
      single_torque_adim(planet) = (corotation_torque + lindblad_torque)
      
    end do
    
    
  end subroutine calculate_torques

end program resultant_torque
