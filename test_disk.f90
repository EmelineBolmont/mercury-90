! Program that test different functions implemented in the module user_module.
program test_disk
  use types_numeriques
  use disk
  use disk_properties
  use iso_fortran_env, only : error_unit
    
  implicit none
  
  call unitary_tests()
  

  
  contains

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
    
    use opacity_hure, only : test_opacity_interpolation
    
    implicit none
    
    real(double_precision) :: stellar_mass ! [Msun * K2]
    
    logical :: isDefined

    stellar_mass = 1.d0 * K2 ! [Msun * K2]
    
    
    
    inquire(file='unitary_tests', exist=isDefined)
    
    ! We create the folder 'dissipation' if he doesn't exists.
    if (.not.isDefined) then
      call system("mkdir unitary_tests")
    end if
    
    inquire(file='disk.out', exist=isDefined)
    
    ! This file is used to say if we continue an integration or not. For test_disk, 
    ! we DO NOT continue an integration, so to avoid it, we delete the file if it exist.
    if (isDefined) then
      call system("rm disk.out")
    end if
    
    write(*,*) 'Initialisation'
    call init_globals(stellar_mass=stellar_mass, time=0.d0)
    ! Note that the initial density profile and temperature profile are calculated inside the 'init_globals' routine.

    
    ! We want to show the torque profile. It is important to check which value has been declared in 'TORQUE_TYPE'
    call study_torques(stellar_mass)
    
    ! we store in a .dat file the temperature profile
    call store_temperature_profile(filename='temperature_profile.dat')
    call store_density_profile(filename='density_profile.dat')
    call store_scaleheight_profile()
    
    ! Unitary tests
    call test_functions_FGK()
    call test_function_zero_temperature(stellar_mass=stellar_mass)
    call test_temperature_interpolation()
    call test_manual_torque_interpolation()
    call test_density_interpolation()
    call test_retrieval_of_orbital_elements(stellar_mass=stellar_mass)
!     call test_turbulence_torque(stellar_mass=stellar_mass)
    call test_turbulence_mode()
    call test_opacity_interpolation() ! subroutine inside the module opacity_hure.f90

    
    ! Physical values and plots
    call study_opacity_profile()
    call study_viscosity(stellar_mass=stellar_mass)
    call study_torques_fixed_a(stellar_mass=stellar_mass)
    call study_torques_fixed_m(stellar_mass=stellar_mass)
    call study_ecc_corot(stellar_mass=stellar_mass)
    call study_eccentricity_effect_on_corotation(stellar_mass=stellar_mass)
    call study_temperature_profile()
    call study_optical_depth_profile()
    call study_thermal_diffusivity_profile()
    call study_scaleheight_profile()
!~     call study_dissipation_at_one_location()
    
    ! Test dissipation
    ! EVERYTHING ABOVE MUST BE COMMENTED BEFORE DECOMMENTING 'ONE' AND ONE ALONE OF THESES ONES
!~     call test_viscous_dissipation(stellar_mass)
!~     call test_disk_dissipation(stellar_mass)
!~     call study_influence_of_dissipation_on_torque(stellar_mass)

    
  end subroutine unitary_tests

! %%% unitary tests of some routines %%%
  subroutine test_functions_FGK
  ! subroutine that test the functions 'get_F', 'get_G' and 'get_K' and 
  
  ! Return:
  !  a data file 'functions_FGK.dat' 
  ! and an associated gnuplot file 'functions_FGK.gnuplot' that display values for get_F, get_G and get_K for a range of p values.
    implicit none
    
    real(double_precision) :: p, f_p, g_p, k_p
    
    real(double_precision), parameter :: p_min = 0.1d0
    real(double_precision), parameter :: p_max = 100.d0
    integer, parameter :: nb_points = 100
    real(double_precision), parameter :: p_step = (p_max/p_min) ** (1/(nb_points-1.d0))
    
    integer :: i ! for loops
    
    write(*,*) 'test of functions F, G, K from (paardekooper, 2010)'
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/functions_FGK.dat')
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
    write(10,*) 'set terminal pdfcairo enhanced'
    write(10,*) 'set output "functions_FGK.pdf"'
    write(10,*) 'set xlabel "p"'
    write(10,*) 'set ylabel "function"'
    write(10,*) 'set logscale x'
    write(10,*) 'set mxtics 10'
    write(10,*) 'set grid xtics ytics mxtics'
    write(10,*) 'set xrange [', p_min, ':', p_max, ']'
    write(10,*) 'plot "functions_FGK.dat" using 1:2 with lines title "F(p)",\'
    write(10,*) "     '' using 1:3 with lines title 'G(p)',\"
    write(10,*) "     '' using 1:4 with lines title 'K(p)'"
    
    close(10)
      
  end subroutine test_functions_FGK

  subroutine test_function_zero_temperature(stellar_mass)
  ! subroutine that test the function 'zero_finding_temperature'
  
  ! Return:
  !  a data file and an associated gnuplot file.
    implicit none
    real(double_precision), intent(in) :: stellar_mass
    
    real(double_precision) :: temperature
    
    real(double_precision), parameter :: T_min = 0.1d0
    real(double_precision), parameter :: T_max = 10000.d0
    integer, parameter :: nb_points = 2000
!~     real(double_precision), parameter :: T_step = (T_max/T_min) ** (1/(nb_points-1.d0))
    real(double_precision), parameter :: T_step = (T_max - T_min) / (nb_points - 1.d0)
    
    real(double_precision), parameter :: mass = 20. * EARTH_MASS * K2
    
    real(double_precision) :: zero_function, tmp, tmp2 ! value that we want to output and a dummy argument 'tmp'
    
    integer :: i ! for loops
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    real(double_precision) :: orbital_position, temperature_old, scaleheight_old, distance_old

  !------------------------------------------------------------------------------
    orbital_position = 100.d0
  
  
    position(1:3) = 0.d0
    velocity(1:3) = 0.d0
    
    ! we define a fantom point
    distance_old = orbital_position * 1.1
    temperature_old = 10.d0 ! we force the temperature to be 10K for this fantom point
    position(1) = distance_old
    velocity(2) = sqrt((stellar_mass + mass) / position(1))
    call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
     p_prop=p_prop) ! Output
    scaleheight_old = get_scaleheight(temperature=temperature_old, angular_speed=p_prop%omega)
    
    write(*,*) 'Test of the zero function used to calculate the temperature at a given radius'
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/function_zero_temperature.dat')
    
    ! We generate cartesian coordinate for the given semi major axis
    position(1) = orbital_position
    
    ! We generate cartesian coordinate for the given mass and semi major axis
    velocity(2) = sqrt((stellar_mass + mass) / position(1))
    
    ! we store in global parameters various properties of the planet
    call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
     p_prop=p_prop) ! Output
    
    
    write(10,*) '# properties of the disk at the location of the planet that influence the value of the temperature'
    write(10,*) '# radial position of the planet (in AU) :', p_prop%radius
    write(10,*) '# viscosity :', p_prop%nu
    write(10,*) '# surface density :', p_prop%sigma
    write(10,*) '# angular velocity :', P_prop%omega
    write(10,*) '# Temperature (K) ; value of the function. The right temperature is when the function is 0'
    
    do i=1, nb_points
!~       temperature = T_min * T_step ** (i-1)
      temperature = (T_min + T_step * (i - 1.d0))
      
      call zero_finding_temperature(temperature=temperature, sigma=p_prop%sigma, omega=p_prop%omega, distance_new=p_prop%radius, & ! Input
                              scaleheight_old=scaleheight_old, distance_old=distance_old,& ! Input
                              funcv=zero_function, optical_depth=tmp, nu=tmp2) ! Output
      
      write(10,*) temperature, zero_function
    end do
    close(10)
    
    open(10, file="unitary_tests/function_zero_temperature.gnuplot")
    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) "set output 'function_zero_temperature.pdf'"
    write(10,*) 'set xlabel "Temperature T"'
    write(10,*) 'set ylabel "zero function"'
    write(10,*) 'set grid'
    write(10,*) 'set xrange [', T_min, ':', T_max, ']'
    write(10,*) "plot 'function_zero_temperature.dat' using 1:2 with lines notitle"
    
    close(10)
    
  end subroutine test_function_zero_temperature

  subroutine test_temperature_interpolation()
  
    implicit none
    
    integer, parameter :: nb_a = 1000
    real(double_precision), parameter :: a_min = 0.d0 ! in AU
    real(double_precision), parameter :: a_max = 100.d0! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision) :: a, temperature, temperature_index, chi, nu
    
    integer :: j ! for loops
    
    write(*,*) 'Test of the temperature interpolation'
    
    open(10, file='unitary_tests/temperature_interpolation.dat')
    do j=1,nb_a
      a = (a_min + a_step * (j - 1.d0))
      ! We generate cartesian coordinate for the given semi major axis
      
      call get_temperature(radius=a, & ! Input
                           temperature=temperature, temperature_index=temperature_index, chi=chi, nu=nu) ! Output
      
      
      write(10,*) a, temperature, temperature_index, chi, nu
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="unitary_tests/temperature_interpolation.gnuplot")
    

    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) '!rm "temperature_interpolation.pdf"'
    write(10,*) "set output 'temperature_interpolation.pdf'"
    
    write(10,*) 'set xlabel "semi major axis a (in AU)"'
    write(10,*) 'set ylabel "temperature [K]"'
      
    write(10,*) 'set grid'


    write(10,*) 'plot "temperature_interpolation.dat" using 1:2 with lines title "Interpolation",\'
    write(10,*) '     "../temperature_profile.dat" using 1:2 with lines title "Profile"'
    
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
    
    open(10, file='unitary_tests/density_interpolation.dat')
    do j=1,nb_a
      a = (a_min + a_step * (j - 1.d0))
      ! We generate cartesian coordinate for the given semi major axis
      
      call get_surface_density(radius=a, sigma=sigma, sigma_index=sigma_index)
      
      
      write(10,*) a, sigma, sigma_index
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="unitary_tests/density_interpolation.gnuplot")
    

    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) '!rm "density_interpolation.pdf"'
    write(10,*) "set output 'density_interpolation.pdf'"
    
    write(10,*) 'set xlabel "semi major axis a (in AU)"'

    
    write(10,*) 'set ylabel "density [g/cm^2]"'
      
    write(10,*) 'set grid'


    write(10,*) 'plot "density_interpolation.dat" using 1:2 with lines title "Interpolation",\'
    write(10,*) '     "../density_profile.dat" using 1:2 with lines title "Profile"'
    
    
    close(10)
  
  end subroutine test_density_interpolation

  subroutine test_manual_torque_interpolation()
  
    implicit none
    
    integer :: i,j ! for loops
    
    logical :: isDefined
    character(len=80) :: filename
    
    write(*,*) 'Test of the manual torque profile interpolation'
    
    filename = 'torque_profile.dat'
    inquire(file=filename, exist=isDefined)
    
    ! If no manual torque profile exists, we create one
    if (.not.isDefined) then
      open(10, file=filename)
      do i=1, 10
        write(10, *) i * 1.d0, 2.d0 - i * 0.45d0
      end do
      do i=1, 10
        write(10, *) 10 + i * 1.d0, -2.d0 + i * 0.45d0
      end do
      close(10)
    end if
    
    call read_torque_profile()
    
    open(10, file='unitary_tests/torque_interpolation.dat')
    do j=1,NB_SAMPLE_PROFILES
      
      
      write(10,*) distance_sample(j), torque_profile(j)
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="unitary_tests/torque_interpolation.gnuplot")
    

    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) '!rm "torque_interpolation.pdf"'
    write(10,*) "set output 'torque_interpolation.pdf'"
    
    write(10,*) 'set xlabel "semi major axis a (in AU)"'

    
    write(10,*) 'set ylabel "torque [?]"'
      
    write(10,*) 'set grid'
    write(10,*) 'set yrange [',minval(torque_profile)-1.d0, ':', maxval(torque_profile)+1.d0,']'

    write(10,*) 'plot "torque_interpolation.dat" using 1:2 with lines linetype -1 title "Interpolated profile",\'
    write(10,*) '     "../torque_profile.dat" using 1:2 with points linetype 1 pointtype 2 linewidth 3 title "Read Profile"'
    
    
    close(10)
  
  end subroutine test_manual_torque_interpolation

  subroutine test_retrieval_of_orbital_elements(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'retrieval_aeI.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
  
    use contour
    
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass ! in [msun * K2]
    
    real(double_precision), parameter :: mass = 10.d0 * EARTH_MASS * K2
    
    integer, parameter :: nb_points = 100
    real(double_precision), parameter :: a_min = 0.01d0
    real(double_precision), parameter :: a_max = 50.d0
    
    ! step for log sampling
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_points-1.d0)
    real(double_precision), dimension(nb_points) :: a
    
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    
    write(*,*) 'Test of the orbital retrieval'
    
    position(1:3) = 0.d0
    velocity(1:3) = 0.d0
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/retrieval_aeI.dat')
    
    write(10,*) '# initial a (AU) ; final a (AU) ; initial e ; final e ; initial I (degrees) ; final I (degrees)'
    
    
    do i=1, nb_points ! loop on the position
      a(i) = a_min + a_step * (i-1)
      
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a(i)
      
      ! We generate cartesian coordinate for the given mass and semi major axis
      velocity(2) = sqrt((stellar_mass + mass) / position(1))
      
      ! we store in global parameters various properties of the planet
      call get_planet_properties(stellar_mass=stellar_mass, & ! Input
       mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
       p_prop=p_prop) ! Output
       
      
              
      write(10,*) a(i), p_prop%semi_major_axis, 0.0, p_prop%eccentricity, 0.0, p_prop%inclination
      
    end do
    close(10)
    close(11)
    close(12)
    
    
    open(10, file="unitary_tests/retrieval_a.gnuplot")
    open(11, file="unitary_tests/retrieval_e.gnuplot")
    open(12, file="unitary_tests/retrieval_I.gnuplot")
    
    do j=10,12
      write(j,*) 'set terminal pdfcairo enhanced'
    end do
    
    write(10,*) 'set output "retrieval_a.pdf"'
    write(11,*) 'set output "retrieval_e.pdf"'
    write(12,*) 'set output "retrieval_I.pdf"'
    
    write(10,*) 'set xlabel "Initial a"'
    write(11,*) 'set xlabel "Initial a (e=0)"'
    write(12,*) 'set xlabel "Initial a (I=0)"'
    
    write(10,*) 'set ylabel "Retrieved a"'
    write(11,*) 'set ylabel "Retrieved e"'
    write(12,*) 'set ylabel "Retrieved I"'
    
    do j=10,12
      write(j,*) 'set mxtics 10'
      write(j,*) 'set grid xtics ytics mxtics'
    end do
    
    write(10,*) 'plot "retrieval_aeI.dat" using 1:2 with lines notitle'
    write(11,*) 'plot "retrieval_aeI.dat" using 1:3 with lines notitle'
    write(12,*) 'plot "retrieval_aeI.dat" using 1:5 with lines notitle'
    
    do j=10, 12
      close(j)
    end do
    
  end subroutine test_retrieval_of_orbital_elements
  
  subroutine test_disk_dissipation(stellar_mass)
  ! function to test the viscous dissipation with a dirac function. 
  
  ! Global parameters
  ! dissipation_timestep : the timestep between two computation of the disk [in days]
  ! INNER_BOUNDARY_RADIUS : the inner radius of the various profiles (all based on the radius profile)
  ! OUTER_BOUNDARY_RADIUS : the outer radius of the various profiles (all based on the radius profile)
  ! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
  ! surface_density_profile : values of the density in MSUN/AU^2 for each value of the 'a' sample
  ! distance_sample : values of 'a' in AU
  ! viscosity : the viscosity of the disk in [cm^2/s] /!\ This parameter is modified here (which must not be done otherwise)
  
    ! Return: (in the folder 'unitary_tests/dissipation' of the current working directory where the tests are launched)
  !  data files with the following paterns :
  !    surface_density*.dat : each file correspond to a data file of the surface density at one timestep (the number correspond to the index)
  ! IN ADDITION : There are 1 gnuplot scripts to generate the output files for the surface density.
  
  use bessel, only : bessik
  
  implicit none
    
    real(double_precision), intent(in) :: stellar_mass
    
    ! time sample
    real(double_precision), parameter :: t_min = 0.d0 ! time in years
    real(double_precision), parameter :: t_max = 1.d7 ! time in years
    real(double_precision), dimension(:), allocatable :: time, time_temp ! time in days
    integer :: time_size ! the size of the array 'time'. 
    
    real(double_precision) :: density_min, density_max
    character(len=80) :: filename_density, filename_density_ref
    character(len=80) :: output_density, output_time, time_format, purcent_format
    integer :: time_length ! the length of the displayed time, usefull for a nice display
    real(double_precision) :: next_dissipation_step
    
    integer :: k ! for loops
    integer :: nb_time ! The total number of 't' values. 
    integer :: error ! to retrieve error, especially during allocations
    logical :: isDefined
    !------------------------------------------------------------------------------
    write(*,*) 'Test dissipation of the disk'
    
    write(*,*) '  Force Initialisation again'
    FIRST_CALL = .true.
    call init_globals(stellar_mass=stellar_mass, time=0.d0)
    
    inquire(file='unitary_tests/dissipation', exist=isDefined)
    
    ! We create the folder 'dissipation' if he doesn't exists.
    if (.not.isDefined) then
      call system("mkdir unitary_tests/dissipation")
    end if
    
    call system("rm unitary_tests/dissipation/*")
    
    ! We want the extremum of the surface density during the dissipation of the disk in order to have nice plots
    density_min = minval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
    density_max = maxval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
    
    ! We want to know the max size of the time display in order to have a nice display, with filled spaces in the final plots
    write(output_time, '(i0)') int(t_max, kind=8)
    time_length = len(trim(output_time))
    write(time_format, *) '(f',time_length,'.0)'
    write(purcent_format, *) '(i',time_length,'"/",i',time_length,'," years")'
    
    !------------------------------------------------------------------------------
    k = 1
    time_size = 512 ! the size of the array. 
    allocate(time(time_size), stat=error)
    time(1) = t_min * 365.25d0 ! days
    
    
    do while (time(k).lt.(t_max*365.d0))
      ! We open the file where we want to write the outputs
      write(filename_density, '(a,i0.5,a)') 'unitary_tests/dissipation/surface_density',k,'.dat'    
      
      ! we store in a .dat file the surface density profile
      call store_density_profile(filename=filename_density)
      
      if (k.eq.1) then
        ! We want the extremum of the surface density during the dissipation of the disk in order to have nice plots
        density_min = minval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
        density_max = maxval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
      end if
      
      ! We calculate the temperature profile for the current time (because the surface density change in function of time)
      
      
      ! we get the new dissipated surface density profile. For that goal, we dissipate as many times as needed to reach the required time for the next frame.
      call dissipate_disk(time(k), next_dissipation_step)
      
      
      ! we expand the 'time' array if the limit is reached
      if (k.eq.time_size) then
        ! If the limit of the 'time' array is reach, we copy the values in a temporary array, allocate with a double size, et paste the 
        ! old values in the new bigger array
        allocate(time_temp(time_size), stat=error)
        time_temp(1:time_size) = time(1:time_size)
        deallocate(time, stat=error)
        time_size = time_size * 2
        allocate(time(time_size), stat=error)
        time(1:time_size/2) = time_temp(1:time_size/2)
        deallocate(time_temp, stat=error)
      end if
      
      write(*,purcent_format) int(time(k)/365.25d0, kind=8), int(t_max, kind=8) ! We display on the screen how far we are from the end of the integration.
      
      
      time(k+1) = next_dissipation_step ! days
      
      
      k = k + 1 ! We increment the integer that point the time in the array (since it's a 'while' and not a 'do' loop)
    end do

    nb_time = k - 1 ! since for the last step with incremented 'k' by one step that is beyond the limit.
    
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the density
    write(filename_density_ref, '(a)') 'theoritical_dissipation.dat'
    open(13, file="unitary_tests/dissipation/density.gnuplot")
    write(13,*) "set terminal pngcairo enhanced size 800, 600 font ',20'"
    write(13,*) 'set xlabel "distance (AU)"'
    write(13,*) 'set ylabel "Surface density (g/cm^2)"'
    write(13,*) 'set grid'
    write(13,*) 'set xrange [', INNER_BOUNDARY_RADIUS, ':', OUTER_BOUNDARY_RADIUS, ']'
    write(13,*) 'set yrange [', density_min, ':', density_max, ']'
    
    do k=1, nb_time
      write(filename_density, '(a,i0.5,a)') 'surface_density',k,'.dat'
      write(output_density, '(a,i0.5,a)') 'surface_density',k,'.png'
      write(output_time, time_format) (time(k)/365.25d0)
      
      write(13,*) "set output '",trim(output_density),"'"
      write(13,*) 'set title "T=', trim(output_time),' years"'
      write(13,*) "plot '",trim(filename_density),"' using 1:2 with lines linetype 0 linewidth 4 notitle"
      write(13,*) ""
    end do
    close(13)
  
  end subroutine test_disk_dissipation
  
  subroutine test_turbulence_torque(stellar_mass)
  
  use turbulence
  use utilities, only : get_polar_coordinates, get_mean, get_stdev, get_histogram
  
  implicit none
  
  ! Inputs
  real(double_precision), intent(in) :: stellar_mass ! in [Msun * K2]
  
  integer, parameter :: nb_points = 10000 ! the time through which we compute the turbulence
  integer, parameter :: nb_bins = 100 ! the number of bins for the histogram of the turbulence torque
  
  real(double_precision) :: initial_time = 0.d0 ! in days
  type(PlanetProperties) :: p_prop ! various properties of a planet
  real(double_precision) :: radius_planet
  real(double_precision) :: theta_planet
  real(double_precision), dimension(3) :: turbulence_acceleration ! in AU/DAY^2
  real(double_precision), dimension(3) :: position
  real(double_precision), dimension(3) :: velocity
  real(double_precision) :: time ! in days
  real(double_precision), dimension(nb_points) :: turbulence_torque ! in AU^2/DAY^2
  
  ! planet parameters
  real(double_precision), parameter :: a = 6.d0 ! in AU
  real(double_precision), parameter :: mass = 1. * K2 * EARTH_MASS ! in [Msun * K2]
  real(double_precision) :: delta_t = 365.25d0!365.25d0 * a**1.5d0 ! the timestep in days between two calculation of the turbulence torque. Must be greater than the coherence time of the turbulence to make the test of the turbulence usefull

  
  ! histogram temp values
  real(double_precision) :: delta_bin
  real(double_precision), dimension(nb_bins) :: bin_x_values, bin_y_values, gauss_fit
  real(double_precision) :: mean, stdev, y_max

  integer :: i ! For loops
  
  !------------------------------------------------------------------------------
  write(*,*) 'Test of turbulence Torque'
  
  call init_turbulence(initial_time)
  
  position(1:3) = 0.d0
  velocity(1:3) = 0.d0
  
  position(1) = a

  ! We generate cartesian coordinate for the given mass and semi major axis
  velocity(2) = sqrt((stellar_mass + mass) / position(1))
  
  ! we store in global parameters various properties of the planet
  call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
   p_prop=p_prop) ! Output
  
  call get_polar_coordinates(position(1), position(2), position(3), radius_planet, theta_planet) 
  ! This radius must be used instead of p_prop%radius because p_prop is the same during the calculation of the derivative
  
  time = initial_time
  open(10, file="unitary_tests/turbulence_torque.dat")
  do i=1, nb_points
    time = time + delta_t
    
    call get_turbulence_acceleration(time, p_prop, position, turbulence_acceleration)
    turbulence_torque(i) = (-sin(theta_planet) * turbulence_acceleration(1) + cos(theta_planet) * turbulence_acceleration(2)) * a 
    
    write(*,*) turbulence_acceleration(1:3)
    write(*,*) theta_planet, radius_planet, position(1:3), velocity(1:3)
    call print_planet_properties(p_prop)
    
    write(10,*) time, turbulence_torque(i)
  end do
  close(10)
  
  call get_histogram(turbulence_torque(1:nb_points), bin_x_values(1:nb_bins), bin_y_values(1:nb_bins)) 
  
  ! We calculate the mean and stdev value of the data set and then generate a supposed gaussian to see if this function fit the datas
  ! the mean of the data set must be 0. So in order to check that, the mean is fixed to 0, to see if the gaussian looks nice.
  delta_bin = (bin_x_values(2) - bin_x_values(1))
  mean = get_mean(turbulence_torque(1:nb_points))
  stdev = get_stdev(turbulence_torque(1:nb_points))
  y_max = 1.d0 / (stdev * sqrt(TWOPI)) * delta_bin
  do i=1,nb_bins
    gauss_fit(i) = y_max * exp(-(bin_x_values(i))**2 / (2. * stdev**2))
  end do
  
  open(10, file="unitary_tests/turbulence_torque.hist")
  do i=1, nb_bins
    write(10,*) bin_x_values(i), bin_y_values(i), gauss_fit(i)
  end do
  close(10)
    
  open(10, file="unitary_tests/turbulence_torque.gnuplot")
  write(10,'(a)') 'set terminal pdfcairo enhanced'
  write(10,'(a)') 'set output "turbulence_torque.pdf"'
  write(10,'(a)') '!rm "turbulence_torque.pdf"'
  write(10,'(a)') 'set xlabel "torque [AU^2/DAY^2]"'
  write(10,'(a)') 'set ylabel "density of probability"'
  write(10,'(a)') 'set grid'
  write(10,'(5(a,es10.2e2))') 'set label " data : {/Symbol m}=',mean,', {/Symbol s}=',stdev,'\n&
                                         & fit  : {/Symbol m}=0, {/Symbol s}=',stdev,'" at graph 0, graph 0.9'
  write(10,'(a)') 'plot "turbulence_torque.hist" using 1:2 with boxes linestyle 3 title "turbulence torque", \'
  write(10,'(a)') '"turbulence_torque.hist" using 1:3 with lines linestyle 1 title "gaussian fit"'
  close(10)
  
  end subroutine test_turbulence_torque
  
  subroutine test_turbulence_mode()
  
  use turbulence
  
  implicit none
    
  integer, parameter :: nb_points = 10000 ! the time through which we compute the turbulence
  
  real(double_precision) :: initial_time = 0.d0
  type(TurbulenceMode) :: turb_mode ! a turbulence mode
  integer :: i ! For loops
  !------------------------------------------------------------------------------
  write(*,*) 'Test of turbulence modes'
  
  call init_turbulence(initial_time)
  
  open(10, file="unitary_tests/turbulence_mode.dat")
  write(10,*) '# wavenumber, r, phi, lifetime, radial_extent, chi'
  
  do i=1, nb_points
    call init_mode(initial_time, turb_mode)
    write(10,*) turb_mode%wavenumber, turb_mode%r, turb_mode%phi, turb_mode%lifetime, &
                turb_mode%radial_extent, turb_mode%chi
  end do
  close(10)
  
!~   open(10, file="unitary_tests/turbulence_mode.gnuplot")
!~   write(10,'(a)') 'set terminal pdfcairo enhanced'
!~   write(10,'(a)') '!rm "turbulence_mode.pdf"'
!~   write(10,'(a)') 'set output "turbulence_mode.pdf"'
!~   write(10,'(a)') 'set xlabel "torque [AU^2/DAY^2]"'
!~   write(10,'(a)') 'set ylabel "density of probability"'
!~   write(10,'(a)') 'set grid'
!~   write(10,'(a)') 'plot "turbulence_mode.dat" using 1:2 with lines '
!~   close(10)
  
  end subroutine test_turbulence_mode

! %%% Physical behaviour %%%
  subroutine study_temperature_profile()
! Subroutine that test the finding of the temperature profile and store a plot of the temperature profile of the disk
! A gnuplot file and a data file are created to display the temperature profile.

    implicit none
    
    integer :: j ! for loops
    
    write(*,*) 'Test of the temperature profile'
    
    open(10, file="unitary_tests/temperature_profile.gnuplot")
    open(11, file="unitary_tests/temperature_index.gnuplot")

    do j=10,11
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) '!rm "temperature_profile.pdf"'
    write(10,*) "set output 'temperature_profile.pdf'"
    write(11,*) '!rm "temperature_index.pdf"'
    write(11,*) "set output 'temperature_index.pdf'"

    do j=10,11
      write(j,*) 'set xlabel "semi major axis a (in AU)"'
      write(j,*) 'set nokey'
    end do
    
    write(10,*) 'set logscale y'
    write(10,*) 'set logscale x'
    write(10,*) 'set ylabel "Temperature [K]"'
    
    
    write(11,*) 'set ylabel "Temperature law index"'
    
    do j=10,11
      write(j,*) 'set grid'
    end do

    write(10,*) 'plot "../temperature_profile.dat" using 1:2 with lines notitle'
    
    write(11,*) 'plot "../temperature_profile.dat" using 1:3 with lines notitle'
    
    close(10)
    close(11)
  
  end subroutine study_temperature_profile
  
  subroutine study_scaleheight_profile
    implicit none
    
    integer :: j
    
    write(*,*) 'Test of the scaleheight'
    
  open(10, file="unitary_tests/scaleheight_profile.gnuplot")
  open(11, file="unitary_tests/aspect_ratio.gnuplot")
  
  do j=10,11
    write(j,*) "set terminal pdfcairo enhanced"
  end do
  
  write(10,*) '!rm "scaleheight_profile.pdf"'
  write(10,*) "set output 'scaleheight_profile.pdf'"
  write(11,*) '!rm "aspect_ratio_profile.pdf"'
  write(11,*) "set output 'aspect_ratio_profile.pdf'"
  
  do j=10,11
    write(j,*) 'set xlabel "semi major axis a (in AU)"'
    write(j,*) 'set nokey'
  end do
  
  write(10,*) 'set ylabel "Scaleheight H [AU]"'
  
  write(11,*) 'set ylabel "Aspect ratio h=H/R"'
  
  do j=10,11
    write(j,*) 'set grid'
  end do

  write(10,*) "plot '../scaleheight_profile.dat' using 1:2 with line linetype -1 notitle, \"
  write(10,*) "     '' using 1:(-$2) with line linetype -1 notitle"
  
  write(11,*) "plot '../scaleheight_profile.dat' using 1:3 with line linetype -1 notitle, \"
  write(11,*) "     '' using 1:(-$3) with line linetype -1 notitle"
  
  close(10)
  close(11)
  
  
  end subroutine study_scaleheight_profile
  
  subroutine study_optical_depth_profile()
! Subroutine that test the finding of the optical depth profile and store a plot of the temperature profile of the disk
! A gnuplot file and a data file are created to display the temperature profile.

    implicit none
    
    open(10, file="unitary_tests/optical_depth_profile.gnuplot")
    
    write(*,*) 'Test of the optical depth'

    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) '!rm "optical_depth_profile.pdf"'
    write(10,*) "set output 'optical_depth_profile.pdf'"
    
    write(10,*) 'set xlabel "semi major axis a (in AU)"'
    write(10,*) 'set nokey'

    
    write(10,*) 'set ylabel "Optical depth {/Symbol t}"'
        
    write(10,*) 'set grid'


    write(10,*) 'plot "../temperature_profile.dat" using 1:5 with lines notitle'
    
    
    close(10)
  
  end subroutine study_optical_depth_profile
  
  subroutine study_thermal_diffusivity_profile()
! Subroutine that test the finding of the optical depth profile and store a plot of the temperature profile of the disk
! A gnuplot file and a data file are created to display the temperature profile.

    implicit none
    
    write(*,*) 'Test of the thermal diffusivity'
    
    open(10, file="unitary_tests/thermal_diffusivity_profile.gnuplot")
    

    write(10,*) "set terminal pdfcairo enhanced"
    
    write(10,*) '!rm "thermal_diffusivity_profile.pdf"'
    write(10,*) "set output 'thermal_diffusivity_profile.pdf'"
    write(10,*) 'set xlabel "semi major axis a (in AU)"'
    write(10,*) 'set nokey'
    write(10,*) 'set ylabel "Thermal diffusivity {/Symbol c} [AU^2/day]"'
        
    write(10,*) 'set grid'

    write(10,*) 'plot "../temperature_profile.dat" using 1:4 with lines notitle'
    
    close(10)
  
  end subroutine study_thermal_diffusivity_profile

  subroutine study_opacity_profile
  ! subroutine that test the function 'get_opacity'
  
  ! Return:
  ! a data file 'opacity.dat' 
  ! and an associated gnuplot file 'opacity.gnuplot' that display values for get_opacity for a range of p values.
    implicit none
    
    real(double_precision), dimension(5) :: bulk_density = (/ 1.d-5, 1.d-6, 1.d-7, 1.d-8, 1.d-9/)
!~     real(double_precision), dimension(5) :: bulk_density = (/1.d-9, 1.d-10, 1.d-11, 1.d-12, 1.d-13/)
    real(double_precision), dimension(5) :: opacity, opacity_zhu, opacity_bell
    real(double_precision) :: temperature
    
    real(double_precision), parameter :: T_min = 10.d0
    real(double_precision), parameter :: T_max = 10000.d0
    integer, parameter :: nb_points = 1000
    real(double_precision), parameter :: T_step = (T_max/T_min) ** (1/(nb_points-1.d0))
    
    integer :: i,j ! for loops
    
    write(*,*) 'Test of the opacity'
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/opacity.dat')
    open(11, file='unitary_tests/opacity_comparison.dat')
    write(10,*) "# Correspond to the opacity type :", OPACITY_TYPE
    write(10,*) '# Temperature (K) ; Opacity for bulk density from 1e-5 to 1e-9 by power of ten'
    write(11,*) '# Temperature (K) ; Opacity bell ; zhu, bell; zhu'
    
    do i=1, nb_points
      temperature = T_min * T_step ** (i-1)
      
      do j=1,5
        opacity(j) = get_opacity(temperature, bulk_density(j)) * MSUN / AU**3
        opacity_bell(j) = get_opacity_bell_lin_1994(temperature, bulk_density(j)) * MSUN / AU**3
        opacity_zhu(j) = get_opacity_zhu_2009(temperature, bulk_density(j)) * MSUN / AU**3
      end do
      
      write(10,*) temperature, opacity(1), opacity(2), opacity(3), opacity(4), opacity(5)
      write(11,*) temperature, opacity_bell(1), opacity_bell(2), opacity_bell(3), &
      opacity_bell(4), opacity_bell(5), opacity_zhu(1), opacity_zhu(2), opacity_zhu(3), opacity_zhu(4), opacity_zhu(5)

    end do
    close(10)
    close(11)
    
    open(10, file="unitary_tests/opacity.gnuplot")
    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) "set output 'opacity.pdf'"
    write(10,*) 'set xlabel "Temperature T"'
    write(10,*) 'set ylabel "Opacity {/Symbol k}"'
    write(10,*) 'set logscale x'
    write(10,*) 'set logscale y'
    write(10,*) 'set grid'
    write(10,*) 'set xrange [', T_min, ':', T_max, ']'
    write(10,*) "plot 'opacity.dat' using 1:2 with lines title '{/Symbol r}=10^{-5}',\"
    write(10,*) "     '' using 1:3 with lines title '{/Symbol r}=10^{-6}',\"
    write(10,*) "     '' using 1:4 with lines title '{/Symbol r}=10^{-7}',\"
    write(10,*) "     '' using 1:5 with lines title '{/Symbol r}=10^{-8}',\"
    write(10,*) "     '' using 1:6 with lines title '{/Symbol r}=10^{-9}'"
    
    close(10)
    
    open(11, file="unitary_tests/opacity_comparison.gnuplot")
    write(11,*) "set terminal pdfcairo enhanced"
    write(11,*) "set output 'opacity_comparison.pdf'"
    write(11,*) 'set xlabel "Temperature T"'
    write(11,*) 'set ylabel "Opacity {/Symbol k}"'
    write(11,*) 'set logscale x'
    write(11,*) 'set logscale y'
    write(11,*) 'set grid'
    write(11,*) 'set termoption dashed'
    write(11,*) 'set title "Dotted line : Bell \& Lin (1994) ; Solid line : Zhu \& Hartmann (2009)"'
    write(11,*) 'set xrange [', T_min, ':', T_max, ']'
    write(11,*) "plot 'opacity_comparison.dat' using 1:2  with lines linetype 3 linecolor 1 notitle,\"
    write(11,*) "'' using 1:3  with lines linetype 3 linecolor 2 notitle,\"
    write(11,*) "'' using 1:4  with lines linetype 3 linecolor 3 notitle,\"
    write(11,*) "'' using 1:5  with lines linetype 3 linecolor 4 notitle,\"
    write(11,*) "'' using 1:6  with lines linetype 3 linecolor 5 notitle,\"
    write(11,*) "'' using 1:7  with lines linetype 1 linecolor 1 title '{/Symbol r}=10^{-5}',\"
    write(11,*) "'' using 1:8  with lines linetype 1 linecolor 2 title '{/Symbol r}=10^{-6}',\"
    write(11,*) "'' using 1:9  with lines linetype 1 linecolor 3 title '{/Symbol r}=10^{-7}',\"
    write(11,*) "'' using 1:10 with lines linetype 1 linecolor 4 title '{/Symbol r}=10^{-8}',\"
    write(11,*) "'' using 1:11 with lines linetype 1 linecolor 5 title '{/Symbol r}=10^{-9}'"
    
    close(11)
    
  end subroutine study_opacity_profile

  subroutine study_torques_fixed_a(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file '.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass ! in [msun * K2]
    
    integer, parameter :: nb_mass = 100
    real(double_precision), parameter :: mass_min = 1.d0 ! in earth mass
    real(double_precision), parameter :: mass_max = 75.d0 ! in earth mass
    real(double_precision), parameter :: mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    
    real(double_precision), parameter :: a = 5.2d0
    
    real(double_precision) :: mass, total_torque, corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
    real(double_precision) :: lindblad_torque_units, corotation_torque_units, total_torque_units
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: j ! for loops
    
    write(*,*) 'Evolution of the torque for a fixed distance "a"'
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/torques_fixed_a.dat')
    open(11, file='unitary_tests/ref_torque_fixed_a.dat')
    open(12, file='unitary_tests/torques_fixed_a_units.dat')
    open(13, file='unitary_tests/specific_torque_fixed_a.dat')
    
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
       
      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to migration
      call get_torques(stellar_mass, mass, p_prop, corotation_torque, &
          lindblad_torque, torque_ref, ecc_corot=ecc_corot)
      
      
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
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) "set output 'torques_fixed_a.pdf'"
    write(11,*) "set output 'ref_torque_fixed_a.pdf'"
    write(12,*) "set output 'torques_fixed_a_units.pdf'"
    write(13,*) "set output 'specific_torque_fixed_a.pdf'"
    
    do j=10,13
      write(j,*) 'set xlabel "Planet mass (m_{earth})"'
      write(j,'(a,f5.1,a)') 'set title "for a semi major axis a= ',a,' AU"'
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

    write(10,*) "plot 'torques_fixed_a.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(10,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(10,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"
    
    write(11,*) "plot 'ref_torque_fixed_a.dat' using 1:2 with lines"
    
    write(12,*) "plot 'torques_fixed_a_units.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(12,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(12,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"
    
    write(13,*) "plot 'specific_torque_fixed_a.dat' using 1:2 with lines"
    
    close(10)
    close(11)
    close(12)
    close(13)

    
  end subroutine study_torques_fixed_a

  subroutine study_torques_fixed_m(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'torques_fixed_m.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass ! in [msun * K2]
    
    integer, parameter :: nb_a = 400
    real(double_precision), parameter :: a_min = 1d0 ! in AU
    real(double_precision), parameter :: a_max = 9d0 ! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision), parameter :: mass = 10. * EARTH_MASS * K2
    
    real(double_precision) :: a, total_torque, corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
    real(double_precision) :: lindblad_torque_units, corotation_torque_units, total_torque_units
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: j ! for loops
    
    write(*,*) 'Evolution of the torque for a fixed planet mass "m"'

    position(1:3) = 0.d0
    velocity(1:3) = 0.d0
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/torques_fixed_m.dat')
    open(11, file='unitary_tests/ref_torque_fixed_m.dat')
    open(12, file='unitary_tests/torques_fixed_m_units.dat')
    
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
      velocity(2) = sqrt((stellar_mass + mass) / position(1))
      
      ! we store in global parameters various properties of the planet
      call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
       p_prop=p_prop) ! Output
      
      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to migration
      call get_torques(stellar_mass, mass, p_prop, corotation_torque, &
          lindblad_torque, torque_ref, ecc_corot=ecc_corot)
      
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
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) "set output 'torques_fixed_m.pdf'"
    write(11,*) "set output 'ref_torque_fixed_m.pdf'"
    write(12,*) "set output 'torques_fixed_m_units.pdf'"
    
    do j=10,12
      write(j,*) 'set xlabel "semi major axis a (in AU)"'
      write(j,'(a,f4.1,a)') 'set title "for planet mass = ',mass / (EARTH_MASS * K2),' m_{earth}"'
    end do
    
    write(10,*) 'set ylabel "torque [{/Symbol G}_0]"'
    
    write(11,*) 'set ylabel "reference torque {/Symbol G}_0 [M_s.AU^2.day^{-2}]"'
    write(11,*) 'set nokey'
    
    write(12,*) 'set ylabel "torque [M_s.AU^2.day^{-2}]"'
    
    do j=10,12
      write(j,*) 'set grid'
      write(j,*) 'set xrange [', a_min, ':', a_max, ']'
    end do

    write(10,*) "plot 'torques_fixed_m.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(10,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(10,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"
    
    write(11,*) "plot 'ref_torque_fixed_m.dat' using 1:2 with lines"
        
    write(12,*) "plot 'torques_fixed_m_units.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(12,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(12,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"
    
    close(10)
    close(11)
    close(12)

    
  end subroutine study_torques_fixed_m
  
  subroutine study_viscosity(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'torques_fixed_m.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass ! in [msun * K2]
    
    integer, parameter :: nb_a = 100
    real(double_precision), parameter :: a_min = 0.1d0 ! in AU
    real(double_precision), parameter :: a_max = 15d0 ! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision), parameter :: mass = 10.d0 * EARTH_MASS * K2
    real(double_precision), parameter :: num2phys = AU**2 / DAY ! convert numerical viscosity to CGS viscosity
    
    real(double_precision) :: a
    real(double_precision) :: viscosity
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    
    integer :: i ! for loops
    
    write(*,*) 'Evolution of the viscosity in function of the position'

    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/viscosity.dat')

    
    write(10,'(a)') '# a in AU ; viscosity'

    
    do i=1,nb_a
      a = (a_min + a_step * (i - 1.d0))
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a
      
      ! We generate cartesian coordinate for the given mass and semi major axis
      velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
      
      ! we store in global parameters various properties of the planet
      call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
       p_prop=p_prop) ! Output
      
      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to migration
      viscosity = get_temp_viscosity(omega=p_prop%omega, scaleheight=p_prop%scaleheight, radius=a)
      
      write(10,*) a, viscosity * num2phys
    end do
    
    close(10)
    
    
    open(10, file="unitary_tests/viscosity.gnuplot")
    
    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) "set output 'viscosity.pdf'"
    
    write(10,*) 'set xlabel "semi major axis a (in AU)"'
    write(10,*) 'set ylabel "viscosity [cm^2/s]"'
    
    write(10,*) 'set grid'
    write(10,*) 'set xrange [', a_min, ':', a_max, ']'
    
    write(10,'(a,i2,a, f4.2,a)') ' plot "viscosity.dat" using 1:2 with lines notitle'
    
    
    close(10)


    
  end subroutine study_viscosity
  
  subroutine study_eccentricity_effect_on_corotation(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'torques_fixed_m.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass ! in [msun * K2]
    
    integer, parameter :: nb_a = 400
    real(double_precision), parameter :: a_min = 0.1d0 ! in AU
    real(double_precision), parameter :: a_max = 15d0 ! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision), parameter :: mass = 10.d0 * EARTH_MASS * K2
    
    real(double_precision) :: a, total_torque, corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    real(double_precision), dimension(5) :: eccentricities
    real(double_precision), dimension(5) :: corot_damping
    real(double_precision), dimension(10) :: outputs
    real(double_precision) :: x_s
    
    integer :: i,j ! for loops
    
    write(*,*) 'Evolution of the torque for a fixed planet mass "m"'
    x_s = 1.d0
    eccentricities(1) = 0.d0   ! in units of x_s 
    eccentricities(2) = 0.25d0 ! in units of x_s 
    eccentricities(3) = 0.5d0  ! in units of x_s
    eccentricities(4) = 1.d0   ! in units of x_s
    eccentricities(5) = 2.d0   ! in units of x_s
    
    do i=1, 5
      corot_damping(i) = get_corotation_damping(e=eccentricities(i), x_s=x_s)
    end do

    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/eccentricity_effect_on_corotation.dat')

    
    write(10,'(10(a,f4.2),a)') '# a in AU ; G_c(e/x_s=',eccentricities(1),') ; G_t(e/x_s=',eccentricities(1),&
                        ') ; G_c(e/x_s=',eccentricities(2),') ; G_t(e/x_s=',eccentricities(2),&
                        ') ; G_c(e/x_s=',eccentricities(3),') ; G_t(e/x_s=',eccentricities(3),&
                        ') ; G_c(e/x_s=',eccentricities(4),') ; G_t(e/x_s=',eccentricities(4),&
                        ') ; G_c(e/x_s=',eccentricities(5),') ; G_t(e/x_s=',eccentricities(5),')'


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
      
      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to migration
      call get_torques(stellar_mass, mass, p_prop, corotation_torque, &
          lindblad_torque, torque_ref, ecc_corot=ecc_corot)
      
      total_torque = lindblad_torque + corotation_torque
      
      do i=1,5
        outputs(2*i-1) = corot_damping(i) * corotation_torque
        outputs(2*i) = lindblad_torque + corot_damping(i) * corotation_torque
      end do
      write(10,*) a, lindblad_torque, outputs
    end do
    
    close(10)
    
    
    open(10, file="unitary_tests/eccentricity_effect_on_corotation.gnuplot")
    
    write(10,*) "set terminal pdfcairo enhanced"
    
    write(10,*) "set output 'eccentricity_effect_on_corotation.pdf'"

    

    write(10,*) 'set xlabel "semi major axis a (in AU)"'
    write(10,'(a,f4.1,a)') 'set title "for planet mass = ',mass / (EARTH_MASS * K2),' m_{earth}"'

    
    write(10,*) 'set ylabel "torque [{/Symbol G}_0]"'
    
    

    write(10,*) 'set grid'
    write(10,*) 'set xrange [', a_min, ':', a_max, ']'
    
    write(10,*) 'plot \'
    do i=1, 5
      write(10,'(a,i2,a, f4.2,a)') "'eccentricity_effect_on_corotation.dat' using 1:",2*i+1, &
      " with lines title '{/Symbol G}_c ; e/x_s=",eccentricities(i),"',\"
      write(10,'(a,i2,a, f4.2,a)') "'eccentricity_effect_on_corotation.dat' using 1:",2*i+2, &
      " with lines title '{/Symbol G}_{tot} ; e/x_s=",eccentricities(i),"',\"
    end do
    write(10,*) "'eccentricity_effect_on_corotation.dat' using 1:2 with lines title '{/Symbol G}_L'"
    
    
    close(10)


    
  end subroutine study_eccentricity_effect_on_corotation
  
  subroutine study_ecc_corot(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'torques_fixed_m.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
  
    use orbital_elements, only : mco_el2x
    
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass ! in [msun * K2]
    
    integer, parameter :: nb_e = 400
    real(double_precision), parameter :: e_min = 1d-5 ! eccentricity
    real(double_precision), parameter :: e_max = 1.d0-1.d-5 ! eccentricity
    real(double_precision), parameter :: e_step = (e_max / e_min)**(1 / (nb_e - 1.d0))
    
    real(double_precision), parameter :: mass = 1.d0 * EARTH_MASS * K2
    real(double_precision), parameter :: a = 6.d0 ! in AU
    real(double_precision) :: I = 0.d0 ! in radians
    
    real(double_precision) :: corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
    real(double_precision) :: e, q, gm, x_s, gamma_eff, Q_p, n
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: j ! for loops
    
    write(*,*) 'Evolution of the torque in function of the eccentricity'
    
    gm = stellar_mass + mass ! Masses are already in [msun * K2]
    position(1:3) = 0.d0
    velocity(1:3) = 0.d0
    
    n = 0.d0
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/ecc_corot.dat')
    
    write(10,*) '# eccentricity ; ecc correction ; total_torque'

    
    do j=1,nb_e
      e = e_min * e_step**(j - 1)
      q = a * (1 - e)
      
      ! We calculate position and velocities in function of orbital elements
      call mco_el2x(gm,q,e,I,0.d0, n, 0.d0,& ! Inputs
                    position(1),position(2),position(3),velocity(1),velocity(2),velocity(3)) ! Outputs
      
      ! we store in global parameters various properties of the planet
      call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
       p_prop=p_prop) ! Output
      
      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to migration
      call get_torques(stellar_mass=stellar_mass, mass=mass, p_prop=p_prop, corotation_torque=corotation_torque, &
          lindblad_torque=lindblad_torque, Gamma_0=torque_ref, ecc_corot=ecc_corot)
      
      !------------------------------------------------------------------------------
      ! Q is needed by the lindblad torque. We set Q for m ~ 2 /3 h (45): 
      Q_p = TWOTHIRD * p_prop%chi / (p_prop%aspect_ratio * p_prop%scaleheight**2 * p_prop%omega) ! p_prop%aspect_ratio**3 * p_prop%radius**2 = aspect_ratio * scaleheight**2
      !------------------------------------------------------------------------------

      gamma_eff = 2.d0 * Q_p * ADIABATIC_INDEX / (ADIABATIC_INDEX * Q_p + 0.5d0 * &
      sqrt(2.d0 * sqrt((ADIABATIC_INDEX * ADIABATIC_INDEX * Q_p * Q_p + 1.d0)**2 - 16.d0 * Q_p * Q_p * (ADIABATIC_INDEX - 1.d0)) &
      + 2.d0 * ADIABATIC_INDEX * ADIABATIC_INDEX * Q_p * Q_p - 2.d0))

      !------------------------------------------------------------------------------

      x_s = X_S_PREFACTOR / gamma_eff**0.25d0 * sqrt(mass / p_prop%aspect_ratio)
      
      
              
      write(10,*) e / x_s, ecc_corot
    end do
    
    close(10)

    
    
    open(10, file="unitary_tests/ecc_corot.gnuplot")

    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) "set output 'ecc_corot.pdf'"
!~     write(10,'(2(a,f4.1),2a)') 'set title "mass = ',mass / (EARTH_MASS * K2),' m_{earth} ; a = ',a ,'AU"'

    write(10,*) 'set xlabel "e/x_s"'
    write(10,*) 'set grid'
    write(10,*) 'set xrange [0:2]'

    write(10,*) "plot 'ecc_corot.dat' using 1:2 with lines title 'eccentricity correction'"    
    close(10)

    
  end subroutine study_ecc_corot


  subroutine study_torques(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
  
    use contour
    
    
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass ! in [msun * K2]
    
    integer, parameter :: nb_mass = 150
    real(double_precision), parameter :: mass_min = 0.1 * EARTH_MASS
    real(double_precision), parameter :: mass_max = 60. * EARTH_MASS
    real(double_precision), parameter :: mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    real(double_precision), dimension(nb_mass) :: mass
    
    integer, parameter :: nb_points = 100
    real(double_precision) :: a_min
    real(double_precision) :: a_max
    ! step for log sampling
    real(double_precision) :: a_step
    real(double_precision), dimension(nb_points) :: a
    
    real(double_precision), dimension(nb_points, nb_mass) :: total_torque, total_torque_units
    
    real(double_precision) :: corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    integer :: deltai, deltaj ! separation between to outputs for the vector part of the display (we do not want a lot of points)
    real(double_precision) :: vector_limit ! The size limit of the vectors, to prevent chevauching.
    
    !------------------------------------------------------------------------------
    a_min = INNER_BOUNDARY_RADIUS + INNER_SMOOTHING_WIDTH
    a_max = 20.d0
    a_step = (a_max - a_min) / (nb_points-1.d0)
!~     a_step = (a_max / a_min)**(1 / (nb_points - 1.d0))
    
    !------------------------------------------------------------------------------
    write(*,*) 'Evolution of the total, lindblad and corotation torques depending on the planet mass and distance'
    
    position(1:3) = 0.d0
    velocity(1:3) = 0.d0
    
    ! We want to have only 10 points both in x and y for vector outputs
    deltai = nb_points / 15
    deltaj = nb_mass / 10
    
    vector_limit = a_step * (deltai - 1)
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/corotation_torque.dat')
    open(11, file='unitary_tests/total_torque.dat')
    open(12, file='unitary_tests/total_torque_units.dat')
    open(13, file='unitary_tests/lindblad_torque.dat')
    open(14, file='unitary_tests/ref_torque.dat')
    open(15, file='unitary_tests/vector_total_torque.dat')
    
    write(10,*) '# semi major axis (AU) ; mass in earth mass ; corotation torque (no dim)'
    write(11,*) '# semi major axis (AU) ; mass in earth mass ; total torque (no dim)'
    write(12,*) '# semi major axis (AU) ; mass in earth mass ; total torque in M_s.AU^2.day^{-2}'
    write(13,*) '# semi major axis (AU) ; mass in earth mass ; lindblad torque (no dim)'
    write(14,*) '# semi major axis (AU) ; mass in earth mass ; reference torque in M_s.AU^2.day^{-2}'
    write(15,*) '# semi major axis (AU) ; mass in earth mass ; Delta x ; Delta y'
    
    
    do i=1, nb_points ! loop on the position
      a(i) = a_min + a_step * (i-1)
!~       a(i) = a_min * a_step**(i - 1)
      
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a(i)
      
      
      do j=1,nb_mass
        mass(j) = (mass_min + mass_step * (j - 1.d0)) * K2 ! mass in [Msun * K2]
        
        ! We generate cartesian coordinate for the given mass and semi major axis
        velocity(2) = sqrt((stellar_mass + mass(j)) / position(1))
        
        ! we store in global parameters various properties of the planet
        call get_planet_properties(stellar_mass=stellar_mass, & ! Input
         mass=mass(j), position=position(1:3), velocity=velocity(1:3),& ! Input
         p_prop=p_prop) ! Output
         
!~         write(*,*) a(i),p_prop%semi_major_axis 
!~         call print_planet_properties(p_prop) 
!~         stop
        ! If the semi major axis is not well determined, we display a warning and give the values
        if (abs(a(i)-p_prop%semi_major_axis)/a(i).gt.1d-2) then
          write(*,*) 'Warning: for manual a=',a(i), 'we get :'
          call print_planet_properties(p_prop) 
        end if
         
        !------------------------------------------------------------------------------
        ! Calculation of the acceleration due to migration
        
        call get_torques(stellar_mass, mass(j), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref, ecc_corot=ecc_corot) ! Output
        
        total_torque(i,j) = lindblad_torque + corotation_torque
        total_torque_units(i,j) = torque_ref * total_torque(i,j)
        
                
        write(10,*) a(i), mass(j) / (EARTH_MASS*K2), corotation_torque
        write(11,*) a(i), mass(j) / (EARTH_MASS*K2), total_torque(i,j)
        write(12,*) a(i), mass(j) / (EARTH_MASS*K2), total_torque_units(i,j)
        write(13,*) a(i), mass(j) / (EARTH_MASS*K2), lindblad_torque
        write(14,*) a(i), mass(j) / (EARTH_MASS*K2), torque_ref
        
        if ((modulo(i-int(deltai/2.),deltai).eq.0).and.(modulo(j-int(deltaj/2.),deltaj).eq.0)) then
          write(15,*) a(i), mass(j) / (EARTH_MASS*K2), 0, &
                      sign(min(sqrt(abs(total_torque(i,j))),vector_limit), total_torque(i,j)), 0, 0
        end if
        
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
    close(15)
    
    
    ! We want to get the contour for the torque equal to 0 for total torque both in physical dimension of units of gamma_0
    mass(1:nb_mass) = mass(1:nb_mass) / (EARTH_MASS*K2)
    call get_contour(total_torque, a, mass,'unitary_tests/contour_total_torque.dat', 0.d0) ! This line opens a file (unit=10), thus this unit must not be used at that time
    
    open(10, file="unitary_tests/corotation_torque.gnuplot")
    open(11, file="unitary_tests/total_torque.gnuplot")
    open(12, file="unitary_tests/total_torque_units.gnuplot")
    open(13, file="unitary_tests/lindblad_torque.gnuplot")
    open(14, file="unitary_tests/ref_torque.gnuplot")

    do j=10,14
      write(j,*) "set terminal pngcairo crop enhanced size 1200, 1000 font ',20'"
    end do
    
    write(10,*) "set output 'corotation_torque.png'"
    write(11,*) "set output 'total_torque.png'"
    write(12,*) "set output 'total_torque_units.png'"
    write(13,*) "set output 'lindblad_torque.png'"
    write(14,*) "set output 'ref_torque.png'"

    do j=10,14
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
      write(j,*) 'set grid xtics ytics linetype 0'
      write(j,*) 'set xrange [', a_min, ':', a_max, ']'
      write(j,*) 'set yrange [', mass_min / EARTH_MASS, ':', mass_max / EARTH_MASS, ']'
    end do

    write(10,*) "splot 'corotation_torque.dat' with pm3d notitle"
    
    write(11,*) "splot 'total_torque.dat' with pm3d notitle, \"
    write(11,*) "      'contour_total_torque.dat' with line linetype -1 linewidth 1 notitle, \"
    write(11,*) "      'vector_total_torque.dat' with vector notitle head filled linestyle -1"
    
    write(12,*) "splot 'total_torque_units.dat' with pm3d notitle, \"
    write(12,*) "      'contour_total_torque.dat' with line linetype -1 linewidth 1 notitle, \"
    write(12,*) "      'vector_total_torque.dat' with vector notitle head filled linestyle -1"
    
    write(13,*) "splot 'lindblad_torque.dat' with pm3d notitle"
    
    write(14,*) "splot 'ref_torque.dat' with pm3d notitle"
        
    
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    
  end subroutine study_torques
  
  subroutine study_dissipation_at_one_location()
  ! Plot the evolution of the surface density at one location through the dissipation
  
  ! Global parameters
  ! dissipation_timestep : the timestep between two computation of the disk [in days]
  ! INNER_BOUNDARY_RADIUS : the inner radius of the various profiles (all based on the radius profile)
  ! OUTER_BOUNDARY_RADIUS : the outer radius of the various profiles (all based on the radius profile)
  ! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
  ! surface_density_profile : values of the density in MSUN/AU^2 for each value of the 'a' sample
  ! distance_sample : values of 'a' in AU
  ! viscosity : the viscosity of the disk in [cm^2/s] /!\ This parameter is modified here (which must not be done otherwise)
  
    ! Return: evolution of density at one location with time
    
  implicit none
    
    ! time sample
    real(double_precision), parameter :: t_min = 0.d0 ! time in years
    real(double_precision), parameter :: t_max = 1.d7 ! time in years
    real(double_precision), dimension(:), allocatable :: time, time_temp ! time in days
    real(double_precision), dimension(:), allocatable :: density, density_temp ! surface density in MSUN/AU^2
    integer :: time_size ! the size of the array 'time'. 
    
    real(double_precision) :: density_position = 2.d0 ! in AU
    real(double_precision) :: x_radius
    integer :: density_idx ! the closest index to get around 1 AU
    
    real(double_precision) :: next_dissipation_step
    
    integer :: k ! for loops
    integer :: nb_time
    integer :: error ! to retrieve error, especially during allocations
    !------------------------------------------------------------------------------
    write(*,*) 'Evolution of surface density at one location through dissipation'
    
    
    if ((density_position .gt. INNER_BOUNDARY_RADIUS) .and. (density_position .lt. distance_sample(NB_SAMPLE_PROFILES-1))) then
      
      x_radius = 2.d0 * sqrt(density_position)
      ! in the range
      density_idx = 1 + int((x_radius - x_sample(1)) / X_SAMPLE_STEP) ! X_SAMPLE_STEP being a global value, x_sample also
      density_position = distance_sample(density_idx)
      
    else
      write(error_unit,*) 'Error: "density_position" is not in the disk range'
      call exit(2)
    end if
    
    
    !------------------------------------------------------------------------------
    k = 1
    time_size = 512 ! the size of the array. 
    allocate(time(time_size), stat=error)
    allocate(density(time_size), stat=error)
    time(1) = t_min * 365.25d0 ! days
    
    
    do while (time(k).lt.(t_max*365.d0))
      ! We calculate the temperature profile for the current time (because the surface density change in function of time)
      
      
      ! we get the new dissipated surface density profile. For that goal, we dissipate as many times as needed to reach the required time for the next frame.
      
      if (DISSIPATION_TYPE.ne.0) then
        call dissipate_disk(time(k), next_dissipation_step)
      end if
      
      if (.not.disk_effect) then
        next_dissipation_step = t_max*365.d0
      end if
      
      density(k) = surface_density_profile(density_idx)
      
      ! we expand the 'time' array if the limit is reached
      if (k.eq.time_size) then
        ! If the limit of the 'time' array is reach, we copy the values in a temporary array, allocate with a double size, et paste the 
        ! old values in the new bigger array
        allocate(time_temp(time_size), stat=error)
        allocate(density_temp(time_size), stat=error)
        time_temp(1:time_size) = time(1:time_size)
        density_temp(1:time_size) = density(1:time_size)
        
        deallocate(time, stat=error)
        deallocate(density, stat=error)
        time_size = time_size * 2
        
        allocate(time(time_size), stat=error)
        allocate(density(time_size), stat=error)
        time(1:time_size/2) = time_temp(1:time_size/2)
        density(1:time_size/2) = density_temp(1:time_size/2)
        deallocate(time_temp, stat=error)
        deallocate(density_temp, stat=error)
      end if
      
      time(k+1) = next_dissipation_step ! days
      
      
      k = k + 1 ! We increment the integer that point the time in the array (since it's a 'while' and not a 'do' loop)
    end do

    nb_time = k - 1 ! since for the last step with incremented 'k' by one step that is beyond the limit.
    
    open(13, file="unitary_tests/local_density_dissipation.dat")
    write(13,*) '# time (years) ; surface density (g/cm^2)'
    do k=1,nb_time
      write(13,*) time(k)/365.25d0, density(k) * SIGMA_NUM2CGS
    end do
    close(13)
    
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the density
    open(13, file="unitary_tests/local_density_dissipation.gnuplot")
    write(13,*) 'set terminal pdfcairo enhanced'
    write(13,*) 'set output "local_density_dissipation.pdf"'
    write(13,*) 'set xlabel "Time (years)"'
    write(13,*) 'set ylabel "Surface density (g/cm^2)"'
    write(13,*) 'set grid'
!~     write(13,*) 'set xrange [', INNER_BOUNDARY_RADIUS, ':', OUTER_BOUNDARY_RADIUS, ']'
!~     write(13,*) 'set yrange [', density_min, ':', density_max, ']'
    
    write(13,'(a,f5.1,a)') 'set title "a=', density_position,' AU"'
    write(13,*) "plot 'local_density_dissipation.dat' using 1:2 with lines notitle"
    close(13)
  
  end subroutine study_dissipation_at_one_location
  
  subroutine study_influence_of_dissipation_on_torque(stellar_mass)
  ! subroutine that plot several values depending on the properties of the disk during its dissipation
  
  ! Global parameters
  ! dissipation_timestep : the timestep between two computation of the disk [in days]
  ! INNER_BOUNDARY_RADIUS : the inner radius of the various profiles (all based on the radius profile)
  ! OUTER_BOUNDARY_RADIUS : the outer radius of the various profiles (all based on the radius profile)
  ! NB_SAMPLE_PROFILES : number of points for the sample of radius of the temperature profile
  ! surface_density_profile : values of the density in MSUN/AU^2 for each value of the 'a' sample
  ! temperature_profile : values of the temperature in K for each value of the 'a' sample
  
  ! Return: (in the folder 'dissipation' of the current working directory where the tests are launched)
  !  data files with the following paterns :
  !    surface_density*.dat : each file correspond to a data file of the surface density at one timestep (the number correspond to the index)
  !    temperature*.dat : each file correspond to a data file of the temperature at one timestep (the number correspond to the index)
  !    total_torque*.dat : each file correspond to a data file of the total torque exerted by the disk at one timestep (the number correspond to the index)
  !    contour*.dat : each file correspond to the zero torque lines of the corresponding total_torque*.dat file.
  ! IN ADDITION : There are 3 gnuplot scripts to generate the output files for the surface density, the temperature and the total torque. The last one need both total_torque and contour data files.
  
    use contour
    
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass ! in [msun * K2]
    
    ! mass sample
    integer, parameter :: nb_mass = 150
    real(double_precision), parameter :: mass_min = 0.1d0 * EARTH_MASS
    real(double_precision), parameter :: mass_max = 20.d0 * EARTH_MASS
    real(double_precision), parameter :: mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    real(double_precision), dimension(nb_mass) :: mass, mass_earth
    
    ! orbital distance sample
    integer, parameter :: nb_points = 100
    real(double_precision), parameter :: a_min = 0.01d0
    real(double_precision), parameter :: a_max = 5.d0
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_points - 1.d0)
    real(double_precision), dimension(nb_points) :: a
    
    ! time sample
    real(double_precision), parameter :: t_min = 0.d0 ! time in years
    real(double_precision), parameter :: t_max = 2.7d6 ! time in years
    real(double_precision), dimension(:), allocatable :: time, time_temp ! time in days
    integer :: time_size ! the size of the array 'time'. 
    
    ! Array to store values of the torque. Used to find the zero torque line
    real(double_precision), dimension(nb_points, nb_mass) :: total_torque
    
    real(double_precision) :: corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    real(double_precision) :: temp_min, temp_max, density_min, density_max, torque_min, torque_max
    character(len=80) :: filename_torque, filename_density, filename_temperature, filename_contour
    character(len=80) :: filename_density_ref, filename_temperature_ref
    character(len=80) :: output_torque, output_density, output_temperature, output_time, time_format, purcent_format
    integer :: time_length ! the length of the displayed time, usefull for a nice display
    real(double_precision) :: next_dissipation_step

    
    integer :: i,j,k ! for loops
    integer :: nb_time ! The total number of 't' values. 
    integer :: error ! to retrieve error, especially during allocations
    logical :: isDefined
    !------------------------------------------------------------------------------
    write(*,*) 'Evolution of the total torque during the dissipation of the disk'
    
    write(*,*) '  Force Initialisation again'
    FIRST_CALL = .true.
    call init_globals(stellar_mass=stellar_mass, time=0.d0)
    
    inquire(file='dissipation', exist=isDefined)
    
    ! We create the folder 'dissipation' if he doesn't exists.
    if (.not.isDefined) then
      call system("mkdir dissipation")
    end if
    
    call system("rm dissipation/*")
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    write(filename_density_ref, '(a,i0.5,a)') 'dissipation/surface_density',0,'.dat'
    write(filename_temperature_ref, '(a,i0.5,a)') 'dissipation/temperature',0,'.dat'
    
    ! we store in a .dat file the temperature profile
    call store_temperature_profile(filename=filename_temperature_ref)
    call store_density_profile(filename=filename_density_ref)
    

    temp_max = 1d4
    temp_min = 1.d0
    
    ! We want the extremum of the surface density during the dissipation of the disk in order to have nice plots
    density_min = 0.d0
    density_max = maxval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * SIGMA_NUM2CGS

    
!~     call system("rm dissipation/*")
    
    do i=1, nb_points ! loop on the position
      a(i) = a_min + a_step * (i-1)
    end do
    
    do j=1,nb_mass
      mass(j) = (mass_min + mass_step * (j - 1.d0)) * K2
      mass_earth(j) = (mass_min + mass_step * (j - 1.d0)) / EARTH_MASS
    end do
    
    ! We want to know the max size of the time display in order to have a nice display, with filled spaces in the final plots
    write(output_time, '(i0)') int(t_max)
    time_length = len(trim(output_time))
    write(time_format, *) '(i',time_length,'.',time_length,')'
    write(purcent_format, *) '(i',time_length,'"/",i',time_length,'," years ; k = ",i5)'
    
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
      
      ! we get the new dissipated surface density profile. For that goal, we dissipate as many times as needed to reach the required time for the next frame.
      call dissipate_disk(time(k), next_dissipation_step)
      
      if (k.eq.time_size) then
        ! If the limit of the array is reach, we copy the values in a temporary array, allocate with a double size, et paste the 
        ! old values in the new bigger array
        allocate(time_temp(time_size), stat=error)
        time_temp(1:time_size) = time(1:time_size)
        deallocate(time, stat=error)
        time_size = time_size * 2
        allocate(time(time_size), stat=error)
        time(1:time_size/2) = time_temp(1:time_size/2)
        deallocate(time_temp, stat=error)
      end if
      
      write(*,purcent_format) int(time(k)/365.25d0), int(t_max), k
      
      if (.not.disk_effect) then
        next_dissipation_step = t_max * 365.d0
      end if

      time(k+1) = next_dissipation_step ! days
      ! we get the temperature profile.
      call calculate_temperature_profile()
      
      ! we store in a .dat file the temperature profile
      call store_temperature_profile(filename=filename_temperature)
      call store_density_profile(filename=filename_density)
      
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
          call get_corotation_torque(stellar_mass, mass(j), p_prop, corotation_torque, lindblad_torque, torque_ref, ecc_corot)
          
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
        torque_max = 7.5d0 ! maxval(total_torque)
      end if
      
      k = k + 1 ! We increment the integer that point the time in the array (since it's a 'while' and not a 'do' loop)
    end do
    
    nb_time = k - 1 ! since for the last step with incremented 'k' by one step that is beyond the limit.
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the total torque
    open(11, file="dissipation/total_torque.gnuplot")
    write(11,*) "set terminal pngcairo enhanced size 1024, 768 font ',20'"
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
      write(output_time, time_format) int(time(k)/365.25)
      
      write(11,*) "set output '",trim(output_torque),"'"
      write(11,*) 'set title "total torque {/Symbol G}_{tot}/{/Symbol G}_0 T=', trim(output_time),' years"'
      write(11,*) "splot '",trim(filename_torque),"' with pm3d notitle, \"
      write(11,*) "      '",trim(filename_contour),"' with line linetype -1 linewidth 2 title '{/Symbol G}=0'"
      write(11,*) ""
    end do
    close(11)
    
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the temperature profile
    write(filename_temperature_ref, '(a,i0.5,a)') 'temperature',0,'.dat'
    open(12, file="dissipation/temperature.gnuplot")
    write(12,*) "set terminal pngcairo enhanced size 1024, 768 font ',20'"
    write(12,*) 'set xlabel "semi major axis (AU)"'
    write(12,*) 'set ylabel "Temperature (K)"'
    write(12,*) 'set mxtics 10'
    write(12,*) 'set mytics 10'
    write(12,*) 'set grid xtics ytics mxtics mytics linetype -1, 0'
    write(12,*) 'set logscale x'
    write(12,*) 'set logscale y'
    write(12,*) 'set xrange [', INNER_BOUNDARY_RADIUS, ':', OUTER_BOUNDARY_RADIUS, ']'
    write(12,*) 'set yrange [', temp_min, ':', temp_max, ']'
    
    do k=1, nb_time
      write(filename_temperature, '(a,i0.5,a)') 'temperature',k,'.dat'
      write(output_temperature, '(a,i0.5,a)') 'temperature',k,'.png'
      write(output_time, time_format) int(time(k)/365.25)
      
      write(12,*) "set output '",trim(output_temperature),"'"
      write(12,*) 'set title "T=', trim(output_time),' years"'
      write(12,*) "plot '",trim(filename_temperature_ref),"' using 1:2 with lines linetype 0 linewidth 3 notitle, \"
      write(12,*) "     '",trim(filename_temperature),"' using 1:2 with lines  linetype 1 notitle"
      write(12,*) ""
    end do
    close(12)
    
    !------------------------------------------------------------------------------
    write(filename_density_ref, '(a,i0.5,a)') 'surface_density',0,'.dat'
    ! Gnuplot script to output the frames of the density
    open(13, file="dissipation/density.gnuplot")
    write(13,*) "set terminal pngcairo enhanced size 1024, 768 font ',20'"
    write(13,*) 'set xlabel "semi major axis (AU)"'
    write(13,*) 'set ylabel "Surface density (g/cm^2)"'
    write(13,*) 'set grid'
    write(13,*) 'set xrange [', INNER_BOUNDARY_RADIUS, ':', OUTER_BOUNDARY_RADIUS, ']'
    write(13,*) 'set yrange [', density_min, ':', density_max, ']'
    
    do k=1, nb_time
      write(filename_density, '(a,i0.5,a)') 'surface_density',k,'.dat'
      write(output_density, '(a,i0.5,a)') 'surface_density',k,'.png'
      write(output_time, time_format) int(time(k)/365.25)
      
      write(13,*) "set output '",trim(output_density),"'"
      write(13,*) 'set title "T=', trim(output_time),' years"'
      write(13,*) "plot '",trim(filename_density_ref),"' using 1:2 with lines linetype 0 linewidth 3 notitle, \"
      write(13,*) "     '",trim(filename_density),"' using 1:2 with lines linetype 1 notitle"
      write(13,*) ""
    end do
    close(13)
  
  end subroutine study_influence_of_dissipation_on_torque

end program test_disk
