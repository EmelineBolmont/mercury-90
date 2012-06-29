! Program that test different functions implemented in the module user_module.
program test_disk
  use types_numeriques
  use disk
    
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
    
    implicit none
    
    real(double_precision) :: stellar_mass ! [Msun * K2]

    stellar_mass = 1.d0 * K2 ! [Msun * K2]
    
    write(*,*) 'Initialisation'
    call init_globals(stellar_mass=stellar_mass, time=0.d0)
    ! Note that the initial density profile and temperature profile are calculated inside the 'init_globals' routine.
    
!~     ! We force the value to be interesting for our tests
!~     TORQUE_TYPE = 'arctan_indep' ! 'real', 'linear_indep', 'arctan_indep', 'mass_dependant', 'manual'
    
    ! We want to show the torque profile. It is important to check which value has been declared in 'TORQUE_TYPE'
    call study_torques(stellar_mass)
    
!~     ! we store in a .dat file the temperature profile
    call store_temperature_profile(filename='temperature_profile.dat')
    call store_density_profile(filename='density_profile.dat')
    call store_scaleheight_profile()
!~     
!~     ! Unitary tests
    call test_functions_FGK()
    call test_function_zero_temperature(stellar_mass)
    call test_temperature_interpolation()
    call test_manual_torque_interpolation()
    call test_density_interpolation()
    call test_retrieval_of_orbital_elements(stellar_mass)
    call test_turbulence()

    
    ! Physical values and plots
    call study_opacity_profile()
    call study_torques_fixed_a(stellar_mass)
    call study_torques_fixed_m(stellar_mass)
    call study_temperature_profile(stellar_mass)
    call study_optical_depth_profile(stellar_mass)
    call study_thermal_diffusivity_profile(stellar_mass)
    call study_scaleheight_profile()
    
    ! Test dissipation
!~     call test_viscous_dissipation()
!~     call test_exponential_dissipation()
!~     call study_dissipation_of_the_disk(stellar_mass)
    
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
    velocity(2) = sqrt((stellar_mass + mass) / position(1))
    
    ! we store in global parameters various properties of the planet
    call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
     p_prop=p_prop) ! Output
    
    ! We calculate this value outside the function because we only have to do this once per step (per radial position)
    prefactor = - (9.d0 * p_prop%nu * p_prop%sigma * p_prop%omega**2 / 16.d0)
    
    write(10,*) '# properties of the disk at the location of the planet that influence the value of the temperature'
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

  subroutine test_manual_torque_interpolation()
  
    implicit none
    
    integer, parameter :: nb_a = 1000
    real(double_precision), parameter :: a_min = 0.d0 ! in AU
    real(double_precision), parameter :: a_max = 100.d0! in AU
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_a - 1.d0)
    
    real(double_precision) :: a, sigma, sigma_index
    
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
    
    open(10, file='unitary_tests/test_torque_interpolation.dat')
    do j=1,NB_SAMPLE_PROFILES
      
      
      write(10,*) distance_sample(j), torque_profile(j)
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="unitary_tests/torque_interpolation.gnuplot")
    

    write(10,*) 'set terminal wxt enhanced'
    write(10,*) 'set xlabel "semi major axis a (in AU)"'

    
    write(10,*) 'set ylabel "torque [?]"'
      
    write(10,*) 'set grid'
    write(10,*) 'set yrange [',minval(torque_profile)-1., ':', maxval(torque_profile)+1,']'

    write(10,*) 'plot "test_torque_interpolation.dat" using 1:2 with lines linetype -1 title "Interpolated profile",\'
    write(10,*) '     "../torque_profile.dat" using 1:2 with points linetype 1 pointtype 2 linewidth 3 title "Read Profile"'
        

    
    write(10,*) "#pause -1 # wait until a carriage return is hit"
    write(10,*) "set terminal pdfcairo enhanced"

    
    write(10,*) '!rm "torque_interpolation.pdf"'
    write(10,*) "set output 'torque_interpolation.pdf'"

    
    
    write(10,*) "replot # to generate the output file"
    
    close(10)
  
  end subroutine test_manual_torque_interpolation

  subroutine test_retrieval_of_orbital_elements(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
  
    use contour
    
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass
    
    integer, parameter :: nb_mass = 150
    real(double_precision), parameter :: mass = 10.d0 * EARTH_MASS * K2
    
    integer, parameter :: nb_points = 100
    real(double_precision), parameter :: a_min = 0.01
    real(double_precision), parameter :: a_max = 50.
    
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
    open(10, file='unitary_tests/test_retrieval_aeI.dat')
    
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
    
    write(10,*) 'set xlabel "Initial a"'
    write(11,*) 'set xlabel "Initial e"'
    write(12,*) 'set xlabel "Initial I"'
    
    write(10,*) 'set ylabel "Final a"'
    write(11,*) 'set ylabel "Final e"'
    write(12,*) 'set ylabel "Final I"'
    
    do j=10,12
      write(j,*) 'set mxtics 10'
      write(j,*) 'set grid xtics ytics mxtics'
    end do
    
    write(10,*) 'plot "test_retrieval_a.dat" using 1:2 with lines title "a",\'
    write(10,*) '                        " " using 1:3 with lines title "r",\'
    write(10,*) '                        " " using 1:4 with lines title "a2",\'
    write(10,*) '                        " " using 1:5 with lines title "r2"'
    write(11,*) 'plot "test_retrieval_e.dat" using 1:2 with lines title "e",\'
    write(11,*) '"" using 1:3 with lines title "e2"'
    write(12,*) 'plot "test_retrieval_I.dat" using 1:2 with lines title "I",\'
    write(12,*) '"" using 1:3 with lines title "I2"'
    
    do j=10,12
      write(10,*) '#pause -1 # wait until a carriage return is hit'
      write(10,*) 'set terminal pdfcairo enhanced'
    end do
    
    write(10,*) 'set output "retrieval_a.pdf"'
    write(11,*) 'set output "retrieval_e.pdf"'
    write(12,*) 'set output "retrieval_I.pdf"'
    
    do j=10, 12
      write(10,*) 'replot'
      close(j)
    end do
    
  end subroutine test_retrieval_of_orbital_elements

  subroutine test_viscous_dissipation()
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
    
    ! time sample
    real(double_precision), parameter :: t_min = 0. ! time in years
    real(double_precision), parameter :: t_max = 1.d5 ! time in years
    real(double_precision), dimension(:), allocatable :: time, time_temp ! time in days
    integer, parameter :: max_frames = 100 ! Parameter to have a control over the number of output frames. (I had near 80 000 once)
    integer :: nb_dissipation_per_step
    integer :: time_size ! the size of the array 'time'. 
    
    real(double_precision), parameter :: a = 1.
    
    
    
    real(double_precision) :: density_min, density_max
    character(len=80) :: filename_density, filename_density_ref
    character(len=80) :: output_density, output_time, time_format, purcent_format
    integer :: time_length ! the length of the displayed time, usefull for a nice display
    
    ! For the definition of the diffusion function
    real(double_precision) :: r_0 = 50.d0 ! position of the dirac function (in AU)
    integer :: r_0i ! integer that correspond to the closest 'a' value from r_0
    real(double_precision), parameter :: t_0 = 0.016 ! the time at which the diffusion begin. Because we test a dirac function that 
                                               ! cannot be computed so we start after the beginning of the theoritical diffusion.
    real(double_precision), parameter :: nu = 1.e-6 ! the viscosity for the dissipation test with a dirac in AU^2/day. (to compare with kley, 1999)
    real(double_precision) :: t_nu ! the viscous spreading time
    real(double_precision) :: sigma_prefactor ! the prefactor of the theoritical function for the spreading of the dirac function (taken from kley, 1999)
    real(double_precision) :: tau ! dimensionless time
    real(double_precision) :: x ! normalized radius
    real(double_precision) :: Ix, Kx, Ixp, Kxp ! the value of the I_1/4 modified bessel function (Kx, Ixp, Kxp are not used here be the subroutine return them anyway)
    
    real(double_precision), dimension(5), parameter :: ref_tau = (/0.018, 0.044, 0.074, 0.120, 0.184/)
    
    integer :: i,k ! for loops
    integer :: nb_time ! The total number of 't' values. 
    integer :: error ! to retrieve error, especially during allocations
    real(double_precision) :: tmp, tmp2(5) ! temporary value for various calculations
    !------------------------------------------------------------------------------
    write(*,*) 'Evolution of the total torque during the dissipation of the disk'
    
    ! we force the dissipation type 
    DISSIPATION_TYPE = 1
    
    ! First, we calculate the value of the viscosity in cm^2/s in order to modify the global variable 'viscosity'
    viscosity = nu / DAY * AU**2 ! in cm^2/s
    write(*,'(a,es10.2,a,es10.2,a)') 'with nu=', get_viscosity(a), 'AU^2/day (',viscosity,'cm^2/s)'
    
    call system("rm unitary_tests/dissipation/*")
    t_nu = r_0**2 / (12.d0 * get_viscosity(a))
    sigma_prefactor = 1 / (PI * r_0**2)
    
    ! We search for the closest radial value from the desired R_0 value
    do i=1, NB_SAMPLE_PROFILES-1
      if ((distance_sample(i).lt.r_0).and.(distance_sample(i+1).ge.r_0)) then 
        r_0i = i
      end if
    end do
    
    if (.not.(((r_0i.gt.0).and.(r_0i.le.NB_SAMPLE_PROFILES)))) then
      write(*,*) 'Error: the r_0 value is not in the radial profile sample'
      write(*,*) r_0i
      stop
    end if
    
    r_0 = distance_sample(r_0i)
    
    ! We store reference profiles
    write(filename_density_ref, '(a)') 'unitary_tests/dissipation/theoritical_dissipation.dat'
    open(10, file=filename_density_ref)
    write(10,*) '# radial length (no dim) ; profile for t=',ref_tau 
    do i=1,200
      x = 0.01 * i
      
      do k=1,5
        tau = ref_tau(k)
        call bessik(2*x/tau,0.25d0,Ix, Kx, Ixp, Kxp)
    
        tmp2(k) = sigma_prefactor / (tau * x**0.25d0) * exp(-(1.d0 + x**2) / tau) * Ix * MSUN / AU**2
      end do
      write(10,*) x, tmp2
    end do
    close(10)

    ! Before dissipating the disk, we erase the 'normal' density profile with a dirac one. 
    tau = t_0

    do i=1, NB_SAMPLE_PROFILES
      x = distance_sample(i)/r_0
      call bessik(2*x/tau,0.25d0,Ix, Kx, Ixp, Kxp)
            
      tmp = sigma_prefactor / (tau * x**0.25d0) * exp(-(1.d0 + x**2) / tau) * Ix
      surface_density_profile(i) = tmp
    end do
    
    ! In some cases, Ix return Infinity. Thus, the surface_density get NaN. To avoid this, we search for NaN and replace by a 0 surface density.
    where(isnan(surface_density_profile))
      surface_density_profile = 0.d0
    end where
    
!~     surface_density_profile(1:NB_SAMPLE_PROFILES) = 0.d0
!~     surface_density_profile(r_0i) = 1/(2*PI*r_0)
    
!~     ! We store the initial profile of the surface density in a reference file.
!~     write(filename_density_ref, '(a,i0.5,a)') 'unitary_tests/dissipation/surface_density',0,'.dat'
!~     call store_density_profile(filename=filename_density_ref)
    
    ! We want the extremum of the surface density during the dissipation of the disk in order to have nice plots
    density_min = minval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
    density_max = maxval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
    
    ! We want to know the max size of the time display in order to have a nice display, with filled spaces in the final plots
    write(output_time, '(i0)') int(t_max, kind=8)
    time_length = len(trim(output_time))
    write(time_format, *) '(f',time_length,'.3)'
    write(purcent_format, *) '(i',time_length,'"/",i',time_length,'," years")'
    
    !------------------------------------------------------------------------------
    k = 1
    time_size = 512 ! the size of the array. 
    allocate(time(time_size), stat=error)
    time(1) = t_min * 365.25d0 ! days
    write(*,*) X_SAMPLE_STEP, get_viscosity(a)
    dissipation_timestep = 0.9d0 * X_SAMPLE_STEP**2 / (4 * get_viscosity(a)) ! a correction factor of 0.5 has been applied. No physical reason to that, just intuition and safety
    ! TODO if the viscosity is not constant anymore, the formulae for the dissipation timestep must be changed
    nb_dissipation_per_step = int((t_max - t_min)* 365.25d0 / (max_frames * dissipation_timestep))
    
    ! if the effective number of timestep is less than the max allowed, then we force to have at least one dissipation timestep 
    ! at each value of 'k' 
    if (nb_dissipation_per_step.eq.0) then
      nb_dissipation_per_step = 1
    end if
    
    do while (time(k).lt.(t_max*365.d0))
      ! We open the file where we want to write the outputs
      write(filename_density, '(a,i0.5,a)') 'unitary_tests/dissipation/surface_density',k,'.dat'    
      
      ! we store in a .dat file the temperature profile
      call store_density_profile(filename=filename_density)
      
      if (k.eq.1) then
        ! We want the extremum of the surface density during the dissipation of the disk in order to have nice plots
        density_min = minval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
        density_max = maxval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
      end if
      
      ! We calculate the temperature profile for the current time (because the surface density change in function of time)
      
      
      ! we get the new dissipated surface density profile. For that goal, we dissipate as many times as needed to reach the required time for the next frame.
      do i=1,nb_dissipation_per_step
        call dissipate_density_profile() ! global parameter 'dissipation_timestep' must exist !
      end do
      
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
      
      
      time(k+1) = time(k) + dfloat(nb_dissipation_per_step) * dissipation_timestep ! days
      
      
      k = k + 1 ! We increment the integer that point the time in the array (since it's a 'while' and not a 'do' loop)
    end do

    nb_time = k - 1 ! since for the last step with incremented 'k' by one step that is beyond the limit.
    
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the density
    write(filename_density_ref, '(a)') 'theoritical_dissipation.dat'
    open(13, file="unitary_tests/dissipation/density.gnuplot")
    write(13,*) "set terminal pngcairo enhanced size 800, 600"
    write(13,*) 'set xlabel "distance (adim)"'
    write(13,*) 'set ylabel "Surface density (g/cm^2)"'
    write(13,*) 'set grid'
    write(13,*) 'set xrange [', 0, ':', 2, ']'
    write(13,*) 'set yrange [', density_min, ':', density_max, ']'
    
    do k=1, nb_time
      write(filename_density, '(a,i0.5,a)') 'surface_density',k,'.dat'
      write(output_density, '(a,i0.5,a)') 'surface_density',k,'.png'
      write(output_time, time_format) (time(k)/t_nu) + t_0
      
      write(13,*) "set output '",trim(output_density),"'"
      write(13,*) 'set title "T=', trim(output_time),' (adim)"'
      write(13,*) "plot '",trim(filename_density_ref),"' using 1:2 with lines title 'theo t=0.018', \"
      write(13,*) " '",trim(filename_density_ref),"' using 1:3 with lines title 'theo t=0.044', \"
      write(13,*) " '",trim(filename_density_ref),"' using 1:4 with lines title 'theo t=0.074', \"
      write(13,*) " '",trim(filename_density_ref),"' using 1:5 with lines title 'theo t=0.120', \"
      write(13,*) " '",trim(filename_density_ref),"' using 1:6 with lines title 'theo t=0.184', \"
      write(13,*) " '",trim(filename_density),"' using ($1/",r_0,"):2 with lines linetype 0 linewidth 4 notitle"
      write(13,*) ""
    end do
    close(13)
  
  end subroutine test_viscous_dissipation
  
  subroutine test_exponential_dissipation()
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
    
    ! time sample
    real(double_precision), parameter :: t_min = 0. ! time in years
    real(double_precision), parameter :: t_max = 1.d7 ! time in years
    real(double_precision), dimension(:), allocatable :: time, time_temp ! time in days
    integer, parameter :: max_frames = 100 ! Parameter to have a control over the number of output frames. (I had near 80 000 once)
    integer :: nb_dissipation_per_step
    integer :: time_size ! the size of the array 'time'. 
    
    real(double_precision), parameter :: a = 1.
    
    
    
    real(double_precision) :: density_min, density_max
    character(len=80) :: filename_density, filename_density_ref
    character(len=80) :: output_density, output_time, time_format, purcent_format
    integer :: time_length ! the length of the displayed time, usefull for a nice display
    
    integer :: i,k ! for loops
    integer :: nb_time ! The total number of 't' values. 
    integer :: error ! to retrieve error, especially during allocations
    real(double_precision) :: tmp, tmp2(5) ! temporary value for various calculations
    !------------------------------------------------------------------------------
    write(*,*) 'Evolution of the total torque during the dissipation of the disk'
    
    ! we force the dissipation type 
    DISSIPATION_TYPE = 2
    TAU_DISSIPATION = 1e6 ! in years
    
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
    
    dissipation_timestep = (t_max - t_min) * 365.25d0 / dfloat(max_frames)
    
    
    do while (time(k).lt.(t_max*365.d0))
      ! We open the file where we want to write the outputs
      write(filename_density, '(a,i0.5,a)') 'unitary_tests/dissipation/surface_density',k,'.dat'    
      
      ! we store in a .dat file the temperature profile
      call store_density_profile(filename=filename_density)
      
      if (k.eq.1) then
        ! We want the extremum of the surface density during the dissipation of the disk in order to have nice plots
        density_min = minval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
        density_max = maxval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2
      end if
      
      ! We calculate the temperature profile for the current time (because the surface density change in function of time)
      
      
      ! we get the new dissipated surface density profile. For that goal, we dissipate as many times as needed to reach the required time for the next frame.
      call exponential_decay_density_profile()
      
      
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
      
      
      time(k+1) = time(k) + dissipation_timestep ! days
      
      
      k = k + 1 ! We increment the integer that point the time in the array (since it's a 'while' and not a 'do' loop)
    end do

    nb_time = k - 1 ! since for the last step with incremented 'k' by one step that is beyond the limit.
    
    !------------------------------------------------------------------------------
    ! Gnuplot script to output the frames of the density
    write(filename_density_ref, '(a)') 'theoritical_dissipation.dat'
    open(13, file="unitary_tests/dissipation/density.gnuplot")
    write(13,*) "set terminal pngcairo enhanced size 800, 600"
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
  
  end subroutine test_exponential_dissipation
  
  subroutine test_turbulence()
  
  use turbulence
  use utilities, only : get_polar_coordinates
  
  implicit none
  
  integer, parameter :: nb_points = 10000 ! the time through which we compute the turbulence
  integer, parameter :: nb_bins = 1000 ! the number of bins for the histogram of the turbulence torque
  
  real(double_precision) :: initial_time = 0.d0
  type(PlanetProperties) :: p_prop ! various properties of a planet
  real(double_precision) :: radius_planet
  real(double_precision) :: theta_planet
  real(double_precision), dimension(3) :: turbulence_acceleration
  real(double_precision), dimension(3) :: position
  real(double_precision), dimension(3) :: velocity
  real(double_precision) :: time ! in days
  real(double_precision), dimension(nb_points) :: turbulence_torque
  
  ! planet parameters
  real(double_precision), parameter :: a = 4. ! in AU
  real(double_precision), parameter :: mass = 15. * EARTH_MASS * K2 ! in [Msun * K2]
  real(double_precision), parameter :: stellar_mass = 1. * K2 ! in [Msun * K2]
  real(double_precision) :: delta_t = 365.25d0 * a**1.5d0 ! the timestep in days between two calculation of the turbulence torque. Must be greater than the coherence time of the turbulence to make the test of the turbulence usefull

  
  ! histogram temp values
  real(double_precision) :: delta_bin, max_value, min_value
  real(double_precision), dimension(nb_bins) :: bin_x_values, bin_y_values
  integer :: index_bin

  integer :: i ! For loops
  
  call init_turbulence(initial_time)
  
  !------------------------------------------------------------------------------
  position(1:3) = 0.d0
  velocity(1:3) = 0.d0
  
  position(1) = a

  ! We generate cartesian coordinate for the given mass and semi major axis
  velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
  
  ! we store in global parameters various properties of the planet
  call get_planet_properties(stellar_mass=stellar_mass, mass=mass, position=position(1:3), velocity=velocity(1:3),& ! Input
   p_prop=p_prop) ! Output
  
  call get_polar_coordinates(position(1), position(2), position(3), radius_planet, theta_planet) 
  ! This radius must be used instead of p_prop%radius because p_prop is the same during the calculation of the derivative
  
  do i=1, nb_points
    time = initial_time + i * delta_t
    
    call get_turbulence_acceleration(time, p_prop, position, turbulence_acceleration)
    turbulence_torque(i) = - sin(theta_planet) * turbulence_acceleration(1) + cos(theta_planet) * turbulence_acceleration(2)
    
  end do
  
  ! Line only used for TESTS !!! To be erased afterwards
!~   turbulence_torque( 1) = 3.5
!~   turbulence_torque( 2) = 0.
!~   turbulence_torque( 3) = 4.
!~   turbulence_torque( 4) = 2.1
!~   turbulence_torque( 5) = 2.9
!~   turbulence_torque( 6) = 2.5
!~   turbulence_torque( 7) = 3.1
!~   turbulence_torque( 8) = 3.0
!~   turbulence_torque( 9) = 2.2
!~   turbulence_torque(10) = 2.3
  ! Je dois avoir
!~   0.5 0.1
!~   1.5 0
!~   2.5 0.7
!~   3.5 0.2
  

!~   ! We initialize the values of the counting array
!~   bin_y_values(1:nb_bins) = 0
!~   
!~ 
!~   ! From the list of values, we get the values for the histogram
!~   max_value = maxval(turbulence_torque(1:nb_points))
!~   min_value = minval(turbulence_torque(1:nb_points))
!~   
!~   delta_bin = (max_value - min_value) / float(nb_bins)
!~   
!~   do i=1,nb_bins
!~     bin_x_values(i) = min_value + (i-0.5d0) * delta_bin
!~   end do
!~     
!~   do i=1, nb_points
!~     ! if the value is exactly equal to max_value, we force the index to be nb_bins
!~     index_bin = min(floor((turbulence_torque(i) - min_value) / delta_bin)+1, nb_bins)
!~     
!~     ! With floor, we get the immediate integer below the value. Since the index start at 1, we add 1 to the value, because the 
!~     ! calculation will get from 0 to the number of bins. Thus, for the max value, we will get nb_bins +1, which is not possible. 
!~     ! As a consequence, we take the lower value between the index and nb_bins, to ensure that for the max value, we get an index 
!~     ! of nb_bins.
!~     
!~     bin_y_values(index_bin) = bin_y_values(index_bin) + 1
!~   end do
!~ 
!~   ! We normalize the histogram
!~   bin_y_values(1:nb_bins) = bin_y_values(1:nb_bins) / float(nb_points)
  
  open(10, file="unitary_tests/test_turbulence_torque.dat")
  do i=1, nb_points
!~     write(10,*) bin_x_values(i), bin_y_values(i)
    write(10,*) turbulence_torque(i)
  end do
  close(10)
  
  open(10, file="unitary_tests/turbulence_torque.gnuplot")
  write(10,*) 'set xlabel "torque"'
  write(10,*) 'set ylabel "density of probability"'
  write(10,*) 'plot "test_turbulence_torque.dat" using 1:2 with boxes'
  write(10,*) '#pause -1 # wait until a carriage return is hit'
  write(10,*) 'set terminal pdfcairo enhanced'
  write(10,*) 'set output "turbulence_torque.pdf"'
  write(10,*) 'replot'
  close(10)
  
  end subroutine test_turbulence
  

! %%% Physical behaviour %%%
  subroutine study_temperature_profile(stellar_mass)
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
  
  end subroutine study_temperature_profile
  
  subroutine study_scaleheight_profile
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
  
  
  end subroutine study_scaleheight_profile
  
  subroutine study_optical_depth_profile(stellar_mass)
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
  
  end subroutine study_optical_depth_profile
  
  subroutine study_thermal_diffusivity_profile(stellar_mass)
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
  
  end subroutine study_thermal_diffusivity_profile

  subroutine study_opacity_profile
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
    
  end subroutine study_opacity_profile

  subroutine study_torques_fixed_a(stellar_mass)
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
    open(10, file='unitary_tests/study_torques_fixed_a.dat')
    open(11, file='unitary_tests/test_ref_torque_fixed_a.dat')
    open(12, file='unitary_tests/study_torques_fixed_a_units.dat')
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
       
      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to migration
      select case(TORQUE_TYPE)
        case('real') ! The normal torque profile, calculated form properties of the disk
          call get_corotation_torque(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        ! for retrocompatibility, 'mass_independant' has been added and refer to the old way of defining a mass-indep convergence zone
        case('linear_indep', 'mass_independant') ! a defined torque profile to get a mass independant convergence zone
          call get_corotation_torque_linear_indep(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        case('arctan_indep') ! a defined torque profile to get a mass independant convergence zone
          call get_corotation_torque_arctan_indep(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        case('mass_dependant')
          call get_corotation_torque_mass_dep_CZ(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        case('manual')
          call get_corotation_torque_manual(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
          
        case default
          write(*,*) 'Warning: The torque rule cannot be found.'
          write(*,*) 'Given value :', TORQUE_TYPE
      end select
      
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

    write(10,*) "plot 'study_torques_fixed_a.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(10,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(10,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"
    
    write(11,*) "plot 'test_ref_torque_fixed_a.dat' using 1:2 with lines"
    
    write(12,*) "plot 'study_torques_fixed_a_units.dat' using 1:2 with lines title '{/Symbol G}_c',\"
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

    
  end subroutine study_torques_fixed_a

  subroutine study_torques_fixed_m(stellar_mass)
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
    
    real(double_precision), parameter :: mass = 15. * EARTH_MASS * K2
    
    real(double_precision) :: a, total_torque, corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: lindblad_torque_units, corotation_torque_units, total_torque_units
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    
    write(*,*) 'Evolution of the torque for a fixed planet mass "m"'

    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/study_torques_fixed_m.dat')
    open(11, file='unitary_tests/test_ref_torque_fixed_m.dat')
    open(12, file='unitary_tests/study_torques_fixed_m_units.dat')
    
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
      
      !------------------------------------------------------------------------------
      ! Calculation of the acceleration due to migration
      select case(TORQUE_TYPE)
        case('real') ! The normal torque profile, calculated form properties of the disk
          call get_corotation_torque(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        ! for retrocompatibility, 'mass_independant' has been added and refer to the old way of defining a mass-indep convergence zone
        case('linear_indep', 'mass_independant') ! a defined torque profile to get a mass independant convergence zone
          call get_corotation_torque_linear_indep(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        case('arctan_indep') ! a defined torque profile to get a mass independant convergence zone
          call get_corotation_torque_arctan_indep(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        case('mass_dependant')
          call get_corotation_torque_mass_dep_CZ(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
        
        case('manual')
          call get_corotation_torque_manual(stellar_mass, mass, p_prop, corotation_torque, lindblad_torque, torque_ref)
          
        case default
          write(*,*) 'Warning: The torque rule cannot be found.'
          write(*,*) 'Given value :', TORQUE_TYPE
      end select
      
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

    write(10,*) "plot 'study_torques_fixed_m.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(10,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(10,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"
    
    write(11,*) "plot 'test_ref_torque_fixed_m.dat' using 1:2 with lines"
        
    write(12,*) "plot 'study_torques_fixed_m_units.dat' using 1:2 with lines title '{/Symbol G}_c',\"
    write(12,*) "                             '' using 1:3 with lines title '{/Symbol G}_L',\"
    write(12,*) "                             '' using 1:4 with lines title '{/Symbol G}_{tot}'"

    
    do j=10,12
      write(j,*) "#pause -1 # wait until a carriage return is hit"
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) "set output 'torques_fixed_m.pdf'"
    write(11,*) "set output 'ref_torque_fixed_m.pdf'"
    write(12,*) "set output 'torques_fixed_m_units.pdf'"
    
    do j=10,12
      write(j,*) "replot # to generate the output file"
    end do
    
    close(10)
    close(11)
    close(12)

    
  end subroutine study_torques_fixed_m
  

  subroutine study_torques(stellar_mass)
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
    integer :: deltai, deltaj ! separation between to outputs for the vector part of the display (we do not want a lot of points)
    real(double_precision) :: vector_limit ! The size limit of the vectors, to prevent chevauching.
    
    write(*,*) 'Evolution of the total, lindblad and corotation torques depending on the planet mass and distance'
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! We want to have only 10 points both in x and y for vector outputs
    deltai = nb_points / 20
    deltaj = nb_mass / 20
    
    vector_limit = a_step * (deltai - 1)
    
    ! We open the file where we want to write the outputs
    open(10, file='unitary_tests/test_corotation_torque.dat')
    open(11, file='unitary_tests/test_total_torque.dat')
    open(12, file='unitary_tests/test_total_torque_units.dat')
    open(13, file='unitary_tests/test_lindblad_torque.dat')
    open(14, file='unitary_tests/test_ref_torque.dat')
    open(15, file='unitary_tests/test_vector_total_torque.dat')
    
    write(10,*) '# semi major axis (AU) ; mass in earth mass ; corotation torque (no dim)'
    write(11,*) '# semi major axis (AU) ; mass in earth mass ; total torque (no dim)'
    write(12,*) '# semi major axis (AU) ; mass in earth mass ; total torque in M_s.AU^2.day^{-2}'
    write(13,*) '# semi major axis (AU) ; mass in earth mass ; lindblad torque (no dim)'
    write(14,*) '# semi major axis (AU) ; mass in earth mass ; reference torque in M_s.AU^2.day^{-2}'
    write(15,*) '# semi major axis (AU) ; mass in earth mass ; Delta x ; Delta y'
    
    
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
         
!~         write(*,*) a(i),p_prop%semi_major_axis 
!~         call print_planet_properties(p_prop) 
!~         stop
        ! If the semi major axis is not well determined, we display a warning and give the values
        if (abs(a(i)-p_prop%semi_major_axis)/a(i).gt.1e-2) then
          write(*,*) 'Warning: for manual a=',a(i), 'we get :'
          call print_planet_properties(p_prop) 
        end if
         
        !------------------------------------------------------------------------------
        ! Calculation of the acceleration due to migration
        select case(TORQUE_TYPE)
          case('real') ! The normal torque profile, calculated form properties of the disk
            call get_corotation_torque(stellar_mass, mass(j), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
          
          ! for retrocompatibility, 'mass_independant' has been added and refer to the old way of defining a mass-indep convergence zone
          case('linear_indep', 'mass_independant') ! a defined torque profile to get a mass independant convergence zone
            call get_corotation_torque_linear_indep(stellar_mass, mass(j), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
          
          case('arctan_indep') ! a defined torque profile to get a mass independant convergence zone
            call get_corotation_torque_arctan_indep(stellar_mass, mass(j), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
          
          case('mass_dependant')
            call get_corotation_torque_mass_dep_CZ(stellar_mass, mass(j), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
          
          case('manual')
            call get_corotation_torque_manual(stellar_mass, mass(j), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
            
          case default
            write(*,*) 'Warning: The torque rule cannot be found.'
            write(*,*) 'Given value :', TORQUE_TYPE
        end select
        
        total_torque(i,j) = lindblad_torque + corotation_torque
        total_torque_units(i,j) = torque_ref * total_torque(i,j)
        
                
        write(10,*) a(i), mass(j) / (EARTH_MASS*K2), corotation_torque
        write(11,*) a(i), mass(j) / (EARTH_MASS*K2), total_torque(i,j)
        write(12,*) a(i), mass(j) / (EARTH_MASS*K2), total_torque_units(i,j)
        write(13,*) a(i), mass(j) / (EARTH_MASS*K2), lindblad_torque
        write(14,*) a(i), mass(j) / (EARTH_MASS*K2), torque_ref
        
        if ((modulo(i-int(deltai/2.),deltai).eq.0).and.(modulo(j,deltaj).eq.0)) then
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
      write(j,*) 'set grid xtics ytics linetype 0'
      write(j,*) 'set xrange [', a_min, ':', a_max, ']'
      write(j,*) 'set yrange [', mass_min / EARTH_MASS, ':', mass_max / EARTH_MASS, ']'
    end do

    write(10,*) "splot 'test_corotation_torque.dat' with pm3d notitle"
    
    write(11,*) "splot 'test_total_torque.dat' with pm3d notitle, \"
    write(11,*) "      'contour_total_torque.dat' with line linetype -1 linewidth 1 notitle, \"
    write(11,*) "      'test_vector_total_torque.dat' with vector notitle head filled linestyle -1"
    
    write(12,*) "splot 'test_total_torque_units.dat' with pm3d notitle, \"
    write(12,*) "      'contour_total_torque.dat' with line linetype -1 linewidth 1 notitle, \"
    write(12,*) "      'test_vector_total_torque.dat' with vector notitle head filled linestyle -1"
    
    write(13,*) "splot 'test_lindblad_torque.dat' with pm3d notitle"
    
    write(14,*) "splot 'test_ref_torque.dat' with pm3d notitle"

    
    do j=10,14
      write(j,*) "#pause -1 # wait until a carriage return is hit"
      write(j,*) "set terminal pngcairo crop enhanced size 1200, 1000"
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
    
  end subroutine study_torques
  
  subroutine study_dissipation_of_the_disk(stellar_mass)
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
    real(double_precision), parameter :: t_max = 1.d6 ! time in years
    real(double_precision), dimension(:), allocatable :: time, time_temp ! time in days
    integer :: time_size ! the size of the array 'time'. 
    
    ! Array to store values of the torque. Used to find the zero torque line
    real(double_precision), dimension(nb_points, nb_mass) :: total_torque
    
    real(double_precision) :: corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    real(double_precision) :: temp_min, temp_max, density_min, density_max, torque_min, torque_max
    character(len=80) :: filename_torque, filename_density, filename_temperature, filename_contour
    character(len=80) :: filename_density_ref, filename_temperature_ref
    character(len=80) :: output_torque, output_density, output_temperature, output_time, time_format, purcent_format
    integer :: time_length ! the length of the displayed time, usefull for a nice display
    
    integer :: i,j,k ! for loops
    integer :: nb_time ! The total number of 't' values. 
    integer :: error ! to retrieve error, especially during allocations
    !------------------------------------------------------------------------------
    write(*,*) 'Evolution of the total torque during the dissipation of the disk'
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    write(filename_density_ref, '(a,i0.5,a)') 'dissipation/surface_density',0,'.dat'
    write(filename_temperature_ref, '(a,i0.5,a)') 'dissipation/temperature',0,'.dat'
    
    ! we store in a .dat file the temperature profile
    call store_temperature_profile(filename=filename_temperature_ref)
    call store_density_profile(filename=filename_density_ref)
    

    temp_max = temperature_profile(1)
    temp_min = 0.
    
    ! We want the extremum of the surface density during the dissipation of the disk in order to have nice plots
    density_min = 0.
    density_max = maxval(surface_density_profile(1:NB_SAMPLE_PROFILES)) * MSUN / AU**2

    
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
    write(purcent_format, *) '(i',time_length,'"/",i',time_length,'," years")'
    
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
      
      dissipation_timestep = 0.5d0 * X_SAMPLE_STEP**2 / (4 * get_viscosity(1.d0)) ! a correction factor of 0.5 has been applied. No physical reason to that, just intuition and safety
      ! TODO if the viscosity is not constant anymore, the formulae for the dissipation timestep must be changed
      
      ! we get the density profile.
      select case(DISSIPATION_TYPE)
        case(1)
          call dissipate_density_profile() ! global parameter 'dissipation_timestep' must exist !
        
        case(2)
          call exponential_decay_density_profile()
          
        case default
          write(*,*) 'Warning: The dissipation rule cannot be found.'
          write(*,*) 'Given value :',DISSIPATION_TYPE
      end select
      
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
      
      write(*,purcent_format) int(time(k)/365.25d0), int(t_max)
      
      
      time(k+1) = time(k) + dissipation_timestep * 365.25d0 ! days
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
    write(12,*) "set terminal pngcairo enhanced size 1024, 768"
    write(12,*) 'set xlabel "semi major axis (AU)"'
    write(12,*) 'set ylabel "Temperature (K)"'
    write(12,*) 'set grid'
    write(12,*) 'set xrange [', a_min, ':', a_max, ']'
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
    write(13,*) "set terminal pngcairo enhanced size 1024, 768"
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
  
  end subroutine study_dissipation_of_the_disk

end program test_disk
