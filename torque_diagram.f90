! Program that test different functions implemented in the module user_module.
program torque_diagram
  use types_numeriques
  use disk
  use disk_properties
  use iso_fortran_env, only : error_unit
    
  implicit none
  
  call init_and_test()
  

  
  contains

! ##################################
! TESTS OF THE ROUTINES
! ##################################
  subroutine init_and_test()
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
    
    logical :: isDefined

    stellar_mass = 1.d0 * K2 ! [Msun * K2]
    
    
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
    
    
  end subroutine init_and_test

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
    
    real(double_precision), dimension(nb_points, nb_mass) :: total_torque
    
    real(double_precision) :: corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    integer :: i,j ! for loops
    
    !------------------------------------------------------------------------------
    a_min = INNER_BOUNDARY_RADIUS
    a_max = OUTER_BOUNDARY_RADIUS
!~     a_step = (a_max - a_min) / (nb_points-1.d0)
    a_step = (a_max / a_min)**(1 / (nb_points - 1.d0))
    
    !------------------------------------------------------------------------------
    write(*,*) 'Evolution of the total, lindblad and corotation torques depending on the planet mass and distance'
    
    position(1:3) = 0.d0
    velocity(1:3) = 0.d0
    
    ! We open the file where we want to write the outputs
    open(11, file='total_torque.dat')
    
    write(11,*) '# semi major axis (AU) ; mass in earth mass ; total torque (no dim)'
    
    do i=1, nb_points ! loop on the position
!~       a(i) = a_min + a_step * (i-1)
      a(i) = a_min * a_step**(i - 1)
      
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
                
        write(11,*) a(i), mass(j) / (EARTH_MASS*K2), total_torque(i,j)
        
      end do
      
      write(11,*) ""! we write a blank line to separate them in the data file, else, gnuplot doesn't want to make the surface plot
    end do
    close(11)
    
    
    ! We want to get the contour for the torque equal to 0 for total torque both in physical dimension of units of gamma_0
    mass(1:nb_mass) = mass(1:nb_mass) / (EARTH_MASS*K2)
    call get_contour(total_torque, a, mass,'contour_total_torque.dat', 0.d0) ! This line opens a file (unit=10), thus this unit must not be used at that time
    
    open(11, file="total_torque.gnuplot")

    write(11,*) "set terminal wxt enhanced"

    write(11,*) 'set xlabel "semi major axis (AU)"'
    write(11,*) 'set ylabel "Planet mass (m_{earth})" center'
    
    write(11,*) 'set title "Evolution of the total torque {/Symbol G}_{tot}/{/Symbol G}_0 "'
    
    write(11,*) 'set pm3d map'
    write(11,*) 'set pm3d explicit'
    write(11,*) 'set palette rgbformulae 22,13,-31'
    write(11,*) 'set grid xtics ytics linetype 0'
    write(11,*) 'set xrange [', a_min, ':', a_max, ']'
    write(11,*) 'set yrange [', mass_min / EARTH_MASS, ':', mass_max / EARTH_MASS, ']'

    
    write(11,*) "splot 'total_torque.dat' with pm3d notitle, \"
    write(11,*) "      'contour_total_torque.dat' with line linetype -1 linewidth 1 notitle"
    write(11,*) "pause -1"
    
    
    close(11)
    
  end subroutine study_torques

end program torque_diagram
