! Program that test different functions implemented in the module user_module.
program torque_diagram
  use types_numeriques
  use disk
  use disk_properties
  use iso_fortran_env, only : error_unit
    
  implicit none
  
  
  ! global parameters
  integer :: nb_distance = 150
  real(double_precision) :: a_min
  real(double_precision) :: a_max
  real(double_precision), dimension(2) :: a_transition = (/0.5d0, 10.d0/)
  integer :: nb_mass = 100
  real(double_precision) :: m_min_em = 0.1d0 ! in earth mass
  real(double_precision) :: m_max_em = 60.d0 ! in earth mass
  real(double_precision) :: mass_min ! in solar mass
  real(double_precision) :: mass_max ! in solar mass
  real(double_precision) :: torque_min = -5.d0 ! Torque boundaries for the gnuplot display
  real(double_precision) :: torque_max =  5.d0 ! Torque boundaries for the gnuplot display
  !--------------------------------
  
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

    inquire(file='torque.in', exist=isDefined)
    ! This file is used to say if we continue an integration or not. For test_disk, 
    ! we DO NOT continue an integration, so to avoid it, we delete the file if it exist.
    if (isDefined) then
      call read_torquein()
    else
      a_min = INNER_BOUNDARY_RADIUS
      a_max = OUTER_BOUNDARY_RADIUS
    end if
    ! We write the torque.in file even if it existed before, to make sure that new paremeters will appear with their default values in there
    call write_torquein()
    
    mass_min = m_min_em * EARTH_MASS
    mass_max = m_max_em * EARTH_MASS
    
    ! We want to show the torque profile. It is important to check which value has been declared in 'TORQUE_TYPE'
    call study_torques(stellar_mass)
    
    
  end subroutine init_and_test
  
  subroutine read_torquein()
  
  implicit none
  
  character(len=80) :: line
  character(len=1), parameter :: comment_character = '!' ! character that will indicate that the rest of the line is a comment
  integer :: comment_position ! the index of the comment character on the line. If zero, there is none on the current string
  integer :: error ! to store the state of a read instruction
  
  logical :: isParameter, isDefined
  character(len=80) :: identificator, value
  !------------------------------------------------------------------------------
  
  inquire(file='torque.in', exist=isDefined)
  if (isDefined) then
  
    open(10, file='torque.in', status='old')
    
    do
      read(10, '(a80)', iostat=error) line
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
        case('a_min')
          read(value, *) a_min
        
        case('a_max')
          read(value, *) a_max
        
        case('a_transition')
          read(value, *) a_transition
          
        case('nb_distance')
          read(value, *) nb_distance
          
        case('m_min')
          read(value, *) m_min_em
        
        case('m_max')
          read(value, *) m_max_em
          
        case('nb_mass')
          read(value, *) nb_mass
        
        case('torque_min')
          read(value, *) torque_min
        
        case('torque_max')
          read(value, *) torque_max
        
        case default
          write(*,*) 'Warning: An unknown parameter has been found'
          write(*,*) "identificator='", trim(identificator), "' ; value(s)='", trim(value),"'"
        end select
      end if
    end do
    close(10)
    
    if (a_transition(1).gt.a_transition(2)) then
      write(error_unit,*) 'Error: a_transition radii must be ordered which is not the case.'
      write(error_unit,*) 'a_transition = ',a_transition
      call exit(2)
    end if
    
  else
    write (*,*) 'Warning: The file "torque.in" does not exist. Default values have been used'
  end if
  
  end subroutine read_torquein
  
subroutine write_torquein()
! subroutine that write the parameters of the user_module into the file 'disk.out'

! Global Parameters
! a_min
! a_max
! nb_distance
! m_min
! m_max
! nb_mass
! torque_min
! torque_max

  implicit none
  
  open(10, file='torque.in')
  write(10,'(a)') '! Parameters to build the migration map'
  write(10,'(a)') '! SAMPLE TUNNING'
  write(10,'(a)') '! 1) Between a_min and 1st transition radius, we will have 10% of the total '
  write(10,'(a)') '!    number of points, log-spaced'
  write(10,'(a)') '! 2) Between 1st and 2nd transition radii, we will have the expected number '
  write(10,'(a)') '!    of points in this range if everything was log-spaced, except this will '
  write(10,'(a)') '!    be linearly spaced.'
  write(10,'(a)') '! 3) Between 2nd transition radius and a_max, remaining points are log-spaced'
  write(10,'(a,2(f5.1),a)') 'a_transition = ', a_transition, " ! (2 transition radii) in AU"
  write(10,'(a)') '! DISTANCE RANGE'
  write(10,'(a,f5.1,a)') 'a_min = ', a_min, " ! in AU"
  write(10,'(a,f5.1,a)') 'a_max = ', a_max, " ! in AU"
  write(10,'(a,i4)') 'nb_distance = ', nb_distance
  write(10,'(a)') '! MASS RANGE'
  write(10,'(a,f5.1,a)') 'm_min = ', m_min_em, " ! in earth mass"
  write(10,'(a,f5.1,a)') 'm_max = ', m_max_em, " ! in earth mass"
  write(10,'(a,i4)') 'nb_mass = ', nb_mass
  write(10,'(a)') '! TORQUE RANGE'
  write(10,'(a,f5.1,a)') 'torque_min = ', torque_min, " ! in units of Gamma_0"
  write(10,'(a,f5.1,a)') 'torque_max = ', torque_max, " ! in units of Gamma_0"




  write(10,*) ''
  close(10)
  
end subroutine write_torquein

  subroutine study_torques(stellar_mass)
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of Semi-major axis.
  
    use contour
    
    
    implicit none
    
    real(double_precision), intent(in) :: stellar_mass ! in [msun * K2]
    
    ! mass sampling
    real(double_precision) :: mass_step
    real(double_precision), dimension(:), allocatable :: mass
    

    ! step for log sampling
    real(double_precision) :: a_step ! the separation between to radius values (can change)
    real(double_precision) :: a_start ! the radius where the second part of the distance array start
    real(double_precision) :: a_stop ! the radius where the second part of the distance array stop
    real(double_precision), dimension(:), allocatable :: a
    
    real(double_precision), dimension(:, :), allocatable :: total_torque
    
    real(double_precision) :: corotation_torque, lindblad_torque, torque_ref
    real(double_precision) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
    real(double_precision) :: position(3), velocity(3)
    type(PlanetProperties) :: p_prop
    
    ! Small number, to be just after the last point and ensure to have a good rendering of the map (else, the last point in often missing)
    real(double_precision), parameter :: range_shift = 1.d-7
    
    integer :: i,j ! for loops
    integer :: idx, nb_tmp ! to fill the radius array step by step
    
    !------------------------------------------------------------------------------
    
    if (.not.allocated(a)) then
      allocate(a(nb_distance))
      allocate(mass(nb_mass))
      allocate(total_torque(nb_distance, nb_mass))
    end if
    
    
    mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    
    ! first regime
    nb_tmp = nb_distance / 10
    if (a_min.lt.a_transition(1)) then
      a_step = (a_transition(1) / a_min)**(1.d0 / (dfloat(nb_tmp)))
      
      do i=1,nb_tmp
        a(i) = a_min * a_step**(i - 1)
      end do
    end if
    idx = i
    
    ! second regime
    a_start = max(a_min, a_transition(1))
    if (a_max.lt.a_transition(2)) then
      nb_tmp = nb_distance - idx
      a_stop = a_max
    else
      ! We calculate the number of log-spaced points we would have in the second regime if all the remaining space was in log, 
      ! and we will space them linearly in the second regime.
      nb_tmp = int(dfloat((nb_distance - idx)) * (log(a_transition(2))-log(a_start)) / (log(a_max)-log(a_start)))
      a_stop = a_transition(2)
    end if
    
    a_step = (a_stop - a_start) / (dfloat(nb_tmp))
    
    do i=idx, idx + nb_tmp
      a(i) = a_start + a_step * (i-idx)
    end do
    idx = i
    
    ! third regime
    a_start = a(idx-1) + a_step
    a_step = (a_max / a_start)**(1.d0 / (dfloat(nb_distance - idx)))
    
    do i=idx, nb_distance
      a(i) = a_start * a_step**(i - idx)
    end do
    
    ! mass array
    do j=1, nb_mass
      mass(j) = (mass_min + mass_step * (j - 1.d0)) * K2 ! mass in [Msun * K2]
    end do
    
    !------------------------------------------------------------------------------
    write(*,*) 'Evolution of the total, lindblad and corotation torques depending on the planet mass and distance'
    
    position(1:3) = 0.d0
    velocity(1:3) = 0.d0
    
    ! We open the file where we want to write the outputs
    open(11, file='total_torque.dat')
    
    write(11,*) '# Semi-major axis (AU) ; mass in earth mass ; total torque (no dim)'
    
    do i=1, nb_distance ! loop on the position
!~       a(i) = a_min + a_step * (i-1)
!~       a(i) = a_min * a_step**(i - 1)
      
      ! We generate cartesian coordinate for the given Semi-major axis
      position(1) = a(i)
      
      
      do j=1,nb_mass
        ! We generate cartesian coordinate for the given mass and Semi-major axis
        velocity(2) = sqrt((stellar_mass + mass(j)) / position(1))
        
        ! we store in global parameters various properties of the planet
        call get_planet_properties(stellar_mass=stellar_mass, & ! Input
         mass=mass(j), position=position(1:3), velocity=velocity(1:3),& ! Input
         p_prop=p_prop) ! Output
         
        ! If the Semi-major axis is not well determined, we display a warning and give the values
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

    write(11,*) "set terminal wxt enhanced font ',20'"

    write(11,*) 'set xlabel "Semi-major axis (AU)"'
    write(11,*) 'set ylabel "Planet mass (m_{earth})" center'
    write(11,*) 'set title "Evolution of the total torque {/Symbol G}_{tot}/{/Symbol G}_0 "'
    write(11,*) 'set pm3d map'
    write(11,*) 'set pm3d explicit'
    write(11,*) 'set palette rgbformulae 22,13,-31'
    write(11,*) 'set grid xtics ytics linetype 0'
    write(11,*) 'set xrange [', a_min, ':', a_max + range_shift, ']'
    write(11,*) 'set yrange [', mass_min / EARTH_MASS, ':', mass_max / EARTH_MASS, ']'
    write(11,'(a,f5.1,a,f5.1,a)') ' set cbrange [',torque_min,':', torque_max,']'
    write(11,*) "splot 'total_torque.dat' with pm3d notitle, \"
    write(11,*) "      'contour_total_torque.dat' with line linetype -1 linewidth 1 notitle"
    write(11,*) "pause -1"
    write(11,*) ""
    
    ! Only the image
    write(11,*) "set terminal png crop enhanced size 1200, 1000"
    write(11,*) "set output 'total_torque_diagram.png'"
    write(11,*) "set pm3d map"
    write(11,*) "set pm3d explicit"
    write(11,*) "set palette rgbformulae 22,13,-31"
    write(11,*) "unset xtics"
    write(11,*) "unset ytics"
    write(11,*) 'unset xlabel'
    write(11,*) 'unset ylabel'
    write(11,*) 'unset title'
    write(11,*) 'unset border'
    write(11,*) 'unset colorbox'
    write(11,*) 'set xrange [', a_min, ':', a_max+range_shift, ']'
    write(11,*) 'set yrange [', mass_min / EARTH_MASS, ':', mass_max / EARTH_MASS, ']'
    write(11,'(a,f5.1,a,f5.1,a)') ' set cbrange [',torque_min,':', torque_max,']'
    write(11,*) "splot 'total_torque.dat' with pm3d notitle"
    write(11,*) ""

    
    ! The final .pdf file
    write(11,*) "set terminal pdfcairo enhanced font ' ,8'"
    write(11,*) "set output 'total_torque.pdf'"
    write(11,*) "set multiplot"
    write(11,*) "# To display the colorbox (without displaying any map)"
    write(11,*) "unset tics"
    write(11,*) "set cbtics"
    write(11,*) "set colorbox"
    write(11,*) "unset border"
    write(11,*) "set pm3d map"
    write(11,*) "set pm3d explicit"
    write(11,*) "set palette rgbformulae 22,13,-31"
    write(11,'(a,f5.1,a,f5.1,a)') ' set cbrange [',torque_min,':', torque_max,']'
    write(11,*) "set lmargin at screen 0.15"
    write(11,*) "set rmargin at screen 0.85"
    write(11,*) "set bmargin at screen 0.175"
    write(11,*) "set tmargin at screen 0.85"
    write(11,*) "unset surface"
    write(11,*) "splot 0 with pm3d notitle"
    write(11,*) ""
    write(11,*) "# We display the bitmap, that we include in the .pdf file"
    write(11,*) "set xrange [0:*]"
    write(11,*) "set yrange [0:*]"
    write(11,*) "set cbrange[*:*] # To have correct display of bitmap colors"
    write(11,*) "plot 'total_torque_diagram.png' binary filetype=png with rgbimage notitle"
    write(11,*) " "
    write(11,*) "set grid xtics ytics linetype 0"
    write(11,*) "set border"
    write(11,*) "set tics"
    write(11,*) 'set xlabel "Semi-major axis (AU)"'
    write(11,*) 'set ylabel "Planet mass (m_{earth})"'
    write(11,*) 'set title "Evolution of the total torque {/Symbol G}_{tot}/{/Symbol G}_0 "'
    write(11,*) 'set xrange [', a_min, ':', a_max+range_shift, ']'
    write(11,*) 'set yrange [', mass_min / EARTH_MASS, ':', mass_max / EARTH_MASS, ']'
    write(11,*) "set border"
    write(11,*) "plot 'contour_total_torque.dat' using 1:2 with line linetype -1 linewidth 1 notitle"
    write(11,*) "unset multiplot"
    
    
    close(11)
    
  end subroutine study_torques

end program torque_diagram
