! Program that test different functions implemented in the module user_module.
program test_mfo_user
  use types_numeriques
  use user_module
  use physical_constant
  
  !rajout pour getNumberOgBodies
!~   use mercury_globals
!~   use mercury_constant
!~   use mercury_outputs
  use utilities
    
  implicit none
  
  integer :: nbod, nbig
  
  call getNumberOfBodies(nb_big_bodies=nbig, nb_bodies=nbod)
  write (*,*) 'nbig=',nbig, 'nbod=',nbod
!~   call unitary_tests()
!~   call test_gamma_eff()
  

  
  contains

subroutine getNumberOfBodies(nb_big_bodies, nb_bodies)
! subroutine that return the number of bodies and the number of big bodies
!
! Return
! nb_bodies : the total number of bodies, including the central object
! nb_big_bodies : the number of big bodies

  implicit none

  integer, intent(out) :: nb_bodies, nb_big_bodies
  
  real(double_precision), dimension(9) :: dummy ! variable to read the value and have the right position in the file without storing the values
  
  ! Local
  integer :: j,k,lim(2,10),nsub, error
  logical test
  character(len=80) :: infile(3),filename,c80
  character(len=150) :: string


  do j = 1, 80
     filename(j:j) = ' '
  end do
  do j = 1, 3
     infile(j)   = filename
  end do  
  
  ! Read in filenames and check for duplicate filenames
  inquire (file='files.in', exist=test)
  if (.not.test) write(*,*) 'Error: the file "files.in" does not exist'
  open (15, file='files.in', status='old')
  
  ! Input files
  do j = 1, 3
     read (15,'(a150)') string
     call mio_spl (150,string,nsub,lim)
     infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
     do k = 1, j - 1
        if (infile(j).eq.infile(k)) write(*,*) 'Error: dans "files.in", certains fichiers sont identiques'
     end do
  end do
  close (15)
  
  !--------------------------------------------------------

  !  READ  IN  DATA  FOR  BIG  AND  SMALL  BODIES
  
  nb_bodies = 1
  do j = 1, 2
     if (j.eq.2) nb_big_bodies = nb_bodies
     
     ! Check if the file containing data for Big bodies exists, and open it
     filename = infile(j)
     inquire (file=filename, exist=test)
     if (.not.test) write (*,*) filename, 'does not exist'
     open (11, file=filename, status='old', iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(filename)
        stop
     end if
     
     ! Read data style
     do
      read (11,'(a150)') string
      if (string(1:1).ne.')') exit
    end do
     
     ! Read epoch of Big bodies
     if (j.eq.1) then
      do
        read (11,'(a150)') string
        if (string(1:1).ne.')') exit
      end do
     end if
     
     ! Read information for each object
     do
      read (11,'(a)',iostat=error) string
      if (error /= 0) exit
      
      if (string(1:1).eq.')') cycle
    
     
     nb_bodies = nb_bodies + 1
     
     
    ! we skip the line(s) that contains informations of the current planet
     do
       read (11,'(a150)') string
       if (string(1:1).ne.')') exit
     end do 
     backspace(11)
     read (11,*) dummy
     
     
     end do
      close (11)
  end do

end subroutine getNumberOfBodies
  
  subroutine test_gamma_eff()
  
  ! subroutine that test the function 'get_opacity'
  
  ! Return:
  !  a data file 'test_opacity.dat' 
  ! and an associated gnuplot file 'opacity.gnuplot' that display values for get_opacity for a range of p values.
    implicit none
    
    real(double_precision) :: adiabatic_index = 5.d0/3.d0
    real(double_precision), dimension(3) :: h = (/ 0.025, 0.05, 0.1/)
    real(double_precision), dimension(3) :: gamma_eff
    real(double_precision) :: alpha
    
    real(double_precision), parameter :: alpha_min = 0.0004
    real(double_precision), parameter :: alpha_max = 80.
    integer, parameter :: nb_points = 400
    real(double_precision), parameter :: alpha_step = (alpha_max/alpha_min) ** (1/(nb_points-1.d0))
    
    integer :: i,j ! for loops
    real(double_precision) :: Q_p
    
    ! We open the file where we want to write the outputs
    open(10, file='test_gamma_eff.dat')
    write(10,*) "Correspond to the figure 8 of Bell & Lin 1994"
    write(10,*) 'alpha_chi ; gamma_eff for aspect ratio equal to 0.025, 0.05 and 0.1'
    
    do i=1, nb_points
      alpha = alpha_min * alpha_step ** (i-1)
      
      do j=1,3
        
        !------------------------------------------------------------------------------
        ! Q is needed by the lindblad torque. We set Q for m ~ 2 /3 h (45): 
!~         Q_p = TWOTHIRD * chi_p / (aspect_ratio_p**3 * radius_p**2 * omega_p)
        Q_p = 2.d0 * alpha / (3.d0 * h(j))
        !------------------------------------------------------------------------------
        
        gamma_eff(j) = 2.d0 * Q_p * adiabatic_index / (adiabatic_index * Q_p + 0.5d0 * &
        sqrt(2.d0 * sqrt((adiabatic_index * adiabatic_index * Q_p * Q_p + 1.d0)**2 - 16.d0 * Q_p * Q_p * (adiabatic_index - 1.d0)) &
        + 2.d0 * adiabatic_index * adiabatic_index * Q_p * Q_p - 2.d0))
          !gamma_eff_t = 2. * Q * gamma / (gamma * Q + 0.5 * sqrt(2. * sqrt((gamma**2 * Q**2 + 1.)**2 -16. * Q**2 * (gamma - 1.)) + 2. * gamma**2 * Q**2 - 2.))

      end do
      
      write(10,*) alpha, gamma_eff(1), gamma_eff(2), gamma_eff(3)
    end do
    close(10)
    write(*,*) "alpha_step=",alpha_step
    
    open(10, file="gamma_eff.gnuplot")
    write(10,*) 'set terminal x11 enhanced'
    write(10,*) 'set xlabel "{/Symbol a}_{/Symbol c}={/Symbol c}/(c_s H)"'
    write(10,*) 'set ylabel "{/Symbol g}_{eff}"'
    write(10,*) 'set logscale x'
    write(10,*) 'set grid'
    write(10,*) 'set xrange [', alpha_min, ':', alpha_max, ']'
    write(10,*) "plot 'test_gamma_eff.dat' using 1:2 with lines title 'h=0.025',\"
    write(10,*) "     '' using 1:3 with lines title 'h=0.05',\"
    write(10,*) "     '' using 1:4 with lines title 'h=0.1'"
    write(10,*) "pause -1 # wait until a carriage return is hit"
    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) "set output 'gamma_eff.pdf'"
    write(10,*) "replot # pour générer le fichier d'output"  
    
    close(10)
  

  end subroutine test_gamma_eff
end program test_mfo_user
