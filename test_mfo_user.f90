! Program that test different functions implemented in the module user_module.
program test_mfo_user
  use types_numeriques
  use user_module
  use physical_constant
  
  implicit none
  
  call unitary_tests()
!~   call test_gamma_eff()
  

  
  contains
  
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
