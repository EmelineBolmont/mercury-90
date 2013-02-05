module user_module

!*************************************************************
!** Modules that contains user defined modules. 
!** Only mfo_user will be public.
!**
!** Version 1.0 - june 2011
!*************************************************************
  use types_numeriques
  use physical_constant

  implicit none
  
  private
  
  public :: mfo_user, unitary_tests ! unitary_tests must be used only to test function, you should not use it out of the debug phase
  
  real(double_precision), parameter :: b_over_h = 0.4 ! the smoothing length for the planet's potential
  real(double_precision), parameter :: adiabatic_index = 1.4 ! the adiabatic index for the gas equation of state
  real(double_precision), parameter :: mean_molecular_weight = 2.35 ! the mean molecular weight in mass of a proton
  
  ! Here we define the power law for surface density sigma(R) = sigma_0 * R^(-sigma_index)
  real(double_precision), parameter :: sigma_0 = 1700 ! the surface density at (R=1AU) [g/cm^2]
  real(double_precision), parameter :: sigma_index = 0.5! the negative slope of the surface density power law (alpha in the paper)
  real(double_precision), parameter :: sigma_0_num = sigma_0 * AU**2 / MSUN ! the surface density at (R=1AU) [Msun/AU^2]
  
  ! Here we define the power law for viscosity viscosity(R) = alpha * c_s * H
  real(double_precision), parameter :: alpha = 0.005 ! alpha parameter for an alpha prescription of the viscosity [No dim]
  
  ! Here we define the power law for temperature T(R) = temperature_0 * R^(-temperature_index)
  real(double_precision), parameter :: temperature_0 = 150. ! the temperature at (R=1AU) [K]
  real(double_precision), parameter :: temperature_index = 1.6! the negativeslope of the temperature power law (beta in the paper)
  
  !prefactors
  real(double_precision) :: x_s_prefactor
  real(double_precision) :: chi_p_prefactor
  real(double_precision) :: scaleheight_prefactor
  
  real(double_precision) :: torque_lindblad ! prefactor for the lindblad torque
  real(double_precision) :: torque_hs_baro ! barotropic part of the horseshoe drag
  real(double_precision) :: torque_c_lin_baro ! barotropic part of the linear corotation torque
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MFO_USER.FOR    (ErikSoft   2 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Applies an arbitrary force, defined by the user.

! If using with the symplectic algorithm MAL_MVS, the force should be
! small compared with the force from the central object.
! If using with the conservative Bulirsch-Stoer algorithm MAL_BS2, the
! force should not be a function of the velocities.

! Code Units are in AU, days and solar mass.

! N.B. All coordinates and velocities must be with respect to central body
! ===
!------------------------------------------------------------------------------

subroutine mfo_user (time,jcen,n_bodies,n_big_bodies,mass,position,velocity,acceleration)
!  mass      = mass (in solar masses)
!  position     = coordinates (x,y,z) with respect to the central body [AU]
!  velocity     = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  n_bodies  = current number of bodies (INCLUDING the central object)
!  n_big_bodies  =    "       "    " big bodies (ones that perturb everything else)
!  time  = current epoch [days]

  use physical_constant
  use mercury_constant  

  implicit none

  
  ! Input
  integer, intent(in) :: n_bodies,n_big_bodies
  real(double_precision),intent(in) :: time,jcen(3),mass(n_bodies),position(3,n_bodies),velocity(3,n_bodies)
  
  ! Output
  real(double_precision),intent(out) :: acceleration(3,n_bodies)
  
  !------------------------------------------------------------------------------ 
  !------Local-------
  
  ! loop integers
  integer :: planet
  
  
  real(double_precision) :: torque ! the torque exerted by the disk on the planet [Ms.AU^2]
  real(double_precision) :: corotation_torque ! the corotation torque exerted by the disk on the planet in unit of torque_ref [No dim]
  real(double_precision) :: torque_ref ! a ref torque that depends on the properties of the planet [Ms.AU^2]
  real(double_precision) :: gamma_eff ! effective adiabatic coefficient [No dim]
  real(double_precision) :: time_mig ! The migration timescale for the planet [day]
  real(double_precision) :: angular_momentum ! the angular momentum of the planet [Ms.AU^2.day^-1]

  !------------------------------------------------------------------------------
  ! Setup
  acceleration(:,:) = 0.d0
  !------------------------------------------------------------------------------
  call init_globals(mass(1))
  
  
  do planet=2,n_big_bodies
    call get_corotation_torque(mass(1), mass(planet), position(3,planet), velocity(3,planet),&
         corotation_torque, torque_ref, gamma_eff)
    
    torque = (torque_ref / gamma_eff) * (torque_lindblad + corotation_torque)
      ! We calculate the angular momentum, that is, the z-composant. We assume that the x and y composant are negligible (to be tested)
    angular_momentum = mass(planet) * (position(1,planet) * velocity(2,planet) - position(2,planet) * velocity(1,planet))
  !~   angular_momentum_x = mass(planet) * (position(2,planet) * velocity(3,planet) - position(3,planet) * velocity(2,planet))
  !~   angular_momentum_y = mass(planet) * (position(3,planet) * velocity(1,planet) - position(1,planet) * velocity(3,planet))
    
    time_mig = 0.5d0 * angular_momentum / torque
    
    acceleration(1,planet) = - velocity(1,planet) / time_mig
    acceleration(2,planet) = - velocity(2,planet) / time_mig
    acceleration(3,planet) = - velocity(3,planet) / time_mig
    
  end do
  
  
  !------------------------------------------------------------------------------
  

  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_user

subroutine get_corotation_torque(stellar_mass, mass, position, velocity, corotation_torque, Gamma_0, gamma_eff)
! function that return the total torque exerted by the disk on the planet 
!
! Parameter:
! mass : an array will all the masses (including the central one as the first element [solar mass]
! radius_p : the distance between the central body and each planet
  use orbital_elements, only : mco_x2a
  implicit none
  real(double_precision), intent(in) :: stellar_mass, mass, position(3), velocity(3)
  
  real(double_precision), intent(out) :: corotation_torque
  real(double_precision), intent(out) :: Gamma_0 ! Normalization factor for all the torques [Ms.AU^2](equation (8) of Paardekooper, Baruteau, 2009)
  real(double_precision), intent(out) :: gamma_eff ! effective adiabatic index depending on several parameters [no dim]

  !Local
  
  ! Meaningless parameters (intermediate constant and so on that doesn't mean anything physically)
  real(double_precision) :: Q_p ! parameter for gamma_eff (equation (45) of Paardekooper, Baruteau, 2010 II. Effects of diffusion)
  real(double_precision) :: k_p ! parameter for p_nu and p_chi for example  !!! This is not 'k' from equation (15)!
  real(double_precision) :: gm ! parameter to calculate the semi major axis. passed to mco_x2a
  
  ! Properties of the planet
  real(double_precision) :: radius_p ! the radial position of the planet [AU]
  real(double_precision) :: velocity_p, vel_squared ! the norm of the speed [AU/day] and the velocity squared.
  real(double_precision) :: omega_p ! the angular rotation [day-1]
  real(double_precision) :: semi_major_axis ! semi major axis of the planet [AU]
  
  !Properties of the disk at the location of the planet

  real(double_precision) :: sigma_p ! the surface density of the gas disk at the planet location [MSUN.AU^-2]
  real(double_precision) :: bulk_density_p ! the bulk_density of the disk at the location of the planet [MSUN.AU^-3]
  real(double_precision) :: scaleheight_p ! the scaleheight of the disk at the location of the planet [AU]
  real(double_precision) :: aspect_ratio_p ! the aspect_ratio of the gas disk at the location of the planet [no dim]
  real(double_precision) :: x_s ! semi-width of the horseshoe region [radius_p (in unity of position of the planet)]
  real(double_precision) :: zeta_eff ! effective entropy index depending on gamma_eff [no dim]
  real(double_precision) :: chi_p ! the thermal diffusion coefficient at the location of the planet [AU^2.day^-1]
  real(double_precision) :: nu_p ! the viscosity of the disk at the location of the planet [AU^2.day^-1]
  real(double_precision) :: opacity_p ! the opacity of the disk at the location of the planet [?]
  real(double_precision) :: temperature_p ! the temperature of the disk at the location of the planet [K]
  real(double_precision) :: p_nu ! parameter for saturation due to viscosity at the location of the planet [no dim]
  real(double_precision) :: p_chi ! parameter for saturation due to thermal diffusion at the location of the planet [no dim]
  
  !Torques (some depends of the planet)
  real(double_precision) :: torque_hs_ent ! entropy related part of the horseshoe drag
  real(double_precision) :: torque_c_lin_ent ! entropy related part of the linear corotation torque
  
  
  ! WE CALCULATE PROPERTIES OF THE PLANETS
  gm = K2 * (stellar_mass + mass)
  ! We get semi_major_axis, radius_p and vel_squared
  call mco_x2a (gm,position(1), position(2), position(3), velocity(1),velocity(2), velocity(3),semi_major_axis,radius_p,vel_squared)
  
  
  !------------------------------------------------------------------------------
  sigma_p = sigma_0_num / radius_p ** sigma_index ! [N.U.]
  temperature_p = temperature_0 / radius_p**temperature_index ! [K]
  
  velocity_p = sqrt(vel_squared)
  omega_p = sqrt(gm / (semi_major_axis * semi_major_axis * semi_major_axis)) ! [day-1]
  
  !------------------------------------------------------------------------------
  ! H = sqrt(k_B * T / (omega^2 * mu * m_H))
  scaleheight_p = scaleheight_prefactor * sqrt(temperature_p) / omega_p
  !------------------------------------------------------------------------------
  nu_p = alpha * omega_p * scaleheight_p**2 ! [AU^2.day-1]
  aspect_ratio_p = scaleheight_p / radius_p
    write(*,*) nu_p * AU**2 / DAY 
  !------------------------------------------------------------------------------
  bulk_density_p = 0.5d0 * sigma_p / scaleheight_p
  !------------------------------------------------------------------------------
  opacity_p = get_opacity(temperature_p, bulk_density_p)
  
  !------------------------------------------------------------------------------
  ! WE CALCULATE TOTAL TORQUE EXERTED BY THE DISK ON THE PLANET
  Gamma_0 = (mass / (stellar_mass * aspect_ratio_p))**2 * sigma_p * radius_p**4 * omega_p**2
  
  chi_p = chi_p_prefactor * temperature_p**4 / (opacity_p * sigma_p**2 * omega_p**2)
  
  !------------------------------------------------------------------------------
  ! Q is needed by the lindblad torque. We set Q for m ~ 2 /3 h (45): 
  Q_p = TWOTHIRD * chi_p / (aspect_ratio_p * scaleheight_p**2 * omega_p) ! aspect_ratio_p**3 * radius_p**2 = aspect_ratio * scaleheight**2
  !------------------------------------------------------------------------------
  
  gamma_eff = 2.d0 * Q_p * adiabatic_index / (adiabatic_index * Q_p + 0.5d0 * &
  sqrt(2.d0 * sqrt((adiabatic_index * adiabatic_index * Q_p * Q_p + 1.d0)**2 - 16.d0 * Q_p * Q_p * (adiabatic_index - 1.d0)) &
  + 2.d0 * adiabatic_index * adiabatic_index * Q_p * Q_p - 2.d0))
  
  !------------------------------------------------------------------------------
  zeta_eff = temperature_index - (gamma_eff - 1.d0) * sigma_index
  
  x_s = x_s_prefactor / gamma_eff**0.25d0 * sqrt(mass / aspect_ratio_p)
  
  !------------------------------------------------------------------------------
  ! k_p is defined to limit the number of operation and to have a value independant from chi_p or nu_p
  k_p = radius_p * radius_p * omega_p * x_s * x_s * x_s / (2.d0 * PI)
  
  !------------------------------------------------------------------------------
  p_nu = TWOTHIRD * sqrt(k_p / nu_p)
  
  p_chi = sqrt(k_p / chi_p)
  
  torque_hs_ent = 7.9d0 * zeta_eff / gamma_eff
  torque_c_lin_ent = (2.2d0 - 1.4d0 / gamma_eff) * zeta_eff
  
  !------------------------------------------------------------------------------

  corotation_torque =  torque_hs_baro * get_F(p_nu) * get_G(p_nu) + torque_c_lin_baro * (1 - get_K(p_nu)) &
    + torque_hs_ent * get_F(p_nu) * get_F(p_chi) * sqrt(get_G(p_nu) * get_G(p_chi)) &
    + torque_c_lin_ent * sqrt((1 - get_K(p_nu)) * (1 - get_K(p_chi)))
      
!~     get_corotation_torque = (Gamma_0 / gamma_eff) * (torque_lindblad &
!~     + torque_hs_baro * get_F(p_nu) * get_G(p_nu) + torque_c_lin_baro * (1 - get_K(p_nu)) &
!~     + torque_hs_ent * get_F(p_nu) * get_F(p_chi) * sqrt(get_G(p_nu) * get_G(p_chi)) &
!~     + torque_c_lin_ent * sqrt((1 - get_K(p_nu)) * (1 - get_K(p_chi))))

  return
end subroutine get_corotation_torque

subroutine init_globals(stellar_mass)
! subroutine that initialize global values that define prefactors or values for torque that does not depend on the planet properties
! Parameters
! stellar_mass : the mass of the central object in solar mass
  
  implicit none
  real(double_precision), intent(in) :: stellar_mass
  
  x_s_prefactor = 1.1d0 * (b_over_h / 0.4d0)**0.25d0 / sqrt(stellar_mass) ! mass(1) is here for the ratio of mass q

  chi_p_prefactor = (16.d0  / 3.d0) * adiabatic_index * (adiabatic_index - 1.d0) * SIGMA_STEFAN
  
  ! AU is in cm, so we must turn into meter before doing the conversion
  ! division of k_B by m_H is done separately for exponant and value to have more precision
  ! sqrt(k_B/m_H) in numerical units, knowing that [k_B]=[m^2.kg.s^-2K^-1] and [m_H]=[kg]. 
  scaleheight_prefactor = sqrt(1.3806503d0/(1.67262158d0 * mean_molecular_weight) * 1.d4) * DAY / (AU * 1.d-2) 
  
  torque_lindblad = -2.5d0 - 1.7d0 * temperature_index + 0.1d0 * sigma_index
  torque_hs_baro = 1.1d0 * (1.5d0 - sigma_index)
  torque_c_lin_baro = 0.7d0 * (1.5d0 - sigma_index)
  
end subroutine init_globals

  function get_F(p)
  ! F function (22) of the paper : "A torque formula for non-isothermal Type I planetary migration - II Effects of diffusion"
  ! By Paardekooper, baruteau and Kley, (2010)
  ! Equation (22) can be approximated within 5% by a much more simpler equation which will be the one we use. 
  !
  ! Parameter:
  ! p : parameter
    implicit none
    real(double_precision), intent(in) :: p
    
    real(double_precision) :: get_F

    !Local
    !------------------------------------------------------------------------------

    
    get_F = 1.d0 / (1.d0 + (p / 1.3d0) * (p / 1.3d0))
    
    return
  end function get_F

  function get_G(p)
  ! G function (30) of the paper : "A torque formula for non-isothermal Type I planetary migration - II Effects of diffusion"
  ! By Paardekooper, baruteau and Kley, (2010)
  ! Equation (30) is a particuliar case of the (29) equation with l=3/4 and tau_0=2
  !
  ! Parameter:
  ! p : parameter
    implicit none
    real(double_precision), intent(in) :: p
    
    real(double_precision) :: get_G

    !Local
    real(double_precision), parameter :: p_0 = sqrt(8.d0 / (45.d0 * PI)) ! With tau_0 = 2
    !------------------------------------------------------------------------------

      
    if (p.lt.p_0) then
      get_G = (16.d0 / 25.d0) * (p / p_0)**(1.5d0)
    else
      get_G = 1.d0 - (9.d0 / 25.d0) * (p_0 / p)**(8.d0/3.d0)
    end if
    
    return
  end function get_G

  function get_K(p)
  ! G function (30) of the paper : "A torque formula for non-isothermal Type I planetary migration - II Effects of diffusion"
  ! By Paardekooper, baruteau and Kley, (2010)
  ! Equation (30) is a particuliar case of the (29) equation with l=3/4 and tau_0=2
  !
  ! Parameter:
  ! p : parameter
    implicit none
    real(double_precision), intent(in) :: p
    
    real(double_precision) :: get_K

    !Local
    real(double_precision), parameter :: p_0 = sqrt(28.d0 / (45.d0 * PI)) ! With tau_0 = 7
    !------------------------------------------------------------------------------

    
    if (p.lt.p_0) then
      get_K = (16.d0 / 25.d0) * (p / p_0)**(3.d0/2.d0)
    else
      get_K = 1.d0 - (9.d0 / 25.d0) * (p_0 / p)**(8.d0/3.d0)
    end if
    
    return
  end function get_K

  function get_opacity(temperature, num_bulk_density)
  ! subroutine that return the opacity of the disk at the location of the planet given various parameters
    implicit none
    
    ! Inputs 
    real(double_precision), intent(in) :: temperature & ! temperature of the disk [K]
                                          , num_bulk_density ! bulk density of the gas disk [MSUN/AU^3] (in numerical units)
    
    ! Output
    real(double_precision) :: get_opacity
    
    ! Local
    real(double_precision) :: temp34, temp45, temp56, temp67
    real(double_precision) :: bulk_density ! [g/cm^3] (in physical units needed for the expression of the opacity)
    real(double_precision), parameter :: num_to_phys_bulk_density = MSUN / AU**3
    real(double_precision), parameter :: phys_to_num_opacity = MSUN / AU**2

    ! we convert the bulk_density from numerical units(AU, MS, DAY) to physical units (CGS)
    bulk_density = num_to_phys_bulk_density * num_bulk_density
    
    
    ! We get the transition point between the various regimes
    temp34 = 2286.787d0 * bulk_density**(2.d0/49.d0)
    temp45 = 2029.764d0 * bulk_density**(1.d0/81.d0)
    temp56 = 10000.d0 * bulk_density**(1.d0/21.d0)
    temp67 = 30000.d0 * bulk_density**(4.d0/75.d0)
    
    if (temperature.le.166.81d0) then ! regime 1
      get_opacity = 2d-4 * temperature * temperature
    elseif ((temperature.gt.166.81d0).and.(temperature.le.202.677d0)) then ! regime 2
      get_opacity = 2.d16 / temperature**7.d0
    elseif ((temperature.gt.202.677d0).and.(temperature.le.temp34)) then ! regime 3
      get_opacity = 0.1d0 * sqrt(temperature)
    elseif ((temperature.gt.temp34).and.(temperature.le.temp45)) then ! regime 4
      get_opacity = 2.d81 * bulk_density / temperature**24.d0
    elseif ((temperature.gt.temp45).and.(temperature.le.temp56)) then ! regime 5
      get_opacity = 1.d-8 * bulk_density**TWOTHIRD * temperature**3
    elseif ((temperature.gt.temp56).and.(temperature.le.temp67)) then ! regime 6
      get_opacity = 1.d-36 * bulk_density**THIRD * temperature**10
    else ! regime 7
      get_opacity = 1.5d20 * bulk_density / temperature**2.5d0
    endif
    
    ! we change the opacity from physical units to numerical units
    get_opacity = phys_to_num_opacity * get_opacity
    
    !------------------------------------------------------------------------------

    return
  end function get_opacity

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
    
    call test_functions_FGK()
    call test_get_opacity()
    call test_torques()
    
  end subroutine unitary_tests

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
    
    
    ! We open the file where we want to write the outputs
    open(10, file='test_functions_FGK.dat')
    write(10,*) "Correspond to the figure 2 of II. Effects of diffusion"
    write(10,*) 'p ; F(p) ; G(p) ; K(p)'
    
    
    do i=1,nb_points
      p = p_min * p_step ** (i-1)
      f_p = get_F(p)
      g_p = get_G(p)
      k_p = get_K(p)
      
      write(10,*) p, f_p, g_p, k_p
    end do
    
    close(10)
    
    open(10, file="functions_FGK.gnuplot")

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

  subroutine test_get_opacity
  ! subroutine that test the function 'get_opacity'
  
  ! Return:
  !  a data file 'test_opacity.dat' 
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
    
    ! We open the file where we want to write the outputs
    open(10, file='test_opacity.dat')
    write(10,*) "Correspond to the figure 9a of Bell & Lin 1994"
    write(10,*) 'Temperature (K) ; Opacity for bulk density from 1e-5 to 1e-9 by power of ten'
    
    do i=1, nb_points
      temperature = T_min * T_step ** (i-1)
      
      do j=1,5
        opacity(j) = get_opacity(temperature, bulk_density(j)) * MSUN / AU**3
      end do
      
      write(10,*) temperature, opacity(1), opacity(2), opacity(3), opacity(4), opacity(5)
    end do
    close(10)
    
    open(10, file="opacity.gnuplot")
    write(10,*) 'set terminal x11 enhanced'
    write(10,*) 'set xlabel "Temperature T"'
    write(10,*) 'set ylabel "Opacity {/Symbol k}"'
    write(10,*) 'set logscale x'
    write(10,*) 'set logscale y'
    write(10,*) 'set grid'
    write(10,*) 'set xrange [', T_min, ':', T_max, ']'
    write(10,*) "plot 'test_opacity.dat' using 1:2 with lines title '{/Symbol r}=10^{-5'},\"
    write(10,*) "     '' using 1:3 with lines title '{/Symbol r}=10^{-6}',\"
    write(10,*) "     '' using 1:4 with lines title '{/Symbol r}=10^{-7}',\"
    write(10,*) "     '' using 1:5 with lines title '{/Symbol r}=10^{-8}',\"
    write(10,*) "     '' using 1:6 with lines title '{/Symbol r}=10^{-9}'"
    write(10,*) "pause -1 # wait until a carriage return is hit"
    write(10,*) "set terminal pdfcairo enhanced"
    write(10,*) "set output 'opacity.pdf'"
    write(10,*) "replot # pour générer le fichier d'output"  
    
    close(10)
    
  end subroutine test_get_opacity
  
  subroutine test_torques
  ! subroutine that test the function 'get_corotation_torque'
  
  ! Return:
  !  a data file 'test_total_torque.dat' 
  ! and an associated gnuplot file 'total_torque.gnuplot' that display values for get_corotation_torque for a range of semi major axis.
    implicit none
    
    integer, parameter :: nb_mass = 100
    real(double_precision), parameter :: mass_min = 1. * EARTH_MASS
    real(double_precision), parameter :: mass_max = 60. * EARTH_MASS
    real(double_precision), parameter :: mass_step = (mass_max - mass_min) / (nb_mass - 1.d0)
    
    integer, parameter :: nb_points = 200
    real(double_precision), parameter :: a_min = 0.1
    real(double_precision), parameter :: a_max = 30.
    ! step for log sampling
    real(double_precision), parameter :: a_step = (a_max - a_min) / (nb_points-1.d0)
    
    real(double_precision) :: a, mass, total_torque, total_torque_units, corotation_torque, gamma_eff, torque_ref
    real(double_precision) :: stellar_mass, position(3), velocity(3)
    
    integer :: i,j ! for loops
    
    position(:) = 0.d0
    velocity(:) = 0.d0
    
    ! stellar mass
    stellar_mass = 1.d0
    
    call init_globals(stellar_mass)
    
    ! We open the file where we want to write the outputs
    open(10, file='test_corotation_torque.dat')
    open(11, file='test_total_torque.dat')
    open(12, file='test_lindblad_torque.dat')
    open(13, file='test_ref_torque.dat')
    open(14, file='test_gamma_eff.dat')
    open(15, file='test_total_torque_units.dat')
    
    write(10,*) 'semi major axis (AU) ; mass in earth mass ; corotation torque (no dim)'
    write(11,*) 'semi major axis (AU) ; mass in earth mass ; total torque (no dim)'
    write(12,*) 'semi major axis (AU) ; mass in earth mass ; lindblad torque (no dim)'
    write(13,*) 'semi major axis (AU) ; mass in earth mass ; reference torque in M_s.AU^2.day^{-2}'
    write(14,*) 'semi major axis (AU) ; mass in earth mass ; gamma effective'
    write(15,*) 'semi major axis (AU) ; mass in earth mass ; total torque in M_s.AU^2.day^{-2}'
    
    do i=1, nb_points ! loop on the position
      a = a_min + a_step * (i-1)
      
      ! We generate cartesian coordinate for the given semi major axis
      position(1) = a
      
      
      do j=1,nb_mass
        mass = mass_min + mass_step * (j - 1.d0)
        
        ! We generate cartesian coordinate for the given mass and semi major axis
        velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
        
        call get_corotation_torque(stellar_mass, mass, position, velocity, corotation_torque, torque_ref, gamma_eff)
        
        total_torque = torque_lindblad + corotation_torque
        total_torque_units = (torque_ref / gamma_eff) * total_torque
        
!~         corotation_torque = torque_ref / gamma_eff * corotation_torque
!~         torque_lindblad = torque_ref / gamma_eff * torque_lindblad
!~         torque = torque_ref / gamma_eff * torque
                
        write(10,*) a, mass / EARTH_MASS, corotation_torque
        write(11,*) a, mass / EARTH_MASS, total_torque
        write(12,*) a, mass / EARTH_MASS, torque_lindblad
        write(13,*) a, mass / EARTH_MASS, torque_ref
        write(14,*) a, mass / EARTH_MASS, gamma_eff
        write(15,*) a, mass / EARTH_MASS, total_torque_units
      end do
      
      do j=10,15
        write(j,*) ""! we write a blank line to separate them in the data file, else, gnuplot doesn't want to make the surface plot
      end do
    end do
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    
    
    open(10, file="corotation_torque.gnuplot")
    open(11, file="total_torque.gnuplot")
    open(12, file="lindblad_torque.gnuplot")
    open(13, file="ref_torque.gnuplot")
    open(14, file="gamma_eff.gnuplot")
    open(15, file="total_torque_units.gnuplot")
    
    do j=10,15
      write(j,*) 'set terminal x11 enhanced'
      write(j,*) 'set xlabel "semi major axis (AU)"'
      write(j,*) 'set ylabel "Planet mass (m_{earth})" center'
    end do
    
    write(10,*) 'set title "Evolution of the corotation torque {/Symbol g}_{eff}{/Symbol G}_c/{/Symbol G}_0"'
    write(11,*) 'set title "Evolution of the total torque {/Symbol g}_{eff}{/Symbol G}_{tot}/{/Symbol G}_0 "'
    write(12,*) 'set title "Evolution of the lindblad torque {/Symbol g}_{eff}{/Symbol G}_L/{/Symbol G}_0 "'
    write(13,*) 'set title "Evolution of the reference torque {/Symbol G}_0 [M_s.AU^2.day^{-2}]"'
    write(14,*) 'set title "Evolution of {/Symbol g}_{eff}"'
    write(15,*) 'set title "Evolution of the total torque {/Symbol G}_{tot} [M_s.AU^2.day^{-2}]"'
    
    do j=10,15
      write(j,*) 'set pm3d map'
!~     write(j,*) 'set xrange [', a_min, ':', a_max, '] noreverse'
!~     write(j,*) 'set yrange [', mass_min, ':', mass_max, '] noreverse'
    end do

    write(10,*) "splot 'test_corotation_torque.dat' title ''"
    write(11,*) "splot 'test_total_torque.dat' title ''"
    write(12,*) "splot 'test_lindblad_torque.dat' title ''"
    write(13,*) "splot 'test_ref_torque.dat' title ''"
    write(14,*) "splot 'test_gamma_eff.dat' title ''"
    write(15,*) "splot 'test_total_torque_units.dat' title ''"
    
    do j=10,15
      write(j,*) "#pause -1 # wait until a carriage return is hit"
      write(j,*) "set terminal pngcairo enhanced"
    end do
    
    write(10,*) "set output 'corotation_torque.png'"
    write(11,*) "set output 'total_torque.png'"
    write(12,*) "set output 'lindblad_torque.png'"
    write(13,*) "set output 'ref_torque.png'"
    write(14,*) "set output 'gamma_eff.png'"
    write(15,*) "set output 'total_torque_units.png'"
    
    do j=10,15
      write(j,*) "replot # pour générer le fichier d'output"
    end do
    
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    
  end subroutine test_torques

end module user_module

! TODO utiliser la masse des objets pour ne pas faire le calcul si trop massif, il faut respecter le domaine de validité des formules des couples
! TODO routine générale de conversion des couples en accélération afin de pouvoir réutiliser ailleurs

!TODO : 
!_ amortissement de l'eccentricité
!_ amortissement de l'inclinaison
!_rajouter le spin comme variable d'entrée de la routine mfo_user
