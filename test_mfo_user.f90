! Program that test different functions implemented in the module user_module.
program test_mfo_user
  use types_numeriques
  use user_module
  use physical_constant
  
  implicit none
  
  real(double_precision), parameter :: b_over_h = 0.4 ! the smoothing length for the planet's potential
  real(double_precision), parameter :: adiabatic_index = 1.4 ! the adiabatic index for the gas equation of state
  real(double_precision), parameter :: aspect_ratio = 0.05 ! the aspect_ratio of the disk. Is used by the function get_aspect_ratio
  
  ! Here we define the power law for surface density sigma(R) = sigma_0 * R^sigma_index
  real(double_precision), parameter :: sigma_0 = 1700 ! the surface density at (R=1AU) [g/cm^2]
  real(double_precision), parameter :: sigma_index = 0.5! the slope of the surface density power law (alpha in the paper)
  real(double_precision), parameter :: sigma_0_num = sigma_0 * AU**2 / MSUN ! the surface density at (R=1AU) [Numerical Units]
  
  ! Here we define the power law for temperature T(R) = temperature_0 * R^temperature_index
  real(double_precision), parameter :: temperature_0 = 150 ! the temperature at (R=1AU) [K]
  real(double_precision), parameter :: temperature_index = 1.! the slope of the temperature power law (beta in the paper)
  
  !prefactors
  real(double_precision) :: x_s_prefactor
  real(double_precision) :: chi_p_prefactor
  
  real(double_precision) :: torque_lindblad ! prefactor for the lindblad torque
  real(double_precision) :: torque_hs_baro ! barotropic part of the horseshoe drag
  real(double_precision) :: torque_c_lin_baro ! barotropic part of the linear corotation torque
  
  ! Units are mass in solar mass, length in AU and time in day
  
  real(double_precision) :: stellar_mass, mass, position(3), velocity(3)
  
  position(:) = 0.d0
  velocity(:) = 0.d0
  
  ! stellar mass
  stellar_mass = 1.d0
  
  ! planet mass
  mass = 1.26d-5
  
  ! planet semi major axis
  position(1) = 1.d0
  
  velocity(2) = sqrt(K2 * (stellar_mass + mass) / position(1))
  
  ! we define the planet mass, radius and velocity, considering that we set a semi major axis and a mass. 
  ! Position will be only on one component (x), the same goes for velocity (y)
  
  call unitary_tests()
  
  x_s_prefactor = 1.1d0 * (b_over_h / 0.4d0)**0.25d0 / sqrt(stellar_mass) ! mass(1) is here for the ratio of mass q

  chi_p_prefactor = (16.d0  / 3.d0) * adiabatic_index * (adiabatic_index - 1.d0) * SIGMA_STEFAN
  
  torque_lindblad = -2.5d0 - 1.7d0 * temperature_index + 0.1d0 * sigma_index
  torque_hs_baro = 1.1d0 * (3.d0/2.d0 - sigma_index)
  torque_c_lin_baro = 0.7d0 * (3.d0/2.d0 - sigma_index)
  
  call test_torque(stellar_mass, mass, position, velocity)
  
  contains

  
  subroutine test_torque(stellar_mass, mass, position, velocity)
! function that return the total torque exerted by the disk on the planet 
!
! Parameter:
! mass : an array will all the masses (including the central one as the first element [solar mass]
! radius_p : the distance between the central body and each planet
  use orbital_elements, only : mco_x2a
  implicit none
  real(double_precision), intent(in) :: stellar_mass, mass, position(3), velocity(3)
  
  real(double_precision) :: get_total_torque

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

  real(double_precision) :: Gamma_0 ! Normalization factor for all the torques. [?](equation (8) of Paardekooper, Baruteau, 2009)
  real(double_precision) :: sigma_p ! the surface density of the gas disk at the planet location [?]
  real(double_precision) :: bulk_density_p ! the bulk_density of the disk at the location of the planet [?]
  real(double_precision) :: aspect_ratio_p ! the aspect_ratio of the gas disk at the location of the planet [no dim]
  real(double_precision) :: x_s ! semi-width of the horseshoe region [radius_p (in unity of position of the planet)]
  real(double_precision) :: gamma_eff ! effective adiabatic index depending on several parameters
  real(double_precision) :: zeta_eff ! effective entropy index depending on gamma_eff
  real(double_precision) :: chi_p ! the thermal diffusion coefficient at the location of the planet [?]
  real(double_precision) :: nu_p ! the viscosity of the disk at the location of the planet [?]
  real(double_precision) :: opacity_p ! the opacity of the disk at the location of the planet [?]
  real(double_precision) :: temperature_p ! the temperature of the disk at the location of the planet [K]
  real(double_precision) :: p_nu ! parameter for saturation due to viscosity at the location of the planet
  real(double_precision) :: p_chi ! parameter for saturation due to thermal diffusion at the location of the planet
  
  !Torques (some depends of the planet)
  real(double_precision) :: torque_hs_ent ! entropy related part of the horseshoe drag
  real(double_precision) :: torque_c_lin_ent ! entropy related part of the linear corotation torque
  
  real(double_precision) :: corotation_torque(7)
  real(double_precision), dimension(4) :: list_chi_p = (/ 1e-5, 2e-6, 1e-6, 1e-7/)
  !------------------------------------------------------------------------------
  ! Parameters for the test
  real(double_precision), parameter :: p_nu_min = 0.1
  real(double_precision), parameter :: p_nu_max = 11.
  integer, parameter :: nb_points = 100
  real(double_precision), parameter :: p_nu_step = (p_nu_max/p_nu_min) ** (1/(nb_points-1.d0))

  integer :: i,j ! for loops
  !------------------------------------------------------------------------------
  
    ! WE CALCULATE PROPERTIES OF THE PLANETS
  gm = K2 * (stellar_mass + mass)
  ! We get semi_major_axis, radius_p and vel_squared
  call mco_x2a (gm,position(1), position(2), position(3), velocity(1),velocity(2), velocity(3),semi_major_axis,radius_p,vel_squared)
  
  
  !------------------------------------------------------------------------------
  sigma_p = sigma_0_num * radius_p ** sigma_index ! [N.U.]
  temperature_p = temperature_0 * radius_p**temperature_index ! [K]
  aspect_ratio_p = aspect_ratio
  
  velocity_p = sqrt(vel_squared)
  omega_p = sqrt(gm / (semi_major_axis * semi_major_axis * semi_major_axis)) ! [day-1]
  
  !------------------------------------------------------------------------------
  bulk_density_p = 0.5d0 * sigma_p / (aspect_ratio_p * radius_p)
  !------------------------------------------------------------------------------
  opacity_p = get_opacity(temperature_p, bulk_density_p)
  
  !------------------------------------------------------------------------------
  ! WE CALCULATE TOTAL TORQUE EXERTED BY THE DISK ON THE PLANET
  Gamma_0 = (mass / (stellar_mass * aspect_ratio_p))**2 * sigma_p * radius_p**4 * omega_p**2
  open(10,file='test_corotation_torque.dat')
  write(10,*) "Correspond to the figure 3 of II. Effects of diffusion"
  write(10,*) 'p_nu ; coration torque (in units of gamma_0)'
  do i=1,nb_points
    do j=1,4
      chi_p = list_chi_p(j) * radius_p**2 * omega_p
      !chi_p = chi_p_prefactor * temperature_p**4 / (opacity_p * sigma_p**2 * omega_p**2)
      
      !------------------------------------------------------------------------------
      ! Q is needed by the lindblad torque. We set Q for m ~ 2 /3 h (45): 
      Q_p = TWOTHIRD * chi_p / (aspect_ratio_p**3 * radius_p**2 * omega_p)
      !------------------------------------------------------------------------------
      
      gamma_eff = 2.d0 * Q_p * adiabatic_index / (adiabatic_index * Q_p + 0.5d0 * &
      sqrt(2.d0 * sqrt((adiabatic_index * adiabatic_index * Q_p * Q_p + 1.d0)**2 + 16.d0 * Q_p * Q_p * (adiabatic_index - 1.d0)) &
      + 2.d0 * adiabatic_index * adiabatic_index * Q_p * Q_p - 2.d0))
      
      !------------------------------------------------------------------------------
      zeta_eff = temperature_index - (gamma_eff - 1.d0) * sigma_index
      
      x_s = x_s_prefactor / gamma_eff**0.25d0 * sqrt(mass / aspect_ratio_p)
      
      !------------------------------------------------------------------------------
      ! k_p is defined to limit the number of operation and to have a value independant from chi_p or nu_p
      k_p = radius_p * radius_p * omega_p * x_s * x_s * x_s / (2.d0 * PI)
      

      
      !------------------------------------------------------------------------------

      
        p_nu = p_nu_min * p_nu_step ** (i-1)
        
    !~   p_nu = TWOTHIRD * sqrt(k_p / nu_p)
      
      p_chi = sqrt(k_p / chi_p)
      
      torque_hs_ent = 7.9d0 * zeta_eff / gamma_eff
      torque_c_lin_ent = (2.2d0 - 1.4d0 / gamma_eff) * zeta_eff
      
      !------------------------------------------------------------------------------

      corotation_torque(j) = (1 / gamma_eff) * (torque_hs_baro * get_F(p_nu) * get_G(p_nu) + torque_c_lin_baro * (1 - get_K(p_nu)) &
        + torque_hs_ent * get_F(p_nu) * get_F(p_chi) * sqrt(get_G(p_nu) * get_G(p_chi)) &
        + torque_c_lin_ent * sqrt((1 - get_K(p_nu)) * (1 - get_K(p_chi))))
        
      
    end do
    write(10,*) p_nu, corotation_torque(1), corotation_torque(2), corotation_torque(3), corotation_torque(4)
    !------------------------------------------------------------------------------
  end do
  close(10)
  
  open(10, file="corotation.gnuplot")
  write(10,*) 'set terminal x11 enhanced'
  write(10,*) 'set xlabel "p_{/Symbol n}"'
  write(10,*) 'set ylabel "corotation torque {/Symbol G}_c/{/Symbol G}_0"'
  write(10,*) 'set logscale x'
  write(10,*) 'set grid'
  write(10,*) 'set xrange [', p_nu_min, ':', p_nu_max, ']'
  write(10,*) "plot 'test_corotation_torque.dat' using 1:2 with lines title '{/Symbol c}_p=10^{-5}', \"
  write(10,*) "     '' using 1:3 with lines title '{/Symbol c}_p=2* 10^{-6}',\"
  write(10,*) "     '' using 1:4 with lines title '{/Symbol c}_p= 10^{-6}',\"
  write(10,*) "     '' using 1:5 with lines title '{/Symbol c}_p= 10^{-7}'"
!~   write(10,*) "pause -1 # wait until a carriage return is hit"
  write(10,*) "set terminal pdfcairo enhanced"
  write(10,*) "set output 'corotation.pdf'"
  write(10,*) "replot # pour générer le fichier d'output"  
  
  close(10)

    
    
  
  return
end subroutine test_torque

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
      get_G = (16.d0 / 25.d0) * (p / p_0)**(3.d0/2.d0)
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
    
    if (temperature.le.166.81d0) then
      get_opacity = 2d-4 * temperature * temperature
    elseif ((temperature.gt.166.81d0).and.(temperature.le.202.677d0)) then
      get_opacity = 2.d16 / temperature**7.d0
    elseif ((temperature.gt.202.677d0).and.(temperature.le.temp34)) then
      get_opacity = 0.1d0 * sqrt(temperature)
    elseif ((temperature.gt.temp34).and.(temperature.le.temp45)) then
      get_opacity = 2.d81 * bulk_density / temperature**24.d0
    elseif ((temperature.gt.temp45).and.(temperature.le.temp56)) then
      get_opacity = 1.d-8 * bulk_density**TWOTHIRD * temperature**3
    elseif ((temperature.gt.temp56).and.(temperature.le.temp67)) then
      get_opacity = 1.d-36 * bulk_density**THIRD * temperature**10
    else
      get_opacity = 1.5d20 * bulk_density / temperature**2.5d0
    endif
    
    ! we change the opacity from physical units to numerical units
    get_opacity = phys_to_num_opacity * get_opacity
    
    !------------------------------------------------------------------------------

    return
  end function get_opacity
end program test_mfo_user
