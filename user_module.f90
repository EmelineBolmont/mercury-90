module user_module

!*************************************************************
!** Modules that contains user defined modules. 
!** Only mfo_user will be public.
!**
!** Version 1.0 - june 2011
!*************************************************************

! The user_module is divided in several parts. First of all, some parameters are calculated at the beginning of the run, in 
! init_globals. This subroutine must be executed before the rest. Since mfo_user is used by mercury, we can't execute init_globals 
! only once, we must include an 'if' test that  run init_globals only if it's the first time we pass in the routine. Another thing 
! to keep in mind is that the temperature profile of the disk is calculated manually. Since the density profile evolve in time, You
! must re-calculate the temperature profile each time you need it. A mistake can be made because the surface density will evolve 
! in time automatically. BUT, the temperature profile will not correspond to it. Another set of routine are prefixed with 'test_', 
! they are all coded to test something, either a routine or a physical value (and plot it). Each plot consist in a data file (*
! .dat) and a gnuplot file (*.gnuplot). You only have to run the gnuplot file with gnuplot with a command "gnuplot file.gnuplot". 

! Most of the plot are stored in a subdirectory "unitary_tests" of the current directory Theses tests are used ONLY in the 
! subroutine "unitary_tests" that can be used outside the module. For instance, I made a fortran source code that use this module 
! and call the routine 'unitary_tests'. That's how I test my module before running it in mercury. The problem is to prepare 
! variable in the same way, both in mfo_user and my tests. That's why init_globals exists, thus I can run it in my tests also.

  use types_numeriques
  use physical_constant
  use disk ! contains all the subroutines that describe the behaviour of the disk

  implicit none
  
  private
  
  public :: mfo_user
  
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

! Code Units are in AU, days and solar mass * K2 (mass are in solar mass, but multiplied by K2 earlier in the code).

! N.B. All coordinates and velocities must be with respect to central body
! ===
!------------------------------------------------------------------------------

subroutine mfo_user (time,jcen,n_bodies,n_big_bodies,mass,position,velocity,acceleration)
!  mass          = mass (in solar masses * K2)
!  position      = coordinates (x,y,z) with respect to the central body [AU]
!  velocity      = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  n_bodies      = current number of bodies (INCLUDING the central object)
!  n_big_bodies  =    "       "    " big bodies (ones that perturb everything else)
!  time          = current epoch [days]

! Global parameters
! DISSIPATION_TYPE : boolean to tell if there is dissipation of the disk or not.
! dissipation_timestep : the timestep between two computation of the disk [in days]
! X_SAMPLE_STEP : the constant step for the x_sample. Indeed, due to diffusion equation, the sample must be constant in X, and not in r. 
  use physical_constant
  use mercury_constant
  use turbulence
  use utilities, only : vect_product
  use utilities, only : get_mean, get_stdev, get_histogram!, vect_product
  
  implicit none

  
  ! Input
  integer, intent(in) :: n_bodies,n_big_bodies
  real(double_precision),intent(in) :: time,jcen(3),mass(n_bodies),position(3,n_bodies),velocity(3,n_bodies)
  
  ! Output
  real(double_precision),intent(out) :: acceleration(3,n_bodies)
  
  !------------------------------------------------------------------------------ 
  !------Local-------

  real(double_precision) :: migration_acc_prefactor ! prefactor for the migration acceleration
  real(double_precision) :: eccentricity_acc_prefactor ! prefactor for the eccentricity acceleration

  ! loop integers
  integer :: planet
  
  
  real(double_precision) :: torque ! the torque exerted by the disk on the planet [Ms.AU^2]
  real(double_precision) :: corotation_torque ! the corotation torque exerted by the disk on the planet in unit of torque_ref [No dim]
  real(double_precision) :: lindblad_torque ! the lindblad torque exerted by the disk on the planet in unit of torque_ref [No dim]
  real(double_precision) :: torque_ref ! a ref torque that depends on the properties of the planet [Ms.AU^2]
  real(double_precision) :: time_mig ! The migration timescale for the planet [day]
  real(double_precision) :: time_wave ! A timescale for the planet that I don't understand for the moment [day]
  real(double_precision) :: time_ecc ! The eccentricity damping timescale for the planet [day]
  real(double_precision) :: time_inc ! The inclination damping timescale for the planet [day]
  
  !local temporary parameters
  type(PlanetProperties) :: p_prop ! various properties of a planet
  real(double_precision) :: e_h ! the ratio between the eccentricity and the aspect ratio for a given planet [no dim]
  real(double_precision) :: i_h ! the ratio between the inclination and the aspect ratio for a given planet [no dim]
  
  real(double_precision), dimension(3) :: migration_acceleration
  real(double_precision), dimension(3) :: eccentricity_acceleration
  real(double_precision), dimension(3) :: turbulence_acceleration
  real(double_precision) :: inclination_acceleration_z
  real(double_precision), save :: next_dissipation_step = -1.d0 ! next time at which we will compute the thermal properties of the disk?
!~   
  ! Temp
!~   integer, parameter :: nb_points = 100000
!~   integer, parameter :: nb_bins = 100
!~   integer :: i = 0
!~   real(double_precision), dimension(nb_points), save :: turbulence_acc_log
!~   real(double_precision) :: delta_bin
!~   real(double_precision), dimension(nb_bins) :: bin_x_values, bin_y_values, gauss_fit
!~   real(double_precision) :: mean, stdev, y_max
!~   real(double_precision), dimension(3) :: tmp3
  !------------------------------------------------------------------------------
  ! Setup
  
  do planet=1,n_bodies
    acceleration(1,planet) = 0.d0
    acceleration(2,planet) = 0.d0
    acceleration(3,planet) = 0.d0
  end do
  
  call init_globals(stellar_mass=mass(1), time=time)
  
  !------------------------------------------------------------------------------
  ! If it's time (depending on the timestep we want between each calculation of the disk properties)
  ! The first 'next_dissipation_step' is set to '-1' to force the calculation for the first timestep. In fact, the first timestep will be done fornothing, but we need this in order to have a clean code.
  if (DISSIPATION_TYPE.ne.0) then
    if (time.gt.next_dissipation_step) then
      ! we get the density profile.
      select case(DISSIPATION_TYPE)
        case(1) ! viscous dissipation
          dissipation_timestep = 0.5d0 * X_SAMPLE_STEP**2 / (4 * get_viscosity(1.d0)) ! a correction factor of 0.5 has been applied. No physical reason to that, just intuition and safety
          ! TODO if the viscosity is not constant anymore, the formulae for the dissipation timestep must be changed
          next_dissipation_step = time + dissipation_timestep
          call dissipate_density_profile() ! global parameter 'dissipation_timestep' must exist !
        
        case(2) ! exponential decay
          ! we want 1% variation : timestep = - tau * ln(1e-2)
          dissipation_timestep = 4.6 * TAU_DISSIPATION
          next_dissipation_step = time + dissipation_timestep
          
          call exponential_decay_density_profile()
          
        case default
          write(*,*) 'Warning: The dissipation rule cannot be found.'
          write(*,*) 'Given value :', DISSIPATION_TYPE
      end select
      
      ! we get the temperature profile.
      call calculate_temperature_profile()
      
      ! we store in a .dat file the temperature profile
      call store_temperature_profile(filename='temperature_profile.dat')
      call store_density_profile(filename='density_profile.dat')
      call store_scaleheight_profile()
    end if
  end if
  !------------------------------------------------------------------------------

  
  do planet=2,n_big_bodies
    ! because ongoing deletion of planets put their mass to 0 for a few steps, we must check. Else, we will have an error "NaN".
    if (mass(planet).gt.TINY) then
      !------------------------------------------------------------------------------
      ! we store in a structure, various properties of the planet usefull in all the 
      ! subroutine to avoid multiple calculation of the same parameters
      call get_planet_properties(stellar_mass=mass(1), mass=mass(planet), & ! input
      position=position(1:3, planet), velocity=velocity(1:3,planet),& ! input
      p_prop=p_prop) ! Output
      
      if ((p_prop%eccentricity.lt.1.d0).and.(p_prop%semi_major_axis.gt.0.3d0)) then
        !------------------------------------------------------------------------------
        ! prefactor calculation for eccentricity and inclination damping
        e_h = p_prop%eccentricity / p_prop%aspect_ratio
        i_h = p_prop%inclination / p_prop%aspect_ratio
        time_wave = mass(1)**2 * p_prop%aspect_ratio**4 / (mass(planet) * K2 * p_prop%sigma * &
                    p_prop%semi_major_axis**2 * p_prop%omega) ! (cresswell, 2008)

        !------------------------------------------------------------------------------
        ! Calculation of the acceleration due to migration
        select case(TORQUE_TYPE)
          case('real') ! The normal torque profile, calculated form properties of the disk
            call get_corotation_torque(mass(1), mass(planet), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
          
          ! for retrocompatibility, 'mass_independant' has been added and refer to the old way of defining a mass-indep convergence zone
          case('linear_indep', 'mass_independant') ! a defined torque profile to get a mass independant convergence zone
            call get_corotation_torque_linear_indep(mass(1), mass(planet), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
          
          case('arctan_indep') ! a defined torque profile to get a mass independant convergence zone
            call get_corotation_torque_arctan_indep(mass(1), mass(planet), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
          
          case('mass_dependant')
            call get_corotation_torque_mass_dep_CZ(mass(1), mass(planet), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
            
          case('manual')
            call get_corotation_torque_manual(mass(1), mass(planet), p_prop, & ! input
            corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref) ! Output
            
          case default
            write(*,*) 'Warning: The torque rule cannot be found.'
            write(*,*) 'Given value :', TORQUE_TYPE
            write(*,*) 'Values possible : real ; linear_indep ; arctan_indep ; manual'
        end select
        

        torque = torque_ref * (lindblad_torque + corotation_torque)      
        
        time_mig = 0.5d0 * p_prop%angular_momentum / torque
        
        migration_acc_prefactor = 1.d0 / time_mig
        
        migration_acceleration(1) = migration_acc_prefactor * velocity(1,planet)
        migration_acceleration(2) = migration_acc_prefactor * velocity(2,planet)
        migration_acceleration(3) = migration_acc_prefactor * velocity(3,planet)
        
        !------------------------------------------------------------------------------
        ! Calculation of the acceleration due to eccentricity damping
        
        time_ecc = time_wave / 0.780d0 * (1.d0 - 0.14d0 * e_h**2 + 0.06 * e_h**3 + 0.18 * e_h * i_h**2)
        
        eccentricity_acc_prefactor = -2.d0 * (position(1,planet) * velocity(1,planet) + position(2,planet) * velocity(2,planet) + &
        position(3,planet) * velocity(3,planet)) / (p_prop%radius**2 * time_ecc)
        
        eccentricity_acceleration(1) = eccentricity_acc_prefactor * position(1,planet)
        eccentricity_acceleration(2) = eccentricity_acc_prefactor * position(2,planet)
        eccentricity_acceleration(3) = eccentricity_acc_prefactor * position(3,planet)
        
        
        !------------------------------------------------------------------------------
        ! Calculation of the acceleration due to the inclination damping
        
        if (p_prop%inclination.gt.INCLINATION_CUTOFF) then
          time_inc = time_wave / 0.544d0 * (1.d0 - 0.30d0 * i_h**2 + 0.24 * i_h**3 + 0.14 * e_h**2 * i_h)
          
          inclination_acceleration_z = - velocity(3,planet) / time_inc
        else
          time_inc = 0.d0
          inclination_acceleration_z = 0.d0
        end if
        
        !------------------------------------------------------------------------------
        ! Calculation of the acceleration due to turbulence, if needed
        turbulence_acceleration(1:3) = 0.d0
        
        if (IS_TURBULENCE) then
          call get_turbulence_acceleration(time, p_prop, position(1:3, planet), turbulence_acceleration(1:3))
        end if
        
        !------------------------------------------------------------------------------
        ! Calculation of the total acceleration on the planet
        
        acceleration(1,planet) = migration_acceleration(1) + turbulence_acceleration(1) + & 
                                 eccentricity_acceleration(1)
        acceleration(2,planet) = migration_acceleration(2) + turbulence_acceleration(2) + & 
                                 eccentricity_acceleration(2)
        acceleration(3,planet) = migration_acceleration(3) + turbulence_acceleration(3) + & 
                                 eccentricity_acceleration(3) + inclination_acceleration_z

                                 
!~         i = i + 1
!~         if (i.eq.1) then 
!~         tmp3(1:3) = vect_product(position(1:3, planet), turbulence_acceleration(1:3))
!~           
!~           turbulence_acc_log(i) = tmp3(3)
!~         
!~         else if (i.le.nb_points) then
!~           tmp3(1:3) = vect_product(position(1:3, planet), turbulence_acceleration(1:3))
!~           
!~           turbulence_acc_log(i) = tmp3(1)!(turbulence_acc_log(i-1) * float((i-1)) + tmp3(3)) / float(i)
!~           open(10, file="turbulence_torque_mean.dat", access='append')
!~           write(10,*) turbulence_acc_log(i)
!~           close(10)
!~         else
!~           call get_histogram(turbulence_acc_log, bin_x_values, bin_y_values) 
!~           
!~           ! We calculate the mean and stdev value of the data set and then generate a supposed gaussian to see if this function fit the datas
!~           ! the mean of the data set must be 0. So in order to check that, the mean is fixed to 0, to see if the gaussian looks nice.
!~           delta_bin = (bin_x_values(2) - bin_x_values(1))
!~           mean = get_mean(turbulence_acc_log(1:nb_points))
!~           stdev = get_stdev(turbulence_acc_log(1:nb_points))
!~           y_max = 1. / (stdev * sqrt(TWOPI)) * delta_bin
!~           do i=1,nb_bins
!~             gauss_fit(i) = y_max * exp(-(bin_x_values(i))**2 / (2. * stdev**2))
!~           end do
!~           
!~           open(10, file="turbulence_torque.hist")
!~           do i=1, nb_bins
!~             write(10,*) bin_x_values(i), bin_y_values(i), gauss_fit(i)
!~           end do
!~           close(10)
!~           
!~           open(10, file="turbulence_torque.gnuplot")
!~           write(10,'(a)') 'set terminal wxt enhanced'
!~           write(10,'(a)') 'set xlabel "torque [AU^2/DAY^2]"'
!~           write(10,'(a)') 'set ylabel "density of probability"'
!~           write(10,'(a)') 'set grid'
!~           write(10,'(5(a,es10.2e2))') 'set label " data : {/Symbol m}=',mean,', {/Symbol s}=',stdev,'\n&
!~                                                  & fit  : {/Symbol m}=0, {/Symbol s}=',stdev,'" at graph 0, graph 0.9'
!~           write(10,'(a)') 'plot "turbulence_torque.hist" using 1:2 with boxes linestyle 3 title "turbulence torque", \'
!~           write(10,'(a)') '"turbulence_torque.hist" using 1:3 with lines linestyle 1 title "gaussian fit"'
!~           write(10,'(a)') 'pause -1'
!~           close(10)
!~           stop
!~         end if
        
!~         if (p_prop%semi_major_axis.gt.4.) then
!~ !          call print_planet_properties(p_prop)
!~           call debug_infos(time, n_bodies, planet, position, velocity, acceleration, &
!~                        time_mig, migration_acceleration, time_ecc, eccentricity_acceleration, &
!~                        turbulence_acceleration)
!~ !          write (*,*) sqrt(sum(vect_product(position(1:3,planet), acceleration(1:3,planet))))
!~ !          stop
!~         end if
!~         if (time.gt.1.2e3) then
!~           stop
!~         end if
      end if
    end if
  end do
  
  !------------------------------------------------------------------------------
  return
end subroutine mfo_user

subroutine debug_infos(time, n_bodies, planet, position, velocity, acceleration, &
                       time_mig, migration_acceleration, time_ecc, eccentricity_acceleration, &
                       turbulence_acceleration)
! subroutine that print informations helpfull to debug, such as migration timescales, accelerations and so on

implicit none
  integer, intent(in) :: n_bodies
  real(double_precision),intent(in) :: time, position(3,n_bodies),velocity(3,n_bodies)
  real(double_precision),intent(in) :: acceleration(3,n_bodies)
  integer, intent(in) :: planet
  
  real(double_precision), intent(in) :: time_mig ! The migration timescale for the planet [day]
  real(double_precision), intent(in) :: time_ecc ! The eccentricity damping timescale for the planet [day]
  
  
  real(double_precision), dimension(3), intent(in) :: migration_acceleration
  real(double_precision), dimension(3), intent(in) :: eccentricity_acceleration
  real(double_precision), dimension(3), intent(in) :: turbulence_acceleration

! Local
character(len=80) :: output_format
character(len=80) :: single_format = 'es10.3e2'

!------------------------------------------------------------------------------
write(output_format, *) '(a,',trim(single_format),'," ",a,',trim(single_format),'," ",a,',trim(single_format),')'

write(*,'(a)') '################################'
write(*,'(a, f10.2, a)') 'time = ', time, ' days'
write(*,'(a)') '################################'
write(*,output_format) ' x = ', position(1, planet), ' y = ', position(2, planet), ' z = ', position(3, planet)
write(*,output_format) 'vx = ', velocity(1, planet), 'vy = ', velocity(2, planet), 'vz = ', velocity(3, planet)
write(*,output_format) 'ax = ', acceleration(1, planet), 'ay = ', acceleration(2, planet), 'az = ', acceleration(3, planet)
write(*,'(a)') '------------------------------------'
write(*,'(a)') '|            Migration             |'
write(*,'(a)') '------------------------------------'
write(*,'(a,es10.3e2,a)') 'migration timescale =', time_mig, ' days'
write(*,output_format) 'amx = ', migration_acceleration(1), &
                       'amy = ', migration_acceleration(2), &
                       'amz = ', migration_acceleration(3)
write(*,'(a)') '------------------------------------'
write(*,'(a)') '|       Eccentricity damping       |'
write(*,'(a)') '------------------------------------'
write(*,'(a,es10.3e2,a)') 'Eccentricity timescale = ', time_ecc, ' days'
write(*,output_format) 'aex = ', eccentricity_acceleration(1), &
                       'aey = ', eccentricity_acceleration(2), &
                       'aez = ', eccentricity_acceleration(3)
write(*,'(a)') '------------------------------------'
write(*,'(a)') '|           Turbulence             |'
write(*,'(a)') '------------------------------------'
write(*,output_format) 'atx = ', turbulence_acceleration(1), &
                       'aty = ', turbulence_acceleration(2), &
                       'atz = ', turbulence_acceleration(3)
write(*,'(a)') '____________________________________'


end subroutine debug_infos


end module user_module
  
! TODO utiliser la masse des objets pour ne pas faire le calcul si trop massif, il faut respecter le domaine de validité des formules des couples
! TODO routine générale de conversion des couples en accélération afin de pouvoir réutiliser ailleurs

!TODO : 
!_rajouter le spin comme variable d'entrée de la routine mfo_user
