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
!
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
  real(double_precision) :: ecc_corot ! prefactor that turns out the corotation torque if the eccentricity is too high (Bitsch & Kley, 2010)
  
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
  
  !------------------------------------------------------------------------------
  ! Setup

  do planet=1,n_bodies
    acceleration(1:3,planet) = 0.d0
  end do
  
  ! By default, there is disk effects. Be carefull, init_globals is only treated if there is disk effects, 
  ! to increase the speed when the disk is no present anymore.
  if (disk_effect) then
    call init_globals(stellar_mass=mass(1), time=time)
    
    !------------------------------------------------------------------------------
    ! If it's time (depending on the timestep we want between each calculation of the disk properties)
    ! The first 'next_dissipation_step' is set to '-1' to force the calculation for the first timestep. In fact, the first timestep will be done fornothing, but we need this in order to have a clean code.
    if (DISSIPATION_TYPE.ne.0) then
      if (time.gt.next_dissipation_step) then
        ! we get the density profile.
        call dissipate_disk(time, next_dissipation_step)
        
        ! we get the temperature profile.
        call calculate_temperature_profile()
        
        ! we store in a .dat file the temperature profile
        call store_temperature_profile(filename='temperature_profile.dat')
        call store_density_profile(filename='surface_density_profile.dat')
        call store_scaleheight_profile()
      end if
    end if
    !------------------------------------------------------------------------------
  
    do planet=2,n_bodies
      ! because ongoing deletion of planets put their mass to 0 for a few steps, we must check. Else, we will have an error "NaN".
      if (mass(planet).gt.TINY) then
        !------------------------------------------------------------------------------
        ! we store in a structure, various properties of the planet usefull in all the 
        ! subroutine to avoid multiple calculation of the same parameters
        call get_planet_properties(stellar_mass=mass(1), mass=mass(planet), & ! input
        position=position(1:3, planet), velocity=velocity(1:3,planet),& ! input
        p_prop=p_prop) ! Output
        
        if (planet.le.n_big_bodies) then
          if ((p_prop%eccentricity.lt.ECCENTRICITY_CUTOFF).and.(p_prop%radius.gt.INNER_BOUNDARY_RADIUS)) then
            !------------------------------------------------------------------------------
            ! prefactor calculation for eccentricity and inclination damping
            e_h = p_prop%eccentricity / p_prop%aspect_ratio
            i_h = p_prop%inclination / p_prop%aspect_ratio
            time_wave = mass(1)**2 * p_prop%aspect_ratio**4 / (mass(planet) * K2 * p_prop%sigma * &
                        p_prop%semi_major_axis**2 * p_prop%omega) ! (cresswell, 2008)

            !------------------------------------------------------------------------------
            ! Calculation of the acceleration due to migration
            call get_torques(mass(1), mass(planet), p_prop, & ! input
                corotation_torque=corotation_torque, lindblad_torque=lindblad_torque, Gamma_0=torque_ref, ecc_corot=ecc_corot) ! Output
            
            
            torque = torque_ref * (lindblad_torque + ecc_corot * corotation_torque)
    !~         ! use this line instead if you want to cutoff the corotation torque damping due to the eccentricity
    !~         torque = torque_ref * (lindblad_torque + corotation_torque)
            

            time_mig = 0.5d0 * p_prop%angular_momentum / torque
            
            migration_acc_prefactor = 1.d0 / time_mig
            
            migration_acceleration(1:3) = migration_acc_prefactor * velocity(1:3,planet)
            
            acceleration(1:3,planet) = acceleration(1:3,planet) + migration_acceleration(1:3)

            
            !------------------------------------------------------------------------------
            ! Calculation of the acceleration due to eccentricity damping
            
            time_ecc = time_wave / 0.780d0 * (1.d0 - 0.14d0 * e_h**2 + 0.06 * e_h**3 + 0.18 * e_h * i_h**2)
            
            eccentricity_acc_prefactor = -2.d0 * (position(1,planet) * velocity(1,planet) + position(2,planet) * velocity(2,planet)&
             + position(3,planet) * velocity(3,planet)) / (p_prop%radius**2 * time_ecc)
            
            eccentricity_acceleration(1:3) = eccentricity_acc_prefactor * position(1:3,planet)
            
            acceleration(1:3,planet) = acceleration(1:3,planet) + eccentricity_acceleration(1:3)
            
            
            !------------------------------------------------------------------------------
            ! Calculation of the acceleration due to the inclination damping
            
            if (p_prop%inclination.gt.INCLINATION_CUTOFF) then
              time_inc = time_wave / 0.544d0 * (1.d0 - 0.30d0 * i_h**2 + 0.24 * i_h**3 + 0.14 * e_h**2 * i_h)
              
              inclination_acceleration_z = - velocity(3,planet) / time_inc
              acceleration(3,planet) = acceleration(3,planet) + inclination_acceleration_z
            end if
            
          end if
        end if
          
        !------------------------------------------------------------------------------
        ! Calculation of the acceleration due to turbulence, if needed
        
        if (IS_TURBULENCE) then
          turbulence_acceleration(1:3) = 0.d0
          
          call get_turbulence_acceleration(time, p_prop, position(1:3, planet), turbulence_acceleration(1:3))
          
          acceleration(1:3,planet) = acceleration(1:3,planet) + turbulence_acceleration(1:3)
        end if
        
  !~       if (time.gt.365.25) then
  !~         if ((p_prop%eccentricity.lt.ECCENTRICITY_CUTOFF).and.(p_prop%radius.gt.INNER_BOUNDARY_RADIUS)) then
  !~         
  !~           call debug_infos(time, n_bodies, planet, position, velocity, acceleration, &
  !~                        time_mig, migration_acceleration, time_ecc, eccentricity_acceleration, &
  !~                        turbulence_acceleration, corotation_torque, lindblad_torque, torque_ref, ecc_corot)
  !~         end if
  !~         open(12, file="debug.out", access='append')
  !~         write (12,*) "time = ", time, " ; planet", planet
  !~         call print_planet_properties(p_prop, output=12)
  !~         close(12)
  !~         stop
  !~       end if
        
      end if
    end do
  end if
  
  !------------------------------------------------------------------------------
  return
end subroutine mfo_user

subroutine debug_infos(time, n_bodies, planet, position, velocity, acceleration, &
                       time_mig, migration_acceleration, time_ecc, eccentricity_acceleration, &
                       turbulence_acceleration, corotation_torque, lindblad_torque, torque_ref, ecc_corot)
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
  
  real(double_precision), intent(in) :: corotation_torque
  real(double_precision), intent(in) :: lindblad_torque
  real(double_precision), intent(in) :: torque_ref
  real(double_precision), intent(in) :: ecc_corot

! Local
character(len=80) :: output_format
character(len=80) :: single_format = 'es10.3e2'

!------------------------------------------------------------------------------
write(output_format, *) '(a,',trim(single_format),'," ",a,',trim(single_format),'," ",a,',trim(single_format),')'

open(12, file="debug.out", access='append')
write(12,'(a)') '################################'
write(12,'(a, f10.2, a)') 'time = ', time, ' days'
write(12,'(a)') '################################'
write(12,output_format) ' x = ', position(1, planet), ' y = ', position(2, planet), ' z = ', position(3, planet)
write(12,output_format) 'vx = ', velocity(1, planet), 'vy = ', velocity(2, planet), 'vz = ', velocity(3, planet)
write(12,output_format) 'ax = ', acceleration(1, planet), 'ay = ', acceleration(2, planet), 'az = ', acceleration(3, planet)
write(12,'(a)') '------------------------------------'
write(12,'(a)') '|            Migration             |'
write(12,'(a)') '------------------------------------'
write(12,'(a,es10.3e2,a)') 'migration timescale =', time_mig, ' days'
write(12,output_format) 'amx = ', migration_acceleration(1), &
                       'amy = ', migration_acceleration(2), &
                       'amz = ', migration_acceleration(3)
write(12,*) 'lindblad torque = ', lindblad_torque
write(12,*) 'corotation torque = ', corotation_torque
write(12,*) 'Gamma_0 = ', torque_ref, '[Ms.AU^2]'
write(12,*) 'corotation damping = ', ecc_corot

write(12,'(a)') '------------------------------------'
write(12,'(a)') '|       Eccentricity damping       |'
write(12,'(a)') '------------------------------------'
write(12,'(a,es10.3e2,a)') 'Eccentricity timescale = ', time_ecc, ' days'
write(12,output_format) 'aex = ', eccentricity_acceleration(1), &
                       'aey = ', eccentricity_acceleration(2), &
                       'aez = ', eccentricity_acceleration(3)
write(12,'(a)') '------------------------------------'
write(12,'(a)') '|           Turbulence             |'
write(12,'(a)') '------------------------------------'
write(12,output_format) 'atx = ', turbulence_acceleration(1), &
                       'aty = ', turbulence_acceleration(2), &
                       'atz = ', turbulence_acceleration(3)
write(12,'(a)') '____________________________________'
close(12)

end subroutine debug_infos


end module user_module
  
! TODO utiliser la masse des objets pour ne pas faire le calcul si trop massif, il faut respecter le domaine de validité des formules des couples
! TODO routine générale de conversion des couples en accélération afin de pouvoir réutiliser ailleurs

!TODO : 
!_rajouter le spin comme variable d'entrée de la routine mfo_user
