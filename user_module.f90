!******************************************************************************
! MODULE: user_module
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contain user defined function. This module can call other module and subroutine. 
!! The only public routine is mfo_user, that return an acceleration that 
!! mimic a random force that depend on what the user want to model.
!
!******************************************************************************

module user_module

  use types_numeriques

  implicit none

  private

  public :: mfo_user

  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 2 March 2001
!
! DESCRIPTION: 
!> @brief Applies an arbitrary force, defined by the user.
!!\n\n
!! If using with the symplectic algorithm MAL_MVS, the force should be
!! small compared with the force from the central object.
!! If using with the conservative Bulirsch-Stoer algorithm MAL_BS2, the
!! force should not be a function of the velocities.
!! \n\n
!! Code Units are in AU, days and solar mass * K2 (mass are in solar mass, but multiplied by K2 earlier in the code).
!
!> @note All coordinates and velocities must be with respect to central body
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine mfo_user (time,jcen,nbod,nbig,m,x,v,a)
!  m             = mass (in solar masses * K2)
!  x             = coordinates (x,y,z) with respect to the central body [AU]
!  v             = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  nbod          = current number of bodies (INCLUDING the central object)
!  nbig          =    "       "    " big bodies (ones that perturb everything else)
!  time          = current epoch [days]

    use physical_constant
    use mercury_constant  
    use tides_constant_GR
    use orbital_elements, only : mco_x2el
    use spline

    implicit none

    ! Input
    integer, intent(in) :: nbod !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
    integer, intent(in) :: nbig !< [in] current number of big bodies (ones that perturb everything else)
    real(double_precision), intent(in) :: time !< [in] current epoch (days)
    real(double_precision), intent(in) :: jcen(3) !< [in] J2,J4,J6 for central body (units of RCEN^i for Ji)
    real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
    real(double_precision), intent(in) :: x(3,nbod)
    real(double_precision), intent(in) :: v(3,nbod)
    
    ! Output
    real(double_precision),intent(out) :: a(3,nbod)
    
    !------------------------------------------------------------------------------ 
    !------Local-------

    ! Local
    integer :: j,kk, error, iPs0, nptmss
    integer :: flagrg2=0
    integer :: flagtime=0
    integer :: ispin=0
    integer :: iwrite=0
    real(double_precision) :: flagbug=0.d0
    real(double_precision) :: timestep!=3.6525d5!3.65d5 !4.56d6 !
    real(double_precision) :: gm,qq,ee,ii,pp,nn,ll,Pst0,Pst
    real(double_precision) :: acc_tid_x,acc_tid_y,acc_tid_z
    real(double_precision) :: acc_rot_x,acc_rot_y,acc_rot_z
    real(double_precision) :: acc_GR_x,acc_GR_y,acc_GR_z
    real(double_precision) :: rr,r_2,r_4,r_5,r_7,r_8,v_2,norm_v,v_rad
    real(double_precision) :: N_tid_px,N_tid_py,N_tid_pz,N_tid_sx,N_tid_sy,N_tid_sz
    real(double_precision) :: N_rot_px,N_rot_py,N_rot_pz,N_rot_sx,N_rot_sy,N_rot_sz
    real(double_precision) :: k_rk_1x,k_rk_2x,k_rk_3x,k_rk_4x,k_rk_5x,k_rk_6x
    real(double_precision) :: k_rk_1y,k_rk_2y,k_rk_3y,k_rk_4y,k_rk_5y,k_rk_6y
    real(double_precision) :: k_rk_1z,k_rk_2z,k_rk_3z,k_rk_4z,k_rk_5z,k_rk_6z
    real(double_precision) :: xintermediate
    real(double_precision) :: dt,hdt,tstop,tmp,tmp1,tmp2,sigmast
    real(double_precision), dimension(2) :: bobo
    real(double_precision), dimension(3) :: totftides
    real(double_precision), dimension(3,nbig+1) :: a1,a2,a3,xh,vh
    real(double_precision), dimension(3,nbig+1) :: horb
    real(double_precision), dimension(ntid+1) :: qa,ea,ia,pa,na,la
    real(double_precision), dimension(3,nbig+1) :: Nts,Ntp,Nrs,Nrp,Ns,Np
    real(double_precision), dimension(3,10) :: spin,spin_bf,xh_bf,vh_bf,xh_bf2,vh_bf2
    real(double_precision), dimension(3,nbig+1) :: xh_1_rk2,vh_1_rk2,xh_2_rk2,vh_2_rk2
    real(double_precision), dimension(3,nbig+1) :: xh_1_rk3,vh_1_rk3,xh_2_rk3,vh_2_rk3
    real(double_precision), dimension(3,nbig+1) :: xh_1_rk4,vh_1_rk4,xh_2_rk4,vh_2_rk4
    real(double_precision), dimension(3,nbig+1) :: xh_1_rk5,vh_1_rk5,xh_2_rk5,vh_2_rk5
    real(double_precision), dimension(3,nbig+1) :: xh_1_rk6,vh_1_rk6,xh_2_rk6,vh_2_rk6
    real(double_precision), dimension(10) :: Rp,sigmap,Rp5,Rp10,tintin,k2p,k2pdeltap,rg2p
    real(double_precision), dimension(10) :: rscalws,rscalwp,normspin2
    ! don't use after collision
    real(double_precision), dimension(nbig+1) :: r,r2,r4,r5,r7,r8,v2,vv,vrad
    real(double_precision), dimension(nbig+1) :: horbn
    real(double_precision), dimension(4161) :: timeBD,radiusBD,lumiBD,HZinGJ,HZoutGJ,HZinb,HZoutb
    real(double_precision), dimension(1065) :: timestar,radiusstar,d2radiusstar
    real(double_precision), dimension(1065) :: timedM,radiusdM
    real(double_precision), dimension(4755) :: timeJup,radiusJup,k2Jup,rg2Jup,spinJup
    real(double_precision), dimension(37) :: rg2st,trg2,rg1,rg2,rg3,rg4,rg5,rg6,rg7,rg8,rg9,rg10,rg11,rg12
    real(double_precision) :: spin0,spinp0,spinb0
    real(double_precision) :: Rst,Rst_5,Rst_10,Rstb
    real(double_precision) :: Rsth,Rsth5,Rsth10,Rstbh
    real(double_precision) :: Rst0,Rst0_5,Rst0_10,Rstb0
    real(double_precision) :: rg2s,rg2sh,rg2s0,k2s
    
    real(double_precision), dimension(nbig+1) :: dEdt
    
    real(double_precision), parameter, dimension(12) :: Ps0 = (/8.d0,13.d0,19.d0,24.d0,30.d0,36.d0,41.d0, &
         47.d0,53.d0,58.d0,64.d0,70.d0/)
    real(double_precision), parameter, dimension(12) :: k2st = (/0.379d0,0.378d0,0.376d0,0.369d0, &
         0.355d0,0.342d0,0.333d0,0.325d0,0.311d0,0.308d0,0.307d0,0.307d0/)

    character(len=80) :: planet_spin_filename
    character(len=80) :: planet_orbt_filename
    character(len=80) :: planet_dEdt_filename

    ! Data 
    save timestar,radiusstar,d2radiusstar
    save timeBD,radiusBD
    save trg2,rg1,rg2,rg3,rg4,rg5,rg6,rg7,rg8,rg9,rg10,rg11,rg12
    save timedM,radiusdM
    save timeJup,radiusJup,k2Jup,rg2Jup,spinJup
    
    save xh_bf,vh_bf,xh_bf2,vh_bf2
    save sigmast,k2s
    save flagrg2,flagtime,ispin,flagbug
    save timestep,nptmss,Rst0,Rst0_5,Rst0_10
    save spin_bf,dt,hdt
    save Rp,sigmap,Rp5,Rp10,tintin,k2p,k2pdeltap,rg2p

    !------------------------------------------------------------------------------
    ! superzoyotte

    ! Error message
    if (ispin.eq.0) then
       if ((brown_dwarf.eq.1).and.((M_dwarf.eq.1).or.(Sun_like_star.eq.1).or.(Jupiter_host.eq.1).or.(Rscst.eq.1))) then 
          write(*,*) "You're trying to have a host body of two types at the same time, this is not going to work"
          stop
       endif
       if ((M_dwarf.eq.1).and.((Sun_like_star.eq.1).or.(Jupiter_host.eq.1).or.(Rscst.eq.1).or.(brown_dwarf.eq.1))) then
          write(*,*) "You're trying to have a host body of two types at the same time, this is not going to work"
          stop
       endif
       if ((Sun_like_star.eq.1).and.((Jupiter_host.eq.1).or.(Rscst.eq.1).or.(brown_dwarf.eq.1).or.(M_dwarf.eq.1))) then 
          write(*,*) "You're trying to have a host body of two types at the same time, this is not going to work"
          stop
       endif
       if ((Jupiter_host.eq.1).and.((Rscst.eq.1).or.(brown_dwarf.eq.1).or.(M_dwarf.eq.1).or.(Sun_like_star.eq.1))) then 
          write(*,*) "You're trying to have a host body of two types at the same time, this is not going to work"
          stop
       endif
       if ((Rscst.eq.1).and.((brown_dwarf.eq.1).or.(M_dwarf.eq.1).or.(Sun_like_star.eq.1).or.(Jupiter_host.eq.1))) then 
          write(*,*) "You're trying to have a host body of two types at the same time, this is not going to work"
          stop
       endif
    endif

    ! Acceleration initialization
    do j = 1, nbod
       a(1,j) = 0.d0
       a(2,j) = 0.d0
       a(3,j) = 0.d0
    end do
    do j=2,ntid+1
       a1(1,j) = 0.d0
       a1(2,j) = 0.d0
       a1(3,j) = 0.d0
       a2(1,j) = 0.d0
       a2(2,j) = 0.d0
       a2(3,j) = 0.d0
       a3(1,j) = 0.d0
       a3(2,j) = 0.d0
       a3(3,j) = 0.d0
       if (ispin.eq.0) then
          qq = 0.d0
          ee = 0.d0
          pp = 0.d0
          ii = 0.d0
          nn = 0.d0
          ll = 0.d0
          qa(j) = 0.d0
          ea(j) = 0.d0
          pa(j) = 0.d0
          ia(j) = 0.d0
          na(j) = 0.d0
          la(j) = 0.d0
          timestep = time
       endif
    end do
    if (iwrite.eq.0) then 
       call write_simus_properties()
       iwrite = 1
    endif

    ! Timestep calculation
    if (flagtime.eq.0) then 
       bobo = get_initial_timestep()
       tstop = bobo(1)
       dt = bobo(2)
       hdt = 0.5d0*dt
       flagtime = flagtime+1
    endif

    ! Following calculations in heliocentric coordinates   
    call conversion_dh2h(nbod,nbig,m,x,v,xh,vh)    
    if (ispin.eq.0) then
       do j=2,ntid+1
          xh_bf(1,j) = xh(1,j)
          xh_bf(2,j) = xh(2,j)
          xh_bf(3,j) = xh(3,j)
          vh_bf(1,j) = vh(1,j)
          vh_bf(2,j) = vh(2,j)
          vh_bf(3,j) = vh(3,j)
       enddo
    endif      
          
    if ((tides.eq.1).or.(rot_flat.eq.1)) then 
       ! Charge host body data 
       ! + initial condition on host body spin, radius
       if (brown_dwarf.eq.1) then 
          if (flagrg2.eq.0) then    
             sigmast = sigma_BD
             ! Radius of BD's gyration data
             open(1,file='rg2BD.dat')
             do nptmss=1,37
                read(1,*,iostat=error)trg2(nptmss),rg1(nptmss),rg2(nptmss), &
                     rg3(nptmss),rg4(nptmss),rg5(nptmss),rg6(nptmss),rg7(nptmss), &
                     rg8(nptmss),rg9(nptmss),rg10(nptmss),rg11(nptmss),rg12(nptmss)
             end do
             
             ! If BD's mass is equal to one of these values, charge radius... 
             if ((m(1).le.0.0101*K2).and.(m(1).ge.0.0099*K2)) then
                iPs0 = 1
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg1(j)
                enddo
                open(1,file='mass_10.0000.dat')
                do nptmss = 1,715
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0121*K2).and.(m(1).ge.0.0119*K2)) then
                iPs0 = 2
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg2(j)
                enddo
                open(1,file='mass_12.0000.dat')
                do nptmss=1,720
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0151*K2).and.(m(1).ge.0.0149*K2)) then
                iPs0 = 3
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg3(j)
                enddo
                open(1,file='mass_15.0000.dat')
                do nptmss=1,856
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0201*K2).and.(m(1).ge.0.0199*K2)) then
                iPs0 = 4
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg4(j)
                enddo
                open(1,file='mass_20.0000.dat')
                do nptmss =1,864
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0301*K2).and.(m(1).ge.0.0299*K2)) then
                iPs0 = 5
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg5(j)
                enddo
                open(1,file='mass_30.0000.dat')
                do nptmss=1,878
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0401*K2).and.(m(1).ge.0.0399*K2)) then
                iPs0 = 6
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg6(j)
                enddo
                open(1,file='mass_40.0000.dat')
                do nptmss = 1,886
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)  
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0501*K2).and.(m(1).ge.0.0499*K2)) then
                iPs0 = 7
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg7(j)
                enddo
                open(1,file='mass_50.0000.dat')
                do nptmss=1,891
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0601*K2).and.(m(1).ge.0.0599*K2)) then
                iPs0 = 8
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg8(j)
                enddo
                open(1,file='mass_60.0000.dat')
                do nptmss=1,1663
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0701*K2).and.(m(1).ge.0.0699*K2)) then
                iPs0 = 9
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg9(j)
                enddo
                open(1,file='mass_70.0000.dat')
                do nptmss =1,3585
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0721*K2).and.(m(1).ge.0.0719*K2)) then
                iPs0 = 10
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg10(j)
                enddo
                open(1,file='mass_72.0000.dat')
                do nptmss =1,3721
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0751*K2).and.(m(1).ge.0.0749*K2)) then
                iPs0 = 11
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg11(j)
                enddo
                open(1,file='mass_75.0000.dat')
                do nptmss =1,3903
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             else if ((m(1).le.0.0801*K2).and.(m(1).ge.0.0799*K2)) then
                iPs0 = 12
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg12(j)
                enddo
                open(1,file='mass_80.0000.dat')
                do nptmss =1,4161
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
                nptmss = nptmss - 1
             endif

             ! Initialization of stellar spin (day-1)
             if (crash.eq.0) then
                spin0 = 24.d0*TWOPI/Pst0 
                spin(1,1) = 0.d0
                spin(2,1) = 0.d0
                spin(3,1) = spin0
             else
                spin(1,1) = rot_crash(1) 
                spin(2,1) = rot_crash(2) 
                spin(3,1) = rot_crash(3) 
             endif
             flagrg2=1
          endif
       endif
       if (M_dwarf.eq.1) then 
          rg2s = rg2_dM
          if (flagrg2.eq.0) then
             ! Charge file of radius of Mdwarf 
             open(1,file='01Msun.dat')
             do nptmss = 1,1065
                read(1,*,iostat=error)timedM(nptmss),radiusdM(nptmss)
             end do
             nptmss = nptmss - 1
             if (crash.eq.0) then
                Pst = Period_st
                ! Initialization of stellar spin (day-1)
                spin0 = TWOPI/Pst 
                spin(1,1) = 0.d0
                spin(2,1) = 0.d0
                spin(3,1) = spin0
             else
                spin(1,1) = rot_crash(1)
                spin(2,1) = rot_crash(2)
                spin(3,1) = rot_crash(3)
             end if
             k2s   = k2st_dM
             sigmast = sigma_dM
             flagrg2 = 1
          endif
       endif
       if (Sun_like_star.eq.1) then 
          rg2s = rg2_Sun
          if (flagrg2.eq.0) then 
             ! Charge file of radius of Mdwarf 
             open(1,file='SRad_Spli_M-1_0000.dat')
             do nptmss = 1,2003
                read(1,*,iostat=error)timestar(nptmss),radiusstar(nptmss),d2radiusstar(nptmss)
             end do
             nptmss = nptmss - 1
             ! Initialization Sun-like star spin (day-1)
             if (crash.eq.0) then
                Pst = Period_st
                spin0 = TWOPI/Pst !day-1
                spin(1,1) = 0.d0
                spin(2,1) = 0.d0
                spin(3,1) = spin0
             else
                spin(1,1) = rot_crash(1) !day-1
                spin(2,1) = rot_crash(2) !day-1
                spin(3,1) = rot_crash(3) !day-1
             end if
             k2s   = k2st_Sun
             sigmast = sigma_Sun
             flagrg2 = 1
          endif
       endif
       if (Jupiter_host.eq.1) then 
          if (flagrg2.eq.0) then 
             ! Charge file of radius of Jupiter 
             open(1,file='Jupiter.dat')
             do nptmss = 1,4755
                read(1,*,iostat=error) timeJup(nptmss),radiusJup(nptmss) &
                     ,k2Jup(nptmss),rg2Jup(nptmss),spinJup(nptmss)
             end do
             nptmss = nptmss - 1
             if (crash.eq.0) then
                call spline_b_val(nptmss,timeJup*365.25-t_init,radiusJup,time,Rstb0)
                Rst0    =  minau * Rstb0
                ! spinJup is s-1
                call spline_b_val(nptmss,timeJup*365.25-t_init,spinJup,time,spinb0)
                ! Initialization of stellar spin (day-1)
                spin0 = spinb0*86400.d0
                spin(1,1) = 0.d0
                spin(2,1) = 0.d0
                spin(3,1) = spin0
             else
                spin(1,1) = rot_crash(1)
                spin(2,1) = rot_crash(2)
                spin(3,1) = rot_crash(3)
                call spline_b_val(nptmss,timeJup*365.25-t_init,radiusJup,0.0d0,Rstb0)
                Rst0    = minau * Rstb0
             end if
             ! Give good value of sigmast
             sigmast = dissstar*2.d0*K2*k2delta_jup/(3.d0*Rst0*Rst0*Rst0*Rst0*Rst0)
             flagrg2 = 1
          endif
       endif
       if (Rscst.eq.1) then 
          rg2s   = rg2_what
          Rsth   = radius_star*rsun
          Rst    = Rsth
          Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
          Rst_5  = Rsth5
          Rst0_5 = Rsth5
          Rsth10 = Rsth5*Rsth5
          Rst_10 = Rsth10
          Rst0_10= Rsth10
          if (flagrg2.eq.0) then
             Pst = Period_st
             k2s   = k2st_what
             sigmast = sigma_what
             if (crash.eq.0) then
                ! Initialization of stellar spin (day-1)
                spin0 = TWOPI/Pst !day-1
                spin(1,1) = 0.d0
                spin(2,1) = 0.d0
                spin(3,1) = spin0
             else
                spin(1,1) = rot_crash(1) !day-1
                spin(2,1) = rot_crash(2) !day-1
                spin(3,1) = rot_crash(3) !day-1
             endif
             flagrg2 = 1
          endif
       endif
    endif
    
    ! Calculation of orbital angular momentum
    ! horb (without mass) in AU^2/day
    do j=2,ntid+1
       horb(1,j)  = (xh(2,j)*vh(3,j)-xh(3,j)*vh(2,j))
       horb(2,j)  = (xh(3,j)*vh(1,j)-xh(1,j)*vh(3,j))
       horb(3,j)  = (xh(1,j)*vh(2,j)-xh(2,j)*vh(1,j))    
       horbn(j) = sqrt(horb(1,j)*horb(1,j)+horb(2,j)*horb(2,j)+horb(3,j)*horb(3,j))
    end do

    ! Initialization of many things concerning planets
    if (ispin.eq.0) then 
       do j=2,ntid+1
          if ((tides.eq.1).or.(rot_flat.eq.1)) then 
             ! Calculate the planets properties   
             ! Planetary radius in AU (rearth in AU) Rocky planet
             if (jupiter(j-1).eq.0) Rp(j) = rearth*((0.0592d0*0.7d0+0.0975d0) &
                  *(dlog10(m(j))+dlog10(m2earth)-dlog10(K2))**2+(0.2337d0*0.7d0+0.4938d0) &
                  *(dlog10(m(j))+dlog10(m2earth)-dlog10(K2))+0.3102d0*0.7d0+0.7932d0)
             if (jupiter(j-1).ne.0) Rp(j) = radius_p(j-1)*rearth
             Rp5(j)  = Rp(j)*Rp(j)*Rp(j)*Rp(j)*Rp(j)
             Rp10(j) = Rp5(j)*Rp5(j)
             ! k2p, rg2p, k2pdeltap and sigmap
             ! Terrestrial for 0 and 1, gas giant for 2, what you want for 3
             if ((jupiter(j-1).eq.0).or.(jupiter(j-1).eq.1)) then
                k2p(j-1) = k2p_terr
                rg2p(j-1) = rg2p_terr
                if (tides.eq.1) then
                   k2pdeltap(j-1) = k2pdeltap_terr
                   sigmap(j) = dissplan(j-1)*2.d0*K2*k2pdeltap(j-1)/(3.d0*Rp5(j))
                endif
             endif
             if (jupiter(j-1).eq.2) then
                k2p(j-1) = k2p_gg
                rg2p(j-1) = rg2p_gg
                if (tides.eq.1) then
                   k2pdeltap(j-1) = k2pdeltap_gg
                   sigmap(j) = dissplan(j-1)*sigma_gg
                endif
             endif
             if (jupiter(j-1).eq.3) then
                k2p(j-1) = k2p_what
                rg2p(j-1) = rg2p_what
                if (tides.eq.1) then
                   k2pdeltap(j-1) = k2pdeltap_what
                   sigmap(j) = dissplan(j-1)*2.d0*K2*k2pdeltap(j-1)/(3.d0*Rp5(j))
                endif
             endif
          endif
          ! Factor used for GR force calculation
          if (GenRel.eq.1) tintin(j) = m(1)*m(j)/(m(1)+m(j))**2
       enddo
       if ((tides.eq.1).or.(rot_flat.eq.1)) then 
          ! Initialization of spin planets
          if (crash.eq.0) then
             do j=2,ntid+1
                gm = m(1) + m(j)
                ! gm is in AU^3.day^-2
                call mco_x2el(gm,xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),qq,ee,ii,pp,nn,ll)
                qa(j) = qq
                ea(j) = ee
                pa(j) = pp
                ia(j) = ii
                na(j) = nn 
                la(j) = ll
                ! Initialization of planetary spin (day-1)
                if (pseudo_rot(j-1).eq.0) spinp0 = 24.d0*TWOPI/Pp0(j-1)
                if (pseudo_rot(j-1).ne.0) then 
                   spinp0 = pseudo_rot(j-1)*(1.d0+15.d0/2.d0*ea(j)**2+45.d0/8.d0*ea(j)**4+5.d0/16.d0*ea(j)**6) &
                        *1.d0/(1.d0+3.d0*ea(j)**2+3.d0/8.d0*ea(j)**4)*1./(1-ea(j)**2)**1.5d0*sqrt(m(1)+m(j)) &
                        *(qa(j)/(1.d0-ea(j)))**(-1.5d0)
                endif
                if (ia(j).eq.0.0) then
                   spin(1,j) = spinp0*sin(oblp(j-1))
                   spin(2,j) = 0.0d0
                   spin(3,j) = spinp0*cos(oblp(j-1)) 
                endif
                if (ia(j).ne.0.0) then
                   spin(1,j) = spinp0*(horb(1,j)/(horbn(j)*sin(ia(j))))*sin(oblp(j-1)+ia(j))
                   spin(2,j) = spinp0*(horb(2,j)/(horbn(j)*sin(ia(j))))*sin(oblp(j-1)+ia(j))
                   spin(3,j) = spinp0*cos(oblp(j-1)+ia(j)) 
                endif
             enddo
          else
             spin(1,1) = rot_crash(1) !day-1
             spin(2,1) = rot_crash(2) !day-1
             spin(3,1) = rot_crash(3) !day-1
             spin(1,2) = rot_crashp1(1) !day-1
             spin(2,2) = rot_crashp1(2) !day-1
             spin(3,2) = rot_crashp1(3) !day-1 
             if (ntid.ge.2) then
                spin(1,3) = rot_crashp2(1) !day-1
                spin(2,3) = rot_crashp2(2) !day-1
                spin(3,3) = rot_crashp2(3) !day-1
             endif
             if (ntid.ge.3) then
                spin(1,4) = rot_crashp3(1) !day-1
                spin(2,4) = rot_crashp3(2) !day-1
                spin(3,4) = rot_crashp3(3) !day-1
             endif
             if (ntid.ge.4) then
                spin(1,5) = rot_crashp4(1) !day-1
                spin(2,5) = rot_crashp4(2) !day-1
                spin(3,5) = rot_crashp4(3) !day-1
             endif
             if (ntid.ge.5) then
                spin(1,6) = rot_crashp5(1) !day-1
                spin(2,6) = rot_crashp5(2) !day-1
                spin(3,6) = rot_crashp5(3) !day-1
             endif
             if (ntid.ge.6) then
                spin(1,7) = rot_crashp6(1) !day-1
                spin(2,7) = rot_crashp6(2) !day-1
                spin(3,7) = rot_crashp6(3) !day-1
             endif
          endif
       endif
    endif

    if ((tides.eq.1).or.(rot_flat.eq.1)) then 
       ! Interpolation to have the host body radius (and radius of gyration for BDs)
       ! Here Rst in AU, Rsth5; Rsth10
       if (brown_dwarf.eq.1) then
          if (crash.eq.0) then
             if (ispin.eq.0) then 
                call spline_b_val(nptmss,timeBD*365.25d0-t_init,radiusBD,time-dt,Rstb0)
                Rst0    = Rsun * Rstb0
                Rst0_5  = Rst0*Rst0*Rst0*Rst0*Rst0
                Rst0_10 = Rst0_5*Rst0_5
             endif
             call spline_b_val(nptmss,timeBD*365.25d0-t_init,radiusBD,time,Rstb)
             Rst    = Rsun * Rstb
             Rst_5  = Rst*Rst*Rst*Rst*Rst
             Rst_10 = Rst_5*Rst_5
             call spline_b_val(nptmss,timeBD*365.25d0-t_init,radiusBD,time-dt*0.5d0,Rstbh)
             Rsth = Rsun * Rstbh
             Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
             Rsth10 = Rsth5*Rsth5             
             call spline_b_val(37,trg2*365.25d0-t_init,rg2st,time-dt,rg2s0)
             call spline_b_val(37,trg2*365.25d0-t_init,rg2st,time,rg2s)
             call spline_b_val(37,trg2*365.25d0-t_init,rg2st,time-dt*0.5d0,rg2sh)
          else
             if (ispin.eq.0) then 
                call spline_b_val(nptmss,timeBD*365.25d0-t_init-t_crash,radiusBD,time-t_crash-dt,Rstb0)
                Rst0    = Rsun * Rstb0
                Rst0_5  = Rst0*Rst0*Rst0*Rst0*Rst0
                Rst0_10 = Rst0_5*Rst0_5
             endif
             call spline_b_val(nptmss,timeBD*365.25d0-t_init-t_crash,radiusBD,time-t_crash,Rstb)
             Rst    = Rsun * Rstb
             Rst_5  = Rst*Rst*Rst*Rst*Rst
             Rst_10 = Rst_5*Rst_5
             call spline_b_val(nptmss,timeBD*365.25d0-t_init-t_crash,radiusBD,time-t_crash-dt*0.5d0,Rstbh)
             Rsth   = Rsun * Rstbh
             Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
             Rsth10 = Rsth5*Rsth5             
             call spline_b_val(37,trg2*365.25d0-t_init-t_crash,rg2st,time-t_crash-dt,rg2s0)
             call spline_b_val(37,trg2*365.25d0-t_init-t_crash,rg2st,time-t_crash,rg2s)
             call spline_b_val(37,trg2*365.25d0-t_init-t_crash,rg2st,time-t_crash-dt*0.5d0,rg2sh)
          endif
       endif
       if (M_dwarf.eq.1) then
          if (crash.eq.0) then
             if (ispin.eq.0) then 
                call spline_b_val(nptmss,timedM*365.25-t_init,radiusdM,time-dt,Rstb0)
                Rst0    = Rsun * Rstb0
                Rst0_5  = Rst0*Rst0*Rst0*Rst0*Rst0
                Rst0_10 = Rst0_5*Rst0_5
             endif
             call spline_b_val(nptmss,timedM*365.25-t_init,radiusdM,time,Rstb)
             Rst    = Rsun * Rstb
             Rst_5  = Rst*Rst*Rst*Rst*Rst
             Rst_10 = Rst_5*Rst_5
             call spline_b_val(nptmss,timedM*365.25-t_init,radiusdM,time-dt*0.5d0,Rstbh)
             Rsth   = Rsun * Rstbh
             Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
             Rsth10 = Rsth5*Rsth5
          else
             if (ispin.eq.0) then 
                call spline_b_val(nptmss,timedM*365.25-t_init-t_crash,radiusdM,time-t_crash-dt,Rstb0)
                Rst0    = Rsun * Rstb0
                Rst0_5  = Rst0*Rst0*Rst0*Rst0*Rst0
                Rst0_10 = Rst0_5*Rst0_5
             endif
             call spline_b_val(nptmss,timedM*365.25d0-t_init-t_crash,radiusdM,time-t_crash,Rstb)
             Rst    = Rsun * Rstb
             Rst_5  = Rst*Rst*Rst*Rst*Rst
             Rst_10 = Rst_5*Rst_5
             call spline_b_val(nptmss,timedM*365.25d0-t_init-t_crash,radiusdM,time-t_crash-dt*0.5d0,Rstbh)
             Rsth = Rsun * Rstbh
             Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
             Rsth10 = Rsth5*Rsth5
          endif
       endif
       if (Sun_like_star.eq.1) then
          if (crash.eq.0) then
             if (ispin.eq.0) then 
                call spline_b_val(nptmss,timestar*365.25-t_init,radiusstar,time-dt,Rstb0)
                Rst0    = minau * Rstb0
                Rst0_5  = Rst0*Rst0*Rst0*Rst0*Rst0
                Rst0_10 = Rst0_5*Rst0_5
             endif
             call spline_b_val(nptmss,timestar*365.25-t_init,radiusstar,time,Rstb)
             Rst    = minau * Rstb
             Rst_5  = Rst*Rst*Rst*Rst*Rst
             Rst_10 = Rst_5*Rst_5
             call spline_b_val(nptmss,timestar*365.25-t_init,radiusstar,time-dt*0.5d0,Rstbh)
             Rsth   = minau * Rstbh
             Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
             Rsth10 = Rsth5*Rsth5
          else
             if (ispin.eq.0) then 
                call spline_b_val(nptmss,timestar*365.25-t_init-t_crash,radiusstar,time-t_crash-dt,Rstb0)
                Rst0    = minau * Rstb0
                Rst0_5  = Rst0*Rst0*Rst0*Rst0*Rst0
                Rst0_10 = Rst0_5*Rst0_5
             endif
             call spline_b_val(nptmss,timestar*365.25-t_init-t_crash,radiusstar,time-t_crash,Rstb)
             Rst    = minau * Rstb
             Rst_5  = Rst*Rst*Rst*Rst*Rst
             Rst_10 = Rst_5*Rst_5
             call spline_b_val(nptmss,timestar*365.25-t_init-t_crash,radiusstar,time-t_crash-dt*0.5d0,Rstbh)
             Rsth   = minau * Rstbh
             Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
             Rsth10 = Rsth5*Rsth5
          endif
       endif
       if (Jupiter_host.eq.1) then
          if (crash.eq.0) then
             if (ispin.eq.0) then 
                call spline_b_val(nptmss,timeJup*365.25d0-t_init,radiusJup,time-dt,Rstb0)
                Rst0    = minau * Rstb0
                Rst0_5  = Rst0*Rst0*Rst0*Rst0*Rst0
                Rst0_10 = Rst0_5*Rst0_5
             endif
             call spline_b_val(nptmss,timeJup*365.25d0-t_init,radiusJup,time,Rstb)
             Rst    = minau * Rstb
             Rst_5  = Rst*Rst*Rst*Rst*Rst
             Rst_10 = Rst_5*Rst_5
             call spline_b_val(nptmss,timeJup*365.25d0-t_init,radiusJup,time-dt*0.5d0,Rstbh)
             Rsth   = minau * Rstbh
             Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
             Rsth10 = Rsth5*Rsth5
             call spline_b_val(nptmss,timeJup*365.25d0-t_init,rg2Jup,time-dt,rg2s0)
             call spline_b_val(nptmss,timeJup*365.25d0-t_init,rg2Jup,time,rg2s)
             call spline_b_val(nptmss,timeJup*365.25d0-t_init,rg2Jup,time-dt*0.5d0,rg2sh)
             call spline_b_val(nptmss,timeJup*365.25d0-t_init,k2Jup,time-dt*0.5d0,k2s)
          else
             if (ispin.eq.0) then 
                call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,radiusJup,time-t_crash-dt,Rstb0)
                Rst0    = minau * Rstb0
                Rst0_5  = Rst0*Rst0*Rst0*Rst0*Rst0
                Rst0_10 = Rst0_5*Rst0_5
             endif
             call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,radiusJup,time-t_crash,Rstb)
             Rst    = minau * Rstb
             Rst_5  = Rst*Rst*Rst*Rst*Rst
             Rst_10 = Rst_5*Rst_5
             call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,radiusJup,time-t_crash-dt*0.5d0,Rstbh)
             Rsth   = minau * Rstbh
             Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
             Rsth10 = Rsth5*Rsth5             
             call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,rg2Jup,time-t_crash-dt,rg2s0)
             call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,rg2Jup,time-t_crash,rg2s)
             call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,rg2Jup,time-t_crash-dt*0.5d0,rg2sh)
             call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,k2Jup,time-t_crash-dt*0.5d0,k2s)
          endif
       endif
    endif

    if ((time.ne.0).and.(flagbug.ge.2)) then

       ! **************************************************************
       ! **************  Runge Kutta integration  *********************
       ! ********************** of spin *******************************
       ! **************************************************************
       
       if ((tides.eq.1).or.(rot_flat.eq.1)) then       
          do j=2,ntid+1
             do kk=1,3
!~                 ! interpolation of value of x at time-dt/2 
!~                 call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
!~                      ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-dt*0.5d0,xintermediate)
!~                 xh_int(kk,j) = xintermediate
!~                 ! interpolation of value of v at time-dt/2 
!~                 call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
!~                      ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-dt*0.5d0,xintermediate)
!~                 vh_int(kk,j) = xintermediate
                ! Interpolation for Runge-Kutta steps on first half of timestep
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-(2.d0-aa(1))*dt*0.5d0,xintermediate)
                xh_1_rk2(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-(2.d0-aa(1))*dt*0.5d0,xintermediate)
                vh_1_rk2(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-(2.d0-aa(2))*dt*0.5d0,xintermediate)
                xh_1_rk3(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-(2.d0-aa(2))*dt*0.5d0,xintermediate)
                vh_1_rk3(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-(2.d0-aa(3))*dt*0.5d0,xintermediate)
                xh_1_rk4(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-(2.d0-aa(3))*dt*0.5d0,xintermediate)
                vh_1_rk4(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-dt*0.5d0,xintermediate)
                xh_1_rk5(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-dt*0.5d0,xintermediate)
                vh_1_rk5(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-(2.d0-aa(5))*dt*0.5d0,xintermediate)
                xh_1_rk6(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-(2.d0-aa(5))*dt*0.5d0,xintermediate)
                vh_1_rk6(kk,j) = xintermediate
                ! Interpolation for Runge-Kutta steps on second half of timestep
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-(1.d0-aa(1))*dt*0.5d0,xintermediate)
                xh_2_rk2(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-(1.d0-aa(1))*dt*0.5d0,xintermediate)
                vh_2_rk2(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-(1.d0-aa(2))*dt*0.5d0,xintermediate)
                xh_2_rk3(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-(1.d0-aa(2))*dt*0.5d0,xintermediate)
                vh_2_rk3(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-(1.d0-aa(3))*dt*0.5d0,xintermediate)
                xh_2_rk4(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-(1.d0-aa(3))*dt*0.5d0,xintermediate)
                vh_2_rk4(kk,j) = xintermediate
                xh_2_rk5(kk,j) = xh(kk,j)
                vh_2_rk5(kk,j) = vh(kk,j)
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/xh_bf2(kk,j),xh_bf(kk,j),xh(kk,j)/),time-(1.d0-aa(5))*dt*0.5d0,xintermediate)
                xh_2_rk6(kk,j) = xintermediate
                call spline_b_val(3,(/time-2.d0*dt,time-dt,time/) &
                     ,(/vh_bf2(kk,j),vh_bf(kk,j),vh(kk,j)/),time-(1.d0-aa(5))*dt*0.5d0,xintermediate)
                vh_2_rk6(kk,j) = xintermediate
             enddo

             
          enddo   
       
          ! **************************************************************
          ! STAR spin evolution
          ! **************************************************************

          ! First half of timestep (zoyotte)
          ! Torque at last time step t=time-dt
          ! Calculation of Runge-kutta factor k1
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_bf(1,j),xh_bf(2,j),xh_bf(3,j) &
                   ,vh_bf(1,j),vh_bf(2,j),vh_bf(3,j) &
                   ,spin_bf(1,1),spin_bf(2,1),spin_bf(3,1) &
                   ,Rst0_10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_bf(1,j),xh_bf(2,j),xh_bf(3,j) &
                     ,spin_bf(1,1),spin_bf(2,1),spin_bf(3,1) &
                     ,Rst0_5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_1x = tmp*totftides(1)
             k_rk_1y = tmp*totftides(2)
             k_rk_1z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rst0*Rst0)
             k_rk_1x = tmp*totftides(1)
             k_rk_1y = tmp*totftides(2)
             k_rk_1z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2s0*Rst0*Rst0)
             k_rk_1x = tmp*totftides(1)
             k_rk_1y = tmp*totftides(2)
             k_rk_1z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(2-a(2)) 
          ! Calculation of Runge-kutta factor k2
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_1_rk2(1,j),xh_1_rk2(2,j),xh_1_rk2(3,j) &
                     ,vh_1_rk2(1,j),vh_1_rk2(2,j),vh_1_rk2(3,j) &
                     ,spin_bf(1,1) + bb2*k_rk_1x &
                     ,spin_bf(2,1) + bb2*k_rk_1y &
                     ,spin_bf(3,1) + bb2*k_rk_1z &
                     ,Rst0_10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_1_rk2(1,j),xh_1_rk2(2,j),xh_1_rk2(3,j) &
                     ,spin_bf(1,1) + bb2*k_rk_1x &
                     ,spin_bf(2,1) + bb2*k_rk_1y &
                     ,spin_bf(3,1) + bb2*k_rk_1z &
                     ,Rst0_5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_2x = tmp*totftides(1)
             k_rk_2y = tmp*totftides(2)
             k_rk_2z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rst0*Rst0)
             k_rk_2x = tmp*totftides(1)
             k_rk_2y = tmp*totftides(2)
             k_rk_2z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2s0*Rst0*Rst0)
             k_rk_2x = tmp*totftides(1)
             k_rk_2y = tmp*totftides(2)
             k_rk_2z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(2-a(3)) 
          ! Calculation of Runge-kutta factor k3
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_1_rk3(1,j),xh_1_rk3(2,j),xh_1_rk3(3,j) &
                     ,vh_1_rk3(1,j),vh_1_rk3(2,j),vh_1_rk3(3,j) &
                     ,spin_bf(1,1) + bb3(1)*k_rk_1x + bb3(2)*k_rk_2x &
                     ,spin_bf(2,1) + bb3(1)*k_rk_1y + bb3(2)*k_rk_2y &
                     ,spin_bf(3,1) + bb3(1)*k_rk_1z + bb3(2)*k_rk_2z &
                     ,Rst0_10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_1_rk3(1,j),xh_1_rk3(2,j),xh_1_rk3(3,j) &
                     ,spin_bf(1,1) + bb3(1)*k_rk_1x + bb3(2)*k_rk_2x &
                     ,spin_bf(2,1) + bb3(1)*k_rk_1y + bb3(2)*k_rk_2y &
                     ,spin_bf(3,1) + bb3(1)*k_rk_1z + bb3(2)*k_rk_2z &
                     ,Rst0_5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_3x = tmp*totftides(1)
             k_rk_3y = tmp*totftides(2)
             k_rk_3z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rst0*Rst0)
             k_rk_3x = tmp*totftides(1)
             k_rk_3y = tmp*totftides(2)
             k_rk_3z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2s0*Rst0*Rst0)
             k_rk_3x = tmp*totftides(1)
             k_rk_3y = tmp*totftides(2)
             k_rk_3z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(2-a(4)) 
          ! Calculation of Runge-kutta factor k4
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_1_rk4(1,j),xh_1_rk4(2,j),xh_1_rk4(3,j) &
                     ,vh_1_rk4(1,j),vh_1_rk4(2,j),vh_1_rk4(3,j) &
                     ,spin_bf(1,1) + bb4(1)*k_rk_1x + bb4(2)*k_rk_2x+ bb4(3)*k_rk_3x &
                     ,spin_bf(2,1) + bb4(1)*k_rk_1y + bb4(2)*k_rk_2y+ bb4(3)*k_rk_3y &
                     ,spin_bf(3,1) + bb4(1)*k_rk_1z + bb4(2)*k_rk_2z+ bb4(3)*k_rk_3z &
                     ,Rsth10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_1_rk4(1,j),xh_1_rk4(2,j),xh_1_rk4(3,j) &
                     ,spin_bf(1,1) + bb4(1)*k_rk_1x + bb4(2)*k_rk_2x+ bb4(3)*k_rk_3x &
                     ,spin_bf(2,1) + bb4(1)*k_rk_1y + bb4(2)*k_rk_2y+ bb4(3)*k_rk_3y &
                     ,spin_bf(3,1) + bb4(1)*k_rk_1z + bb4(2)*k_rk_2z+ bb4(3)*k_rk_3z &
                     ,Rsth5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_4x = tmp*totftides(1)
             k_rk_4y = tmp*totftides(2)
             k_rk_4z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rsth*Rsth)
             k_rk_4x = tmp*totftides(1)
             k_rk_4y = tmp*totftides(2)
             k_rk_4z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2sh*Rsth*Rsth)
             k_rk_4x = tmp*totftides(1)
             k_rk_4y = tmp*totftides(2)
             k_rk_4z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(2-a(5)) 
          ! Calculation of Runge-kutta factor k5
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_1_rk5(1,j),xh_1_rk5(2,j),xh_1_rk5(3,j) &
                     ,vh_1_rk5(1,j),vh_1_rk5(2,j),vh_1_rk5(3,j) &
                     ,spin_bf(1,1) + bb5(1)*k_rk_1x + bb5(2)*k_rk_2x+ bb5(3)*k_rk_3x+ bb5(4)*k_rk_4x &
                     ,spin_bf(2,1) + bb5(1)*k_rk_1y + bb5(2)*k_rk_2y+ bb5(3)*k_rk_3y+ bb5(4)*k_rk_4y &
                     ,spin_bf(3,1) + bb5(1)*k_rk_1z + bb5(2)*k_rk_2z+ bb5(3)*k_rk_3z+ bb5(4)*k_rk_4z &
                     ,Rsth10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_1_rk5(1,j),xh_1_rk5(2,j),xh_1_rk5(3,j) &
                     ,spin_bf(1,1) + bb5(1)*k_rk_1x + bb5(2)*k_rk_2x+ bb5(3)*k_rk_3x+ bb5(4)*k_rk_4x &
                     ,spin_bf(2,1) + bb5(1)*k_rk_1y + bb5(2)*k_rk_2y+ bb5(3)*k_rk_3y+ bb5(4)*k_rk_4y &
                     ,spin_bf(3,1) + bb5(1)*k_rk_1z + bb5(2)*k_rk_2z+ bb5(3)*k_rk_3z+ bb5(4)*k_rk_4z &
                     ,Rsth5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_5x = tmp*totftides(1)
             k_rk_5y = tmp*totftides(2)
             k_rk_5z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rsth*Rsth)
             k_rk_5x = tmp*totftides(1)
             k_rk_5y = tmp*totftides(2)
             k_rk_5z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2sh*Rsth*Rsth)
             k_rk_5x = tmp*totftides(1)
             k_rk_5y = tmp*totftides(2)
             k_rk_5z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(2-a(6)) 
          ! Calculation of Runge-kutta factor k6
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_1_rk6(1,j),xh_1_rk6(2,j),xh_1_rk6(3,j) &
                     ,vh_1_rk6(1,j),vh_1_rk6(2,j),vh_1_rk6(3,j) &
                     ,spin_bf(1,1) + bb6(1)*k_rk_1x + bb6(2)*k_rk_2x+ bb6(3)*k_rk_3x+ bb6(4)*k_rk_4x + bb6(5)*k_rk_5x &
                     ,spin_bf(2,1) + bb6(1)*k_rk_1y + bb6(2)*k_rk_2y+ bb6(3)*k_rk_3y+ bb6(4)*k_rk_4y + bb6(5)*k_rk_5y &
                     ,spin_bf(3,1) + bb6(1)*k_rk_1z + bb6(2)*k_rk_2z+ bb6(3)*k_rk_3z+ bb6(4)*k_rk_4z + bb6(5)*k_rk_5z &
                     ,Rsth10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_1_rk6(1,j),xh_1_rk6(2,j),xh_1_rk6(3,j) &
                     ,spin_bf(1,1) + bb6(1)*k_rk_1x + bb6(2)*k_rk_2x+ bb6(3)*k_rk_3x+ bb6(4)*k_rk_4x + bb6(5)*k_rk_5x &
                     ,spin_bf(2,1) + bb6(1)*k_rk_1y + bb6(2)*k_rk_2y+ bb6(3)*k_rk_3y+ bb6(4)*k_rk_4y + bb6(5)*k_rk_5y &
                     ,spin_bf(3,1) + bb6(1)*k_rk_1z + bb6(2)*k_rk_2z+ bb6(3)*k_rk_3z+ bb6(4)*k_rk_4z + bb6(5)*k_rk_5z &
                     ,Rsth5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_6x = tmp*totftides(1)
             k_rk_6y = tmp*totftides(2)
             k_rk_6z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rsth*Rsth)
             k_rk_6x = tmp*totftides(1)
             k_rk_6y = tmp*totftides(2)
             k_rk_6z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2sh*Rsth*Rsth)
             k_rk_6x = tmp*totftides(1)
             k_rk_6y = tmp*totftides(2)
             k_rk_6z = tmp*totftides(3)
          endif
          
          ! Integration on first half timestep
          spin(1,1) = spin_bf(1,1) + cc(1)*k_rk_1x + cc(2)*k_rk_2x + cc(3)*k_rk_3x + cc(4)*k_rk_4x + cc(5)*k_rk_5x + cc(6)*k_rk_6x
          spin(2,1) = spin_bf(2,1) + cc(1)*k_rk_1y + cc(2)*k_rk_2y + cc(3)*k_rk_3y + cc(4)*k_rk_4y + cc(5)*k_rk_5y + cc(6)*k_rk_6y
          spin(3,1) = spin_bf(3,1) + cc(1)*k_rk_1z + cc(2)*k_rk_2z + cc(3)*k_rk_3z + cc(4)*k_rk_4z + cc(5)*k_rk_5z + cc(6)*k_rk_6z
          

          ! **************************************************************
          ! Second half of timestep
          ! Torque at time t=time-dt/2
          ! Calculation of Runge-kutta factor k1
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_1_rk5(1,j),xh_1_rk5(2,j),xh_1_rk5(3,j) &
                   ,vh_1_rk5(1,j),vh_1_rk5(2,j),vh_1_rk5(3,j) &
                   ,spin(1,1),spin(2,1),spin(3,1) &
                   ,Rsth10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_1_rk5(1,j),xh_1_rk5(2,j),xh_1_rk5(3,j) &
                     ,spin(1,1),spin(2,1),spin(3,1) &
                     ,Rsth5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_1x = tmp*totftides(1)
             k_rk_1y = tmp*totftides(2)
             k_rk_1z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rsth*Rsth)
             k_rk_1x = tmp*totftides(1)
             k_rk_1y = tmp*totftides(2)
             k_rk_1z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2sh*Rsth*Rsth)
             k_rk_1x = tmp*totftides(1)
             k_rk_1y = tmp*totftides(2)
             k_rk_1z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(1-a(2)) 
          ! Calculation of Runge-kutta factor k2
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_2_rk2(1,j),xh_2_rk2(2,j),xh_2_rk2(3,j) &
                     ,vh_2_rk2(1,j),vh_2_rk2(2,j),vh_2_rk2(3,j) &
                     ,spin(1,1) + bb2*k_rk_1x &
                     ,spin(2,1) + bb2*k_rk_1y &
                     ,spin(3,1) + bb2*k_rk_1z &
                     ,Rsth10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_2_rk2(1,j),xh_2_rk2(2,j),xh_2_rk2(3,j) &
                     ,spin(1,1) + bb2*k_rk_1x &
                     ,spin(2,1) + bb2*k_rk_1y &
                     ,spin(3,1) + bb2*k_rk_1z &
                     ,Rsth5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_2x = tmp*totftides(1)
             k_rk_2y = tmp*totftides(2)
             k_rk_2z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rsth*Rsth)
             k_rk_2x = tmp*totftides(1)
             k_rk_2y = tmp*totftides(2)
             k_rk_2z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2sh*Rsth*Rsth)
             k_rk_2x = tmp*totftides(1)
             k_rk_2y = tmp*totftides(2)
             k_rk_2z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(1-a(3)) 
          ! Calculation of Runge-kutta factor k3
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_2_rk3(1,j),xh_2_rk3(2,j),xh_2_rk3(3,j) &
                     ,vh_2_rk3(1,j),vh_2_rk3(2,j),vh_2_rk3(3,j) &
                     ,spin(1,1) + bb3(1)*k_rk_1x + bb3(2)*k_rk_2x &
                     ,spin(2,1) + bb3(1)*k_rk_1y + bb3(2)*k_rk_2y &
                     ,spin(3,1) + bb3(1)*k_rk_1z + bb3(2)*k_rk_2z &
                     ,Rsth10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_2_rk3(1,j),xh_2_rk3(2,j),xh_2_rk3(3,j) &
                     ,spin(1,1) + bb3(1)*k_rk_1x + bb3(2)*k_rk_2x &
                     ,spin(2,1) + bb3(1)*k_rk_1y + bb3(2)*k_rk_2y &
                     ,spin(3,1) + bb3(1)*k_rk_1z + bb3(2)*k_rk_2z &
                     ,Rsth5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_3x = tmp*totftides(1)
             k_rk_3y = tmp*totftides(2)
             k_rk_3z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rsth*Rsth)
             k_rk_3x = tmp*totftides(1)
             k_rk_3y = tmp*totftides(2)
             k_rk_3z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2sh*Rsth*Rsth)
             k_rk_3x = tmp*totftides(1)
             k_rk_3y = tmp*totftides(2)
             k_rk_3z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(1-a(4)) 
          ! Calculation of Runge-kutta factor k4
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_2_rk4(1,j),xh_2_rk4(2,j),xh_2_rk4(3,j) &
                     ,vh_2_rk4(1,j),vh_2_rk4(2,j),vh_2_rk4(3,j) &
                     ,spin(1,1) + bb4(1)*k_rk_1x + bb4(2)*k_rk_2x+ bb4(3)*k_rk_3x &
                     ,spin(2,1) + bb4(1)*k_rk_1y + bb4(2)*k_rk_2y+ bb4(3)*k_rk_3y &
                     ,spin(3,1) + bb4(1)*k_rk_1z + bb4(2)*k_rk_2z+ bb4(3)*k_rk_3z &
                     ,Rst_10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_2_rk4(1,j),xh_2_rk4(2,j),xh_2_rk4(3,j) &
                     ,spin(1,1) + bb4(1)*k_rk_1x + bb4(2)*k_rk_2x+ bb4(3)*k_rk_3x &
                     ,spin(2,1) + bb4(1)*k_rk_1y + bb4(2)*k_rk_2y+ bb4(3)*k_rk_3y &
                     ,spin(3,1) + bb4(1)*k_rk_1z + bb4(2)*k_rk_2z+ bb4(3)*k_rk_3z &
                     ,Rst_5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_4x = tmp*totftides(1)
             k_rk_4y = tmp*totftides(2)
             k_rk_4z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_4x = tmp*totftides(1)
             k_rk_4y = tmp*totftides(2)
             k_rk_4z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_4x = tmp*totftides(1)
             k_rk_4y = tmp*totftides(2)
             k_rk_4z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(1-a(5)) 
          ! Calculation of Runge-kutta factor k5
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_2_rk5(1,j),xh_2_rk5(2,j),xh_2_rk5(3,j) &
                     ,vh_2_rk5(1,j),vh_2_rk5(2,j),vh_2_rk5(3,j) &
                     ,spin(1,1) + bb5(1)*k_rk_1x + bb5(2)*k_rk_2x+ bb5(3)*k_rk_3x+ bb5(4)*k_rk_4x &
                     ,spin(2,1) + bb5(1)*k_rk_1y + bb5(2)*k_rk_2y+ bb5(3)*k_rk_3y+ bb5(4)*k_rk_4y &
                     ,spin(3,1) + bb5(1)*k_rk_1z + bb5(2)*k_rk_2z+ bb5(3)*k_rk_3z+ bb5(4)*k_rk_4z &
                     ,Rst_10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_2_rk5(1,j),xh_2_rk5(2,j),xh_2_rk5(3,j) &
                     ,spin(1,1) + bb5(1)*k_rk_1x + bb5(2)*k_rk_2x+ bb5(3)*k_rk_3x+ bb5(4)*k_rk_4x &
                     ,spin(2,1) + bb5(1)*k_rk_1y + bb5(2)*k_rk_2y+ bb5(3)*k_rk_3y+ bb5(4)*k_rk_4y &
                     ,spin(3,1) + bb5(1)*k_rk_1z + bb5(2)*k_rk_2z+ bb5(3)*k_rk_3z+ bb5(4)*k_rk_4z &
                     ,Rst_5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_5x = tmp*totftides(1)
             k_rk_5y = tmp*totftides(2)
             k_rk_5z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_5x = tmp*totftides(1)
             k_rk_5y = tmp*totftides(2)
             k_rk_5z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_5x = tmp*totftides(1)
             k_rk_5y = tmp*totftides(2)
             k_rk_5z = tmp*totftides(3)
          endif
          
          ! Torque at t=time-hdt*(1-a(6)) 
          ! Calculation of Runge-kutta factor k6
          do j=2,ntid+1 
             if (tides.eq.1) then
                call Torque_tides_s (nbod,m,xh_2_rk6(1,j),xh_2_rk6(2,j),xh_2_rk6(3,j) &
                     ,vh_2_rk6(1,j),vh_2_rk6(2,j),vh_2_rk6(3,j) &
                     ,spin(1,1) + bb6(1)*k_rk_1x + bb6(2)*k_rk_2x+ bb6(3)*k_rk_3x+ bb6(4)*k_rk_4x + bb6(5)*k_rk_5x &
                     ,spin(2,1) + bb6(1)*k_rk_1y + bb6(2)*k_rk_2y+ bb6(3)*k_rk_3y+ bb6(4)*k_rk_4y + bb6(5)*k_rk_5y &
                     ,spin(3,1) + bb6(1)*k_rk_1z + bb6(2)*k_rk_2z+ bb6(3)*k_rk_3z+ bb6(4)*k_rk_4z + bb6(5)*k_rk_5z &
                     ,Rst_10,dissstar,sigmast,j,N_tid_sx,N_tid_sy,N_tid_sz)
             else
                N_tid_sx = 0.0d0
                N_tid_sy = 0.0d0
                N_tid_sz = 0.0d0
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_s (nbod,m,xh_2_rk6(1,j),xh_2_rk6(2,j),xh_2_rk6(3,j) &
                     ,spin(1,1) + bb6(1)*k_rk_1x + bb6(2)*k_rk_2x+ bb6(3)*k_rk_3x+ bb6(4)*k_rk_4x + bb6(5)*k_rk_5x &
                     ,spin(2,1) + bb6(1)*k_rk_1y + bb6(2)*k_rk_2y+ bb6(3)*k_rk_3y+ bb6(4)*k_rk_4y + bb6(5)*k_rk_5y &
                     ,spin(3,1) + bb6(1)*k_rk_1z + bb6(2)*k_rk_2z+ bb6(3)*k_rk_3z+ bb6(4)*k_rk_4z + bb6(5)*k_rk_5z &
                     ,Rst_5,k2s,j,N_rot_sx,N_rot_sy,N_rot_sz)
             else
                N_rot_sx = 0.0d0
                N_rot_sy = 0.0d0
                N_rot_sz = 0.0d0
             endif    
             Ns(1,j) =  tides*N_tid_sx + rot_flat*N_rot_sx
             Ns(2,j) =  tides*N_tid_sy + rot_flat*N_rot_sy  
             Ns(3,j) =  tides*N_tid_sz + rot_flat*N_rot_sz
          enddo
          totftides(1) = 0.d0
          totftides(2) = 0.d0
          totftides(3) = 0.d0
          do j=2,ntid+1 
             tmp = K2/(m(1)+m(j))
             totftides(1) = totftides(1) + tmp*Ns(1,j)
             totftides(2) = totftides(2) + tmp*Ns(2,j)
             totftides(3) = totftides(3) + tmp*Ns(3,j)
          end do   
          if (Rscst.eq.1) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_6x = tmp*totftides(1)
             k_rk_6y = tmp*totftides(2)
             k_rk_6z = tmp*totftides(3)           
          else if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_6x = tmp*totftides(1)
             k_rk_6y = tmp*totftides(2)
             k_rk_6z = tmp*totftides(3)
          else if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
             tmp = - hdt/(rg2s*Rst*Rst)
             k_rk_6x = tmp*totftides(1)
             k_rk_6y = tmp*totftides(2)
             k_rk_6z = tmp*totftides(3)
          endif
          
          ! Integration on first half timestep
          spin(1,1) = spin(1,1) + cc(1)*k_rk_1x + cc(2)*k_rk_2x + cc(3)*k_rk_3x + cc(4)*k_rk_4x + cc(5)*k_rk_5x + cc(6)*k_rk_6x
          spin(2,1) = spin(2,1) + cc(1)*k_rk_1y + cc(2)*k_rk_2y + cc(3)*k_rk_3y + cc(4)*k_rk_4y + cc(5)*k_rk_5y + cc(6)*k_rk_6y
          spin(3,1) = spin(3,1) + cc(1)*k_rk_1z + cc(2)*k_rk_2z + cc(3)*k_rk_3z + cc(4)*k_rk_4z + cc(5)*k_rk_5z + cc(6)*k_rk_6z
          
          
          
          ! **************************************************************
          ! **************************************************************
          ! **************************************************************
          ! PLANETS spin evolution
          ! **************************************************************
          
          do j=2,ntid+1 
             tmp = - dt*K2*m(1)/(m(j)*(m(j)+m(1))*rg2p(j-1)*Rp(j)*Rp(j))
             
             ! Integration on first half timestep
             ! Torque at last time step t=time-dt
             ! Calculation of Runge-kutta factor k1
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_bf(1,j),xh_bf(2,j),xh_bf(3,j) &
                     ,vh_bf(1,j),vh_bf(2,j),vh_bf(3,j) &
                     ,spin_bf(1,j),spin_bf(2,j),spin_bf(3,j),Rp10(j),sigmap(j),j &
                     ,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_bf(1,j),xh_bf(2,j),xh_bf(3,j) &
                     ,spin_bf(1,j),spin_bf(2,j),spin_bf(3,j) &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif
             k_rk_1x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_1y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_1z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Torque at t=time-hdt*(2-a(2)) 
             ! Calculation of Runge-kutta factor k2
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_1_rk2(1,j),xh_1_rk2(2,j),xh_1_rk2(3,j) &
                     ,vh_1_rk2(1,j),vh_1_rk2(2,j),vh_1_rk2(3,j) &
                     ,spin_bf(1,j)+bb2*k_rk_1x &
                     ,spin_bf(2,j)+bb2*k_rk_1y &
                     ,spin_bf(3,j)+bb2*k_rk_1z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_1_rk2(1,j),xh_1_rk2(2,j),xh_1_rk2(3,j) &
                     ,spin_bf(1,j)+bb2*k_rk_1x &
                     ,spin_bf(2,j)+bb2*k_rk_1y &
                     ,spin_bf(3,j)+bb2*k_rk_1z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif       
             k_rk_2x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_2y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_2z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             
             ! Torque at t=time-hdt*(2-a(3)) 
             ! Calculation of Runge-kutta factor k3
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_1_rk3(1,j),xh_1_rk3(2,j),xh_1_rk3(3,j) &
                     ,vh_1_rk3(1,j),vh_1_rk3(2,j),vh_1_rk3(3,j) &
                     ,spin_bf(1,j)+ bb3(1)*k_rk_1x + bb3(2)*k_rk_2x &
                     ,spin_bf(2,j)+ bb3(1)*k_rk_1y + bb3(2)*k_rk_2y &
                     ,spin_bf(3,j)+ bb3(1)*k_rk_1z + bb3(2)*k_rk_2z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_1_rk3(1,j),xh_1_rk3(2,j),xh_1_rk3(3,j) &
                     ,spin_bf(1,j)+ bb3(1)*k_rk_1x + bb3(2)*k_rk_2x &
                     ,spin_bf(2,j)+ bb3(1)*k_rk_1y + bb3(2)*k_rk_2y &
                     ,spin_bf(3,j)+ bb3(1)*k_rk_1z + bb3(2)*k_rk_2z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif           
             k_rk_3x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_3y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_3z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Torque at t=time-hdt*(2-a(4)) 
             ! Calculation of Runge-kutta factor k4
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_1_rk4(1,j),xh_1_rk4(2,j),xh_1_rk4(3,j) &
                     ,vh_1_rk4(1,j),vh_1_rk4(2,j),vh_1_rk4(3,j) &
                     ,spin_bf(1,j)+ bb4(1)*k_rk_1x + bb4(2)*k_rk_2x + bb4(3)*k_rk_3x &
                     ,spin_bf(2,j)+ bb4(1)*k_rk_1y + bb4(2)*k_rk_2y + bb4(3)*k_rk_3y &
                     ,spin_bf(3,j)+ bb4(1)*k_rk_1z + bb4(2)*k_rk_2z + bb4(3)*k_rk_3z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_1_rk4(1,j),xh_1_rk4(2,j),xh_1_rk4(3,j) &
                     ,spin_bf(1,j)+ bb4(1)*k_rk_1x + bb4(2)*k_rk_2x + bb4(3)*k_rk_3x &
                     ,spin_bf(2,j)+ bb4(1)*k_rk_1y + bb4(2)*k_rk_2y + bb4(3)*k_rk_3y &
                     ,spin_bf(3,j)+ bb4(1)*k_rk_1z + bb4(2)*k_rk_2z + bb4(3)*k_rk_3z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif         
             k_rk_4x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_4y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_4z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Torque at t=time-hdt*(2-a(5)) 
             ! Calculation of Runge-kutta factor k5
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_1_rk5(1,j),xh_1_rk5(2,j),xh_1_rk5(3,j) &
                     ,vh_1_rk5(1,j),vh_1_rk5(2,j),vh_1_rk5(3,j) &
                     ,spin_bf(1,j)+ bb5(1)*k_rk_1x + bb5(2)*k_rk_2x + bb5(3)*k_rk_3x + bb5(4)*k_rk_4x &
                     ,spin_bf(2,j)+ bb5(1)*k_rk_1y + bb5(2)*k_rk_2y + bb5(3)*k_rk_3y + bb5(4)*k_rk_4y &
                     ,spin_bf(3,j)+ bb5(1)*k_rk_1z + bb5(2)*k_rk_2z + bb5(3)*k_rk_3z + bb5(4)*k_rk_4z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_1_rk5(1,j),xh_1_rk5(2,j),xh_1_rk5(3,j) &
                     ,spin_bf(1,j)+ bb5(1)*k_rk_1x + bb5(2)*k_rk_2x + bb5(3)*k_rk_3x + bb5(4)*k_rk_4x &
                     ,spin_bf(2,j)+ bb5(1)*k_rk_1y + bb5(2)*k_rk_2y + bb5(3)*k_rk_3y + bb5(4)*k_rk_4y &
                     ,spin_bf(3,j)+ bb5(1)*k_rk_1z + bb5(2)*k_rk_2z + bb5(3)*k_rk_3z + bb5(4)*k_rk_4z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif         
             k_rk_5x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_5y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_5z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Torque at t=time-hdt*(2-a(6)) 
             ! Calculation of Runge-kutta factor k6
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_1_rk6(1,j),xh_1_rk6(2,j),xh_1_rk6(3,j) &
                     ,vh_1_rk6(1,j),vh_1_rk6(2,j),vh_1_rk6(3,j) &
                     ,spin_bf(1,j)+ bb6(1)*k_rk_1x + bb6(2)*k_rk_2x + bb6(3)*k_rk_3x + bb6(4)*k_rk_4x + bb6(5)*k_rk_5x &
                     ,spin_bf(2,j)+ bb6(1)*k_rk_1y + bb6(2)*k_rk_2y + bb6(3)*k_rk_3y + bb6(4)*k_rk_4y + bb6(5)*k_rk_5y &
                     ,spin_bf(3,j)+ bb6(1)*k_rk_1z + bb6(2)*k_rk_2z + bb6(3)*k_rk_3z + bb6(4)*k_rk_4z + bb6(5)*k_rk_5z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_1_rk6(1,j),xh_1_rk6(2,j),xh_1_rk6(3,j) &
                     ,spin_bf(1,j)+ bb6(1)*k_rk_1x + bb6(2)*k_rk_2x + bb6(3)*k_rk_3x + bb6(4)*k_rk_4x + bb6(5)*k_rk_5x &
                     ,spin_bf(2,j)+ bb6(1)*k_rk_1y + bb6(2)*k_rk_2y + bb6(3)*k_rk_3y + bb6(4)*k_rk_4y + bb6(5)*k_rk_5y &
                     ,spin_bf(3,j)+ bb6(1)*k_rk_1z + bb6(2)*k_rk_2z + bb6(3)*k_rk_3z + bb6(4)*k_rk_4z + bb6(5)*k_rk_5z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif         
             k_rk_6x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_6y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_6z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Integration for first half of timestep
         spin(1,j) = spin_bf(1,j) + cc(1)*k_rk_1x + cc(2)*k_rk_2x + cc(3)*k_rk_3x + cc(4)*k_rk_4x + cc(5)*k_rk_5x + cc(6)*k_rk_6x
         spin(2,j) = spin_bf(2,j) + cc(1)*k_rk_1y + cc(2)*k_rk_2y + cc(3)*k_rk_3y + cc(4)*k_rk_4y + cc(5)*k_rk_5y + cc(6)*k_rk_6y
         spin(3,j) = spin_bf(3,j) + cc(1)*k_rk_1z + cc(2)*k_rk_2z + cc(3)*k_rk_3z + cc(4)*k_rk_4z + cc(5)*k_rk_5z + cc(6)*k_rk_6z
             
             ! *********************************************
             ! Integration on second half timestep
             ! Torque at last time step t=time-dt/2
             ! Calculation of Runge-kutta factor k1
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_1_rk5(1,j),xh_1_rk5(2,j),xh_1_rk5(3,j) &
                     ,vh_1_rk5(1,j),vh_1_rk5(2,j),vh_1_rk5(3,j) &
                     ,spin(1,j),spin(2,j),spin(3,j),Rp10(j),sigmap(j),j &
                     ,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_1_rk5(1,j),xh_1_rk5(2,j),xh_1_rk5(3,j) &
                     ,spin(1,j),spin(2,j),spin(3,j) &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif
             k_rk_1x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_1y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_1z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Torque at t=time-hdt*(1-a(2)) 
             ! Calculation of Runge-kutta factor k2
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_2_rk2(1,j),xh_2_rk2(2,j),xh_2_rk2(3,j) &
                     ,vh_2_rk2(1,j),vh_2_rk2(2,j),vh_2_rk2(3,j) &
                     ,spin(1,j)+bb2*k_rk_1x &
                     ,spin(2,j)+bb2*k_rk_1y &
                     ,spin(3,j)+bb2*k_rk_1z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_2_rk2(1,j),xh_2_rk2(2,j),xh_2_rk2(3,j) &
                     ,spin(1,j)+bb2*k_rk_1x &
                     ,spin(2,j)+bb2*k_rk_1y &
                     ,spin(3,j)+bb2*k_rk_1z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif       
             k_rk_2x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_2y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_2z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             
             ! Torque at t=time-hdt*(1-a(3)) 
             ! Calculation of Runge-kutta factor k3
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_2_rk3(1,j),xh_2_rk3(2,j),xh_2_rk3(3,j) &
                     ,vh_2_rk3(1,j),vh_2_rk3(2,j),vh_2_rk3(3,j) &
                     ,spin(1,j)+ bb3(1)*k_rk_1x + bb3(2)*k_rk_2x &
                     ,spin(2,j)+ bb3(1)*k_rk_1y + bb3(2)*k_rk_2y &
                     ,spin(3,j)+ bb3(1)*k_rk_1z + bb3(2)*k_rk_2z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_2_rk3(1,j),xh_2_rk3(2,j),xh_2_rk3(3,j) &
                     ,spin(1,j)+ bb3(1)*k_rk_1x + bb3(2)*k_rk_2x &
                     ,spin(2,j)+ bb3(1)*k_rk_1y + bb3(2)*k_rk_2y &
                     ,spin(3,j)+ bb3(1)*k_rk_1z + bb3(2)*k_rk_2z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif           
             k_rk_3x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_3y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_3z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Torque at t=time-hdt*(1-a(4)) 
             ! Calculation of Runge-kutta factor k4
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_2_rk4(1,j),xh_2_rk4(2,j),xh_2_rk4(3,j) &
                     ,vh_2_rk4(1,j),vh_2_rk4(2,j),vh_2_rk4(3,j) &
                     ,spin(1,j)+ bb4(1)*k_rk_1x + bb4(2)*k_rk_2x + bb4(3)*k_rk_3x &
                     ,spin(2,j)+ bb4(1)*k_rk_1y + bb4(2)*k_rk_2y + bb4(3)*k_rk_3y &
                     ,spin(3,j)+ bb4(1)*k_rk_1z + bb4(2)*k_rk_2z + bb4(3)*k_rk_3z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_2_rk4(1,j),xh_2_rk4(2,j),xh_2_rk4(3,j) &
                     ,spin(1,j)+ bb4(1)*k_rk_1x + bb4(2)*k_rk_2x + bb4(3)*k_rk_3x &
                     ,spin(2,j)+ bb4(1)*k_rk_1y + bb4(2)*k_rk_2y + bb4(3)*k_rk_3y &
                     ,spin(3,j)+ bb4(1)*k_rk_1z + bb4(2)*k_rk_2z + bb4(3)*k_rk_3z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif         
             k_rk_4x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_4y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_4z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Torque at t=time-hdt*(1-a(5)) 
             ! Calculation of Runge-kutta factor k5
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_2_rk5(1,j),xh_2_rk5(2,j),xh_2_rk5(3,j) &
                     ,vh_2_rk5(1,j),vh_2_rk5(2,j),vh_2_rk5(3,j) &
                     ,spin(1,j)+ bb5(1)*k_rk_1x + bb5(2)*k_rk_2x + bb5(3)*k_rk_3x + bb5(4)*k_rk_4x &
                     ,spin(2,j)+ bb5(1)*k_rk_1y + bb5(2)*k_rk_2y + bb5(3)*k_rk_3y + bb5(4)*k_rk_4y &
                     ,spin(3,j)+ bb5(1)*k_rk_1z + bb5(2)*k_rk_2z + bb5(3)*k_rk_3z + bb5(4)*k_rk_4z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_2_rk5(1,j),xh_2_rk5(2,j),xh_2_rk5(3,j) &
                     ,spin(1,j)+ bb5(1)*k_rk_1x + bb5(2)*k_rk_2x + bb5(3)*k_rk_3x + bb5(4)*k_rk_4x &
                     ,spin(2,j)+ bb5(1)*k_rk_1y + bb5(2)*k_rk_2y + bb5(3)*k_rk_3y + bb5(4)*k_rk_4y &
                     ,spin(3,j)+ bb5(1)*k_rk_1z + bb5(2)*k_rk_2z + bb5(3)*k_rk_3z + bb5(4)*k_rk_4z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif         
             k_rk_5x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_5y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_5z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Torque at t=time-hdt*(1-a(6)) 
             ! Calculation of Runge-kutta factor k6
             if (tides.eq.1) then
                call Torque_tides_p (nbod,m,xh_2_rk6(1,j),xh_2_rk6(2,j),xh_2_rk6(3,j) &
                     ,vh_2_rk6(1,j),vh_2_rk6(2,j),vh_2_rk6(3,j) &
                     ,spin(1,j)+ bb6(1)*k_rk_1x + bb6(2)*k_rk_2x + bb6(3)*k_rk_3x + bb6(4)*k_rk_4x + bb6(5)*k_rk_5x &
                     ,spin(2,j)+ bb6(1)*k_rk_1y + bb6(2)*k_rk_2y + bb6(3)*k_rk_3y + bb6(4)*k_rk_4y + bb6(5)*k_rk_5y &
                     ,spin(3,j)+ bb6(1)*k_rk_1z + bb6(2)*k_rk_2z + bb6(3)*k_rk_3z + bb6(4)*k_rk_4z + bb6(5)*k_rk_5z &
                     ,Rp10(j),sigmap(j),j,N_tid_px,N_tid_py,N_tid_pz)
             else
                N_tid_px = 0.0d0
                N_tid_py = 0.0d0
                N_tid_pz = 0.0d0 
             endif
             if (rot_flat.eq.1) then 
                call Torque_rot_p (nbod,m,xh_2_rk6(1,j),xh_2_rk6(2,j),xh_2_rk6(3,j) &
                     ,spin(1,j)+ bb6(1)*k_rk_1x + bb6(2)*k_rk_2x + bb6(3)*k_rk_3x + bb6(4)*k_rk_4x + bb6(5)*k_rk_5x &
                     ,spin(2,j)+ bb6(1)*k_rk_1y + bb6(2)*k_rk_2y + bb6(3)*k_rk_3y + bb6(4)*k_rk_4y + bb6(5)*k_rk_5y &
                     ,spin(3,j)+ bb6(1)*k_rk_1z + bb6(2)*k_rk_2z + bb6(3)*k_rk_3z + bb6(4)*k_rk_4z + bb6(5)*k_rk_5z &
                     ,Rp5(j),k2p(j-1),j,N_rot_px,N_rot_py,N_rot_pz)
             else
                N_rot_px = 0.0d0
                N_rot_py = 0.0d0
                N_rot_pz = 0.0d0 
             endif         
             k_rk_6x = tmp*(tides*N_tid_px + rot_flat*N_rot_px)
             k_rk_6y = tmp*(tides*N_tid_py + rot_flat*N_rot_py)
             k_rk_6z = tmp*(tides*N_tid_pz + rot_flat*N_rot_pz)
             
             ! Integration for second half of timestep
             spin(1,j) = spin(1,j) + cc(1)*k_rk_1x + cc(2)*k_rk_2x + cc(3)*k_rk_3x + cc(4)*k_rk_4x + cc(5)*k_rk_5x + cc(6)*k_rk_6x
             spin(2,j) = spin(2,j) + cc(1)*k_rk_1y + cc(2)*k_rk_2y + cc(3)*k_rk_3y + cc(4)*k_rk_4y + cc(5)*k_rk_5y + cc(6)*k_rk_6y
             spin(3,j) = spin(3,j) + cc(1)*k_rk_1z + cc(2)*k_rk_2z + cc(3)*k_rk_3z + cc(4)*k_rk_4z + cc(5)*k_rk_5z + cc(6)*k_rk_6z

          enddo
       endif
       
       
       !****************************************************************
       !******************** final accelerations ***********************
       !****************************************************************

       do j=2,ntid+1 
          if (tides.eq.1) then 
             call acc_tides (nbod,m,xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),spin &
                  ,Rsth5,Rsth10,k2s,dissstar,sigmast &
                  ,Rp5(j),Rp10(j),k2p(j-1),sigmap(j) &
                  ,j,acc_tid_x,acc_tid_y,acc_tid_z)
             a1(1,j) = acc_tid_x
             a1(2,j) = acc_tid_y
             a1(3,j) = acc_tid_z
             call dEdt_tides (nbod,m,xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),spin &
                  ,Rp10(j),sigmap(j),j,tmp)
             dEdt(j) = tmp
          else
             a1(1,j) = 0.0d0
             a1(2,j) = 0.0d0
             a1(3,j) = 0.0d0
          endif
          if (rot_flat.eq.1) then 
             call acc_rotation (nbod,m,xh(1,j),xh(2,j),xh(3,j),spin &
                  ,Rsth5,k2s,Rp5(j),k2p(j-1) &
                  ,j,acc_rot_x,acc_rot_y,acc_rot_z)
             a3(1,j) = acc_rot_x
             a3(2,j) = acc_rot_y
             a3(3,j) = acc_rot_z
          else
             a3(1,j) = 0.0d0
             a3(2,j) = 0.0d0
             a3(3,j) = 0.0d0
          endif
          if (GenRel.eq.1) then 
             call acc_GenRel (nbod,m,xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j) &
                  ,tintin(j),C2,j,acc_GR_x,acc_GR_y,acc_GR_z)
             a2(1,j) = acc_GR_x
             a2(2,j) = acc_GR_y
             a2(3,j) = acc_GR_z
          else
             a2(1,j) = 0.0d0
             a2(2,j) = 0.0d0
             a2(3,j) = 0.0d0
          endif
          a(1,j) = tides*a1(1,j)+GenRel*a2(1,j)+rot_flat*a3(1,j)
          a(2,j) = tides*a1(2,j)+GenRel*a2(2,j)+rot_flat*a3(2,j)
          a(3,j) = tides*a1(3,j)+GenRel*a2(3,j)+rot_flat*a3(3,j)
       end do
    endif
    
    ! Write stuff
    if ((tides.eq.1).or.(rot_flat.eq.1)) then          
       if (time.ge.timestep) then
          open(13, file="spins.out", access="append")
          write(13,'(8("  ", es19.9e3))') time/365.25d0,spin(1,1),spin(2,1),spin(3,1),Rst/rsun,rg2s,k2s,sigmast
          close(13)
          do j=2,ntid+1
             write(planet_spin_filename,('(a,i1,a)')) 'spinp',j-1,'.out'
             write(planet_orbt_filename,('(a,i1,a)')) 'horb',j-1,'.out'
             write(planet_dEdt_filename,('(a,i1,a)')) 'dEdt',j-1,'.out'
             open(13, file=planet_spin_filename, access='append')
             write(13,'(6("  ", es19.9e3))') time/365.25d0,spin(1,j),spin(2,j),spin(3,j),Rp(j)/rsun,rg2p(j-1)
             close(13)
             open(13, file=planet_orbt_filename, access='append')
             write(13,'(5("  ", es19.9e3))') time/365.25d0,horb(1,j)/horbn(j),horb(2,j)/horbn(j) &
                  ,horb(3,j)/horbn(j),horbn(j)
             close(13)
             open(13, file=planet_dEdt_filename, access='append')
             write(13,'(4("  ", es19.9e3))') time/365.25d0,dEdt(j)
             close(13)
          enddo
          timestep = timestep + output*365.25d0
       endif
    endif    
    
    do j=2,ntid+1 
       xh_bf2(1,j) = xh_bf(1,j)
       xh_bf2(2,j) = xh_bf(2,j)
       xh_bf2(3,j) = xh_bf(3,j)
       vh_bf2(1,j) = vh_bf(1,j)
       vh_bf2(2,j) = vh_bf(2,j)
       vh_bf2(3,j) = vh_bf(3,j)
       xh_bf(1,j) = xh(1,j)
       xh_bf(2,j) = xh(2,j)
       xh_bf(3,j) = xh(3,j)
       vh_bf(1,j) = vh(1,j)
       vh_bf(2,j) = vh(2,j)
       vh_bf(3,j) = vh(3,j)
       spin_bf(1,j) = spin(1,j)
       spin_bf(2,j) = spin(2,j)
       spin_bf(3,j) = spin(3,j)
    enddo
    spin_bf(1,1) = spin(1,1)
    spin_bf(2,1) = spin(2,1)
    spin_bf(3,1) = spin(3,1)
    Rst0    = Rst
    Rst0_5  = Rst_5
    Rst0_10 = Rst_10

    flagbug = flagbug+1
    ispin=1

    return
  end subroutine mfo_user

  subroutine conversion_dh2h (nbod,nbig,m,x,v,xh,vh)

    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,nbig
    real(double_precision),intent(in) :: m(nbod),x(3,nbod),v(3,nbod)

    real(double_precision), intent(out) :: xh(3,nbod),vh(3,nbod)

    ! Local
    integer :: j
    real(double_precision) :: mvsum(3),temp

    !------------------------------------------------------------------------------

    mvsum(1) = 0.d0
    mvsum(2) = 0.d0
    mvsum(3) = 0.d0

    do j = 2, nbod
       xh(1,j) = x(1,j)
       xh(2,j) = x(2,j)
       xh(3,j) = x(3,j)
       mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
       mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
       mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
    end do

    temp = 1.d0 / m(1)
    mvsum(1) = temp * mvsum(1)
    mvsum(2) = temp * mvsum(2)
    mvsum(3) = temp * mvsum(3)
    !~ 
    do j = 2, nbod
       vh(1,j) = v(1,j) + mvsum(1)
       vh(2,j) = v(2,j) + mvsum(2)
       vh(3,j) = v(3,j) + mvsum(3)
    end do

    !------------------------------------------------------------------------------

    return
  end subroutine conversion_dh2h
 
  ! Calculation of r(j), powers of r(j)
  ! Distances in AU
  subroutine rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)

    implicit none

    ! Input/Output
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision), intent(out) :: r_2,rr,r_4,r_5,r_7,r_8
    !------------------------------------------------------------------------------
    r_2 = xhx*xhx+xhy*xhy+xhz*xhz
    rr  = sqrt(r_2)
    r_4 = r_2*r_2
    r_5 = r_4*rr
    r_7 = r_4*r_2*rr
    r_8 = r_4*r_4
    !------------------------------------------------------------------------------
    return
  end subroutine rad_power
  
  ! Calculation of velocity vv(j), radial velocity vrad(j)
  ! velocities in AU/day
  subroutine velocities (xhx,xhy,xhz,vhx,vhy,vhz,v_2,norm_v,v_rad)

    implicit none

    ! Input/Output
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision), intent(out) :: v_2,norm_v,v_rad
    
    !------------------------------------------------------------------------------
    v_2    = vhx*vhx+vhy*vhy+vhz*vhz
    norm_v = sqrt(v_2)         
    ! Radial velocity
    v_rad  = (xhx*vhx+xhy*vhy+xhz*vhz)/sqrt(xhx*xhx+xhy*xhy+xhz*xhz)
    !------------------------------------------------------------------------------
    return
  end subroutine velocities
   
  ! Calculation of r scalar spin
  subroutine r_scal_spin (xhx,xhy,xhz,spinx,spiny,spinz,rscalspin)

    implicit none

    ! Input/Output
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision),intent(in) :: spinx,spiny,spinz

    real(double_precision), intent(out) :: rscalspin

    !------------------------------------------------------------------------------
    rscalspin = xhx*spinx+xhy*spiny+xhz*spinz
    !------------------------------------------------------------------------------

    return
  end subroutine r_scal_spin
  
  ! Calculation of the norm square of the spin
  subroutine norm_spin_2 (spinx,spiny,spinz,normspin_2)

    implicit none

    ! Input/Output
    real(double_precision),intent(in) :: spinx,spiny,spinz

    real(double_precision), intent(out) :: normspin_2

    ! Local
    ! none

    !------------------------------------------------------------------------------
    normspin_2 = spinx*spinx+spiny*spiny+spinz*spinz
    !------------------------------------------------------------------------------

    return
  end subroutine norm_spin_2
  
  !*******************************************
  !**************** TIDES ********************
  !*******************************************
  ! ****************** tidal force *********************
  ! Ftidr in Msun.AU.day-2
  ! Ftidos and Ftidop in Msun.AU.day-1
  ! K2 = G in AU^3.Msun-1.day-2
  
  ! Conservative part of the radial tidal force
  subroutine F_tides_rad_cons (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz &
       ,R_star5,k2_star &
       ,R_plan5,k2_plan,j,Ftidr_cons)

    use physical_constant
    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: R_star5,k2_star
    real(double_precision),intent(in) :: R_plan5,k2_plan
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: Ftidr_cons

    ! Local
!~     integer :: j
    real(double_precision) :: tmp1,tmp2
    real(double_precision) :: r_2,rr,r_4,r_5,r_7,r_8

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    tmp1 = m(1)*m(1)
    tmp2 = m(j)*m(j)
    Ftidr_cons =  -3.0d0/(r_7*K2) &
              *(tmp2*R_star5*k2_star+tmp1*R_plan5*k2_plan)
                          
    !------------------------------------------------------------------------------
    return
  end subroutine F_tides_rad_cons
  
  ! Dissipative part of the radial tidal force
  subroutine F_tides_rad_diss (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz &
       ,R_star10,diss_star,sigma_star &
       ,R_plan10,sigma_plan,j,Ftidr_diss)

    use physical_constant
    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: R_star10,diss_star,sigma_star
    real(double_precision),intent(in) :: R_plan10,sigma_plan
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: Ftidr_diss

    ! Local
!~     integer :: j
    real(double_precision) :: tmp,tmp1,tmp2
    real(double_precision) :: r_2,rr,r_4,r_5,r_7,r_8,v_2,norm_v,v_rad

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    call velocities (xhx,xhy,xhz,vhx,vhy,vhz,v_2,norm_v,v_rad)
    tmp  = K2*K2
    tmp1 = m(1)*m(1)
    tmp2 = m(j)*m(j)
    Ftidr_diss =  - 13.5d0*v_rad/(r_8*tmp) &
              *(tmp2*R_star10*diss_star*sigma_star &
              +tmp1*R_plan10*sigma_plan)              
    !------------------------------------------------------------------------------
    return
  end subroutine F_tides_rad_diss
  
  ! Sum of the dissipative and conservative part of the radial force
  subroutine F_tides_rad (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz &
       ,R_star5,R_star10,k2_star,diss_star,sigma_star &
       ,R_plan5,R_plan10,k2_plan,sigma_plan,j,Ftidr)

    use physical_constant
    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: R_star5,R_star10,k2_star,diss_star,sigma_star
    real(double_precision),intent(in) :: R_plan5,R_plan10,k2_plan,sigma_plan
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: Ftidr

    ! Local
!~     integer :: j
    real(double_precision) :: Ftidr_cons,Ftidr_diss

    !------------------------------------------------------------------------------
    call F_tides_rad_cons (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz &
       ,R_star5,k2_star &
       ,R_plan5,k2_plan,j,Ftidr_cons)
    call F_tides_rad_diss (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz &
       ,R_star10,diss_star,sigma_star &
       ,R_plan10,sigma_plan,j,Ftidr_diss)

    Ftidr = Ftidr_cons + Ftidr_diss            
    !------------------------------------------------------------------------------
    return
  end subroutine F_tides_rad
  
  subroutine F_tides_ortho_star (nbod,m,xhx,xhy,xhz,R_star10,diss_star,sigma_star,j,Ftidos)
  
  use physical_constant
  implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision),intent(in) :: R_star10,diss_star,sigma_star
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: Ftidos

    ! Local
    real(double_precision) :: r_2,rr,r_4,r_5,r_7,r_8

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    Ftidos = 4.5d0*m(j)*m(j)*R_star10*diss_star*sigma_star/(K2*K2*r_7)
    !------------------------------------------------------------------------------
    return
  end subroutine F_tides_ortho_star 
    
  subroutine F_tides_ortho_plan (nbod,m,xhx,xhy,xhz,R_plan10,sigma_plan,j,Ftidop)
  
  use physical_constant
  implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision),intent(in) :: R_plan10,sigma_plan
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: Ftidop

    ! Local
    real(double_precision) :: r_2,rr,r_4,r_5,r_7,r_8

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    Ftidop = 4.5d0*m(1)*m(1)*R_plan10*sigma_plan/(K2*K2*r_7)
    !------------------------------------------------------------------------------
    return
  end subroutine F_tides_ortho_plan   
  
  subroutine Torque_tides_p (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz,spinx,spiny,spinz &
       ,R_plan10,sigma_plan,j,N_tid_px,N_tid_py,N_tid_pz)
  
  implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: spinx,spiny,spinz
    real(double_precision),intent(in) :: R_plan10,sigma_plan
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: N_tid_px,N_tid_py,N_tid_pz

    ! Local
    real(double_precision) :: r_2,rr,r_4,r_5,r_7,r_8,rscalspin,Ftidop

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    call r_scal_spin (xhx,xhy,xhz,spinx,spiny,spinz,rscalspin)
    call F_tides_ortho_plan (nbod,m,xhx,xhy,xhz,R_plan10,sigma_plan,j,Ftidop)        
    N_tid_px = Ftidop*(rr*spinx-rscalspin*xhx/rr-1.0d0/rr*(xhy*vhz-xhz*vhy))
    N_tid_py = Ftidop*(rr*spiny-rscalspin*xhy/rr-1.0d0/rr*(xhz*vhx-xhx*vhz))
    N_tid_pz = Ftidop*(rr*spinz-rscalspin*xhz/rr-1.0d0/rr*(xhx*vhy-xhy*vhx))              
    !------------------------------------------------------------------------------
    return
  end subroutine Torque_tides_p 
  
  subroutine Torque_tides_s (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz,spinx,spiny,spinz &
       ,R_star10,diss_star,sigma_star,j,N_tid_sx,N_tid_sy,N_tid_sz)
  
  implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: spinx,spiny,spinz
    real(double_precision),intent(in) :: R_star10,diss_star,sigma_star
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: N_tid_sx,N_tid_sy,N_tid_sz

    ! Local
    real(double_precision) :: r_2,rr,r_4,r_5,r_7,r_8,rscalspin,Ftidos

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    call r_scal_spin (xhx,xhy,xhz,spinx,spiny,spinz,rscalspin)
    call F_tides_ortho_star (nbod,m,xhx,xhy,xhz,R_star10,diss_star,sigma_star,j,Ftidos)          
    N_tid_sx = Ftidos*(rr*spinx-rscalspin*xhx/rr-1.0d0/rr*(xhy*vhz-xhz*vhy))
    N_tid_sy = Ftidos*(rr*spiny-rscalspin*xhy/rr-1.0d0/rr*(xhz*vhx-xhx*vhz))
    N_tid_sz = Ftidos*(rr*spinz-rscalspin*xhz/rr-1.0d0/rr*(xhx*vhy-xhy*vhx))              
    !------------------------------------------------------------------------------
    return
  end subroutine Torque_tides_s 
  
  subroutine acc_tides (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz,spin,R_star5,R_star10,k2_star,diss_star,sigma_star &
       ,R_plan5,R_plan10,k2_plan,sigma_plan,j,acc_tid_x,acc_tid_y,acc_tid_z)

    use physical_constant
    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: R_star5,R_star10,k2_star,diss_star,sigma_star
    real(double_precision),intent(in) :: R_plan5,R_plan10,k2_plan,sigma_plan
    real(double_precision),intent(in) :: m(nbod),spin(3,10)
    
    real(double_precision), intent(out) :: acc_tid_x,acc_tid_y,acc_tid_z

    ! Local
!~     integer :: j
    real(double_precision) :: tmp,tmp1,r_2,rr,r_4,r_5,r_7,r_8,v_2,norm_v,v_rad
    real(double_precision) :: Ftidr,Ftidos,Ftidop

    !------------------------------------------------------------------------------
    
    call F_tides_rad (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz &
       ,R_star5,R_star10,k2_star,diss_star,sigma_star &
       ,R_plan5,R_plan10,k2_plan,sigma_plan,j,Ftidr)
    call F_tides_ortho_star (nbod,m,xhx,xhy,xhz,R_star10,diss_star,sigma_star,j,Ftidos)
    call F_tides_ortho_plan (nbod,m,xhx,xhy,xhz,R_plan10,sigma_plan,j,Ftidop)
    call velocities (xhx,xhy,xhz,vhx,vhy,vhz,v_2,norm_v,v_rad)
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    
    tmp  = K2/m(j)
    tmp1 = Ftidr+(Ftidos+Ftidop)*v_rad/rr
    acc_tid_x = tmp*(tmp1*xhx/rr &
         + Ftidos/rr*(spin(2,1)*xhz-spin(3,1)*xhy-vhx) &
         + Ftidop/rr*(spin(2,j)*xhz-spin(3,j)*xhy-vhx))
    acc_tid_y = tmp*(tmp1*xhy/rr &
         + Ftidos/rr*(spin(3,1)*xhx-spin(1,1)*xhz-vhy) &
         + Ftidop/rr*(spin(3,j)*xhx-spin(1,j)*xhz-vhy))
    acc_tid_z = tmp*(tmp1*xhz/rr &
         + Ftidos/rr*(spin(1,1)*xhy-spin(2,1)*xhx-vhz) &
         + Ftidop/rr*(spin(1,j)*xhy-spin(2,j)*xhx-vhz)) 
         
    !------------------------------------------------------------------------------
    return
  end subroutine acc_tides
  
  ! Instantaneous energy loss dE/dt due to tides
  ! in Msun.AU^2.day^(-3)
  subroutine dEdt_tides (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz,spin &
       ,R_plan10,sigma_plan,j,dEdt)

    use physical_constant
    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: R_plan10,sigma_plan
    real(double_precision),intent(in) :: m(nbod),spin(3,10)
    
    real(double_precision), intent(out) :: dEdt

    ! Local
!~     integer :: j
    real(double_precision) :: tmp,tmp1
    real(double_precision) :: Ftidr_diss,Ftidop
	real(double_precision) :: r_2,rr,r_4,r_5,r_7,r_8,v_2,norm_v,v_rad

    !------------------------------------------------------------------------------
    call F_tides_rad_diss (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz &
       ,0.0d0,0.0d0,0.0d0 &
       ,R_plan10,sigma_plan,j,Ftidr_diss)
    call F_tides_ortho_plan (nbod,m,xhx,xhy,xhz,R_plan10,sigma_plan,j,Ftidop)
    call velocities (xhx,xhy,xhz,vhx,vhy,vhz,v_2,norm_v,v_rad)
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    
    tmp = Ftidop/rr
    tmp1 = 1.d0/rr*(Ftidr_diss + tmp*v_rad)

    dEdt = tmp1*(xhx*vhx+xhy*vhy+xhz*vhz) &
			+ tmp*((spin(2,j)*xhz-spin(3,j)*xhy-vhx)*vhx &
				  +(spin(3,j)*xhx-spin(1,j)*xhz-vhy)*vhy &
				  +(spin(1,j)*xhy-spin(2,j)*xhx-vhz)*vhz)
         
    !------------------------------------------------------------------------------
    return
  end subroutine dEdt_tides
  
  !********************************************
  !************** ROTATION ********************
  !********************************************
  ! ****************** force due to rotational deformation ********************* 
  ! J2 of planet and star: Jpi and Jsi  (no unit)
  ! Cpi in Msun.AU^5.day-2
  ! Frot_r in Msun.day-2
  ! Frot_os and Frot_op in Msun.AU.day-1
  
  subroutine F_rot_rad (nbod,m,xhx,xhy,xhz,spin,R_star5,k2_star,R_plan5,k2_plan,j,Frot_r)

    use physical_constant
    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: R_star5,k2_star,R_plan5,k2_plan
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision),intent(in) :: spin(3,10)
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: Frot_r

    ! Local
    real(double_precision) :: Cpi,Csi,rscalspinp,rscalspins
    real(double_precision) :: normspin_2p,normspin_2s,r_2,rr,r_4,r_5,r_7,r_8

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    call r_scal_spin (xhx,xhy,xhz,spin(1,j),spin(2,j),spin(3,j),rscalspinp)
    call r_scal_spin (xhx,xhy,xhz,spin(1,1),spin(2,1),spin(3,1),rscalspins)
    call norm_spin_2 (spin(1,j),spin(2,j),spin(3,j),normspin_2p)
    call norm_spin_2 (spin(1,1),spin(2,1),spin(3,1),normspin_2s)
    
    Cpi = m(1)*k2_plan*normspin_2p*R_plan5/(6.d0*K2)
    Csi = m(j)*k2_star*normspin_2s*R_star5/(6.d0*K2)
    
    Frot_r = -3.d0/r_5*(Csi+Cpi) &
         + 15.d0/r_7*(Csi*rscalspins*rscalspins/normspin_2s &
         +Cpi*rscalspinp*rscalspinp/normspin_2p)    
    !------------------------------------------------------------------------------
    return
  end subroutine F_rot_rad
  
  subroutine F_rot_ortho_s (nbod,m,xhx,xhy,xhz,spinx,spiny,spinz,R_star5,k2_star,j,Frot_os)

    use physical_constant
    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: R_star5,k2_star
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision),intent(in) :: spinx,spiny,spinz
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: Frot_os

    ! Local
    real(double_precision) :: Cpi,Csi,rscalspins,normspin_2s,r_2,rr,r_4,r_5,r_7,r_8

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    call r_scal_spin (xhx,xhy,xhz,spinx,spiny,spinz,rscalspins)
    call norm_spin_2 (spinx,spiny,spinz,normspin_2s)
    Csi = m(j)*k2_star*normspin_2s*R_star5/(6.d0*K2)
    Frot_os =  -6.d0*Csi*rscalspins/(normspin_2s*r_5)  
    !------------------------------------------------------------------------------
    return
  end subroutine F_rot_ortho_s
  
  subroutine F_rot_ortho_p (nbod,m,xhx,xhy,xhz,spinx,spiny,spinz,R_plan5,k2_plan,j,Frot_op)

    use physical_constant
    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: R_plan5,k2_plan
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision),intent(in) :: spinx,spiny,spinz
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: Frot_op

    ! Local
    real(double_precision) :: Cpi,Csi,rscalspinp,normspin_2p,r_2,rr,r_4,r_5,r_7,r_8

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    call r_scal_spin (xhx,xhy,xhz,spinx,spiny,spinz,rscalspinp)
    call norm_spin_2 (spinx,spiny,spinz,normspin_2p)
    Cpi = m(1)*k2_plan*normspin_2p*R_plan5/(6.d0*K2)
    Frot_op =  -6.d0*Cpi*rscalspinp/(normspin_2p*r_5)  
    !------------------------------------------------------------------------------
    return
  end subroutine F_rot_ortho_p
  
  subroutine Torque_rot_p (nbod,m,xhx,xhy,xhz,spinx,spiny,spinz &
       ,R_plan5,k2_plan,j,N_rot_px,N_rot_py,N_rot_pz)
  
  implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision),intent(in) :: spinx,spiny,spinz
    real(double_precision),intent(in) :: R_plan5,k2_plan
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: N_rot_px,N_rot_py,N_rot_pz

    ! Local
    real(double_precision) :: Frot_op

    !------------------------------------------------------------------------------
    call F_rot_ortho_p (nbod,m,xhx,xhy,xhz,spinx,spiny,spinz,R_plan5,k2_plan,j,Frot_op)               
    N_rot_px = Frot_op*(xhy*spinz-xhz*spiny)
    N_rot_py = Frot_op*(xhz*spinx-xhx*spinz)
    N_rot_pz = Frot_op*(xhx*spiny-xhy*spinx) 
    !------------------------------------------------------------------------------
    return
  end subroutine Torque_rot_p 
  
  subroutine Torque_rot_s (nbod,m,xhx,xhy,xhz,spinx,spiny,spinz &
       ,R_star5,k2_star,j,N_rot_sx,N_rot_sy,N_rot_sz)
  
  implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision),intent(in) :: spinx,spiny,spinz
    real(double_precision),intent(in) :: R_star5,k2_star
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: N_rot_sx,N_rot_sy,N_rot_sz

    ! Local
    real(double_precision) :: Frot_os

    !------------------------------------------------------------------------------
    call F_rot_ortho_s (nbod,m,xhx,xhy,xhz,spinx,spiny,spinz,R_star5,k2_star,j,Frot_os)               
    N_rot_sx = Frot_os*(xhy*spinz-xhz*spiny)
    N_rot_sy = Frot_os*(xhz*spinx-xhx*spinz)
    N_rot_sz = Frot_os*(xhx*spiny-xhy*spinx) 
    !------------------------------------------------------------------------------
    return
  end subroutine Torque_rot_s 
  
  subroutine acc_rotation (nbod,m,xhx,xhy,xhz,spin,R_star5,k2_star,R_plan5,k2_plan &
       ,j,acc_rot_x,acc_rot_y,acc_rot_z)

    use physical_constant
    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz
    real(double_precision),intent(in) :: R_star5,k2_star
    real(double_precision),intent(in) :: R_plan5,k2_plan
    real(double_precision),intent(in) :: m(nbod),spin(3,10)
    
    real(double_precision), intent(out) :: acc_rot_x,acc_rot_y,acc_rot_z

    ! Local
!~     integer :: j
    real(double_precision) :: tmp
    real(double_precision) :: Frot_r,Frot_os,Frot_op

    !------------------------------------------------------------------------------
    
    call F_rot_rad (nbod,m,xhx,xhy,xhz,spin,R_star5,k2_star,R_plan5,k2_plan,j,Frot_r)
    call F_rot_ortho_s (nbod,m,xhx,xhy,xhz,spin(1,1),spin(2,1),spin(3,1),R_star5,k2_star,j,Frot_os)
    call F_rot_ortho_p (nbod,m,xhx,xhy,xhz,spin(1,j),spin(2,j),spin(3,j),R_plan5,k2_plan,j,Frot_op)
    
    tmp = K2/m(j)
    acc_rot_x = tmp*(Frot_r*xhx + Frot_op*spin(1,j) + Frot_os*spin(1,1))
    acc_rot_y = tmp*(Frot_r*xhy + Frot_op*spin(2,j) + Frot_os*spin(2,1))
    acc_rot_z = tmp*(Frot_r*xhz + Frot_op*spin(3,j) + Frot_os*spin(3,1))
    !------------------------------------------------------------------------------
    return
  end subroutine acc_rotation
  
  !********************************************
  !*********** GENERAL RELATIVITY *************
  !********************************************
  ! FGR_rad in AU.day-2 and FGR_ort in day-1
  
  subroutine F_GR_rad (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz,tintin,C2,j,FGR_rad)

    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: tintin,C2
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: FGR_rad

    ! Local
    real(double_precision) :: r_2,rr,r_4,r_5,r_7,r_8,v_2,norm_v,v_rad

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    call velocities (xhx,xhy,xhz,vhx,vhy,vhz,v_2,norm_v,v_rad)
    
    FGR_rad = -(m(1)+m(j))/(r_2*C2*C2) &
         *((1.0d0+3.0d0*tintin)*v_2 &  
         -2.d0*(2.d0+tintin)*(m(1)+m(j))/rr &
         -1.5d0*tintin*v_rad*v_rad) 
    !------------------------------------------------------------------------------
    return
  end subroutine F_GR_rad
  
  subroutine F_GR_ortho (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz,tintin,C2,j,FGR_ort)

    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: tintin,C2
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: FGR_ort

    ! Local
    real(double_precision) :: r_2,rr,r_4,r_5,r_7,r_8,v_2,norm_v,v_rad

    !------------------------------------------------------------------------------
    call rad_power (xhx,xhy,xhz,r_2,rr,r_4,r_5,r_7,r_8)
    call velocities (xhx,xhy,xhz,vhx,vhy,vhz,v_2,norm_v,v_rad)
    
    FGR_ort = (m(1)+m(j))/(r_2*C2*C2) &
                  *2.0d0*(2.0d0-tintin)*v_rad*norm_v 
    !------------------------------------------------------------------------------
    return
  end subroutine F_GR_ortho
  
  subroutine acc_GenRel (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz,tintin,C2,j,acc_GR_x,acc_GR_y,acc_GR_z)

    implicit none

    ! Input/Output
    integer,intent(in) :: nbod,j
    real(double_precision),intent(in) :: xhx,xhy,xhz,vhx,vhy,vhz
    real(double_precision),intent(in) :: tintin,C2
    real(double_precision),intent(in) :: m(nbod)
    
    real(double_precision), intent(out) :: acc_GR_x,acc_GR_y,acc_GR_z

    ! Local
    real(double_precision) :: tmp,tmp1,FGR_rad,FGR_ort

    !------------------------------------------------------------------------------
    call F_GR_rad (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz,tintin,C2,j,FGR_rad)
    call F_GR_ortho (nbod,m,xhx,xhy,xhz,vhx,vhy,vhz,tintin,C2,j,FGR_ort)
    tmp  = sqrt(xhx*xhx+xhy*xhy+xhz*xhz)
    tmp1 = sqrt(vhx*vhx+vhy*vhy+vhz*vhz)
    acc_GR_x = (FGR_rad*xhx/tmp+FGR_ort*vhx/tmp1)
    acc_GR_y = (FGR_rad*xhy/tmp+FGR_ort*vhy/tmp1)
    acc_GR_z = (FGR_rad*xhz/tmp+FGR_ort*vhz/tmp1)
    !------------------------------------------------------------------------------
    return
  end subroutine acc_GenRel
  
  
  
end module user_module
