module user_module

  !*************************************************************
  !** modules that contains user defined modules. 
  !** Only mfo_user will be public.
  !**
  !** Version 1.0 - june 2011
  !*************************************************************
  use types_numeriques

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

  ! N.B. All coordinates and velocities must be with respect to central body
  ! mercury gives x,v in democratic heliocentric !
  ! ===
  !------------------------------------------------------------------------------

  subroutine mfo_user (time,jcen,nbod,nbig,m,x,v,a)

    use physical_constant
    use mercury_constant  
    use tides_constant_GR
    use orbital_elements, only : mco_x2el
    use spline

    implicit none


    ! Input/Output
    integer, intent(in) :: nbod, nbig
    real(double_precision),intent(in) :: time,jcen(3),m(nbod),x(3,nbod),v(3,nbod)
    real(double_precision),intent(out) :: a(3,nbod)
    !real(double_precision),intent(inout) entrée sortie pour spin

    ! Local
    integer :: j, error, nptmss, iPs0
    integer :: flagrg2=0
    integer :: flagtime=0
    integer :: ispin=0
    integer :: iwrite=0
    real(double_precision) :: flagbug=0.d0
    real(double_precision) :: timestep!=3.6525d5!3.65d5 !4.56d6 !
    real(double_precision) :: gm,qq,ee,ii,pp,nn,ll,Pst0,Pst
    real(double_precision) :: dt,tstop,tmp,tmp1,tmp2,sigmast,Jpi,Jsi,Cpi,Csi
    real(double_precision), dimension(2) :: bobo
    real(double_precision), dimension(3) :: totftides
    real(double_precision), dimension(3,nbig+1) :: a1,a2,a3,xh,vh
    real(double_precision), dimension(3,nbig+1) :: horb,trueanom
    real(double_precision), dimension(ntid+1) :: qa,ea,ia,pa,na,la
    real(double_precision), dimension(3,nbig+1) :: Nts,Ntp
    real(double_precision), dimension(3,8) :: spin
    real(double_precision), dimension(8) :: Rp,sigmap,Rp5,Rp10,tintin,k2p,k2pdeltap,rg2p
    real(double_precision), dimension(8) :: rscalws2,rscalwp2,normspin2
    ! don't use after collision
    real(double_precision), dimension(nbig+1) :: r,r2,r4,r5,r7,r8,v2,vv,vrad
    real(double_precision), dimension(nbig+1) :: horbn
    real(double_precision), dimension(4161) :: timeBD,radiusBD,lumiBD,HZinGJ,HZoutGJ,HZinb,HZoutb
    real(double_precision), dimension(1065) :: timestar,radiusstar,d2radiusstar
    real(double_precision), dimension(1065) :: timedM,radiusdM
    real(double_precision), dimension(4755) :: timeJup,radiusJup,k2Jup,rg2Jup,spinJup
    real(double_precision), dimension(37) :: rg2st,trg2,rg1,rg2,rg3,rg4,rg5,rg6,rg7,rg8,rg9,rg10,rg11,rg12
    real(double_precision) :: Rst,Rst5,Rst10,Rstb,Rsth,Rsth5,Rsth10,Rstbh,Rst0,Rstj0,Rstb0,spin0,spinp0,spinb0
    real(double_precision) :: rg2s,rg2sh,rg2s0,k2s
    real(double_precision), dimension(nbig+1) :: Ftr,Ftso,Ftpo,Frotr,Frotos,Frotop
    real(double_precision), dimension(nbig+1) :: FGRo,FGRr
    real(double_precision), parameter, dimension(12) :: Ps0 = (/8.d0,13.d0,19.d0,24.d0,30.d0,36.d0,41.d0, &
         47.d0,53.d0,58.d0,64.d0,70.d0/)
    real(double_precision), parameter, dimension(12) :: k2st = (/0.379d0,0.378d0,0.376d0,0.369d0, &
         0.355d0,0.342d0,0.333d0,0.325d0,0.311d0,0.308d0,0.307d0,0.307d0/)

    ! Data 
    save timestar,radiusstar,d2radiusstar
    save timeBD,radiusBD
    save trg2,rg1,rg2,rg3,rg4,rg5,rg6,rg7,rg8,rg9,rg10,rg11,rg12
    save timedM,radiusdM
    save timeJup,radiusJup,k2Jup,rg2Jup,spinJup
    
    save sigmast,k2s
    save flagrg2,flagtime,ispin,flagbug
    save timestep,nptmss
    save spin,dt
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
       flagtime = flagtime+1
    endif 

    ! Following calculations in heliocentric coordinates   
    call conversion_dh2h(nbod,nbig,m,x,v,xh,vh)    

    if (tides.eq.1) then 
       ! Charge host body data 
       ! + initial condition on host body spin, radius
      
       if (brown_dwarf.eq.1) then 
          if (flagrg2.eq.0) then    
             sigmast = sigma_BD
             ! Radius of BD's gyration data
             open(1,file='rg2BD.dat')
             nptmss = 0
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
                nptmss = 0
                do nptmss = 1,715
                read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0121*K2).and.(m(1).ge.0.0119*K2)) then
                iPs0 = 2
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg2(j)
                enddo
                open(1,file='mass_12.0000.dat')
                nptmss = 0
                do nptmss=1,720
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0151*K2).and.(m(1).ge.0.0149*K2)) then
                iPs0 = 3
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg3(j)
                enddo     
                open(1,file='mass_15.0000.dat')
                nptmss = 0
                do nptmss=1,856
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0201*K2).and.(m(1).ge.0.0199*K2)) then
                iPs0 = 4
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg4(j)
                enddo
                open(1,file='mass_20.0000.dat')
                nptmss = 0
                do nptmss =1,864
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0301*K2).and.(m(1).ge.0.0299*K2)) then
                iPs0 = 5
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg5(j)
                enddo
                open(1,file='mass_30.0000.dat')
                nptmss = 0
                do nptmss=1,878
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0401*K2).and.(m(1).ge.0.0399*K2)) then
                iPs0 = 6
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg6(j)
                enddo
                open(1,file='mass_40.0000.dat')
                nptmss = 0
                do nptmss = 1,886
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)  
                end do
             endif
             if ((m(1).le.0.0501*K2).and.(m(1).ge.0.0499*K2)) then
                iPs0 = 7
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg7(j)
                enddo
                open(1,file='mass_50.0000.dat')
                nptmss = 0
                do nptmss=1,891
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0601*K2).and.(m(1).ge.0.0599*K2)) then
                iPs0 = 8
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg8(j)
                enddo
                open(1,file='mass_60.0000.dat')
                nptmss = 0
                do nptmss=1,1663
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0701*K2).and.(m(1).ge.0.0699*K2)) then
                iPs0 = 9
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg9(j)
                enddo
                open(1,file='mass_70.0000.dat')
                nptmss = 0
                do nptmss =1,3585
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0721*K2).and.(m(1).ge.0.0719*K2)) then
                iPs0 = 10
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg10(j)
                enddo
                open(1,file='mass_72.0000.dat')
                nptmss = 0
                do nptmss =1,3721
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0751*K2).and.(m(1).ge.0.0749*K2)) then
                iPs0 = 11
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg11(j)
                enddo
                open(1,file='mass_75.0000.dat')
                nptmss = 0
                do nptmss =1,3903
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             if ((m(1).le.0.0801*K2).and.(m(1).ge.0.0799*K2)) then
                iPs0 = 12
                k2s  = k2st(iPs0)
                Pst0 = Ps0(iPs0)
                do j = 1,37       
                   rg2st(j) = rg12(j)
                enddo
                open(1,file='mass_80.0000.dat')
                nptmss = 0
                do nptmss =1,4161
                   read(1,*,iostat=error)timeBD(nptmss),radiusBD(nptmss),lumiBD(nptmss), &
                        HZinGJ(nptmss),HZoutGJ(nptmss),HZinb(nptmss),HZoutb(nptmss)
                end do
             endif
             
             ! Initialization of stellar spin (day-1)
             if (crash.eq.0) then
                spin0 = 24.d0*TWOPI/Pst0 
                spin(1,1) = 0.d0
		        spin(2,1) = 0.d0
		        spin(3,1) = spin0
             endif
             if (crash.eq.1) then
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
             nptmss = 0
             do nptmss = 1,1065
                read(1,*,iostat=error)timedM(nptmss),radiusdM(nptmss)
             end do
       
             if (crash.eq.0) then
                Pst = Period_st
                ! Initialization of stellar spin (day-1)
                spin0 = TWOPI/Pst 
                spin(1,1) = 0.d0
	            spin(2,1) = 0.d0
	            spin(3,1) = spin0
             end if
             if (crash.eq.1) then
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
             nptmss = 0
             do nptmss = 1,2003
                read(1,*,iostat=error)timestar(nptmss),radiusstar(nptmss),d2radiusstar(nptmss)
             end do
            ! Initialization Sun-like star spin (day-1)
             if (crash.eq.0) then
                Pst = Period_st
                spin0 = TWOPI/Pst !day-1
                spin(1,1) = 0.d0
	            spin(2,1) = 0.d0
	            spin(3,1) = spin0
             end if
             if (crash.eq.1) then
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
             nptmss = 0
             do nptmss = 1,4755
                read(1,*,iostat=error) timeJup(nptmss),radiusJup(nptmss) &
                     ,k2Jup(nptmss),rg2Jup(nptmss),spinJup(nptmss)
             end do
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
             end if
             if (crash.eq.1) then
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
          Rsth10 = Rsth5*Rsth5
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
		     end if
             if (crash.eq.1) then
                spin(1,1) = rot_crash(1) !day-1
			    spin(2,1) = rot_crash(2) !day-1
                spin(3,1) = rot_crash(3) !day-1
             endif
             flagrg2 = 1
          endif
       endif
       
     endif
       
!~       ! Calculation of r(j), powers of r(j), 
!~       ! velocity vv(j), radial velocity vrad(j) and orbital angular momentum
!~       ! theta_point x er = trueanom (cf textbook)
!~       ! Distances in AU, velocities in AU/day, horb in AU^2/day
          do j=2,ntid+1

             r2(j) = xh(1,j)*xh(1,j)+xh(2,j)*xh(2,j)+xh(3,j)*xh(3,j)
             r(j)  = sqrt(r2(j))
             r4(j) = r2(j)*r2(j)
             r5(j) = r2(j)*r2(j)*r(j)
             r8(j) = r2(j)*r2(j)*r2(j)*r2(j)
             r7(j) = r2(j)*r2(j)*r2(j)*r(j) 
             v2(j) = vh(1,j)*vh(1,j)+vh(2,j)*vh(2,j)+vh(3,j)*vh(3,j)
             vv(j) = sqrt(v2(j))         
             ! Radial velocity
             vrad(j)    = (xh(1,j)*vh(1,j)+xh(2,j)*vh(2,j)+xh(3,j)*vh(3,j))/r(j)
             ! Orbital angular momentum (without mass) in AU^2.day-1
			 horb(1,j)  = (xh(2,j)*vh(3,j)-xh(3,j)*vh(2,j))
		     horb(2,j)  = (xh(3,j)*vh(1,j)-xh(1,j)*vh(3,j))
		     horb(3,j)  = (xh(1,j)*vh(2,j)-xh(2,j)*vh(1,j))    
		     horbn(j) = sqrt(horb(1,j)*horb(1,j)+horb(2,j)*horb(2,j)+horb(3,j)*horb(3,j))
		     
		     ! Term called truenanom in textbook = theta_point x er
		     trueanom(1,j) = 1.d0/r2(j)*(horb(2,j)*xh(3,j)-horb(3,j)*xh(2,j))
		     trueanom(2,j) = 1.d0/r2(j)*(horb(3,j)*xh(1,j)-horb(1,j)*xh(3,j))
		     trueanom(3,j) = 1.d0/r2(j)*(horb(1,j)*xh(2,j)-horb(2,j)*xh(1,j))

          end do


!~       ! Initialization of many things concerning planets
          if (ispin.eq.0) then 
          
             do j=2,ntid+1
                if (tides.eq.1) then 
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
                      k2pdeltap(j-1) = k2pdeltap_terr
                      sigmap(j) = dissplan(j-1)*2.d0*K2*k2pdeltap(j-1)/(3.d0*Rp5(j))
                   endif
                   if (jupiter(j-1).eq.2) then
                      k2p(j-1) = k2p_gg
                      rg2p(j-1) = rg2p_gg
                      k2pdeltap(j-1) = k2pdeltap_gg
                      sigmap(j) = dissplan(j-1)*sigma_gg
                   endif
                   if (jupiter(j-1).eq.3) then
                      k2p(j-1) = k2p_what
                      rg2p(j-1) = rg2p_what
                      k2pdeltap(j-1) = k2pdeltap_what
                      sigmap(j) = dissplan(j-1)*2.d0*K2*k2pdeltap(j-1)/(3.d0*Rp5(j))
                   endif
                endif
                
                ! Factor used for GR force calculation
                if (GenRel.eq.1) tintin(j) = m(1)*m(j)/(m(1)+m(j))**2
             enddo
          
             if (tides.eq.1) then 
!~              ! Initialization of spin planets
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
                end if   
			    if (crash.eq.1) then
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
          
          if (tides.eq.1) then 
             ! Interpolation to have the host body radius (and radius of gyration for BDs)
             ! Here Rst in AU, Rsth5; Rsth10
             if (brown_dwarf.eq.1) then
                if (crash.eq.0) then
                   call spline_b_val(nptmss,timeBD*365.25d0-t_init,radiusBD,time,Rstb0)
                   Rst0= Rsun * Rstb0
                   call spline_b_val(nptmss,timeBD*365.25d0-t_init,radiusBD,time+dt,Rstb)
                   Rst = Rsun * Rstb
                   call spline_b_val(nptmss,timeBD*365.25d0-t_init,radiusBD,time+dt/2.d0,Rstbh)
                   Rsth = Rsun * Rstbh
                   call spline_b_val(37,trg2*365.25d0-t_init,rg2st,time,rg2s0)
                   call spline_b_val(37,trg2*365.25d0-t_init,rg2st,time+dt,rg2s)
                   call spline_b_val(37,trg2*365.25d0-t_init,rg2st,time+dt/2.d0,rg2sh)
                   Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
                   Rsth10 = Rsth5*Rsth5
                endif
                if (crash.eq.1) then
                   call spline_b_val(nptmss,timeBD*365.25d0-t_init-t_crash,radiusBD,time-t_crash,Rstb0)
                   Rst0= Rsun * Rstb0
                   call spline_b_val(nptmss,timeBD*365.25d0-t_init-t_crash,radiusBD,time-t_crash+dt,Rstb)
                   Rst = Rsun * Rstb
                   call spline_b_val(nptmss,timeBD*365.25d0-t_init-t_crash,radiusBD,time-t_crash+dt/2.d0,Rstbh)
                   Rsth = Rsun * Rstbh
                   call spline_b_val(37,trg2*365.25d0-t_init-t_crash,rg2st,time-t_crash,rg2s0)
                   call spline_b_val(37,trg2*365.25d0-t_init-t_crash,rg2st,time-t_crash+dt,rg2s)
                   call spline_b_val(37,trg2*365.25d0-t_init-t_crash,rg2st,time-t_crash+dt/2.d0,rg2sh)
                   Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
                   Rsth10 = Rsth5*Rsth5
                endif
             endif
             if (M_dwarf.eq.1) then
                if (crash.eq.0) then
                   call spline_b_val(nptmss,timedM*365.25-t_init,radiusdM,time,Rstb0)
                   Rst0= Rsun * Rstb0
                   call spline_b_val(nptmss,timedM*365.25-t_init,radiusdM,time+dt,Rstb)
                   Rst = Rsun * Rstb
                   call spline_b_val(nptmss,timedM*365.25-t_init,radiusdM,time+dt/2.d0,Rstbh)
                   Rsth = Rsun * Rstbh
                   Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
                   Rsth10 = Rsth5*Rsth5
                endif
                if (crash.eq.1) then
                   call spline_b_val(nptmss,timedM*365.25-t_init-t_crash,radiusdM,time-t_crash,Rstb0)
                   Rst0= Rsun * Rstb0
                   call spline_b_val(nptmss,timedM*365.25d0-t_init-t_crash,radiusdM,time-t_crash+dt,Rstb)
                   Rst = Rsun * Rstb
                   call spline_b_val(nptmss,timedM*365.25d0-t_init-t_crash,radiusdM,time-t_crash+dt/2.d0,Rstbh)
                   Rsth = Rsun * Rstbh
                   Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
                   Rsth10 = Rsth5*Rsth5
                endif
             endif
             if (Sun_like_star.eq.1) then
                if (crash.eq.0) then
                   call spline_b_val(nptmss,timestar*365.25-t_init,radiusstar,time,Rstb0)
                   Rst0= minau * Rstb0
                   call spline_b_val(nptmss,timestar*365.25-t_init,radiusstar,time+dt,Rstb)
                   Rst = minau * Rstb
                   call spline_b_val(nptmss,timestar*365.25-t_init,radiusstar,time+dt/2.d0,Rstbh)
                   Rsth = minau * Rstbh
                   Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
                   Rsth10 = Rsth5*Rsth5
                endif
                if (crash.eq.1) then
                   call spline_b_val(nptmss,timestar*365.25-t_init-t_crash,radiusstar,time-t_crash,Rstb0)
                   Rst0= minau * Rstb0
                   call spline_b_val(nptmss,timestar*365.25-t_init-t_crash,radiusstar,time-t_crash+dt,Rstb)
                   Rst = minau * Rstb
                   call spline_b_val(nptmss,timestar*365.25-t_init-t_crash,radiusstar,time-t_crash+dt/2.d0,Rstbh)
                   Rsth = minau * Rstbh
                   Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
                   Rsth10 = Rsth5*Rsth5
                endif
             endif
             if (Jupiter_host.eq.1) then
                if (crash.eq.0) then
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init,radiusJup,time,Rstb0)
                   Rst0= minau * Rstb0
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init,radiusJup,time+dt,Rstb)
                   Rst = minau * Rstb
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init,radiusJup,time+dt/2.d0,Rstbh)
                   Rsth = minau * Rstbh
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init,rg2Jup,time,rg2s0)
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init,rg2Jup,time+dt,rg2s)
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init,rg2Jup,time+dt/2.d0,rg2sh)
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init,k2Jup,time+dt/2.d0,k2s)
                   Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
                   Rsth10 = Rsth5*Rsth5
                endif
                if (crash.eq.1) then
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,radiusJup,time-t_crash,Rstb0)
                   Rst0= minau * Rstb0
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,radiusJup,time-t_crash+dt,Rstb)
                   Rst = minau * Rstb
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,radiusJup,time-t_crash+dt/2.d0,Rstbh)
                   Rsth = minau * Rstbh
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,rg2Jup,time-t_crash,rg2s0)
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,rg2Jup,time-t_crash+dt,rg2s)
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,rg2Jup,time-t_crash+dt/2.d0,rg2sh)
                   call spline_b_val(nptmss,timeJup*365.25d0-t_init-t_crash,k2Jup,time-t_crash+dt/2.d0,k2s)
                   Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
                   Rsth10 = Rsth5*Rsth5
                endif
             endif
          endif


          ! **************************************************************
          ! **************************************************************

          ! Radial and orthoradial force of BD on planets
          ! GR forces
          ! See note book for details !!!
          do j=2,ntid+1
             if (tides.eq.1) then 
                ! ****************** tidal force *********************
                ! Ftr in Msun.AU.day-2
                ! Ftso and Ftpo in Msun.AU.day-1
                ! K2 = G in AU^3.Msun-1.day-2
			    tmp  = K2*K2
			    tmp1 = m(1)*m(1)
			    tmp2 = m(j)*m(j)
                Ftr(j) = -3.0d0/(r7(j)*K2) &
                     *(tmp2*Rsth5*k2s+tmp1*Rp5(j)*dissplan(j-1)*k2p(j-1)) & 
                     - 13.5d0*vrad(j)/(r8(j)*tmp) &
                     *(tmp2*Rsth10*dissstar*sigmast &
                     +tmp1*Rp10(j)*sigmap(j))           
                Ftso(j) = 4.5d0*tmp2*Rsth10*dissstar*sigmast/(tmp*r7(j))
                Ftpo(j) = 4.5d0*tmp1*Rp10(j)*sigmap(j)/(tmp*r7(j))
             endif
             if (tides.eq.0) then 
                Ftr(j)  = 0.0d0
                Ftso(j) = 0.0d0
                Ftpo(j) = 0.0d0
             endif
             if (rot_flat.eq.1) then   
                ! ****************** force due to rotational deformation ********************* 
                ! J2 of planet and star: Jpi and Jsi  (no unit)
                ! Cpi in Msun.AU^5.day-2
                ! Frotr in Msun.AU.day-2
                ! Froto in Msun.AU.day-2 
                
                ! Square of the norm of the spin
                normspin2(1) = spin(1,1)*spin(1,1)+spin(2,1)*spin(2,1)+spin(3,1)*spin(3,1)
                normspin2(j) = spin(1,j)*spin(1,j)+spin(2,j)*spin(2,j)+spin(3,j)*spin(3,j)
                ! (r scalar w/w)^2
                rscalws2(j) =(xh(1,j)*spin(1,1)+xh(2,j)*spin(2,1)+xh(3,j)*spin(3,1))**2/normspin2(1)
                rscalwp2(j) =(xh(1,j)*spin(1,j)+xh(2,j)*spin(2,j)+xh(3,j)*spin(3,j))**2/normspin2(j)
                
                Jpi = k2p(j-1)*normspin2(j)*Rp(j)*Rp(j)*Rp(j)/(3.d0*m(j))
                Jsi = k2s*normspin2(1)*Rsth*Rsth*Rsth/(3.d0*m(1))
                Cpi = (m(j)*m(1))/(2.d0*K2)*Jpi*Rp(j)*Rp(j)
                Csi = (m(j)*m(1))/(2.d0*K2)*Jsi*Rsth*Rsth 
                
                Frotr(j) = -3.d0/r4(j)*(Csi+Cpi) &
                     + 15.d0*r(j)/r7(j)*(Csi*rscalws2(j)+Cpi*rscalwp2(j))
                Frotos(j) = -6.d0/r5(j)*Csi*sqrt(rscalws2(j))
                Frotop(j) = -6.d0/r5(j)*Cpi*sqrt(rscalwp2(j))
             endif
             if (rot_flat.eq.0) then 
                Frotr(j)  = 0.0d0
                Frotos(j) = 0.0d0
                Frotop(j) = 0.0d0
             endif  

             !****************** GR forces ****************************
             ! FGRr in AU.day-2 and FGRo in day-1
             if (GenRel.eq.1) then
                FGRr(j) = -(m(1)+m(j))/(r2(j)*C2*C2) &
                     *((1.0d0+3.0d0*tintin(j))*v2(j) &  
                       -2.d0*(2.d0+tintin(j))*(m(1)+m(j))/r(j) &
                       -1.5d0*tintin(j)*vrad(j)*vrad(j))                          
                FGRo(j) = (m(1)+m(j))/(r2(j)*C2*C2) &
                     *2.0d0*(2.0d0-tintin(j))*vrad(j)*vv(j)
             endif
          enddo   
          
          if (ispin.eq.1) then
             if ((tides.eq.1).or.(rot_flat.eq.1)) then 
                ! Calculation of tidal torque !AU,Msun,day
                do j=2,ntid+1  
                   ! star
                   Nts(1,j)  =  Ftso(j)*1.d0/r(j) &
                        *(xh(2,j)*(spin(1,1)*xh(2,j)-spin(2,1)*xh(1,j)-trueanom(3,j)) &
                        - xh(3,j)*(spin(3,1)*xh(1,j)-spin(1,1)*xh(3,j)-trueanom(2,j))) &
                        + Frotos(j)*(spin(3,1)*xh(2,j)-spin(2,1)*xh(3,j))/sqrt(normspin2(1))
                   Nts(2,j)  =  Ftso(j)*1.d0/r(j) &
                        *(xh(3,j)*(spin(2,1)*xh(3,j)-spin(3,1)*xh(2,j)-trueanom(1,j)) &
                        - xh(1,j)*(spin(1,1)*xh(2,j)-spin(2,1)*xh(1,j)-trueanom(3,j))) &  
                        + Frotos(j)*(spin(1,1)*xh(3,j)-spin(3,1)*xh(1,j))/sqrt(normspin2(1))
                   Nts(3,j)  =  Ftso(j)*1.d0/r(j) &
                        *(xh(1,j)*(spin(3,1)*xh(1,j)-spin(1,1)*xh(3,j)-trueanom(2,j)) &
                        - xh(2,j)*(spin(2,1)*xh(3,j)-spin(3,1)*xh(2,j)-trueanom(1,j))) &
                        + Frotos(j)*(spin(2,1)*xh(1,j)-spin(1,1)*xh(2,j))/sqrt(normspin2(1))  
                
                   ! planet
                   Ntp(1,j)  =  Ftpo(j)*1.d0/r(j) &
                        *(xh(2,j)*(spin(1,j)*xh(2,j)-spin(2,j)*xh(1,j)-trueanom(3,j)) &
                        - xh(3,j)*(spin(3,j)*xh(1,j)-spin(1,j)*xh(3,j)-trueanom(2,j))) &
                        + Frotop(j)*(spin(3,j)*xh(2,j)-spin(2,j)*xh(3,j))/sqrt(normspin2(j))
                   Ntp(2,j)  =  Ftpo(j)*1.d0/r(j) &
                        *(xh(3,j)*(spin(2,j)*xh(3,j)-spin(3,j)*xh(2,j)-trueanom(1,j)) &
                        - xh(1,j)*(spin(1,j)*xh(2,j)-spin(2,j)*xh(1,j)-trueanom(3,j))) & 
                        + Frotop(j)*(spin(1,j)*xh(3,j)-spin(3,j)*xh(1,j))/sqrt(normspin2(j))
                   Ntp(3,j)  =  Ftpo(j)*1.d0/r(j) &
                        *(xh(1,j)*(spin(3,j)*xh(1,j)-spin(1,j)*xh(3,j)-trueanom(2,j)) &
                        - xh(2,j)*(spin(2,j)*xh(3,j)-spin(3,j)*xh(2,j)-trueanom(1,j))) & 
                        + Frotop(j)*(spin(2,j)*xh(1,j)-spin(1,j)*xh(2,j))/sqrt(normspin2(j))
                end do
			        
		        ! **************************************************************
                ! **************************************************************
                
                ! Spin evolution
                
                ! STAR
                ! Sum of the different contribution from the different planets        
                totftides(1) = 0.d0
                totftides(2) = 0.d0
                totftides(3) = 0.d0
                do j=2,ntid+1 
                   tmp = K2/(m(1)+m(j))
                   totftides(1) = totftides(1) + tmp*Nts(1,j)
                   totftides(2) = totftides(2) + tmp*Nts(2,j)
                   totftides(3) = totftides(3) + tmp*Nts(3,j)
                end do
                ! d/dt(I.Omega) = - (r x F)
                if (Rscst.eq.1) then 
                   tmp = - dt/(rg2s*Rst*Rst)
                   spin(1,1) = spin(1,1) + tmp*totftides(1)
                   spin(2,1) = spin(2,1) + tmp*totftides(2)
                   spin(3,1) = spin(3,1) + tmp*totftides(3)
                endif
                if ((M_dwarf.eq.1).or.(Sun_like_star.eq.1)) then
                   tmp  = Rst0*Rst0/(Rst*Rst)
                   tmp1 = - dt/(rg2s*Rsth*Rsth)
                   tmp2 = - dt/(rg2s*Rst*Rst)
                   spin(1,1) = tmp*spin(1,1)+tmp2*totftides(1)
                   spin(2,1) = tmp*spin(2,1)+tmp2*totftides(2)
                   spin(3,1) = tmp*spin(3,1)*exp(tmp1*totftides(3)/spin(3,1))
                endif
                if ((brown_dwarf.eq.1).or.(Jupiter_host.eq.1)) then 
                   tmp  = rg2s0/rg2s*Rst0*Rst0/(Rst*Rst)
                   tmp1 = - dt/(rg2sh*Rsth*Rsth)
                   tmp2 = - dt/(rg2s*Rst*Rst)
                   spin(1,1) = tmp*spin(1,1)+tmp2*totftides(1)
                   spin(2,1) = tmp*spin(2,1)+tmp2*totftides(2)
                   spin(3,1) = tmp*spin(3,1)*exp(tmp1*totftides(3)/spin(3,1))
                endif
                   
		        !PLANETS
		        do j=2,ntid+1 
		           tmp = - dt*K2*m(1)/(m(j)*(m(j)+m(1))*rg2p(j-1)*Rp(j)*Rp(j))
		           spin(1,j) = spin(1,j) + tmp*Ntp(1,j)
		           spin(2,j) = spin(2,j) + tmp*Ntp(2,j)
		           spin(3,j) = spin(3,j) + tmp*Ntp(3,j)
		        enddo
             endif
          endif   
             ! Write stuff
          if (tides.eq.1) then   
             if (flagbug.eq.0.0d0) then 
                write(*,*) "time(yr)    spin x,y,z(day-1)     R(Rsun)   rg2  k2s sigmast"
             endif         
             if (time.ge.timestep) then
                write(*,*) "s",time/365.25d0,spin(1,1),spin(2,1),spin(3,1),Rst/rsun,rg2s,k2s,sigmast
                write(*,*) "p1",time/365.25d0,spin(1,2),spin(2,2),spin(3,2),Rp(2)/rsun,rg2p(1)
                write(*,*) "h1",time/365.25d0,horb(1,2)/horbn(2),horb(2,2)/horbn(2) &
                     ,horb(3,2)/horbn(2),horbn(2)
                if (ntid.ge.2) then
                   write(*,*) "p2",time/365.25d0,spin(1,3),spin(2,3),spin(3,3),Rp(3)/rsun,rg2p(2)
                   write(*,*) "h2",time/365.25d0,horb(1,3)/horbn(3),horb(2,3)/horbn(3) &
                     ,horb(3,3)/horbn(3),horbn(3)
                endif
                if (ntid.ge.3) then
                   write(*,*) "p3",time/365.25d0,spin(1,4),spin(2,4),spin(3,4),Rp(4)/rsun,rg2p(3)
                   write(*,*) "h3",time/365.25d0,horb(1,4)/horbn(4),horb(2,4)/horbn(4) &
                     ,horb(3,4)/horbn(4),horbn(4)
                endif
                if (ntid.ge.4) then
                   write(*,*) "p4",time/365.25d0,spin(1,5),spin(2,5),spin(3,5),Rp(5)/rsun,rg2p(4)
                   write(*,*) "h4",time/365.25d0,horb(1,5)/horbn(5),horb(2,5)/horbn(5) &
                     ,horb(3,5)/horbn(5),horbn(5)
                endif
                if (ntid.ge.5) then
                   write(*,*) "p5",time/365.25d0,spin(1,6),spin(2,6),spin(3,6),Rp(6)/rsun,rg2p(5)
                   write(*,*) "h5",time/365.25d0,horb(1,6)/horbn(6),horb(2,6)/horbn(6) &
                     ,horb(3,6)/horbn(6),horbn(6)
                endif
                if (ntid.ge.6) then
                   write(*,*) "p6",time/365.25d0,spin(1,7),spin(2,7),spin(3,7),Rp(7)/rsun,rg2p(6)
                   write(*,*) "h6",time/365.25d0,horb(1,7)/horbn(7),horb(2,7)/horbn(7) &
                     ,horb(3,7)/horbn(7),horbn(7)
                endif
                timestep = timestep + output*365.25d0
                
             endif
          endif

          !****************************************************************
          !******************** final accelerations ***********************

          ! We add acceleration due to tides and GR
          do j=2,ntid+1 
             if ((tides.eq.1).or.(rot_flat.eq.1)) then 
                a1(1,j) = K2/m(j)*((Ftr(j)+Frotr(j))*xh(1,j)/r(j) &
                     +Ftso(j)/r(j)*(spin(2,1)*xh(3,j)-spin(3,1)*xh(2,j)-trueanom(1,j)) &
                     +Ftpo(j)/r(j)*(spin(2,j)*xh(3,j)-spin(3,j)*xh(2,j)-trueanom(1,j)) &
                     +Frotop(j)*spin(1,j)/sqrt(normspin2(j)) &
                     +Frotos(j)*spin(1,1)/sqrt(normspin2(1)))
                a1(2,j) = K2/m(j)*((Ftr(j)+Frotr(j))*xh(2,j)/r(j) &
                     +Ftso(j)/r(j)*(spin(3,1)*xh(1,j)-spin(1,1)*xh(3,j)-trueanom(2,j)) &
                     +Ftpo(j)/r(j)*(spin(3,j)*xh(1,j)-spin(1,j)*xh(3,j)-trueanom(2,j)) &
                     +Frotop(j)*spin(2,j)/sqrt(normspin2(j)) &
                     +Frotos(j)*spin(2,1)/sqrt(normspin2(1)))
                a1(3,j) = K2/m(j)*((Ftr(j)+Frotr(j))*xh(3,j)/r(j) &
                     +Ftso(j)/r(j)*(spin(1,1)*xh(2,j)-spin(2,1)*xh(1,j)-trueanom(3,j)) &
                     +Ftpo(j)/r(j)*(spin(1,j)*xh(2,j)-spin(2,j)*xh(1,j)-trueanom(3,j)) &
                     +Frotop(j)*spin(3,j)/sqrt(normspin2(j)) &
                     +Frotos(j)*spin(3,1)/sqrt(normspin2(1)))
             endif
		     if (GenRel.eq.1) then 
                a2(1,j) = (FGRr(j)*xh(1,j)/r(j)+FGRo(j)*vh(1,j)/vv(j))
                a2(2,j) = (FGRr(j)*xh(2,j)/r(j)+FGRo(j)*vh(2,j)/vv(j))
                a2(3,j) = (FGRr(j)*xh(3,j)/r(j)+FGRo(j)*vh(3,j)/vv(j))
             endif

             a(1,j) = a1(1,j)+a2(1,j)
             a(2,j) = a1(2,j)+a2(2,j)
             a(3,j) = a1(3,j)+a2(3,j)
          end do

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

end module user_module
