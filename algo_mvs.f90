module algo_mvs

!*************************************************************
!** Modules that gather various functions about the mvs
!** algorithm.
!**
!** Version 1.0 - june 2011
!*************************************************************
  use user_module
  
  private
  
  public :: mdt_mvs, mco_h2mvs, mco_mvs2h
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2MVS.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Applies an inverse symplectic corrector, which converts coordinates with
! respect to the central body to integrator coordinates for a second-order
! mixed-variable symplectic integrator.
!
!------------------------------------------------------------------------------
!
subroutine mco_h2mvs (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use drift
  
  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
  real(double_precision) :: v(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j,k,iflag,stat(NMAX)
  real(double_precision) :: minside,msofar,gm(NMAX),a(3,NMAX),xj(3,NMAX),vj(3,NMAX)
  real(double_precision) :: ha(2),hb(2),rt10,angf(3,NMAX),ausr(3,NMAX)
  !
  !------------------------------------------------------------------------------
  !
  rt10 = sqrt(10.d0)
  ha(1) = -h * rt10 / 5.d0
  hb(1) = -h * rt10 / 24.d0
  ha(2) = -h * rt10 * 3.d0 / 10.d0
  hb(2) =  h * rt10 / 72.d0
  do j = 2, nbod
     angf(1,j) = 0.d0
     angf(2,j) = 0.d0
     angf(3,j) = 0.d0
     ausr(1,j) = 0.d0
     ausr(2,j) = 0.d0
     ausr(3,j) = 0.d0
  end do
  call mco_iden (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  ! Calculate effective central masses for Kepler drifts
  minside = m(1)
  do j = 2, nbig
     msofar = minside + m(j)
     gm(j) = m(1) * msofar / minside
     minside = msofar
  end do
  !
  do k = 1, 2
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j (time,jcen,nbig,nbig,h,m,x,v,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),ha(k),iflag)
     end do
     !
     ! Advance Interaction Hamiltonian
     call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
     !
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
     !
     do j = 2, nbod
        v(1,j) = v(1,j)  +  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  +  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  +  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j (time,jcen,nbig,nbig,h,m,x,v,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),-2.d0*ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),-2.d0*ha(k),iflag)
     end do
     !
     ! Advance Interaction Hamiltonian
     call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
     !
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
     !
     do j = 2, nbod
        v(1,j) = v(1,j)  -  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  -  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  -  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j (time,jcen,nbig,nbig,h,m,x,v,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),ha(k),iflag)
     end do
     call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag,opt)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_h2mvs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_MVS2H.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Applies a symplectic corrector, which converts coordinates for a second-
! order mixed-variable symplectic integrator to ones with respect to the
! central body.
!
!------------------------------------------------------------------------------
!
subroutine mco_mvs2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag,opt)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use drift
  
  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
  real(double_precision) :: vh(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j,k,iflag,stat(NMAX)
  real(double_precision) :: minside,msofar,gm(NMAX),a(3,NMAX),xj(3,NMAX),vj(3,NMAX)
  real(double_precision) :: ha(2),hb(2),rt10,angf(3,NMAX),ausr(3,NMAX)
  !
  !------------------------------------------------------------------------------
  !
  rt10 = sqrt(10.d0)
  ha(1) =  h * rt10 * 3.d0 / 10.d0
  hb(1) = -h * rt10 / 72.d0
  ha(2) =  h * rt10 / 5.d0
  hb(2) =  h * rt10 / 24.d0
  do j = 2, nbod
     angf(1,j) = 0.d0
     angf(2,j) = 0.d0
     angf(3,j) = 0.d0
     ausr(1,j) = 0.d0
     ausr(2,j) = 0.d0
     ausr(3,j) = 0.d0
  end do
  call mco_iden (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag,opt)
  !
  ! Calculate effective central masses for Kepler drifts
  minside = m(1)
  do j = 2, nbig
     msofar = minside + m(j)
     gm(j) = m(1) * msofar / minside
     minside = msofar
  end do
  !
  ! Two step corrector
  do k = 1, 2
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j(time,jcen,nbig,nbig,h,m,xh,vh,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),ha(k),iflag)
     end do
     !
     ! Advance Interaction Hamiltonian
     call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,xh,xj,a,stat)
     !
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,xh,vh,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,xh,vh,angf,ngf)
     !
     do j = 2, nbod
        vh(1,j) = vh(1,j)  -  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        vh(2,j) = vh(2,j)  -  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        vh(3,j) = vh(3,j)  -  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j(time,jcen,nbig,nbig,h,m,xh,vh,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),-2.d0*ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),-2.d0*ha(k),iflag)
     end do
     !
     ! Advance Interaction Hamiltonian
     call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,xh,xj,a,stat)
     !
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,xh,vh,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,xh,vh,angf,ngf)
     !
     do j = 2, nbod
        vh(1,j) = vh(1,j)  +  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        vh(2,j) = vh(2,j)  +  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        vh(3,j) = vh(3,j)  +  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j(time,jcen,nbig,nbig,h,m,xh,vh,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),ha(k),iflag)
     end do
     call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag,opt)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_mvs2h

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_MVS.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H
! using a second-order mixed-variable symplectic integrator.
!
! DTFLAG = 0 implies first ever call to this subroutine, 
!        = 1 implies first call since number/masses of objects changed.
!        = 2 normal call
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
subroutine mdt_mvs (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag,ngflag,&
     opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile,mem,lmem)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use drift
  use dynamic
  
  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,stat(nbod),algor,opt(8),dtflag,ngflag,opflag
  integer :: colflag,lmem(NMESS),nclo,iclo(CMAX),jclo(CMAX)
  real(double_precision) :: time,tstart,h0,tol,rmax,en(3),am(3),jcen(3),rcen
  real(double_precision) :: m(nbod),x(3,nbod),v(3,nbod),s(3,nbod),rphys(nbod)
  real(double_precision) :: rce(nbod),rcrit(nbod),ngf(4,nbod),tclo(CMAX),dclo(CMAX)
  real(double_precision) :: ixvclo(6,CMAX),jxvclo(6,CMAX)
  character*80 outfile(3),mem(NMESS)
  character*8 id(nbod)
  !
  ! Local
  integer :: j,iflag,nhit,ihit(CMAX),jhit(CMAX),chit(CMAX),nowflag
  real(double_precision) :: xj(3,NMAX),vj(3,NMAX),a(3,NMAX),gm(NMAX),hby2,thit1,temp
  real(double_precision) :: msofar,minside,x0(3,NMAX),v0(3,NMAX),dhit(CMAX),thit(CMAX)
  real(double_precision) :: angf(3,NMAX),ausr(3,NMAX)
  !
  !------------------------------------------------------------------------------
  !
  save a, xj, gm, angf, ausr
  hby2 = .5d0 * h0
  nclo = 0
  !
  ! If accelerations from previous call are not valid, calculate them now,
  ! and also the Jacobi coordinates XJ, and effective central masses GM.
  if (dtflag.ne.2) then
     dtflag = 2
     call mco_h2j (time,jcen,nbig,nbig,h0,m,x,v,xj,vj,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
     !
     minside = m(1)
     do j = 2, nbig
        msofar = minside + m(j)
        gm(j) = m(1) * msofar / minside
        minside = msofar
        angf(1,j) = 0.d0
        angf(2,j) = 0.d0
        angf(3,j) = 0.d0
        ausr(1,j) = 0.d0
        ausr(2,j) = 0.d0
        ausr(3,j) = 0.d0
     end do
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
  end if
  !
  ! Advance interaction Hamiltonian for H/2
  do j = 2, nbod
     v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
     v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
     v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
  end do
  !
  ! Save current coordinates and velocities
  call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,x0,v0,ngf,ngflag,opt)
  !
  ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
  call mco_h2j (time,jcen,nbig,nbig,h0,m,x,v,xj,vj,ngf,ngflag,opt)
  do j = 2, nbig
     iflag = 0
     call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),h0,iflag)
  end do
  do j = nbig + 1, nbod
     iflag = 0
     call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),h0,iflag)
  end do
  call mco_j2h (time,jcen,nbig,nbig,h0,m,xj,vj,x,v,ngf,ngflag,opt)
  !
  ! Check for close-encounter minima during drift step
  temp = time + h0
  call mce_stat (temp,h0,rcen,nbod,nbig,m,x0,v0,x,v,rce,rphys,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit,chit,dhit,thit,&
       thit1,nowflag,stat,outfile(3),mem,lmem)
  !
  ! Advance interaction Hamiltonian for H/2
  call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
  if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
  if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
  !
  do j = 2, nbod
     v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
     v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
     v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mdt_mvs

end module algo_mvs