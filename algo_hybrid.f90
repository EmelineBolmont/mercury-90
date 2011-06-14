module algo_hybrid

!*************************************************************
!** Modules that gather various functions about the hybrid
!** algorithm.
!**
!** Version 1.0 - june 2011
!*************************************************************
  
  private
  
  public :: mdt_hy, mco_h2dh, mco_dh2h
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_HY.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H
! using a second-order hybrid-symplectic integrator algorithm
!
! DTFLAG = 0 implies first ever call to this subroutine, 
!        = 1 implies first call since number/masses of objects changed.
!        = 2 normal call
!
! N.B. Input/output must be in democratic heliocentric coordinates.
! ===
!
!------------------------------------------------------------------------------
!
subroutine mdt_hy (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag,ngflag,&
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
  integer :: j,nce,ice(NMAX),jce(NMAX),ce(NMAX),iflag
  real(double_precision) :: a(3,NMAX),hby2,hrec,x0(3,NMAX),v0(3,NMAX),mvsum(3),temp
  real(double_precision) :: angf(3,NMAX),ausr(3,NMAX)
  external mfo_hkce
  !
  !------------------------------------------------------------------------------
  !
  save a, hrec, angf, ausr
  hby2 = h0 * .5d0
  nclo = 0
  colflag = 0
  !
  ! If accelerations from previous call are not valid, calculate them now
  if (dtflag.ne.2) then
     if (dtflag.eq.0) hrec = h0
     call mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
     dtflag = 2
     do j = 2, nbod
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
  ! Advance solar Hamiltonian for H/2
  mvsum(1) = 0.d0
  mvsum(2) = 0.d0
  mvsum(3) = 0.d0
  do j = 2, nbod
     mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
     mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
     mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
  end do
  !
  temp = hby2 / m(1)
  mvsum(1) = temp * mvsum(1)
  mvsum(2) = temp * mvsum(2)
  mvsum(3) = temp * mvsum(3)
  do j = 2, nbod
     x(1,j) = x(1,j)  +  mvsum(1)
     x(2,j) = x(2,j)  +  mvsum(2)
     x(3,j) = x(3,j)  +  mvsum(3)
  end do
  !
  ! Save the current coordinates and velocities
  call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,x0,v0,ngf,ngflag,opt)
  !
  ! Advance H_K for H
  do j = 2, nbod
     iflag = 0
     call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),h0,iflag)
  end do
  !
  ! Check whether any object separations were < R_CRIT whilst advancing H_K
  call mce_snif (h0,2,nbod,nbig,x0,v0,x,v,rcrit,ce,nce,ice,jce)
  !
  ! If objects had close encounters, advance H_K using Bulirsch-Stoer instead
  if (nce.gt.0) then
     do j = 2, nbod
        if (ce(j).ne.0) then
           x(1,j) = x0(1,j)
           x(2,j) = x0(2,j)
           x(3,j) = x0(3,j)
           v(1,j) = v0(1,j)
           v(2,j) = v0(2,j)
           v(3,j) = v0(3,j)
        end if
     end do
     call mdt_hkce (time,tstart,h0,hrec,tol,rmax,en(3),jcen,rcen,nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,ngflag,&
          colflag,ce,nce,ice,jce,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile,mem,lmem,mfo_hkce)
  end if
  !
  ! Advance solar Hamiltonian for H/2
  mvsum(1) = 0.d0
  mvsum(2) = 0.d0
  mvsum(3) = 0.d0
  do j = 2, nbod
     mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
     mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
     mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
  end do
  !
  temp = hby2 / m(1)
  mvsum(1) = temp * mvsum(1)
  mvsum(2) = temp * mvsum(2)
  mvsum(3) = temp * mvsum(3)
  do j = 2, nbod
     x(1,j) = x(1,j)  +  mvsum(1)
     x(2,j) = x(2,j)  +  mvsum(2)
     x(3,j) = x(3,j)  +  mvsum(3)
  end do
  !
  ! Advance interaction Hamiltonian for H/2
  call mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
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
end subroutine mdt_hy
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_HKCE.FOR    (ErikSoft   1 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H under
! the Hamiltonian H_K, including close-encounter terms.
!
!------------------------------------------------------------------------------
!
subroutine mdt_hkce (time,tstart,h0,hrec,tol,rmax,elost,jcen,rcen,nbod,nbig,m,x,v,s,rphy,rcrit,rce,stat,id,ngf,algor,opt,ngflag,&
     colflag,ce,nce,ice,jce,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile,mem,lmem,force)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use dynamic

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,nce,ice(nce),jce(nce),stat(nbod),ngflag,ce(nbod)
  integer :: algor,opt(8),colflag,lmem(NMESS),nclo,iclo(CMAX)
  integer :: jclo(CMAX)
  real(double_precision) :: time,tstart,h0,hrec,tol,rmax,elost,jcen(3),rcen
  real(double_precision) :: m(nbod),x(3,nbod),v(3,nbod),s(3,nbod)
  real(double_precision) :: rce(nbod),rphy(nbod),rcrit(nbod),ngf(4,nbod)
  real(double_precision) :: tclo(CMAX),dclo(CMAX),ixvclo(6,CMAX),jxvclo(6,CMAX)
  character*80 outfile(3),mem(NMESS)
  character*8 id(nbod)
  external force
  !
  ! Local
  integer :: iback(NMAX),index(NMAX),ibs(NMAX),jbs(NMAX),nclo_old
  integer :: i,j,k,nbs,nbsbig,statbs(NMAX)
  integer :: nhit,ihit(CMAX),jhit(CMAX),chit(CMAX),nowflag,dtflag
  real(double_precision) :: tlocal,hlocal,hdid,tmp0
  real(double_precision) :: mbs(NMAX),xbs(3,NMAX),vbs(3,NMAX),sbs(3,NMAX)
  real(double_precision) :: rcritbs(NMAX),rcebs(NMAX),rphybs(NMAX)
  real(double_precision) :: ngfbs(4,NMAX),x0(3,NMAX),v0(3,NMAX)
  real(double_precision) :: thit(CMAX),dhit(CMAX),thit1,temp
  character*8 idbs(NMAX)
  !
  !------------------------------------------------------------------------------
  !
  ! N.B. Don't set nclo to zero!!
  nbs = 1
  nbsbig = 0
  mbs(1) = m(1)
  if (algor.eq.11) mbs(1) = m(1) + m(2)
  sbs(1,1) = s(1,1)
  sbs(2,1) = s(2,1)
  sbs(3,1) = s(3,1)
  !
  ! Put data for close-encounter bodies into local arrays for use with BS routine
  do j = 2, nbod
     if (ce(j).ne.0) then
        nbs = nbs + 1
        if (j.le.nbig) nbsbig = nbs
        mbs(nbs)   = m(j)
        xbs(1,nbs) = x(1,j)
        xbs(2,nbs) = x(2,j)
        xbs(3,nbs) = x(3,j)
        vbs(1,nbs) = v(1,j)
        vbs(2,nbs) = v(2,j)
        vbs(3,nbs) = v(3,j)
        sbs(1,nbs) = s(1,j)
        sbs(2,nbs) = s(2,j)
        sbs(3,nbs) = s(3,j)
        rcebs(nbs) = rce(j)
        rphybs(nbs) = rphy(j)
        statbs(nbs) = stat(j)
        rcritbs(nbs) = rcrit(j)
        idbs(nbs) = id(j)
        index(nbs) = j
        iback(j) = nbs
     end if
  end do
  !
  do k = 1, nce
     ibs(k) = iback(ice(k))
     jbs(k) = iback(jce(k))
  end do
  !
  tlocal = 0.d0
  hlocal = sign(hrec,h0)
  !
  ! Begin the Bulirsch-Stoer integration
50 continue
  tmp0 = abs(h0) - abs(tlocal)
  hrec = hlocal
  if (abs(hlocal).gt.tmp0) hlocal = sign (tmp0, h0)
  !
  ! Save old coordinates and integrate
  call mco_iden (time,jcen,nbs,0,h0,mbs,xbs,vbs,x0,v0,ngf,ngflag,opt)
  call mdt_bs2 (time,hlocal,hdid,tol,jcen,nbs,nbsbig,mbs,xbs,vbs,sbs,rphybs,rcritbs,ngfbs,statbs,dtflag,ngflag,opt,nce,ibs,jbs,&
       force)
  tlocal = tlocal + hdid
  !
  ! Check for close-encounter minima
  nclo_old = nclo
  temp = time + tlocal
  call mce_stat (temp,hdid,rcen,nbs,nbsbig,mbs,x0,v0,xbs,vbs,rcebs,rphybs,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit,&
       chit,dhit,thit,thit1,nowflag,statbs,outfile(3),mem,lmem)
  !
  ! If collisions occurred, resolve the collision and return a flag
  if ((nhit.gt.0).and.(opt(2).ne.0)) then
     do k = 1, nhit
        if (chit(k).eq.1) then
           i = ihit(k)
           j = jhit(k)
           call mce_coll (thit(k),tstart,elost,jcen,i,j,nbs,nbsbig,mbs,xbs,vbs,sbs,rphybs,statbs,idbs,opt,mem,lmem,outfile(3))
           colflag = colflag + 1
        end if
     end do
  end if
  !
  ! If necessary, continue integrating objects undergoing close encounters
  if ((tlocal - h0)*h0.lt.0) goto 50
  !
  ! Return data for the close-encounter objects to global arrays
  do k = 2, nbs
     j = index(k)
     m(j)   = mbs(k)
     x(1,j) = xbs(1,k)
     x(2,j) = xbs(2,k)
     x(3,j) = xbs(3,k)
     v(1,j) = vbs(1,k)
     v(2,j) = vbs(2,k)
     v(3,j) = vbs(3,k)
     s(1,j) = sbs(1,k)
     s(2,j) = sbs(2,k)
     s(3,j) = sbs(3,k)
     stat(j) = statbs(k)
  end do
  do k = 1, nclo
     iclo(k) = index(iclo(k))
     jclo(k) = index(jclo(k))
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mdt_hkce

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2DH.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Convert coordinates with respect to the central body to democratic
! heliocentric coordinates.
!
!------------------------------------------------------------------------------
!
subroutine mco_h2dh (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
  real(double_precision) :: v(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: mtot,temp,mvsum(3)
  !
  !------------------------------------------------------------------------------
  !
  mtot = 0.d0
  mvsum(1) = 0.d0
  mvsum(2) = 0.d0
  mvsum(3) = 0.d0
  !
  do j = 2, nbod
     x(1,j) = xh(1,j)
     x(2,j) = xh(2,j)
     x(3,j) = xh(3,j)
     mtot = mtot + m(j)
     mvsum(1) = mvsum(1)  +  m(j) * vh(1,j)
     mvsum(2) = mvsum(2)  +  m(j) * vh(2,j)
     mvsum(3) = mvsum(3)  +  m(j) * vh(3,j)
  end do
  !
  temp = 1.d0 / (m(1) + mtot)
  !
  mvsum(1) = temp * mvsum(1)
  mvsum(2) = temp * mvsum(2)
  mvsum(3) = temp * mvsum(3)
  !
  do j = 2, nbod
     v(1,j) = vh(1,j) - mvsum(1)
     v(2,j) = vh(2,j) - mvsum(2)
     v(3,j) = vh(3,j) - mvsum(3)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_h2dh

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_DH2H.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts democratic heliocentric coordinates to coordinates with respect
! to the central body.
!
!------------------------------------------------------------------------------
!
subroutine mco_dh2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag,opt)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
  real(double_precision) :: vh(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: mvsum(3),temp
  !
  !------------------------------------------------------------------------------
  !
  mvsum(1) = 0.d0
  mvsum(2) = 0.d0
  mvsum(3) = 0.d0
  !
  do j = 2, nbod
     xh(1,j) = x(1,j)
     xh(2,j) = x(2,j)
     xh(3,j) = x(3,j)
     mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
     mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
     mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
  end do
  !
  temp = 1.d0 / m(1)
  mvsum(1) = temp * mvsum(1)
  mvsum(2) = temp * mvsum(2)
  mvsum(3) = temp * mvsum(3)
  !
  do j = 2, nbod
     vh(1,j) = v(1,j) + mvsum(1)
     vh(2,j) = v(2,j) + mvsum(2)
     vh(3,j) = v(3,j) + mvsum(3)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_dh2h

end module algo_hybrid