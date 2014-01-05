!******************************************************************************
! MODULE: algo_mvs
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that gather various functions about the mvs
!! algorithm.\n\n
!! The standard symplectic (MVS) algorithm is described in J.Widsom and
!! M.Holman (1991) Astronomical Journal, vol 102, pp1528.
!
!******************************************************************************

module algo_mvs

!*************************************************************
!** Modules that gather various functions about the mvs
!** algorithm.
!**
!** Version 1.0 - june 2011
!*************************************************************
  use types_numeriques
  use mercury_globals
  use system_properties
  use user_module
  use forces, only : mfo_ngf
  
  implicit none
  
  private
  
  ! values that need to be saved in mdt_mvs
  real(double_precision), dimension(:,:), allocatable :: a, xj, angf, ausr ! (3,Number of bodies)
  real(double_precision), dimension(:), allocatable :: gm ! (Number of bodies)
  
  public :: mdt_mvs, mco_h2mvs, mco_mvs2h, mco_h2j, allocate_mvs
  
  contains
  
subroutine allocate_mvs(nb_bodies)
! subroutine that allocate various private variable of the module to avoid testing at each timestep
!
! Parameters :
! nb_bodies : number of bodies, that is, the size of the arrays

  implicit none
  
  integer, intent(in) :: nb_bodies
  
  if (.not. allocated(a)) then
    allocate(a(3,nb_bodies))
    a(1:3,1:nb_bodies) = 0.d0
  end if
  
  if (.not. allocated(angf)) then
    allocate(angf(3,nb_bodies))
    angf(1:3,1:nb_bodies) = 0.d0
  end if
  
  if (.not. allocated(ausr)) then
    allocate(ausr(3,nb_bodies))
    ausr(1:3,1:nb_bodies) = 0.d0
  end if
  
  if (.not. allocated(xj)) then
    allocate(xj(3,nb_bodies))
    xj(1:3,1:nb_bodies) = 0.d0
  end if
  
  if (.not. allocated(gm)) then
    allocate(gm(nb_bodies))
    gm(1:nb_bodies) = 0.d0
  end if
  
end subroutine allocate_mvs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCO_H2MVS.FOR    (ErikSoft   28 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Applies an inverse symplectic corrector, which converts coordinates with
! respect to the central body to integrator coordinates for a second-order
! mixed-variable symplectic integrator.

!------------------------------------------------------------------------------

subroutine mco_h2mvs (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag)
  
  use physical_constant
  use mercury_constant
  use drift
  
  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  integer, intent(in) :: ngflag
  
  real(double_precision), intent(in) :: time
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: h
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: xh(3,nbod) !< [in] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision), intent(in) :: vh(3,nbod) !< [in] velocities (vx,vy,vz) with respect to the central body (AU/day)
  real(double_precision), intent(in) :: ngf(4,nbod)
  
  ! Input/Output
  real(double_precision), intent(inout) :: x(3,nbod)
  real(double_precision), intent(inout) :: v(3,nbod)
  
  ! Local
  integer :: j,k,iflag,stat(nb_bodies_initial)
  real(double_precision) :: minside,msofar,gm(nb_bodies_initial)
  real(double_precision) :: a(3,nb_bodies_initial),xj(3,nb_bodies_initial),vj(3,nb_bodies_initial)
  real(double_precision) :: ha(2),hb(2),rt10,angf(3,nb_bodies_initial),ausr(3,nb_bodies_initial)
  
  !------------------------------------------------------------------------------
  
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

  call mco_iden (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag)
  
  ! Calculate effective central masses for Kepler drifts
  minside = m(1)
  do j = 2, nbig
     msofar = minside + m(j)
     gm(j) = m(1) * msofar / minside
     minside = msofar
  end do
  
  do k = 1, 2
     
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j (jcen,nbig,nbig,h,m,x,v,xj,vj)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),ha(k),iflag)
     end do
     
     ! Advance Interaction Hamiltonian
     call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag)
     call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
     
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
     if ((ngflag.eq.1).or.(ngflag.eq.3)) call mfo_ngf (nbod,x,v,angf,ngf)
     
     do j = 2, nbod
        v(1,j) = v(1,j)  +  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  +  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  +  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j (jcen,nbig,nbig,h,m,x,v,xj,vj)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),-2.d0*ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),-2.d0*ha(k),iflag)
     end do
     
     ! Advance Interaction Hamiltonian
     call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag)
     call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
     
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
     if ((ngflag.eq.1).or.(ngflag.eq.3)) call mfo_ngf (nbod,x,v,angf,ngf)
     
     do j = 2, nbod
        v(1,j) = v(1,j)  -  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  -  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  -  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j (jcen,nbig,nbig,h,m,x,v,xj,vj)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),ha(k),iflag)
     end do
     call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag)
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mco_h2mvs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCO_MVS2H.FOR    (ErikSoft   28 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Applies a symplectic corrector, which converts coordinates for a second-
! order mixed-variable symplectic integrator to ones with respect to the
! central body.

!------------------------------------------------------------------------------

subroutine mco_mvs2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag)
  
  use physical_constant
  use mercury_constant
  use drift
  
  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  integer, intent(in) :: ngflag
  real(double_precision), intent(in) :: time
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: h
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  real(double_precision), intent(in) :: ngf(4,nbod)
  
   ! Input/Output
  real(double_precision), intent(inout) :: xh(3,nbod) !< [in,out] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision), intent(inout) :: vh(3,nbod) !< [in,out] velocities (vx,vy,vz) with respect to the central body (AU/day)
  
  ! Local
  integer :: j,k,iflag,stat(nb_bodies_initial)
  real(double_precision) :: minside,msofar,gm(nb_bodies_initial)
  real(double_precision) :: a(3,nb_bodies_initial),xj(3,nb_bodies_initial),vj(3,nb_bodies_initial)
  real(double_precision) :: ha(2),hb(2),rt10,angf(3,nb_bodies_initial),ausr(3,nb_bodies_initial)
  
  !------------------------------------------------------------------------------
  
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

  call mco_iden (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag)
  
  ! Calculate effective central masses for Kepler drifts
  minside = m(1)
  do j = 2, nbig
     msofar = minside + m(j)
     gm(j) = m(1) * msofar / minside
     minside = msofar
  end do
  
  ! Two step corrector
  do k = 1, 2
     
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j(jcen,nbig,nbig,h,m,xh,vh,xj,vj)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),ha(k),iflag)
     end do
     
     ! Advance Interaction Hamiltonian
     call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag)
     call mfo_mvs (jcen,nbod,nbig,m,xh,xj,a,stat)
     
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,xh,vh,ausr)
     if ((ngflag.eq.1).or.(ngflag.eq.3)) call mfo_ngf (nbod,xh,vh,angf,ngf)
     
     do j = 2, nbod
        vh(1,j) = vh(1,j)  -  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        vh(2,j) = vh(2,j)  -  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        vh(3,j) = vh(3,j)  -  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j(jcen,nbig,nbig,h,m,xh,vh,xj,vj)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),-2.d0*ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),-2.d0*ha(k),iflag)
     end do
     
     ! Advance Interaction Hamiltonian
     call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag)
     call mfo_mvs (jcen,nbod,nbig,m,xh,xj,a,stat)
     
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,xh,vh,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,xh,vh,angf,ngf)
     
     do j = 2, nbod
        vh(1,j) = vh(1,j)  +  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        vh(2,j) = vh(2,j)  +  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        vh(3,j) = vh(3,j)  +  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j(jcen,nbig,nbig,h,m,xh,vh,xj,vj)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),ha(k),iflag)
     end do
     call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag)
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mco_mvs2h

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MDT_MVS.FOR    (ErikSoft   28 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Integrates NBOD bodies (of which NBIG are Big) for one timestep H
! using a second-order mixed-variable symplectic integrator.

! DTFLAG = 0 implies first ever call to this subroutine, 
!        = 1 implies first call since number/masses of objects changed.
!        = 2 normal call

! N.B. Input/output must be in coordinates with respect to the central body.
! ===

!------------------------------------------------------------------------------

subroutine mdt_mvs (time,h0,tol,en,am,jcen,rcen,nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,dtflag,ngflag,&
     opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo)
  
  use physical_constant
  use mercury_constant
  use drift
  use dynamic
  
  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  integer, intent(in) :: ngflag
  integer, intent(in) :: opflag
  real(double_precision), intent(in) :: time
  real(double_precision), intent(in) :: h0
  real(double_precision), intent(in) :: tol
  real(double_precision), intent(in) :: am(3)
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: rcen
  real(double_precision), intent(in) :: rphys(nbod)
  real(double_precision), intent(in) :: rce(nbod)
  real(double_precision), intent(in) :: rcrit(nbod)
  real(double_precision), intent(in) :: ngf(4,nbod)
  character(len=8), intent(in) :: id(nbod)
  
  ! Output
  integer, intent(out) :: colflag
  integer, intent(out) :: nclo
  integer, intent(out) :: iclo(CMAX)
  integer, intent(out) :: jclo(CMAX)
  real(double_precision), intent(out) :: tclo(CMAX)
  real(double_precision), intent(out) :: dclo(CMAX)
  real(double_precision), intent(out) :: ixvclo(6,CMAX)
  real(double_precision), intent(out) :: jxvclo(6,CMAX)
  
  ! Input/Output
  integer, intent(inout) :: stat(nbod)
  integer, intent(inout) :: dtflag
  real(double_precision), intent(inout) :: m(nbod) !< [in,out] mass (in solar masses * K2)
  real(double_precision), intent(inout) :: x(3,nbod)
  real(double_precision), intent(inout) :: v(3,nbod)
  real(double_precision), intent(inout) :: s(3,nbod) !< [in,out] spin angular momentum (solar masses AU^2/day)
  real(double_precision), intent(inout) :: en(3)
  
  ! Local
  integer :: j,iflag,nhit,ihit(CMAX),jhit(CMAX),chit(CMAX),nowflag
  real(double_precision) :: vj(3,nb_bodies_initial)
  real(double_precision) :: hby2,thit1,temp
  real(double_precision) :: msofar,minside,x0(3,nb_bodies_initial),v0(3,nb_bodies_initial),dhit(CMAX),thit(CMAX)
  
  !------------------------------------------------------------------------------
  
  hby2 = .5d0 * h0
  nclo = 0
  
  ! If accelerations from previous call are not valid, calculate them now,
  ! and also the Jacobi coordinates XJ, and effective central masses GM.
  if (dtflag.ne.2) then
     dtflag = 2
     call mco_h2j (jcen,nbig,nbig,h0,m,x,v,xj,vj)
     call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
     
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
  
  ! Advance interaction Hamiltonian for H/2
  do j = 2, nbod
     v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
     v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
     v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
  end do
  
  ! Save current coordinates and velocities

  call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,x0,v0,ngf,ngflag)
  
  ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
  call mco_h2j (jcen,nbig,nbig,h0,m,x,v,xj,vj)
  do j = 2, nbig
     iflag = 0
     call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),h0,iflag)
  end do
  do j = nbig + 1, nbod
     iflag = 0
     call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),h0,iflag)
  end do
  call mco_j2h (time,jcen,nbig,nbig,h0,m,xj,vj,x,v,ngf,ngflag)
  
  ! Check for close-encounter minima during drift step
  temp = time + h0
  call mce_stat (temp,h0,rcen,nbod,nbig,m,x0,v0,x,v,rce,rphys,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit,chit,dhit,&
       thit,thit1,nowflag,stat,outfile(3))
  
  ! Advance interaction Hamiltonian for H/2
  call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
  if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
  if ((ngflag.eq.1).or.(ngflag.eq.3)) call mfo_ngf (nbod,x,v,angf,ngf)
  
  do j = 2, nbod
     v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
     v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
     v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mdt_mvs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!     MFO_MVS.FOR    (ErikSoft   2 October 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates accelerations on a set of NBOD bodies (of which NBIG are Big)
! due to gravitational perturbations by all the other bodies.
! This routine is designed for use with a mixed-variable symplectic
! integrator using Jacobi coordinates.

! Based upon routines from Levison and Duncan's SWIFT integrator.

!------------------------------------------------------------------------------

subroutine mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
  
  use physical_constant
  use mercury_constant
  use forces, only : mfo_obl

  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  integer, intent(in) :: stat(nbod)
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: xj(3,nbod)
  
  ! Output
  real(double_precision), intent(out) :: a(3,nbod)
  
  ! Local
  integer :: i,j,k,k1
  real(double_precision) :: fac0,fac1,fac12,fac2,minside,dx,dy,dz,s_1,s2,s_3,faci,facj
  real(double_precision), dimension(3) :: a0,a0tp
  real(double_precision), dimension(3,nb_bodies_initial) :: a1,a2,a3,aobl
  real(double_precision) :: r,r2,r3,rj,rj2,rj3,q,q2,q3,q4,q5,q6,q7,acen(3)
  
  !------------------------------------------------------------------------------
  
  ! Initialize variables
  a0(1) = 0.d0
  a0(2) = 0.d0
  a0(3) = 0.d0
  a1(1,2) = 0.d0
  a1(2,2) = 0.d0
  a1(3,2) = 0.d0
  a2(1,2) = 0.d0
  a2(2,2) = 0.d0
  a2(3,2) = 0.d0
  minside = 0.d0
  
  ! Calculate acceleration terms
  do k = 3, nbig
     k1 = k - 1
     minside = minside + m(k1)
     r2   = x(1,k)  * x(1,k)  +  x(2,k) * x(2,k)  +  x(3,k) * x(3,k)
     rj2  = xj(1,k) * xj(1,k) + xj(2,k) * xj(2,k) + xj(3,k) * xj(3,k)
     r  = 1.d0 / sqrt(r2)
     rj = 1.d0 / sqrt(rj2)
     r3  = r  * r  * r
     rj3 = rj * rj * rj
     
     fac0 = m(k) * r3
     fac12 = m(1) * rj3
     fac2 = m(k) * fac12 / (minside + m(1))
     q = (r2 - rj2) * .5d0 / rj2
     q2 = q  * q
     q3 = q  * q2
     q4 = q2 * q2
     q5 = q2 * q3
     q6 = q3 * q3
     q7 = q3 * q4
     fac1 = 402.1875d0*q7 - 187.6875d0*q6 + 86.625d0*q5   - 39.375d0*q4 + 17.5d0*q3 - 7.5d0*q2 + 3.d0*q - 1.d0
     
     ! Add to A0 term
     a0(1) = a0(1)  -  fac0 * x(1,k)
     a0(2) = a0(2)  -  fac0 * x(2,k)
     a0(3) = a0(3)  -  fac0 * x(3,k)
     
     ! Calculate A1 for this body
     a1(1,k) = fac12 * (xj(1,k) + fac1*x(1,k))
     a1(2,k) = fac12 * (xj(2,k) + fac1*x(2,k))
     a1(3,k) = fac12 * (xj(3,k) + fac1*x(3,k))
     
     ! Calculate A2 for this body
     a2(1,k) = a2(1,k1)  +  fac2 * xj(1,k)
     a2(2,k) = a2(2,k1)  +  fac2 * xj(2,k)
     a2(3,k) = a2(3,k1)  +  fac2 * xj(3,k)
  end do
  
  r2   = x(1,2)  * x(1,2)  +  x(2,2) * x(2,2)  +  x(3,2) * x(3,2)
  r  = 1.d0 / sqrt(r2)
  r3  = r  * r  * r
  fac0 = m(2) * r3
  a0tp(1) = a0(1)  -  fac0 * x(1,2)
  a0tp(2) = a0(2)  -  fac0 * x(2,2)
  a0tp(3) = a0(3)  -  fac0 * x(3,2)
  
  ! Calculate A3 (direct terms)
  do k = 2, nbod
     a3(1,k) = 0.d0
     a3(2,k) = 0.d0
     a3(3,k) = 0.d0
  end do
  do i = 2, nbig
     do j = i + 1, nbig
        dx = x(1,j) - x(1,i)
        dy = x(2,j) - x(2,i)
        dz = x(3,j) - x(3,i)
        s2 = dx*dx + dy*dy + dz*dz
        s_1 = 1.d0 / sqrt(s2)
        s_3 = s_1 * s_1 * s_1
        faci = m(i) * s_3
        facj = m(j) * s_3
        a3(1,j) = a3(1,j)  -  faci * dx
        a3(2,j) = a3(2,j)  -  faci * dy
        a3(3,j) = a3(3,j)  -  faci * dz
        a3(1,i) = a3(1,i)  +  facj * dx
        a3(2,i) = a3(2,i)  +  facj * dy
        a3(3,i) = a3(3,i)  +  facj * dz
     end do
     
     do j = nbig + 1, nbod
        dx = x(1,j) - x(1,i)
        dy = x(2,j) - x(2,i)
        dz = x(3,j) - x(3,i)
        s2 = dx*dx + dy*dy + dz*dz
        s_1 = 1.d0 / sqrt(s2)
        s_3 = s_1 * s_1 * s_1
        faci = m(i) * s_3
        a3(1,j) = a3(1,j)  -  faci * dx
        a3(2,j) = a3(2,j)  -  faci * dy
        a3(3,j) = a3(3,j)  -  faci * dz
     end do
  end do
  
  ! Big-body accelerations
  do k = 2, nbig
     a(1,k) = a0(1) + a1(1,k) + a2(1,k) + a3(1,k)
     a(2,k) = a0(2) + a1(2,k) + a2(2,k) + a3(2,k)
     a(3,k) = a0(3) + a1(3,k) + a2(3,k) + a3(3,k)
  end do
  
  ! Small-body accelerations
  do k = nbig + 1, nbod
     a(1,k) = a0tp(1) + a3(1,k)
     a(2,k) = a0tp(2) + a3(2,k)
     a(3,k) = a0tp(3) + a3(3,k)
  end do
  
  ! Correct for oblateness of the central body
  if ((jcen(1).ne.0).or.(jcen(2).ne.0).or.(jcen(3).ne.0)) then
     call mfo_obl (jcen,nbod,m,x,aobl,acen)
     do k = 2, nbod
        a(1,k) = a(1,k) + (aobl(1,k) - acen(1))
        a(2,k) = a(2,k) + (aobl(2,k) - acen(2))
        a(3,k) = a(3,k) + (aobl(3,k) - acen(3))
     end do
  end if
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_mvs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCO_J2H.FOR    (ErikSoft   2 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Converts Jacobi coordinates to coordinates with respect to the central
! body.

! N.B. The Jacobi coordinates of the small bodies are assumed to be equal
! ===  to their coordinates with respect to the central body.

!------------------------------------------------------------------------------

subroutine mco_j2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag)
  

  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  integer, intent(in) :: ngflag
  real(double_precision), intent(in) :: time
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: h
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  real(double_precision), intent(in) :: ngf(4,nbod)
  
  ! Output
  real(double_precision), intent(out) :: xh(3,nbod) !< [out] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision), intent(out) :: vh(3,nbod) !< [out] velocities (vx,vy,vz) with respect to the central body (AU/day)
  
  ! Local
  integer :: j
  real(double_precision) :: mtot, mx, my, mz, mu, mv, mw, temp
  
  !------------------------------------------------------------------------------
  
  xh(1,2) = x(1,2)
  xh(2,2) = x(2,2)
  xh(3,2) = x(3,2)
  vh(1,2) = v(1,2)
  vh(2,2) = v(2,2)
  vh(3,2) = v(3,2)
  mtot = m(2)
  temp = m(2) / (mtot + m(1))
  mx = temp * x(1,2)
  my = temp * x(2,2)
  mz = temp * x(3,2)
  mu = temp * v(1,2)
  mv = temp * v(2,2)
  mw = temp * v(3,2)
  
  do j = 3, nbig - 1
     xh(1,j) = x(1,j) + mx
     xh(2,j) = x(2,j) + my
     xh(3,j) = x(3,j) + mz
     vh(1,j) = v(1,j) + mu
     vh(2,j) = v(2,j) + mv
     vh(3,j) = v(3,j) + mw
     mtot = mtot + m(j)
     temp = m(j) / (mtot + m(1))
     mx = mx  +  temp * x(1,j)
     my = my  +  temp * x(2,j)
     mz = mz  +  temp * x(3,j)
     mu = mu  +  temp * v(1,j)
     mv = mv  +  temp * v(2,j)
     mw = mw  +  temp * v(3,j)
  enddo
  
  if (nbig.gt.2) then
     xh(1,nbig) = x(1,nbig) + mx
     xh(2,nbig) = x(2,nbig) + my
     xh(3,nbig) = x(3,nbig) + mz
     vh(1,nbig) = v(1,nbig) + mu
     vh(2,nbig) = v(2,nbig) + mv
     vh(3,nbig) = v(3,nbig) + mw
  end if
  
  do j = nbig + 1, nbod
     xh(1,j) = x(1,j)
     xh(2,j) = x(2,j)
     xh(3,j) = x(3,j)
     vh(1,j) = v(1,j)
     vh(2,j) = v(2,j)
     vh(3,j) = v(3,j)
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mco_j2h

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCO_H2J.FOR    (ErikSoft   2 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Converts coordinates with respect to the central body to Jacobi coordinates.

! N.B. The coordinates respect to the central body for the small bodies
! ===  are assumed to be equal to their Jacobi coordinates.

!------------------------------------------------------------------------------

subroutine mco_h2j (jcen,nbod,nbig,h,m,xh,vh,x,v)

  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: h
  real(double_precision), intent(in) :: m(nbig) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: xh(3,nbig) !< [in] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision), intent(in) :: vh(3,nbig) !< [in] velocities (vx,vy,vz) with respect to the central body (AU/day)
  
  ! Output
  real(double_precision), intent(out) :: x(3,nbig)
  real(double_precision), intent(out) :: v(3,nbig)
  
  ! Local
  integer :: j
  real(double_precision) :: mtot, mx, my, mz, mu, mv, mw, temp
  
  !------------------------------------------------------------------------------c
  mtot = m(2)
  x(1,2) = xh(1,2)
  x(2,2) = xh(2,2)
  x(3,2) = xh(3,2)
  v(1,2) = vh(1,2)
  v(2,2) = vh(2,2)
  v(3,2) = vh(3,2)
  mx = m(2) * xh(1,2)
  my = m(2) * xh(2,2)
  mz = m(2) * xh(3,2)
  mu = m(2) * vh(1,2)
  mv = m(2) * vh(2,2)
  mw = m(2) * vh(3,2)
  
  do j = 3, nbig - 1
     temp = 1.d0 / (mtot + m(1))
     mtot = mtot + m(j)
     x(1,j) = xh(1,j)  -  temp * mx
     x(2,j) = xh(2,j)  -  temp * my
     x(3,j) = xh(3,j)  -  temp * mz
     v(1,j) = vh(1,j)  -  temp * mu
     v(2,j) = vh(2,j)  -  temp * mv
     v(3,j) = vh(3,j)  -  temp * mw
     mx = mx  +  m(j) * xh(1,j)
     my = my  +  m(j) * xh(2,j)
     mz = mz  +  m(j) * xh(3,j)
     mu = mu  +  m(j) * vh(1,j)
     mv = mv  +  m(j) * vh(2,j)
     mw = mw  +  m(j) * vh(3,j)
  enddo
  
  if (nbig.gt.2) then
     temp = 1.d0 / (mtot + m(1))
     x(1,nbig) = xh(1,nbig)  -  temp * mx
     x(2,nbig) = xh(2,nbig)  -  temp * my
     x(3,nbig) = xh(3,nbig)  -  temp * mz
     v(1,nbig) = vh(1,nbig)  -  temp * mu
     v(2,nbig) = vh(2,nbig)  -  temp * mv
     v(3,nbig) = vh(3,nbig)  -  temp * mw
  end if
  
  do j = nbig + 1, nbod
     x(1,j) = xh(1,j)
     x(2,j) = xh(2,j)
     x(3,j) = xh(3,j)
     v(1,j) = vh(1,j)
     v(2,j) = vh(2,j)
     v(3,j) = vh(3,j)
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mco_h2j

end module algo_mvs
