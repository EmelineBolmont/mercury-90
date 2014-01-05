!******************************************************************************
! MODULE: forces
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that contains all common forces
!
!******************************************************************************

module forces

  use types_numeriques
  use mercury_globals

  implicit none
  
  private
  
  public :: mfo_all
  public :: mfo_ngf ! Needed by HYBRID and MVS
  public :: mfo_obl ! Needed by HYBRID and MVS on respectively mfo_hy and mfo_mvs
  
  contains
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MFO_ALL.FOR    (ErikSoft   2 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates accelerations on a set of NBOD bodies (of which NBIG are Big)
! due to Newtonian gravitational perturbations, post-Newtonian
! corrections (if required), cometary non-gravitational forces (if required)
! and user-defined forces (if required).

! N.B. Input/output must be in coordinates with respect to the central body.
! ===

!------------------------------------------------------------------------------

subroutine mfo_all (time,jcen,nbod,nbig,m,x,v,s,rcrit,a,stat,ngf,ngflag,nce,ice,jce)
  
  use physical_constant
  use mercury_constant
  use user_module

  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  integer, intent(in) :: ngflag
  integer, intent(in) :: stat(nbod) !< [in] status (0 => alive, <>0 => to be removed)
  integer, intent(in) :: nce
  integer, intent(in) :: ice(nce)
  integer, intent(in) :: jce(nce)
  real(double_precision), intent(in) :: time,jcen(3)
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  real(double_precision), intent(in) :: s(3,nbod) !< [in] spin angular momentum (solar masses AU^2/day)
  real(double_precision), intent(in) :: ngf(4,nbod) !< [in] non gravitational forces parameters
  !! \n(1-3) cometary non-gravitational (jet) force parameters
  !! \n(4)  beta parameter for radiation pressure and P-R drag
  real(double_precision), intent(in) :: rcrit(nbod)
  
  ! Output
  real(double_precision), intent(out) :: a(3,nbod)
  
  ! Local
  integer :: j
  real(double_precision) :: acor(3,nb_bodies_initial),acen(3)
  
  !------------------------------------------------------------------------------
  
  ! Newtonian gravitational forces
  call mfo_grav (nbod,nbig,m,x,v,a,stat)
  
  ! Correct for oblateness of the central body
  if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
     call mfo_obl (jcen,nbod,m,x,acor,acen)
     do j = 2, nbod
        a(1,j) = a(1,j) + (acor(1,j) - acen(1))
        a(2,j) = a(2,j) + (acor(2,j) - acen(2))
        a(3,j) = a(3,j) + (acor(3,j) - acen(3))
     end do
  end if
  
  ! Include non-gravitational (cometary jet) accelerations if necessary
  if ((ngflag.eq.1).or.(ngflag.eq.3)) then
     call mfo_ngf (nbod,x,v,acor,ngf)
     do j = 2, nbod
        a(1,j) = a(1,j) + acor(1,j)
        a(2,j) = a(2,j) + acor(2,j)
        a(3,j) = a(3,j) + acor(3,j)
     end do
  end if
  
  ! Include radiation pressure/Poynting-Robertson drag if necessary
  if (ngflag.eq.2.or.ngflag.eq.3) then
     call mfo_pr (nbod,nbig,m,x,v,acor,ngf)
     do j = 2, nbod
        a(1,j) = a(1,j) + acor(1,j)
        a(2,j) = a(2,j) + acor(2,j)
        a(3,j) = a(3,j) + acor(3,j)
     end do
  end if
  
  ! Include post-Newtonian corrections if required
  if (opt(7).eq.1) then
     call mfo_pn (nbod,nbig,m,x,v,acor)
     do j = 2, nbod
        a(1,j) = a(1,j) + acor(1,j)
        a(2,j) = a(2,j) + acor(2,j)
        a(3,j) = a(3,j) + acor(3,j)
     end do
  end if
  
  ! Include user-defined accelerations if required
  if (opt(8).eq.1) then
     call mfo_user (time,jcen,nbod,nbig,m,x,v,acor)
     do j = 2, nbod
        a(1,j) = a(1,j) + acor(1,j)
        a(2,j) = a(2,j) + acor(2,j)
        a(3,j) = a(3,j) + acor(3,j)
     end do
  end if
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_all

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MFO_GRAV.FOR    (ErikSoft   3 October 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates accelerations on a set of NBOD bodies (NBIG of which are Big)
! due to gravitational perturbations by all the other bodies, except that
! Small bodies do not interact with one another.

! The positions and velocities are stored in arrays X, V with the format
! (x,y,z) and (vx,vy,vz) for each object in succession. The accelerations 
! are stored in the array A (ax,ay,az).

! N.B. All coordinates and velocities must be with respect to central body!!!
! ===
!------------------------------------------------------------------------------

subroutine mfo_grav (nbod,nbig,m,x,v,a,stat)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  integer, intent(in) :: stat(nbod) !< [in] status (0 => alive, <>0 => to be removed)
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  
  ! Output
  real(double_precision), intent(out) :: a(3,nbod)
  
  ! Local
  integer :: i, j
  real(double_precision) :: sx, sy, sz, dx, dy, dz, tmp1, tmp2, s_1, s2, s_3, r3(nb_bodies_initial)
  
  !------------------------------------------------------------------------------
  
  sx = 0.d0
  sy = 0.d0
  sz = 0.d0
  do i = 2, nbod
     a(1,i) = 0.d0
     a(2,i) = 0.d0
     a(3,i) = 0.d0
     s2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
     s_1  = 1.d0 / sqrt(s2)
     r3(i) = s_1 * s_1 * s_1
  end do
  
  do i = 2, nbod
     tmp1 = m(i) * r3(i)
     sx = sx  -  tmp1 * x(1,i)
     sy = sy  -  tmp1 * x(2,i)
     sz = sz  -  tmp1 * x(3,i)
  end do
  
  ! Direct terms
  do i = 2, nbig
     do j = i + 1, nbod
        dx = x(1,j) - x(1,i)
        dy = x(2,j) - x(2,i)
        dz = x(3,j) - x(3,i)
        s2 = dx*dx + dy*dy + dz*dz
        s_1 = 1.d0 / sqrt(s2)
        s_3 = s_1 * s_1 * s_1
        tmp1 = s_3 * m(i)
        tmp2 = s_3 * m(j)
        a(1,j) = a(1,j)  -  tmp1 * dx
        a(2,j) = a(2,j)  -  tmp1 * dy
        a(3,j) = a(3,j)  -  tmp1 * dz
        a(1,i) = a(1,i)  +  tmp2 * dx
        a(2,i) = a(2,i)  +  tmp2 * dy
        a(3,i) = a(3,i)  +  tmp2 * dz
     end do
  end do
  
  ! Indirect terms (add these on last to reduce roundoff error)
  do i = 2, nbod
     tmp1 = m(1) * r3(i)
     a(1,i) = a(1,i)  +  sx  -  tmp1 * x(1,i)
     a(2,i) = a(2,i)  +  sy  -  tmp1 * x(2,i)
     a(3,i) = a(3,i)  +  sz  -  tmp1 * x(3,i)
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_grav

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MFO_NGF.FOR    (ErikSoft  29 November 1999)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates accelerations on a set of NBOD bodies due to cometary
! non-gravitational jet forces. The positions and velocities are stored in
! arrays X, V with the format (x,y,z) and (vx,vy,vz) for each object in
! succession. The accelerations are stored in the array A (ax,ay,az). The
! non-gravitational accelerations follow a force law described by Marsden
! et al. (1973) Astron. J. 211-225, with magnitude determined by the
! parameters NGF(1,2,3) for each object.

! N.B. All coordinates and velocities must be with respect to central body!!!
! ===
!------------------------------------------------------------------------------

subroutine mfo_ngf (nbod,x,v,a,ngf)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  real(double_precision), intent(in) :: ngf(4,nbod) !< [in] non gravitational forces parameters
  !! \n(1-3) cometary non-gravitational (jet) force parameters
  !! \n(4)  beta parameter for radiation pressure and P-R drag
  
  ! Output
  real(double_precision), intent(out) :: a(3,nbod)
  
  ! Local
  integer :: j
  real(double_precision) :: r2,r,rv,q,g,tx,ty,tz,nx,ny,nz,a1,a2,a3
  
  !------------------------------------------------------------------------------
  
  do j = 2, nbod
     r2 = x(1,j)*x(1,j) + x(2,j)*x(2,j) +x(3,j)*x(3,j)
     
     ! Only calculate accelerations if body is close to the Sun (R < 9.36 AU), 
     ! or if the non-gravitational force parameters are exceptionally large.
     if (r2.lt.88.d0.or.abs(ngf(1,j)).gt.1d-7.or.abs(ngf(2,j)).gt.1d-7.or.abs(ngf(3,j)).gt.1d-7) then
        r = sqrt(r2)
        rv = x(1,j)*v(1,j) + x(2,j)*v(2,j) + x(3,j)*v(3,j)
        
        ! Calculate Q = R / R0, where R0 = 2.808 AU
        q = r * .3561253561253561d0
        g = .111262d0 * q**(-2.15d0) * (1.d0+q**5.093d0)**(-4.6142d0)
        
        ! Within-orbital-plane transverse vector components
        tx = r2*v(1,j) - rv*x(1,j)
        ty = r2*v(2,j) - rv*x(2,j)
        tz = r2*v(3,j) - rv*x(3,j)
        
        ! Orbit-normal vector components
        nx = x(2,j)*v(3,j) - x(3,j)*v(2,j)
        ny = x(3,j)*v(1,j) - x(1,j)*v(3,j)
        nz = x(1,j)*v(2,j) - x(2,j)*v(1,j)
        
        ! Multiplication factors
        a1 = ngf(1,j) * g / r
        a2 = ngf(2,j) * g / sqrt(tx*tx + ty*ty + tz*tz)
        a3 = ngf(3,j) * g / sqrt(nx*nx + ny*ny + nz*nz)
        
        ! X,Y and Z components of non-gravitational acceleration
        a(1,j) = a1*x(1,j) + a2*tx + a3*nx
        a(2,j) = a1*x(2,j) + a2*ty + a3*ny
        a(3,j) = a1*x(3,j) + a2*tz + a3*nz
     else
        a(1,j) = 0.d0
        a(2,j) = 0.d0
        a(3,j) = 0.d0
     end if
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_ngf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MFO_PN.FOR    (ErikSoft   3 October 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers
! TODO
! ****** To be completed at a later date ******

! Calculates post-Newtonian relativistic corrective accelerations for a set
! of NBOD bodies (NBIG of which are Big).

! This routine should not be called from the symplectic algorithm MAL_MVS
! or the conservative Bulirsch-Stoer algorithm MAL_BS2.

! N.B. All coordinates and velocities must be with respect to central body!!!
! ===
!------------------------------------------------------------------------------

subroutine mfo_pn (nbod,nbig,m,x,v,a)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  
  ! Output
  real(double_precision), intent(out) :: a(3,nbod)
  
  ! Local
  integer :: j
  
  !------------------------------------------------------------------------------
  
  do j = 1, nbod
     a(1,j) = 0.d0
     a(2,j) = 0.d0
     a(3,j) = 0.d0
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_pn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MFO_PR.FOR    (ErikSoft   3 October 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers
! TODO
! ****** To be completed at a later date ******

! Calculates radiation pressure and Poynting-Robertson drag for a set
! of NBOD bodies (NBIG of which are Big).

! This routine should not be called from the symplectic algorithm MAL_MVS
! or the conservative Bulirsch-Stoer algorithm MAL_BS2.

! N.B. All coordinates and velocities must be with respect to central body!!!
! ===
!------------------------------------------------------------------------------

subroutine mfo_pr (nbod,nbig,m,x,v,a,ngf)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  real(double_precision), intent(in) :: ngf(4,nbod) !< [in] non gravitational forces parameters
  !! \n(1-3) cometary non-gravitational (jet) force parameters
  !! \n(4)  beta parameter for radiation pressure and P-R drag
  
  ! Output
  real(double_precision), intent(out) :: a(3,nbod)
  
  ! Local
  integer :: j
  
  !------------------------------------------------------------------------------
  
  do j = 1, nbod
     a(1,j) = 0.d0
     a(2,j) = 0.d0
     a(3,j) = 0.d0
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_pr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!     MFO_OBL.FOR    (ErikSoft   2 October 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates barycentric accelerations of NBOD bodies due to oblateness of
! the central body. Also returns the corresponding barycentric acceleration
! of the central body.

! N.B. All coordinates must be with respect to the central body!!!
! ===
!------------------------------------------------------------------------------

subroutine mfo_obl (jcen,nbod,m,x,a,acen)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input
  integer, intent(in) :: nbod
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  
  ! Output
  real(double_precision), intent(out) :: acen(3)
  real(double_precision), intent(out) :: a(3,nbod)
  
  ! Local
  integer :: i
  real(double_precision) :: jr2,jr4,jr6,r2,r_1,r_2,r_3,u2,u4,u6,tmp1,tmp2,tmp3,tmp4
  
  !------------------------------------------------------------------------------
  
  acen(1) = 0.d0
  acen(2) = 0.d0
  acen(3) = 0.d0
  
  do i = 2, nbod
     
     ! Calculate barycentric accelerations on the objects
     r2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
     r_1 = 1.d0 / sqrt(r2)
     r_2 = r_1 * r_1
     r_3 = r_2 * r_1
     jr2 = jcen(1) * r_2
     jr4 = jcen(2) * r_2 * r_2
     jr6 = jcen(3) * r_2 * r_2 * r_2
     u2 = x(3,i) * x(3,i) * r_2
     u4 = u2 * u2
     u6 = u4 * u2
     
     tmp1 = m(1) * r_3
     tmp2 = jr2*(7.5d0*u2 - 1.5d0) + jr4*(39.375d0*u4 - 26.25d0*u2 + 1.875d0) + &
          jr6*(187.6875d0*u6 -216.5625d0*u4 +59.0625d0*u2 -2.1875d0)
     tmp3 = jr2*3.d0 + jr4*(17.5d0*u2 - 7.5d0) + jr6*(86.625d0*u4 - 78.75d0*u2 + 13.125d0)
     
     a(1,i) = x(1,i) * tmp1 * tmp2
     a(2,i) = x(2,i) * tmp1 * tmp2
     a(3,i) = x(3,i) * tmp1 * (tmp2 - tmp3)
     
     ! Calculate barycentric accelerations on the central body
     tmp4 = m(i) / m(1)
     acen(1) = acen(1)  -  tmp4 * a(1,i)
     acen(2) = acen(2)  -  tmp4 * a(2,i)
     acen(3) = acen(3)  -  tmp4 * a(3,i)
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mfo_obl
  
end module forces
