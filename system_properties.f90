module system_properties

!*************************************************************
!** Modules that compute various properties of the system like
!** the total energy and angular momentum, jocobi constants, 
!** hill radii, if there are ejections and so on.
!** Version 1.0 - june 2011
!*************************************************************
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCE_HILL.FOR    (ErikSoft   4 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates the Hill radii for all objects given their masses, M,
! coordinates, X, and velocities, V; plus the mass of the central body, M(1)
! Where HILL = a * (m/3*m(1))^(1/3)
!
! If the orbit is hyperbolic or parabolic, the Hill radius is calculated using:
!       HILL = r * (m/3*m(1))^(1/3)
! where R is the current distance from the central body.
!
! The routine also gives the semi-major axis, A, of each object's orbit.
!
! N.B. Designed to use heliocentric coordinates, but should be adequate using
! ===  barycentric coordinates.
!
!------------------------------------------------------------------------------
!
subroutine mce_hill (nbod,m,x,v,hill,a)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use orbital_elements

  implicit none

  real(double_precision) :: THIRD
  parameter (THIRD = .3333333333333333d0)
  !
  ! Input/Output
  integer :: nbod
  real(double_precision) :: m(nbod),x(3,nbod),v(3,nbod),hill(nbod),a(nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: r, v2, gm
  !
  !------------------------------------------------------------------------------
  !
  do j = 2, nbod
     gm = m(1) + m(j)
     call mco_x2a (gm,x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),a(j),r,v2)
     ! If orbit is hyperbolic, use the distance rather than the semi-major axis
     if (a(j).le.0) a(j) = r
     hill(j) = a(j) * (THIRD * m(j) / m(1))**THIRD
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mce_hill
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCE_INIT.FOR    (ErikSoft   28 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates close-approach limits RCE (in AU) and physical radii RPHYS
! (in AU) for all objects, given their masses M, coordinates X, velocities
! V, densities RHO, and close-approach limits RCEH (in Hill radii).
!
! Also calculates the changeover distance RCRIT, used by the hybrid
! symplectic integrator. RCRIT is defined to be the larger of N1*HILL and
! N2*H*VMAX, where HILL is the Hill radius, H is the timestep, VMAX is the
! largest expected velocity of any body, and N1, N2 are parameters (see
! section 4.2 of Chambers 1999, Monthly Notices, vol 304, p793-799).
!
! N.B. Designed to use heliocentric coordinates, but should be adequate using
! ===  barycentric coordinates.
!
!------------------------------------------------------------------------------
!
subroutine mce_init (tstart,algor,h,jcen,rcen,rmax,cefac,nbod,nbig,m,x,v,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile,rcritflag)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use ascii_conversion

  implicit none

  !
  real(double_precision) :: N2,THIRD
  parameter (N2=.4d0,THIRD=.3333333333333333d0)
  !
  ! Input/Output
  integer :: nbod,nbig,algor,opt(8),rcritflag
  real(double_precision) :: tstart,h,jcen(3),rcen,rmax,cefac,m(nbod),x(3,nbod)
  real(double_precision) :: v(3,nbod),s(3,nbod),rho(nbod),rceh(nbod),rphys(nbod)
  real(double_precision) :: rce(nbod),rcrit(nbod)
  character*8 id(nbod)
  character*80 outfile
  !
  ! Local
  integer :: j, error
  real(double_precision) :: a(NMAX),hill(NMAX),temp,amin,vmax,k_2,rhocgs,rcen_2
  character*80 header,c(NMAX)
  !
  !------------------------------------------------------------------------------
  !
  rhocgs = AU * AU * AU * K2 / MSUN
  k_2 = 1.d0 / K2
  rcen_2 = 1.d0 / (rcen * rcen)
  amin = HUGE
  !
  ! Calculate the Hill radii
  call mce_hill (nbod,m,x,v,hill,a)
  !
  ! Determine the maximum close-encounter distances, and the physical radii
  temp = 2.25d0 * m(1) / PI
  do j = 2, nbod
     rce(j)   = hill(j) * rceh(j)
     rphys(j) = hill(j) / a(j) * (temp/rho(j))**THIRD
     amin = min (a(j), amin)
  end do
  !
  ! If required, calculate the changeover distance used by hybrid algorithm
  if (rcritflag.eq.1) then
     vmax = sqrt (m(1) / amin)
     temp = N2 * h * vmax
     do j = 2, nbod
        rcrit(j) = max(hill(j) * cefac, temp)
     end do
  end if
  !
  ! Write list of object's identities to close-encounter output file
  header(1:8)   = mio_fl2c (tstart)
  header(9:16)  = mio_re2c (dble(nbig - 1),   0.d0, 11239423.99d0)
  header(12:19) = mio_re2c (dble(nbod - nbig),0.d0, 11239423.99d0)
  header(15:22) = mio_fl2c (m(1) * k_2)
  header(23:30) = mio_fl2c (jcen(1) * rcen_2)
  header(31:38) = mio_fl2c (jcen(2) * rcen_2 * rcen_2)
  header(39:46) = mio_fl2c (jcen(3) * rcen_2 * rcen_2 * rcen_2)
  header(47:54) = mio_fl2c (rcen)
  header(55:62) = mio_fl2c (rmax)
  !
  do j = 2, nbod
     c(j)(1:8) = mio_re2c (dble(j - 1), 0.d0, 11239423.99d0)
     c(j)(4:11) = id(j)
     c(j)(12:19) = mio_fl2c (m(j) * k_2)
     c(j)(20:27) = mio_fl2c (s(1,j) * k_2)
     c(j)(28:35) = mio_fl2c (s(2,j) * k_2)
     c(j)(36:43) = mio_fl2c (s(3,j) * k_2)
     c(j)(44:51) = mio_fl2c (rho(j) / rhocgs)
  end do
  !
  ! Write compressed output to file
  open (22, file=outfile, status='old', access='append', iostat=error)
  if (error /= 0) then
     write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile)
     stop
  end if
  write (22,'(a1,a2,i2,a62,i1)') char(12),'6a',algor,header(1:62),opt(4)
  do j = 2, nbod
     write (22,'(a51)') c(j)(1:51)
  end do
  close (22)
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mce_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MXX_EN.FOR    (ErikSoft   21 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates the total energy and angular-momentum for a system of objects
! with masses M, coordinates X, velocities V and spin angular momenta S.
!
! N.B. All coordinates and velocities must be with respect to the central
! ===  body.
!
!------------------------------------------------------------------------------
!
subroutine mxx_en  (jcen,nbod,nbig,m,xh,vh,s,e,l2)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig
  real(double_precision) :: jcen(3),m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),e,l2
  !
  ! Local
  integer :: j,k,iflag,itmp(8)
  real(double_precision) :: x(3,NMAX),v(3,NMAX),temp,dx,dy,dz,r2,tmp,ke,pe,l(3)
  real(double_precision) :: r_1,r_2,r_4,r_6,u2,u4,u6,tmp2(4,NMAX)
  !
  !------------------------------------------------------------------------------
  !
  ke = 0.d0
  pe = 0.d0
  l(1) = 0.d0
  l(2) = 0.d0
  l(3) = 0.d0
  !
  ! Convert to barycentric coordinates and velocities
  call mco_h2b(temp,jcen,nbod,nbig,temp,m,xh,vh,x,v,tmp2,iflag,itmp)
  !
  ! Do the spin angular momenta first (probably the smallest terms)
  do j = 1, nbod
     l(1) = l(1) + s(1,j)
     l(2) = l(2) + s(2,j)
     l(3) = l(3) + s(3,j)
  end do
  !
  ! Orbital angular momentum and kinetic energy terms
  do j = 1, nbod
     l(1) = l(1)  +  m(j)*(x(2,j) * v(3,j)  -  x(3,j) * v(2,j))
     l(2) = l(2)  +  m(j)*(x(3,j) * v(1,j)  -  x(1,j) * v(3,j))
     l(3) = l(3)  +  m(j)*(x(1,j) * v(2,j)  -  x(2,j) * v(1,j))
     ke = ke + m(j)*(v(1,j)*v(1,j)+v(2,j)*v(2,j)+v(3,j)*v(3,j))
  end do
  !
  ! Potential energy terms due to pairs of bodies
  do j = 2, nbod
     tmp = 0.d0
     do k = j + 1, nbod
        dx = x(1,k) - x(1,j)
        dy = x(2,k) - x(2,j)
        dz = x(3,k) - x(3,j)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2.ne.0) tmp = tmp + m(k) / sqrt(r2)
     end do
     pe = pe  -  tmp * m(j)
  end do
  !
  ! Potential energy terms involving the central body
  do j = 2, nbod
     dx = x(1,j) - x(1,1)
     dy = x(2,j) - x(2,1)
     dz = x(3,j) - x(3,1)
     r2 = dx*dx + dy*dy + dz*dz
     if (r2.ne.0) pe = pe  -  m(1) * m(j) / sqrt(r2)
  end do
  !
  ! Corrections for oblateness
  if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
     do j = 2, nbod
        r2 = xh(1,j)*xh(1,j) + xh(2,j)*xh(2,j) + xh(3,j)*xh(3,j)
        r_1 = 1.d0 / sqrt(r2)
        r_2 = r_1 * r_1
        r_4 = r_2 * r_2
        r_6 = r_4 * r_2
        u2 = xh(3,j) * xh(3,j) * r_2
        u4 = u2 * u2
        u6 = u4 * u2
        pe = pe + m(1) * m(j) * r_1 * (jcen(1) * r_2 * (1.5d0*u2 - 0.5d0) +  jcen(2) * r_4 * (4.375d0*u4 - 3.75d0*u2 + .375d0)&
             +  jcen(3) * r_6 *(14.4375d0*u6 - 19.6875d0*u4 + 6.5625d0*u2 - .3125d0))
     end do
  end if
  !
  e = .5d0 * ke  +  pe
  l2 = sqrt(l(1)*l(1) + l(2)*l(2) + l(3)*l(3))
  !
  !------------------------------------------------------------------------------
  !
  return	
end subroutine mxx_en
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MXX_JAC.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates the Jacobi constant for massless particles. This assumes that
! there are only 2 massive bodies (including the central body) moving on
! circular orbits.
!
! N.B. All coordinates and velocities must be heliocentric!!
! ===
!
!------------------------------------------------------------------------------
!
subroutine mxx_jac (jcen,nbod,nbig,m,xh,vh,jac)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig
  real(double_precision) :: jcen(3),m(nbod),xh(3,nbod),vh(3,nbod)
  !
  ! Local
  integer :: j,itmp(8),iflag
  real(double_precision) :: x(3,NMAX),v(3,NMAX),temp,dx,dy,dz,r,d,a2,n,jac(NMAX)
  real(double_precision) :: tmp2(4,NMAX)
  !
  !------------------------------------------------------------------------------
  !
  call mco_h2b(temp,jcen,nbod,nbig,temp,m,xh,vh,x,v,tmp2,iflag,itmp)
  dx = x(1,2) - x(1,1)
  dy = x(2,2) - x(2,1)
  dz = x(3,2) - x(3,1)
  a2 = dx*dx + dy*dy + dz*dz
  n = sqrt((m(1)+m(2)) / (a2*sqrt(a2)))
  !
  do j = nbig + 1, nbod
     dx = x(1,j) - x(1,1)
     dy = x(2,j) - x(2,1)
     dz = x(3,j) - x(3,1)
     r = sqrt(dx*dx + dy*dy + dz*dz)
     dx = x(1,j) - x(1,2)
     dy = x(2,j) - x(2,2)
     dz = x(3,j) - x(3,2)
     d = sqrt(dx*dx + dy*dy + dz*dz)
     !
     jac(j) = .5d0*(v(1,j)*v(1,j) + v(2,j)*v(2,j) + v(3,j)*v(3,j))     - m(1)/r - m(2)/d - n*(x(1,j)*v(2,j) - x(2,j)*v(1,j))
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mxx_jac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MXX_EJEC.FOR    (ErikSoft   2 November 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates the distance from the central body of each object with index
! I >= I0. If this distance exceeds RMAX, the object is flagged for ejection 
! (STAT set to -3). If any object is to be ejected, EJFLAG = 1 on exit,
! otherwise EJFLAG = 0.
!
! Also updates the values of EN(3) and AM(3)---the change in energy and
! angular momentum due to collisions and ejections.
!
!
! N.B. All coordinates must be with respect to the central body!!
! ===
!
!------------------------------------------------------------------------------
!
subroutine mxx_ejec (time,tstart,rmax,en,am,jcen,i0,nbod,nbig,m,x,v,s,stat,id,opt,ejflag,outfile,mem,lmem)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: i0, nbod, nbig, stat(nbod), opt(8), ejflag, lmem(NMESS)
  real(double_precision) :: time, tstart, rmax, en(3), am(3), jcen(3)
  real(double_precision) :: m(nbod), x(3,nbod), v(3,nbod), s(3,nbod)
  character*80 outfile, mem(NMESS)
  character*8 id(nbod)
  !
  ! Local
  integer :: j, year, month
  real(double_precision) :: r2,rmax2,t1,e,l
  character*38 flost
  character*6 tstring
  integer :: error
  !
  !------------------------------------------------------------------------------
  !
  if (i0.le.0) i0 = 2
  ejflag = 0
  rmax2 = rmax * rmax
  !
  ! Calculate initial energy and angular momentum
  call mxx_en (jcen,nbod,nbig,m,x,v,s,e,l)
  !
  ! Flag each object which is ejected, and set its mass to zero
  do j = i0, nbod
     r2 = x(1,j)*x(1,j) + x(2,j)*x(2,j) + x(3,j)*x(3,j)
     if (r2.gt.rmax2) then
        ejflag = 1
        stat(j) = -3
        m(j) = 0.d0
        s(1,j) = 0.d0
        s(2,j) = 0.d0
        s(3,j) = 0.d0
        !
        ! Write message to information file
        open  (23,file=outfile,status='old',access='append',iostat=error)
        if (error /= 0) then
           write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile)
           stop
        end if
        if (opt(3).eq.1) then
           call mio_jd2y (time,year,month,t1)
           flost = '(1x,a8,a,i10,1x,i2,1x,f8.5)'
           write (23,flost) id(j),mem(68)(1:lmem(68)),year,month,t1
        else
           if (opt(3).eq.3) then
              t1 = (time - tstart) / 365.25d0
              tstring = mem(2)
              flost = '(1x,a8,a,f18.7,a)'
           else
              if (opt(3).eq.0) t1 = time
              if (opt(3).eq.2) t1 = time - tstart
              tstring = mem(1)
              flost = '(1x,a8,a,f18.5,a)'
           end if
           write (23,flost) id(j),mem(68)(1:lmem(68)),t1,tstring
        end if
        close (23)
     end if
  end do
  !
  ! If ejections occurred, update ELOST and LLOST
  if (ejflag.ne.0) then
     call mxx_en (jcen,nbod,nbig,m,x,v,s,en(2),am(2))
     en(3) = en(3) + (e - en(2))
     am(3) = am(3) + (l - am(2))
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mxx_ejec


end module system_properties
