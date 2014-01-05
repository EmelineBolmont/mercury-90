!******************************************************************************
! MODULE: system_properties
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that compute various properties of the system like
!! the total energy and angular momentum, jacobi constants, 
!! hill radii, if there are ejections and so on.
!
!******************************************************************************

module system_properties

  use types_numeriques
  use mercury_globals

  implicit none
  
  private m_sfunc
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCE_HILL.FOR    (ErikSoft   4 October 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates the Hill radii for all objects given their masses, M,
! coordinates, X, and velocities, V; plus the mass of the central body, M(1)
! Where HILL = a * (m/3*m(1))^(1/3)

! If the orbit is hyperbolic or parabolic, the Hill radius is calculated using:
!       HILL = r * (m/3*m(1))^(1/3)
! where R is the current distance from the central body.

! The routine also gives the semi-major axis, A, of each object's orbit.

! N.B. Designed to use heliocentric coordinates, but should be adequate using
! ===  barycentric coordinates.

!------------------------------------------------------------------------------

subroutine mce_hill (nbod,m,x,v,hill,a)
  
  use physical_constant
  use mercury_constant
  use orbital_elements

  implicit none
  
  ! Input/Output
  integer, intent(in) :: nbod
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  real(double_precision), intent(out) :: hill(nbod)
  real(double_precision), intent(out) :: a(nbod)
  
  ! Local
  integer :: j
  real(double_precision) :: r, v2, gm
  
  !------------------------------------------------------------------------------
  
  do j = 2, nbod
     gm = m(1) + m(j)
     call mco_x2a (gm,x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),a(j),r,v2)
     ! If orbit is hyperbolic, use the distance rather than the semi-major axis
     if (a(j).le.0) a(j) = r
     hill(j) = a(j) * (THIRD * m(j) / m(1))**THIRD
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_hill

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCE_INIT.FOR    (ErikSoft   28 February 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates close-approach limits RCE (in AU) and physical radii RPHYS
! (in AU) for all objects, given their masses M, coordinates X, velocities
! V, densities RHO, and close-approach limits RCEH (in Hill radii).

! Also calculates the changeover distance RCRIT, used by the hybrid
! symplectic integrator. RCRIT is defined to be the larger of N1*HILL and
! N2*H*VMAX, where HILL is the Hill radius, H is the timestep, VMAX is the
! largest expected velocity of any body, and N1, N2 are parameters (see
! section 4.2 of Chambers 1999, Monthly Notices, vol 304, p793-799).

! N.B. Designed to use heliocentric coordinates, but should be adequate using
! ===  barycentric coordinates.

!------------------------------------------------------------------------------

subroutine mce_init (h,jcen,rcen,cefac,nbod,nbig,m,x,v,s,rho,rceh,rphys,rce,rcrit,id,outfile,rcritflag)
  
  use physical_constant
  use mercury_constant
  use ascii_conversion

  implicit none

  
  real(double_precision), parameter :: N2=.4d0
  
  ! Input/Output
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  integer, intent(in) :: rcritflag
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  real(double_precision), intent(in) :: rho(nbod) !< [in] physical density (g/cm^3)
  real(double_precision), intent(in) :: rceh(nbod) !< [in] close-encounter limit (Hill radii)
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: s(3,nbod) !< [in] spin angular momentum (solar masses AU^2/day)
  real(double_precision), intent(in) :: h
  real(double_precision), intent(in) :: rcen
  real(double_precision), intent(in) :: cefac
  character(len=8), intent(in) :: id(nbod)
  character(len=80), intent(in) :: outfile
  
  real(double_precision), intent(out) :: rce(nbod)
  real(double_precision), intent(out) :: rphys(nbod)
  real(double_precision), intent(out) :: rcrit(nbod)

  
  ! Local
  integer :: j, error
  real(double_precision) :: a(nb_bodies_initial),hill(nb_bodies_initial),temp,amin,vmax,k_2,rhocgs,rcen_2
  character(len=80) :: header,c(nb_bodies_initial)
  
  !------------------------------------------------------------------------------
  
  rhocgs = AU * AU * AU * K2 / MSUN
  k_2 = 1.d0 / K2
  rcen_2 = 1.d0 / (rcen * rcen)
  amin = HUGE
  
  ! Calculate the Hill radii
  call mce_hill (nbod,m,x,v,hill,a)
  
  ! Determine the maximum close-encounter distances, and the physical radii
  temp = 2.25d0 * m(1) / PI
  do j = 2, nbod
     rce(j)   = hill(j) * rceh(j)
     rphys(j) = hill(j) / a(j) * (temp/rho(j))**THIRD
     amin = min (a(j), amin)
  end do
  
  ! If required, calculate the changeover distance used by hybrid algorithm
  if (rcritflag.eq.1) then
     vmax = sqrt (m(1) / amin)
     temp = N2 * h * vmax
     do j = 2, nbod
        rcrit(j) = max(hill(j) * cefac, temp)
     end do
  end if
  
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
  
  do j = 2, nbod
     c(j)(1:8) = mio_re2c (dble(j - 1), 0.d0, 11239423.99d0)
     c(j)(4:11) = id(j)
     c(j)(12:19) = mio_fl2c (m(j) * k_2)
     c(j)(20:27) = mio_fl2c (s(1,j) * k_2)
     c(j)(28:35) = mio_fl2c (s(2,j) * k_2)
     c(j)(36:43) = mio_fl2c (s(3,j) * k_2)
     c(j)(44:51) = mio_fl2c (rho(j) / rhocgs)
  end do
  
  ! Write compressed output to file
  open (22, file=outfile, status='old', position='append', iostat=error)
  if (error /= 0) then
     write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile)
     stop
  end if
  write (22,'(a1,a2,i2,a62,i1)') char(12),'6a',algor,header(1:62),opt(4)
  do j = 2, nbod
     write (22,'(a51)') c(j)(1:51)
  end do
  close (22)
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!     MXX_EN.FOR    (ErikSoft   21 February 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates the total energy and angular-momentum for a system of objects
! with masses M, coordinates X, velocities V and spin angular momenta S.

! N.B. All coordinates and velocities must be with respect to the central
! ===  body.

!------------------------------------------------------------------------------

subroutine mxx_en  (jcen,nbod,nbig,m,xh,vh,s,e,l2)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input/Output
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: xh(3,nbod) !< [in] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision), intent(in) :: vh(3,nbod) !< [in] velocities (vx,vy,vz) with respect to the central body (AU/day)
  real(double_precision), intent(in) :: s(3,nbod) !< [in] spin angular momentum (solar masses AU^2/day)
  
  real(double_precision), intent(out) :: e ! energy
  real(double_precision), intent(out) :: l2 ! angular momentum
  
  ! Local
  integer :: j,k,iflag,itmp(8)
  real(double_precision) :: x(3,nb_bodies_initial),v(3,nb_bodies_initial),temp,dx,dy,dz,r2,tmp,ke,pe,l(3)
  real(double_precision) :: r_1,r_2,r_4,r_6,u2,u4,u6,tmp2(4,nb_bodies_initial)
  
  !------------------------------------------------------------------------------
  
  ke = 0.d0
  pe = 0.d0
  l(1) = 0.d0
  l(2) = 0.d0
  l(3) = 0.d0
  
  ! Convert to barycentric coordinates and velocities
  call mco_h2b(jcen,nbod,nbig,temp,m,xh,vh,x,v)
  
  ! Do the spin angular momenta first (probably the smallest terms)
  do j = 1, nbod
     l(1) = l(1) + s(1,j)
     l(2) = l(2) + s(2,j)
     l(3) = l(3) + s(3,j)
  end do
  
  ! Orbital angular momentum and kinetic energy terms
  do j = 1, nbod
     l(1) = l(1)  +  m(j)*(x(2,j) * v(3,j)  -  x(3,j) * v(2,j))
     l(2) = l(2)  +  m(j)*(x(3,j) * v(1,j)  -  x(1,j) * v(3,j))
     l(3) = l(3)  +  m(j)*(x(1,j) * v(2,j)  -  x(2,j) * v(1,j))
     ke = ke + m(j)*(v(1,j)*v(1,j)+v(2,j)*v(2,j)+v(3,j)*v(3,j))
  end do
  
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
  
  ! Potential energy terms involving the central body
  do j = 2, nbod
     dx = x(1,j) - x(1,1)
     dy = x(2,j) - x(2,1)
     dz = x(3,j) - x(3,1)
     r2 = dx*dx + dy*dy + dz*dz
     if (r2.ne.0) pe = pe  -  m(1) * m(j) / sqrt(r2)
  end do
  
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
  
  e = .5d0 * ke  +  pe
  l2 = sqrt(l(1)*l(1) + l(2)*l(2) + l(3)*l(3))
  
  !------------------------------------------------------------------------------
  
  return  
end subroutine mxx_en

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!     MXX_JAC.FOR    (ErikSoft   2 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates the Jacobi constant for massless particles. This assumes that
! there are only 2 massive bodies (including the central body) moving on
! circular orbits.

! N.B. All coordinates and velocities must be heliocentric!
! ===

!------------------------------------------------------------------------------

subroutine mxx_jac (jcen,nbod,nbig,m,xh,vh,jac)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input/Output
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: xh(3,nbod) !< [in] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision), intent(in) :: vh(3,nbod) !< [in] velocities (vx,vy,vz) with respect to the central body (AU/day)
  
  real(double_precision), intent(out) :: jac(nbod)
  
  ! Local
  integer :: j,itmp(8),iflag
  real(double_precision) :: x(3,nb_bodies_initial),v(3,nb_bodies_initial),temp,dx,dy,dz,r,d,a2,n
  real(double_precision) :: tmp2(4,nb_bodies_initial)
  
  !------------------------------------------------------------------------------
  
  call mco_h2b(jcen,nbod,nbig,temp,m,xh,vh,x,v)
  dx = x(1,2) - x(1,1)
  dy = x(2,2) - x(2,1)
  dz = x(3,2) - x(3,1)
  a2 = dx*dx + dy*dy + dz*dz
  n = sqrt((m(1)+m(2)) / (a2*sqrt(a2)))
  
  do j = nbig + 1, nbod
     dx = x(1,j) - x(1,1)
     dy = x(2,j) - x(2,1)
     dz = x(3,j) - x(3,1)
     r = sqrt(dx*dx + dy*dy + dz*dz)
     dx = x(1,j) - x(1,2)
     dy = x(2,j) - x(2,2)
     dz = x(3,j) - x(3,2)
     d = sqrt(dx*dx + dy*dy + dz*dz)
     
     jac(j) = .5d0*(v(1,j)*v(1,j) + v(2,j)*v(2,j) + v(3,j)*v(3,j))     - m(1)/r - m(2)/d - n*(x(1,j)*v(2,j) - x(2,j)*v(1,j))
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mxx_jac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MXX_EJEC.FOR    (ErikSoft   2 November 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates the distance from the central body of each object with index
! I >= I0. If this distance exceeds RMAX, the object is flagged for ejection 
! (STAT set to -3). If any object is to be ejected, EJFLAG = 1 on exit,
! otherwise EJFLAG = 0.

! Also updates the values of EN(3) and AM(3)---the change in energy and
! angular momentum due to collisions and ejections.


! N.B. All coordinates must be with respect to the central body!
! ===

!------------------------------------------------------------------------------

subroutine mxx_ejec (time,en,am,jcen,i0,nbod,nbig,m,x,v,s,stat,id,ejflag,outfile)
  
  use physical_constant
  use mercury_constant
  use utilities

  implicit none

  
  ! Input/Output
  integer, intent(in) :: i0, nbod, nbig

  real(double_precision), intent(in) :: time
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: x(3,nbod)
  real(double_precision), intent(in) :: v(3,nbod)
  character(len=80), intent(in) :: outfile
  character(len=8), intent(in) :: id(nbod)
  
  real(double_precision), intent(out) :: en(3)
  real(double_precision), intent(out) :: am(3)
  integer, intent(out) :: ejflag
  integer, intent(out) :: stat(nbod)
  
  real(double_precision), intent(inout) :: m(nbod) !< [in,out] mass (in solar masses * K2)
  real(double_precision), intent(inout) :: s(3,nbod) !< [in,out] spin angular momentum (solar masses AU^2/day)

  ! Local
  integer :: j,j0, year, month
  real(double_precision) :: r2,rmax2,t1,e,l
  character(len=38) :: flost
  character(len=6) :: tstring
  integer :: error
  
  !------------------------------------------------------------------------------
  
  if (i0.le.0) then
    j0 = 2
  else
    j0 = i0
  end if
  
  ejflag = 0
  rmax2 = rmax * rmax
  
  ! Calculate initial energy and angular momentum
  call mxx_en (jcen,nbod,nbig,m,x,v,s,e,l)
  
  ! Flag each object which is ejected, and set its mass to zero
  do j = j0, nbod
     r2 = x(1,j)*x(1,j) + x(2,j)*x(2,j) + x(3,j)*x(3,j)
     if (r2.gt.rmax2) then
        ejflag = 1
        stat(j) = -3
        m(j) = 0.d0
        s(1,j) = 0.d0
        s(2,j) = 0.d0
        s(3,j) = 0.d0
        
        ! Write message to information file
        open  (23,file=outfile,status='old',position='append',iostat=error)
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
  
  ! If ejections occurred, update ELOST and LLOST
  if (ejflag.ne.0) then
     call mxx_en (jcen,nbod,nbig,m,x,v,s,en(2),am(2))
     en(3) = en(3) + (e - en(2))
     am(3) = am(3) + (l - am(2))
  end if
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mxx_ejec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCO_B2H.FOR    (ErikSoft   2 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Converts barycentric coordinates to coordinates with respect to the central
! body.

!------------------------------------------------------------------------------

subroutine mco_b2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag)
  

  implicit none

  
  ! Input/Output
  integer,intent(in) :: nbod
  integer,intent(in) :: nbig
  integer,intent(in) :: ngflag
  real(double_precision),intent(in) :: time
  real(double_precision),intent(in) :: h
  real(double_precision),intent(in) :: jcen(3)
  real(double_precision),intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision),intent(in) :: x(3,nbod)
  real(double_precision),intent(in) :: v(3,nbod)
  real(double_precision),intent(in) :: ngf(4,nbod)
  
  real(double_precision),intent(out) :: xh(3,nbod) !< [out] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision),intent(out) :: vh(3,nbod) !< [out] velocities (vx,vy,vz) with respect to the central body (AU/day)
  
  ! Local
  integer :: j
  
  !------------------------------------------------------------------------------
  
  do j = 2, nbod
     xh(1,j) = x(1,j) - x(1,1)
     xh(2,j) = x(2,j) - x(2,1)
     xh(3,j) = x(3,j) - x(3,1)
     vh(1,j) = v(1,j) - v(1,1)
     vh(2,j) = v(2,j) - v(2,1)
     vh(3,j) = v(3,j) - v(3,1)
  enddo
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mco_b2h

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCO_H2B.FOR    (ErikSoft   2 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Converts coordinates with respect to the central body to barycentric
! coordinates.

!------------------------------------------------------------------------------

subroutine mco_h2b (jcen,nbod,nbig,h,m,xh,vh,x,v)
  

  implicit none

  
  ! Input/Output
  integer,intent(in) :: nbod
  integer,intent(in) :: nbig
  real(double_precision),intent(in) :: h
  real(double_precision),intent(in) :: jcen(3)
  real(double_precision),intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision),intent(in) :: xh(3,nbod) !< [in] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision),intent(in) :: vh(3,nbod) !< [in] velocities (vx,vy,vz) with respect to the central body (AU/day)
  real(double_precision),intent(out) :: x(3,nbod)
  real(double_precision),intent(out) :: v(3,nbod)
  
  ! Local
  integer :: j
  real(double_precision) :: mtot,temp
  
  !------------------------------------------------------------------------------
  
  mtot = 0.d0
  x(1,1) = 0.d0
  x(2,1) = 0.d0
  x(3,1) = 0.d0
  v(1,1) = 0.d0
  v(2,1) = 0.d0
  v(3,1) = 0.d0
  
  ! Calculate coordinates and velocities of the central body
  do j = 2, nbod
     mtot = mtot  +  m(j)
     x(1,1) = x(1,1)  +  m(j) * xh(1,j)
     x(2,1) = x(2,1)  +  m(j) * xh(2,j)
     x(3,1) = x(3,1)  +  m(j) * xh(3,j)
     v(1,1) = v(1,1)  +  m(j) * vh(1,j)
     v(2,1) = v(2,1)  +  m(j) * vh(2,j)
     v(3,1) = v(3,1)  +  m(j) * vh(3,j)
  enddo
  
  temp = -1.d0 / (mtot + m(1))
  x(1,1) = temp * x(1,1)
  x(2,1) = temp * x(2,1)
  x(3,1) = temp * x(3,1)
  v(1,1) = temp * v(1,1)
  v(2,1) = temp * v(2,1)
  v(3,1) = temp * v(3,1)
  
  ! Calculate the barycentric coordinates and velocities
  do j = 2, nbod
     x(1,j) = xh(1,j) + x(1,1)
     x(2,j) = xh(2,j) + x(2,1)
     x(3,j) = xh(3,j) + x(3,1)
     v(1,j) = vh(1,j) + v(1,1)
     v(2,j) = vh(2,j) + v(2,1)
     v(3,j) = vh(3,j) + v(3,1)
  enddo
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mco_h2b

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCO_H2CB.FOR    (ErikSoft   2 November 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Convert coordinates with respect to the central body to close-binary
! coordinates.

!------------------------------------------------------------------------------

subroutine mco_h2cb (jcen,nbod,nbig,h,m,xh,vh,x,v)
  

  implicit none

  
  ! Input/Output
  integer,intent(in) :: nbod
  integer,intent(in) :: nbig
  real(double_precision),intent(in) :: jcen(3)
  real(double_precision),intent(in) :: h
  real(double_precision),intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision),intent(in) :: xh(3,nbod) !< [in] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision),intent(in) :: vh(3,nbod) !< [in] velocities (vx,vy,vz) with respect to the central body (AU/day)
  
  real(double_precision),intent(out) :: x(3,nbod)
  real(double_precision),intent(out) :: v(3,nbod)
  
  ! Local
  integer :: j
  real(double_precision) :: msum,mvsum(3),temp,mbin,mbin_1,mtot_1
  
  !------------------------------------------------------------------------------
  
  msum = 0.d0
  mvsum(1) = 0.d0
  mvsum(2) = 0.d0
  mvsum(3) = 0.d0
  mbin = m(1) + m(2)
  mbin_1 = 1.d0 / mbin
  
  x(1,2) = xh(1,2)
  x(2,2) = xh(2,2)
  x(3,2) = xh(3,2)
  temp = m(1) * mbin_1
  v(1,2) = temp * vh(1,2)
  v(2,2) = temp * vh(2,2)
  v(3,2) = temp * vh(3,2)
  
  do j = 3, nbod
     msum = msum + m(j)
     mvsum(1) = mvsum(1)  +  m(j) * vh(1,j)
     mvsum(2) = mvsum(2)  +  m(j) * vh(2,j)
     mvsum(3) = mvsum(3)  +  m(j) * vh(3,j)
  end do
  mtot_1 = 1.d0 / (msum + mbin)
  mvsum(1) = mtot_1 * (mvsum(1) + m(2)*vh(1,2))
  mvsum(2) = mtot_1 * (mvsum(2) + m(2)*vh(2,2))
  mvsum(3) = mtot_1 * (mvsum(3) + m(2)*vh(3,2))
  
  temp = m(2) * mbin_1
  do j = 3, nbod
     x(1,j) = xh(1,j)  -  temp * xh(1,2)
     x(2,j) = xh(2,j)  -  temp * xh(2,2)
     x(3,j) = xh(3,j)  -  temp * xh(3,2)
     v(1,j) = vh(1,j)  -  mvsum(1)
     v(2,j) = vh(2,j)  -  mvsum(2)
     v(3,j) = vh(3,j)  -  mvsum(3)
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mco_h2cb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCO_IDEN.FOR    (ErikSoft   2 November 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Fake subroutine that simulate a change of coordinate system. In fact, 
! it only duplicate the existing coordinates. So it's the identity. 
! It is used to test the program with a fake algorithm (I think)

!------------------------------------------------------------------------------

subroutine mco_iden (time,jcen,nbod,nbig,h,m,x_in,v_in,x_out,v_out,ngf,ngflag)
  implicit none

  
  ! Input/Output
  integer,intent(in) :: nbod
  integer,intent(in) :: nbig
  integer,intent(in) :: ngflag
  real(double_precision),intent(in) :: time
  real(double_precision),intent(in) :: jcen(3)
  real(double_precision),intent(in) :: h
  real(double_precision),intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision),intent(in) :: x_in(3,nbod)
  real(double_precision),intent(in) :: v_in(3,nbod)
  real(double_precision),intent(in) :: ngf(4,nbod)

  real(double_precision), intent(out) :: x_out(3,nbod)
  real(double_precision), intent(out) :: v_out(3,nbod)
  
  ! Local
  integer :: j
  
  !------------------------------------------------------------------------------
  
  do j = 1, nbod
    x_out(1,j) = x_in(1,j)
    x_out(2,j) = x_in(2,j)
    x_out(3,j) = x_in(3,j)
    v_out(1,j) = v_in(1,j)
    v_out(2,j) = v_in(2,j)
    v_out(3,j) = v_in(3,j)
  enddo
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mco_iden

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCE_SPIN.FOR    (ErikSoft  2 December 1999)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Calculates the spin rate (in rotations per day) for a fluid body given
! its mass, spin angular momentum and density. The routine assumes the
! body is a MacClaurin ellipsoid, whose axis ratio is defined by the
! quantity SS = SQRT(A^2/C^2 - 1), where A and C are the
! major and minor axes.

!------------------------------------------------------------------------------

subroutine mce_spin (g,mass,spin,rho,rote)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input/Output
  real(double_precision),intent(in) :: g
  real(double_precision),intent(in) :: mass !< [in] mass (in solar masses * K2)
  real(double_precision),intent(in) :: spin !< [in] spin angular momentum (solar masses AU^2/day)
  real(double_precision),intent(in) :: rho !< [in] physical density (g/cm^3)
  
  real(double_precision),intent(out) :: rote !< [out] spin rate (in rotations per day)
  
  ! Local
  integer :: k
  real(double_precision) :: ss,s2,f,df,z,dz,tmp0,tmp1,t23
  
  !------------------------------------------------------------------------------
  
  t23 = 2.d0 / 3.d0
  tmp1 = spin * spin / (2.d0 * PI * rho * g)    * ( 250.d0*PI*PI*rho*rho / (9.d0*mass**5) )**t23
  
  ! Calculate SS using Newton's method
  ss = 1.d0
  do k = 1, 20
     s2 = ss * ss
     tmp0 = (1.d0 + s2)**t23
     call m_sfunc (ss,z,dz)
     f = z * tmp0  -  tmp1
     df = tmp0 * ( dz  +  4.d0 * ss * z / (3.d0*(1.d0 + s2)) )
     ss = ss - f/df
  end do
  
  rote = sqrt(TWOPI * g * rho * z) / TWOPI
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_spin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      M_SFUNC.FOR     (ErikSoft  14 November 1998)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Calculates Z = [ (3 + S^2)arctan(S) - 3S ] / S^3 and its derivative DZ,
! for S > 0.

!------------------------------------------------------------------------------

subroutine m_sfunc (s,z,dz)
  

  implicit none

  
  ! Input/Output
  real(double_precision),intent(in) :: s
  real(double_precision),intent(out) :: z
  real(double_precision),intent(out) :: dz
  
  ! Local
  real(double_precision) :: s2,s4,s6,s8,a
  
  !------------------------------------------------------------------------------
  
  s2 = s * s
  
  if (s.gt.1.d-2) then
     a  = atan(s)
     z  = ((3.d0 + s2)*a - 3.d0*s) / (s * s2)
     dz = (2.d0*s*a - 3.d0 + (3.d0+s2)/(1.d0+s2)) / (s * s2) - 3.d0 * z / s
  else
     s4 = s2 * s2
     s6 = s2 * s4
     s8 = s4 * s4
     z  = - .1616161616161616d0*s8   + .1904761904761905d0*s6   - .2285714285714286d0*s4   + .2666666666666667d0*s2
     dz = s * (- 1.292929292929293d0*s6 + 1.142857142857143d0*s4 - 0.914285714285714d0*s2 + 0.533333333333333d0)
  end if
  
  !------------------------------------------------------------------------------
  
  return
end subroutine m_sfunc

end module system_properties
