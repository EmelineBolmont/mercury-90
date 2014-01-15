!******************************************************************************
! MODULE: dynamic
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that compute various dynamic behaviour such as 
!! collision, close encounter, close approach to central body
!
!******************************************************************************

module dynamic

  use types_numeriques
  use utilities
  use mercury_globals
  
  implicit none
  
  contains
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !> @author 
   !> John E. Chambers
   !
   !> @date 4 March 2001
   !
   ! DESCRIPTION: 
   !> @brief Checks all objects with index I >= I0, to see if they have had a collision
   !! with the central body in a time interval H, when given the initial and 
   !! final coordinates and velocities. The routine uses cubic interpolation
   !! to estimate the minimum separations.
   !
   !> @attention All coordinates & velocities must be with respect to the central body!
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine mce_cent (time,h,rcen,jcen,start_index,nbod,nbig,m,x0,v0,x1,v1,nhit,jhit,thit,dhit,ngf,ngflag)
  
  use physical_constant
  use mercury_constant

  implicit none
  
  ! Input
  integer, intent(in) :: start_index
  integer, intent(in) :: nbod !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
  integer, intent(in) :: nbig !< [in] current number of big bodies (ones that perturb everything else)
  integer, intent(in) :: ngflag !< [in] do any bodies experience non-grav. forces?
!!\n                            ( 0 = no non-grav forces)
!!\n                              1 = cometary jets only
!!\n                              2 = radiation pressure/P-R drag only
!!\n                              3 = both
  real(double_precision), intent(in) :: time !< [in] current epoch (days)
  real(double_precision), intent(in) :: h !< [in] current integration timestep (days)
  real(double_precision), intent(in) :: rcen !< [in] radius of central body (AU)
  real(double_precision), intent(in) :: jcen(3) !< [in] J2,J4,J6 for central body (units of RCEN^i for Ji)
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x0(3,nbod)
  real(double_precision), intent(in) :: v0(3,nbod)
  real(double_precision), intent(in) :: x1(3,nbod)
  real(double_precision), intent(in) :: v1(3,nbod)
  real(double_precision), intent(in) :: ngf(4,nbod) !< [in] non gravitational forces parameters
  !! \n(1-3) cometary non-gravitational (jet) force parameters
  !! \n(4)  beta parameter for radiation pressure and P-R drag
  
  ! Output
  integer, intent(out) :: nhit
  integer, intent(out) :: jhit(CMAX)
  real(double_precision), intent(out) :: thit(CMAX)
  real(double_precision), intent(out) :: dhit(CMAX)
  
  ! Local
  integer :: i0 !< starting index for checking close encounters, mainly equal to start_index
  integer :: j
  real(double_precision) :: rcen2,a,q,u0,uhit,m0,mhit,mm,r0,mcen
  real(double_precision) :: hx,hy,hz,h2,p,rr0,rr1,rv0,rv1,temp,e,v2
  real(double_precision) :: xu0(3,nb_bodies_initial),xu1(3,nb_bodies_initial),vu0(3,nb_bodies_initial),vu1(3,nb_bodies_initial)
  
  !------------------------------------------------------------------------------
  
  if (start_index.le.0) then
    i0 = 2
  else
    i0 = start_index
  endif
  nhit = 0
  rcen2 = rcen * rcen
  mcen = m(1)
  
  ! If using close-binary code, convert to coords with respect to the binary
  !      if (algor.eq.11) then
  !        mcen = m(1) + m(2)
  !        call mco_h2ub (temp,jcen,nbod,nbig,h,m,x0,v0,xu0,vu0,ngf,ngflag)
  !        call mco_h2ub (temp,jcen,nbod,nbig,h,m,x1,v1,xu1,vu1,ngf,ngflag)
  !      end if
  
  ! Check for collisions with the central body
  do j = i0, nbod
     if (algor.eq.11) then
        rr0 = xu0(1,j)*xu0(1,j) + xu0(2,j)*xu0(2,j) +xu0(3,j)*xu0(3,j)
        rr1 = xu1(1,j)*xu1(1,j) + xu1(2,j)*xu1(2,j) +xu1(3,j)*xu1(3,j)
        rv0 = vu0(1,j)*xu0(1,j) + vu0(2,j)*xu0(2,j) +vu0(3,j)*xu0(3,j)
        rv1 = vu1(1,j)*xu1(1,j) + vu1(2,j)*xu1(2,j) +vu1(3,j)*xu1(3,j)
     else
        rr0 = x0(1,j)*x0(1,j) + x0(2,j)*x0(2,j) + x0(3,j)*x0(3,j)
        rr1 = x1(1,j)*x1(1,j) + x1(2,j)*x1(2,j) + x1(3,j)*x1(3,j)
        rv0 = v0(1,j)*x0(1,j) + v0(2,j)*x0(2,j) + v0(3,j)*x0(3,j)
        rv1 = v1(1,j)*x1(1,j) + v1(2,j)*x1(2,j) + v1(3,j)*x1(3,j)
     end if
     
     ! If inside the central body, or passing through pericentre, use 2-body approx.
     if ((rv0*h.le.0.and.rv1*h.ge.0).or.min(rr0,rr1).le.rcen2) then
        if (algor.eq.11) then
           hx = xu0(2,j) * vu0(3,j)  -  xu0(3,j) * vu0(2,j)
           hy = xu0(3,j) * vu0(1,j)  -  xu0(1,j) * vu0(3,j)
           hz = xu0(1,j) * vu0(2,j)  -  xu0(2,j) * vu0(1,j)
           v2 = vu0(1,j)*vu0(1,j) +vu0(2,j)*vu0(2,j) +vu0(3,j)*vu0(3,j)
        else
           hx = x0(2,j) * v0(3,j)  -  x0(3,j) * v0(2,j)
           hy = x0(3,j) * v0(1,j)  -  x0(1,j) * v0(3,j)
           hz = x0(1,j) * v0(2,j)  -  x0(2,j) * v0(1,j)
           v2 = v0(1,j)*v0(1,j) + v0(2,j)*v0(2,j) + v0(3,j)*v0(3,j)
        end if
        h2 = hx*hx + hy*hy + hz*hz
        p = h2 / (mcen + m(j))
        r0 = sqrt(rr0)
        temp = 1.d0 + p*(v2/(mcen + m(j)) - 2.d0/r0)
        e = sqrt( max(temp,0.d0) )
        q = p / (1.d0 + e)
        
        ! If the object hit the central body
        if (q.le.rcen) then
           nhit = nhit + 1
           jhit(nhit) = j
           dhit(nhit) = rcen
           
           ! Time of impact relative to the end of the timestep
           if (e.lt.1) then
              a = q / (1.d0 - e)
              uhit = sign (acos((1.d0 - rcen/a)/e), -h)
              u0   = sign (acos((1.d0 - r0/a  )/e), rv0)
              mhit = mod (uhit - e*sin(uhit) + PI, TWOPI) - PI
              m0   = mod (u0   - e*sin(u0)   + PI, TWOPI) - PI
           else
              a = q / (e - 1.d0)
              uhit = sign (arcosh((1.d0 - rcen/a)/e), -h)
              u0   = sign (arcosh((1.d0 - r0/a  )/e), rv0)
              mhit = mod (uhit - e*sinh(uhit) + PI, TWOPI) - PI
              m0   = mod (u0   - e*sinh(u0)   + PI, TWOPI) - PI
           end if
           mm = sqrt((mcen + m(j)) / (a*a*a))
           thit(nhit) = (mhit - m0) / mm + time
        end if
     end if
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_cent

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCE_COLL.FOR    (ErikSoft   2 October 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Resolves a collision between two objects, using the collision model chosen
! by the user. Also writes a message to the information file, and updates the
! value of ELOST, the change in energy due to collisions and ejections.

! N.B. All coordinates and velocities must be with respect to central body.
! ===

!------------------------------------------------------------------------------

subroutine mce_coll (time,elost,jcen,planet_id_1,planet_id_2,nbod,nbig,m,xh,vh,s,rphys,stat,id,outfile)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input/Output
  integer, intent(in) :: planet_id_1
  integer, intent(in) :: planet_id_2
  integer, intent(in) :: nbod !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
  integer, intent(in) :: nbig !< [in] current number of big bodies (ones that perturb everything else)
  integer, intent(inout) :: stat(nbod) !< [in,out] status (0 => alive, <>0 => to be removed)
  real(double_precision), intent(in) :: time !< [in] current epoch (days)
  real(double_precision), intent(in) :: jcen(3) !< [in] J2,J4,J6 for central body (units of RCEN^i for Ji)
  real(double_precision), intent(inout) :: m(nbod) !< [in,out] mass (in solar masses * K2)
  real(double_precision), intent(inout) :: xh(3,nbod) !< [in,out] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision), intent(inout) :: vh(3,nbod) !< [in,out] velocities (vx,vy,vz) with respect to the central body (AU/day)
  real(double_precision), intent(inout) :: s(3,nbod) !< [in,out] spin angular momentum (solar masses AU^2/day)
  real(double_precision), intent(in) :: rphys(nbod)
  character(len=80), intent(in) :: outfile
  character(len=8), intent(in) :: id(nbod) !< [in] name of the object (8 characters)
  
  real(double_precision), intent(inout) :: elost
  
  ! Local
  integer :: i
  integer :: j
  integer :: year,month,itmp, error
  real(double_precision) :: t1
  character(len=38) :: flost,fcol
  character(len=6) :: tstring
  
  !------------------------------------------------------------------------------
  
  ! If two bodies collided, check that the less massive one is removed
  ! (unless the more massive one is a Small body)
  if (planet_id_1.ne.0) then
     if (m(planet_id_2).gt.m(planet_id_1).and.planet_id_2.le.nbig) then
        i = planet_id_2
        j = planet_id_1
     else
        i = planet_id_1
        j = planet_id_2
     end if
  end if
  
  ! Write message to info file (I=0 implies collision with the central body)
  open (23, file=outfile, status='old', position='append', iostat=error)
  if (error /= 0) then
     write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile)
     stop
  end if
  
  if (opt(3).eq.1) then
     call mio_jd2y (time,year,month,t1)
     if (i.eq.0) then
        flost = '(1x,a8,a,i10,1x,i2,1x,f8.5)'
        write (23,flost) id(j),mem(67)(1:lmem(67)),year,month,t1
     else
        fcol  = '(1x,a8,a,a8,a,i10,1x,i2,1x,f4.1)'
        write (23,fcol) id(i),mem(69)(1:lmem(69)),id(j),mem(71)(1:lmem(71)),year,month,t1
     end if
  else
     if (opt(3).eq.3) then
        t1 = (time - tstart) / 365.25d0
        tstring = mem(2)
        flost = '(1x,a8,a,f18.7,a)'
        fcol  = '(1x,a8,a,a8,a,1x,f14.3,a)'
     else
        if (opt(3).eq.0) t1 = time
        if (opt(3).eq.2) t1 = time - tstart
        tstring = mem(1)
        flost = '(1x,a8,a,f18.5,a)'
        fcol  = '(1x,a8,a,a8,a,1x,f14.1,a)'
     end if
     if (i.eq.0.or.i.eq.1) then
        write (23,flost) id(j),mem(67)(1:lmem(67)),t1,tstring
     else
        write (23,fcol) id(i),mem(69)(1:lmem(69)),id(j),mem(71)(1:lmem(71)),t1,tstring
     end if
  end if
  close (23)
  
  ! Do the collision (inelastic merger)
  call mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost)
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_coll



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MCE_MERG.FOR    (ErikSoft   2 October 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c
! Author: John E. Chambers

! Merges objects I and J inelastically to produce a single new body by 
! conserving mass and linear momentum.
!   If J <= NBIG, then J is a Big body
!   If J >  NBIG, then J is a Small body
!   If I = 0, then I is the central body

! N.B. All coordinates and velocities must be with respect to central body.
! ===

!------------------------------------------------------------------------------

subroutine mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost)
  
  use physical_constant
  use mercury_constant
  use system_properties

  implicit none

  
  ! Input
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: nbod !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
  integer, intent(in) :: nbig !< [in] current number of big bodies (ones that perturb everything else)
  real(double_precision), intent(in) :: jcen(3) !< [in] J2,J4,J6 for central body (units of RCEN^i for Ji)
  
  ! Output
  integer, intent(inout) :: stat(nbod) !< [in,out] status (0 => alive, <>0 => to be removed)
  real(double_precision), intent(inout) :: m(nbod) !< [in,out] mass (in solar masses * K2)
  real(double_precision), intent(inout) :: xh(3,nbod) !< [in,out] coordinates (x,y,z) with respect to the central body (AU)
  real(double_precision), intent(inout) :: vh(3,nbod) !< [in,out] velocities (vx,vy,vz) with respect to the central body (AU/day)
  real(double_precision), intent(inout) :: s(3,nbod) !< [in,out] spin angular momentum (solar masses AU^2/day)
  real(double_precision), intent(inout) :: elost
  
  ! Local
  integer :: k
  real(double_precision) :: tmp1, tmp2, dx, dy, dz, du, dv, dw, msum, mredu, msum_1
  real(double_precision) :: e0, e1, l2
  
  !------------------------------------------------------------------------------
  
  ! If a body hits the central body
  if (i.le.1) then
     call mxx_en (jcen,nbod,nbig,m,xh,vh,s,e0,l2)
     
     ! If a body hit the central body...
     msum   = m(1) + m(j)
     msum_1 = 1.d0 / msum
     mredu  = m(1) * m(j) * msum_1
     dx = xh(1,j)
     dy = xh(2,j)
     dz = xh(3,j)
     du = vh(1,j)
     dv = vh(2,j)
     dw = vh(3,j)
     
     ! Calculate new spin angular momentum of the central body
     s(1,1) = s(1,1)  +  s(1,j)  +  mredu * (dy * dw  -  dz * dv)
     s(2,1) = s(2,1)  +  s(2,j)  +  mredu * (dz * du  -  dx * dw)
     s(3,1) = s(3,1)  +  s(3,j)  +  mredu * (dx * dv  -  dy * du)
     
     ! Calculate shift in barycentric coords and velocities of central body
     tmp2 = m(j) * msum_1
     xh(1,1) = tmp2 * xh(1,j)
     xh(2,1) = tmp2 * xh(2,j)
     xh(3,1) = tmp2 * xh(3,j)
     vh(1,1) = tmp2 * vh(1,j)
     vh(2,1) = tmp2 * vh(2,j)
     vh(3,1) = tmp2 * vh(3,j)
     m(1) = msum
     m(j) = 0.d0
     s(1,j) = 0.d0
     s(2,j) = 0.d0
     s(3,j) = 0.d0
     
     ! Shift the heliocentric coordinates and velocities of all bodies
     do k = 2, nbod
        xh(1,k) = xh(1,k) - xh(1,1)
        xh(2,k) = xh(2,k) - xh(2,1)
        xh(3,k) = xh(3,k) - xh(3,1)
        vh(1,k) = vh(1,k) - vh(1,1)
        vh(2,k) = vh(2,k) - vh(2,1)
        vh(3,k) = vh(3,k) - vh(3,1)
     end do
     
     ! Calculate energy loss due to the collision
     call mxx_en (jcen,nbod,nbig,m,xh,vh,s,e1,l2)
     elost = elost + (e0 - e1)
  else
     
     ! If two bodies collided...
     msum   = m(i) + m(j)
     msum_1 = 1.d0 / msum
     mredu  = m(i) * m(j) * msum_1
     dx = xh(1,i) - xh(1,j)
     dy = xh(2,i) - xh(2,j)
     dz = xh(3,i) - xh(3,j)
     du = vh(1,i) - vh(1,j)
     dv = vh(2,i) - vh(2,j)
     dw = vh(3,i) - vh(3,j)
     
     ! Calculate energy loss due to the collision
     elost = elost  +  .5d0 * mredu * (du*du + dv*dv + dw*dw)    -  m(i) * m(j) / sqrt(dx*dx + dy*dy + dz*dz)
     
     ! Calculate spin angular momentum of the new body
     s(1,i) = s(1,i)  +  s(1,j)  +  mredu * (dy * dw  -  dz * dv)
     s(2,i) = s(2,i)  +  s(2,j)  +  mredu * (dz * du  -  dx * dw)
     s(3,i) = s(3,i)  +  s(3,j)  +  mredu * (dx * dv  -  dy * du)
     
     ! Calculate new coords and velocities by conserving centre of mass & momentum
     tmp1 = m(i) * msum_1
     tmp2 = m(j) * msum_1
     xh(1,i) = xh(1,i) * tmp1  +  xh(1,j) * tmp2
     xh(2,i) = xh(2,i) * tmp1  +  xh(2,j) * tmp2
     xh(3,i) = xh(3,i) * tmp1  +  xh(3,j) * tmp2
     vh(1,i) = vh(1,i) * tmp1  +  vh(1,j) * tmp2
     vh(2,i) = vh(2,i) * tmp1  +  vh(2,j) * tmp2
     vh(3,i) = vh(3,i) * tmp1  +  vh(3,j) * tmp2
     m(i) = msum
  end if
  
  ! Flag the lost body for removal, and move it away from the new body
  stat(j) = -2
  xh(1,j) = -xh(1,j)
  xh(2,j) = -xh(2,j)
  xh(3,j) = -xh(3,j)
  vh(1,j) = -vh(1,j)
  vh(2,j) = -vh(2,j)
  vh(3,j) = -vh(3,j)
  m(j)   = 0.d0
  s(1,j) = 0.d0
  s(2,j) = 0.d0
  s(3,j) = 0.d0
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_merg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MXX_ELIM.FOR    (ErikSoft   13 February 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Removes any objects with STAT < 0 (i.e. those that have been flagged for 
! removal) and reindexes all the appropriate arrays for the remaining objects.

!------------------------------------------------------------------------------

subroutine mxx_elim (nbod,nbig,m,x,v,s,rho,rceh,rcrit,ngf,stat,id,outfile,nelim)
  
  use physical_constant
  use mercury_constant

  implicit none
  
  ! Input
  real(double_precision), intent(in) :: rcrit(nbod)
  character(len=80), intent(in) :: outfile
  
  ! Input/Output
  integer, intent(inout) :: nbod !< [in,out] current number of bodies (INCLUDING the central object)
  integer, intent(inout) :: nbig !< [in,out] current number of big bodies (ones that perturb everything else)
  integer, intent(inout) :: nelim
  integer, intent(inout) :: stat(nbod) !< [in,out] status (0 => alive, <>0 => to be removed)
  real(double_precision), intent(inout) :: m(nbod) !< [in,out] mass (in solar masses * K2)
  real(double_precision), intent(inout) :: x(3,nbod)
  real(double_precision), intent(inout) :: v(3,nbod)
  real(double_precision), intent(inout) :: s(3,nbod) !< [in,out] spin angular momentum (solar masses AU^2/day)
  real(double_precision), intent(inout) :: rho(nbod) !< [in,out] physical density (g/cm^3)
  real(double_precision), intent(inout) :: rceh(nbod) !< [in,out] close-encounter limit (Hill radii)
  real(double_precision), intent(inout) :: ngf(4,nbod) !< [in,out] non gravitational forces parameters
  !! \n(1-3) cometary non-gravitational (jet) force parameters
  !! \n(4)  beta parameter for radiation pressure and P-R drag
  
  character(len=8), intent(inout) :: id(nbod) !< [in,out] name of the object (8 characters)
  
  ! Local
  integer :: j, k, l, nbigelim, elim(nb_bodies_initial+1)
  integer :: error
  
  !------------------------------------------------------------------------------
  
  ! Find out how many objects are to be removed
  nelim = 0
  nbigelim = 0
  do j = 2, nbod
     if (stat(j).lt.0) then
        nelim = nelim + 1
        elim(nelim) = j
        if (j.le.nbig) nbigelim = nbigelim + 1
     end if
  end do
  elim(nelim+1) = nbod + 1
  
  ! Eliminate unwanted objects
  do k = 1, nelim
     do j = elim(k) - k + 1, elim(k+1) - k - 1
        l = j + k
        x(1,j) = x(1,l)
        x(2,j) = x(2,l)
        x(3,j) = x(3,l)
        v(1,j) = v(1,l)
        v(2,j) = v(2,l)
        v(3,j) = v(3,l)
        m(j)   = m(l)
        s(1,j) = s(1,l)
        s(2,j) = s(2,l)
        s(3,j) = s(3,l)
        rho(j) = rho(l)
        rceh(j) = rceh(l)
        stat(j) = stat(l)
        id(j) = id(l)
        ngf(1,j) = ngf(1,l)
        ngf(2,j) = ngf(2,l)
        ngf(3,j) = ngf(3,l)
        ngf(4,j) = ngf(4,l)
     end do
  end do
  
  ! Update total number of bodies and number of Big bodies
  nbod = nbod - nelim
  nbig = nbig - nbigelim
  
  ! If no massive bodies remain, stop the integration
  if (nbig.lt.1) then
     open (23,file=outfile,status='old',position='append',iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile)
        stop
     end if
     write (23,'(2a)') mem(81)(1:lmem(81)),mem(124)(1:lmem(124))
     close (23)
     stop
  end if
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mxx_elim

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    MCE_SNIF.FOR    (ErikSoft   3 October 2000)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Given initial and final coordinates and velocities X and V, and a timestep 
! H, the routine estimates which objects were involved in a close encounter
! during the timestep. The routine examines all objects with indices I >= I0.

! Returns an array CE, which for each object is: 
!                           0 if it will undergo no encounters
!                           2 if it will pass within RCRIT of a Big body

! Also returns arrays ICE and JCE, containing the indices of each pair of
! objects estimated to have undergone an encounter.

! N.B. All coordinates must be with respect to the central body!!!
! ===

!------------------------------------------------------------------------------

subroutine mce_snif (h,start_index,nbod,nbig,x0,v0,x1,v1,rcrit,ce,nce,ice,jce)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input
  integer, intent(in) :: start_index
  integer, intent(in) :: nbod !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
  integer, intent(in) :: nbig !< [in] current number of big bodies (ones that perturb everything else)
  real(double_precision), intent(in) :: x0(3,nbod)
  real(double_precision), intent(in) :: v0(3,nbod)
  real(double_precision), intent(in) :: x1(3,nbod)
  real(double_precision), intent(in) :: v1(3,nbod)
  real(double_precision), intent(in) :: h !< [in] current integration timestep (days)
  real(double_precision), intent(in) :: rcrit(nbod)
  
  ! Output
  integer, intent(out) :: ce(nbod) !< [out] close encounter status
  integer, intent(out) :: nce
  integer, intent(out) :: ice(CMAX)
  integer, intent(out) :: jce(CMAX)
  
  ! Local
  integer :: i0 !< starting index for checking close encounters, mainly equal to start_index
  integer :: i,j
  real(double_precision) :: d0,d1,d0t,d1t,d2min,temp,tmin,rc,rc2
  real(double_precision) :: dx0,dy0,dz0,du0,dv0,dw0,dx1,dy1,dz1,du1,dv1,dw1
  real(double_precision) :: xmin(nb_bodies_initial),xmax(nb_bodies_initial),ymin(nb_bodies_initial),ymax(nb_bodies_initial)
  
  !------------------------------------------------------------------------------
  
  if (start_index.le.0) then
    i0 = 2
  else
    i0 = start_index
  endif
  
  nce = 0
  do j = 2, nbod
     ce(j) = 0
  end do
  
  ! Calculate maximum and minimum values of x and y coordinates
  call mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
  
  ! Adjust values for the Big bodies by symplectic close-encounter distance
  do j = i0, nbig
     xmin(j) = xmin(j) - rcrit(j)
     xmax(j) = xmax(j) + rcrit(j)
     ymin(j) = ymin(j) - rcrit(j)
     ymax(j) = ymax(j) + rcrit(j)
  end do
  
  ! Identify pairs whose X-Y boxes overlap, and calculate minimum separation
  do i = i0, nbig
     do j = i + 1, nbod
        if (xmax(i).ge.xmin(j).and.xmax(j).ge.xmin(i).and.ymax(i).ge.ymin(j).and.ymax(j).ge.ymin(i)) then
           
           ! Determine the maximum separation that would qualify as an encounter
           rc = max(rcrit(i), rcrit(j))
           rc2 = rc * rc
           
           ! Calculate initial and final separations
           dx0 = x0(1,i) - x0(1,j)
           dy0 = x0(2,i) - x0(2,j)
           dz0 = x0(3,i) - x0(3,j)
           dx1 = x1(1,i) - x1(1,j)
           dy1 = x1(2,i) - x1(2,j)
           dz1 = x1(3,i) - x1(3,j)
           d0 = dx0*dx0 + dy0*dy0 + dz0*dz0
           d1 = dx1*dx1 + dy1*dy1 + dz1*dz1
           
           ! Check for a possible minimum in between
           du0 = v0(1,i) - v0(1,j)
           dv0 = v0(2,i) - v0(2,j)
           dw0 = v0(3,i) - v0(3,j)
           du1 = v1(1,i) - v1(1,j)
           dv1 = v1(2,i) - v1(2,j)
           dw1 = v1(3,i) - v1(3,j)
           d0t = (dx0*du0 + dy0*dv0 + dz0*dw0) * 2.d0
           d1t = (dx1*du1 + dy1*dv1 + dz1*dw1) * 2.d0
           
           ! If separation derivative changes sign, find the minimum separation
           d2min = HUGE
           if (d0t*h.le.0.and.d1t*h.ge.0) call mce_min (d0,d1,d0t,d1t,h,d2min,tmin)
           
           ! If minimum separation is small enough, flag this as a possible encounter
           temp = min (d0,d1,d2min)
           if (temp.le.rc2) then
              ce(i) = 2
              ce(j) = 2
              if (nce.lt.CMAX) then
                nce = nce + 1
                ice(nce) = i
                jce(nce) = j
              else
                write(*,*) 'Error: The number of close encounters exceed the limit (CMAX=',CMAX,').'
                write(*,*) '       All remaining close encounters are skipped'
                return
              end if
           end if
        end if
     end do
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_snif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    MCE_STAT.FOR    (ErikSoft   1 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Returns details of all close-encounter minima involving at least one Big
! body during a timestep. The routine estimates minima using the initial
! and final coordinates X(0),X(1) and velocities V(0),V(1) of the step, and
! the stepsize H.

!  ICLO, JCLO contain the indices of the objects
!  DCLO is their minimum separation
!  TCLO is the time of closest approach relative to current time

! The routine also checks for collisions/near misses given the physical radii 
! RPHYS, and returns the time THIT of the collision/near miss closest to the
! start of the timestep, and the identities IHIT and JHIT of the objects
! involved.

!  NHIT = +1 implies a collision
!         -1    "    a near miss

! N.B. All coordinates & velocities must be with respect to the central body!
! ===
!------------------------------------------------------------------------------

subroutine mce_stat (time,h,rcen,nbod,nbig,m,x0,v0,x1,v1,rce,rphys,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit,chit,dhit,&
     thit,thit1,nowflag,stat,outfile)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input/Output
  integer, intent(in) :: nbod !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
  integer, intent(in) :: nbig !< [in] current number of big bodies (ones that perturb everything else)
  integer, intent(in) :: stat(nbod) !< [in] status (0 => alive, <>0 => to be removed)
  integer, intent(out) :: nowflag
  
  integer, intent(out) :: nclo
  integer, intent(out) :: iclo(CMAX)
  integer, intent(out) :: jclo(CMAX)
  
  integer, intent(out) :: nhit
  integer, intent(out) :: ihit(CMAX)
  integer, intent(out) :: jhit(CMAX)
  integer, intent(out) :: chit(CMAX)
  
  real(double_precision), intent(in) :: time !< [in] current epoch (days)
  real(double_precision), intent(in) :: h !< [in] current integration timestep (days)
  real(double_precision), intent(in) :: rcen !< [in] radius of central body (AU)
  real(double_precision), intent(in) :: m(nbod) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: x0(3,nbod)
  real(double_precision), intent(in) :: v0(3,nbod)
  
  real(double_precision), intent(in) :: x1(3,nbod)
  real(double_precision), intent(in) :: v1(3,nbod)
  real(double_precision), intent(in) :: rce(nbod)
  real(double_precision), intent(in) :: rphys(nbod)
  real(double_precision), intent(out) :: dclo(CMAX)
  real(double_precision), intent(out) :: tclo(CMAX)
  real(double_precision), intent(out) :: thit(CMAX)
  real(double_precision), intent(out) :: dhit(CMAX)
  real(double_precision), intent(out) :: thit1
  real(double_precision), intent(out) :: ixvclo(6,CMAX)
  real(double_precision), intent(out) :: jxvclo(6,CMAX)
  character(len=80), intent(in) :: outfile
  
  ! Local
  integer :: i,j, error
  real(double_precision) :: d0,d1,d0t,d1t,hm1,tmp0,tmp1
  real(double_precision) :: dx0,dy0,dz0,du0,dv0,dw0,dx1,dy1,dz1,du1,dv1,dw1
  real(double_precision) :: xmin(nb_bodies_initial),xmax(nb_bodies_initial),ymin(nb_bodies_initial),ymax(nb_bodies_initial)
  real(double_precision) :: d2min,d2ce,d2near,d2hit,temp,tmin
  
  !------------------------------------------------------------------------------
  
  nhit = 0
  thit1 = sign(HUGE, h)
  hm1 = 1.d0 / h
  
  ! Calculate maximum and minimum values of x and y coords for each object
  call mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
  
  ! Adjust values by the maximum close-encounter radius plus a fudge factor
  do j = 2, nbod
     temp = rce(j) * 1.2d0
     xmin(j) = xmin(j)  -  temp
     xmax(j) = xmax(j)  +  temp
     ymin(j) = ymin(j)  -  temp
     ymax(j) = ymax(j)  +  temp
  end do
  
  ! Check for close encounters between each pair of objects
  do i = 2, nbig
     do j = i + 1, nbod
        if (xmax(i).ge.xmin(j).and.xmax(j).ge.xmin(i).and.ymax(i).ge.ymin(j).and.ymax(j).ge.ymin(i).and.&
             stat(i).ge.0.and.stat(j).ge.0) then
           
           ! If the X-Y boxes for this pair overlap, check circumstances more closely
           dx0 = x0(1,i) - x0(1,j)
           dy0 = x0(2,i) - x0(2,j)
           dz0 = x0(3,i) - x0(3,j)
           du0 = v0(1,i) - v0(1,j)
           dv0 = v0(2,i) - v0(2,j)
           dw0 = v0(3,i) - v0(3,j)
           d0t = (dx0*du0 + dy0*dv0 + dz0*dw0) * 2.d0
           
           dx1 = x1(1,i) - x1(1,j)
           dy1 = x1(2,i) - x1(2,j)
           dz1 = x1(3,i) - x1(3,j)
           du1 = v1(1,i) - v1(1,j)
           dv1 = v1(2,i) - v1(2,j)
           dw1 = v1(3,i) - v1(3,j)
           d1t = (dx1*du1 + dy1*dv1 + dz1*dw1) * 2.d0
           
           ! Estimate minimum separation during the time interval, using interpolation
           d0 = dx0*dx0 + dy0*dy0 + dz0*dz0
           d1 = dx1*dx1 + dy1*dy1 + dz1*dz1
           call mce_min (d0,d1,d0t,d1t,h,d2min,tmin)
           d2ce  = max (rce(i), rce(j))
           d2hit = rphys(i) + rphys(j)
           d2ce   = d2ce  * d2ce
           d2hit  = d2hit * d2hit
           d2near = d2hit * 4.d0
           !if (abs(time - 2452983.63D0) <= 10.D0 .and. (i == 3 .or. j == 3)) then
           !  print *, 'clo', time, i, j, d2min,d2ce,d0t,h,d1t,d2hit
           !end if
           
           ! If the minimum separation qualifies as an encounter or if a collision
           ! is in progress, store details
           if ((d2min.le.d2ce.and.d0t*h.le.0.and.d1t*h.ge.0).or.(d2min.le.d2hit)) then
              nclo = nclo + 1
              !if (i == 3 .or. j == 3) then
               ! print *, 'clo', time, nclo, i, j
              !end if
              if (nclo.gt.CMAX) then
                 open (23,file=outfile,status='old',position='append',iostat=error)
                 if (error /= 0) then
                    write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile)
                    stop
                 end if
                 write (23,'(/,2a,/,a)') mem(121)(1:lmem(121)),mem(132)(1:lmem(132)),mem(82)(1:lmem(82))
                 close (23)
              else
                 tclo(nclo) = tmin + time
                 dclo(nclo) = sqrt (max(0.d0,d2min))
                 iclo(nclo) = i
                 jclo(nclo) = j
                 
                 ! Make sure the more massive body is listed first
                 if (m(j).gt.m(i).and.j.le.nbig) then
                    iclo(nclo) = j
                    jclo(nclo) = i
                 end if
                 
                 ! Make linear interpolation to get coordinates at time of closest approach
                 tmp0 = 1.d0 + tmin*hm1
                 tmp1 = -tmin*hm1
                 ixvclo(1,nclo) = tmp0 * x0(1,i)  +  tmp1 * x1(1,i)
                 ixvclo(2,nclo) = tmp0 * x0(2,i)  +  tmp1 * x1(2,i)
                 ixvclo(3,nclo) = tmp0 * x0(3,i)  +  tmp1 * x1(3,i)
                 ixvclo(4,nclo) = tmp0 * v0(1,i)  +  tmp1 * v1(1,i)
                 ixvclo(5,nclo) = tmp0 * v0(2,i)  +  tmp1 * v1(2,i)
                 ixvclo(6,nclo) = tmp0 * v0(3,i)  +  tmp1 * v1(3,i)
                 jxvclo(1,nclo) = tmp0 * x0(1,j)  +  tmp1 * x1(1,j)
                 jxvclo(2,nclo) = tmp0 * x0(2,j)  +  tmp1 * x1(2,j)
                 jxvclo(3,nclo) = tmp0 * x0(3,j)  +  tmp1 * x1(3,j)
                 jxvclo(4,nclo) = tmp0 * v0(1,j)  +  tmp1 * v1(1,j)
                 jxvclo(5,nclo) = tmp0 * v0(2,j)  +  tmp1 * v1(2,j)
                 jxvclo(6,nclo) = tmp0 * v0(3,j)  +  tmp1 * v1(3,j)
              end if
           end if
           
           ! Check for a near miss or collision
           if (d2min.le.d2near) then
              nhit = nhit + 1
              ihit(nhit) = i
              jhit(nhit) = j
              thit(nhit) = tmin + time
              dhit(nhit) = sqrt(d2min)
              chit(nhit) = -1
              if (d2min.le.d2hit) chit(nhit) = 1
              
              ! Make sure the more massive body is listed first
              if (m(jhit(nhit)).gt.m(ihit(nhit)).and.j.le.nbig) then
                 ihit(nhit) = j
                 jhit(nhit) = i
              end if
              
              ! Is this the collision closest to the start of the time step?
              if ((tmin-thit1)*h.lt.0) then
                 thit1 = tmin
                 nowflag = 0
                 if (d1.le.d2hit) nowflag = 1
              end if
           end if
        end if
        
        ! Move on to the next pair of objects
     end do
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_stat

end module dynamic
