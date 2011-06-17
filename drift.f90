module drift

!*************************************************************
!** Modules that group all the subroutine linked to drift  ***
!** (don't know exactly what it is)
!** Version 1.0 - june 2011
!*************************************************************
  use types_numeriques

  implicit none
  
  private
  
  public :: drift_one ! Only drift_one can be used externaly.
  
  contains

!------------------------------------------------------------------

!*************************************************************************
!                        DRIFT_ONE.F
!*************************************************************************
! This subroutine does the danby-type drift for one particle, using 
! appropriate vbles and redoing a drift if the accuracy is too poor 
! (as flagged by the integer iflg).

!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 x,y,z         ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 vx,vy,vz      ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 dt            ==>  time step
!             Output:
!                 x,y,z         ==>  final position in jacobi coord 
!                                       (real scalars)
!                 vx,vy,vz      ==>  final position in jacobi coord 
!                                       (real scalars)
!                 iflg          ==>  integer :: (zero for successful step)

! Authors:  Hal Levison & Martin Duncan 
! Date:    2/10/93
! Last revision: 2/10/93


subroutine drift_one(mu,x,y,z,vx,vy,vz,dt,iflg)

  use mercury_constant

  implicit none


  !...  Inputs Only: 
  real(double_precision) :: mu,dt

  !...  Inputs and Outputs:
  real(double_precision) :: x,y,z
  real(double_precision) :: vx,vy,vz

  !...  Output
  integer :: iflg

  !...  Internals:
  integer :: i
  real(double_precision) :: dttmp

  !----
  !...  Executable code 

  call drift_dan(mu,x,y,z,vx,vy,vz,dt,iflg)

  if(iflg .ne. 0) then

     do i = 1,10
        dttmp = dt/10.d0
        call drift_dan(mu,x,y,z,vx,vy,vz,dttmp,iflg)
        if(iflg .ne. 0) return
     enddo

  endif

  return
end subroutine drift_one    ! drift_one

!*************************************************************************
!                        DRIFT_DAN.F
!*************************************************************************
! This subroutine does the Danby and decides which vbles to use

!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 x0,y0,z0         ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 vx0,vy0,vz0      ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 dt0            ==>  time step
!             Output:
!                 x0,y0,z0         ==>  final position in jacobi coord 
!                                       (real scalars)
!                 vx0,vy0,vz0      ==>  final position in jacobi coord 
!                                       (real scalars)
!                 iflg             ==>  integer :: flag (zero if satisfactory)
!					      (non-zero if nonconvergence)

! Authors:  Hal Levison & Martin Duncan  
! Date:    2/10/93
! Last revision: April 6/93 - MD adds dt and keeps dt0 unchanged

subroutine drift_dan(mu,x0,y0,z0,vx0,vy0,vz0,dt0,iflg)

  use mercury_constant
  use physical_constant

  implicit none


  !...  Inputs Only: 
  real(double_precision) :: mu,dt0

  !...  Inputs and Outputs:
  real(double_precision) :: x0,y0,z0
  real(double_precision) :: vx0,vy0,vz0

  !...  Output
  integer :: iflg

  !...  Internals:
  real(double_precision) :: x,y,z,vx,vy,vz,dt
  real(double_precision) :: f,g,fdot,c1,c2
  real(double_precision) :: c3,gdot
  real(double_precision) :: u,alpha,fp,r0,v0s
  real(double_precision) :: a,asq,en
  real(double_precision) :: dm,ec,es,esq,xkep
  real(double_precision) :: fchk,s,c

  !----
  !...  Executable code 

  !...  Set dt = dt0 to be sure timestep is not altered while solving
  !...  for new coords.
  dt = dt0
  iflg = 0
  r0 = sqrt(x0*x0 + y0*y0 + z0*z0)
  v0s = vx0*vx0 + vy0*vy0 + vz0*vz0
  u = x0*vx0 + y0*vy0 + z0*vz0
  alpha = 2.0*mu/r0 - v0s

  if (alpha.gt.0.d0) then
     a = mu/alpha
     asq = a*a
     en = sqrt(mu/(a*asq))
     ec = 1.0d0 - r0/a
     es = u/(en*asq)
     esq = ec*ec + es*es
     dm = dt*en - int(dt*en/TWOPI)*TWOPI
     dt = dm/en
     if((dm*dm .gt. 0.16d0) .or. (esq.gt.0.36d0)) goto 100

     if(esq*dm*dm .lt. 0.0016) then

        call drift_kepmd(dm,es,ec,xkep,s,c)
        fchk = (xkep - ec*s +es*(1.-c) - dm)

        if(fchk*fchk .gt. DANBYB) then
           iflg = 1
           return
        endif

        fp = 1. - ec*c + es*s
        f = (a/r0) * (c-1.) + 1.
        g = dt + (s-xkep)/en
        fdot = - (a/(r0*fp))*en*s
        gdot = (c-1.)/fp + 1.

        x = x0*f + vx0*g
        y = y0*f + vy0*g
        z = z0*f + vz0*g
        vx = x0*fdot + vx0*gdot
        vy = y0*fdot + vy0*gdot
        vz = z0*fdot + vz0*gdot

        x0 = x
        y0 = y
        z0 = z
        vx0 = vx
        vy0 = vy
        vz0 = vz

        iflg = 0
        return

     endif

  endif

100 call drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

  if(iflg .eq.0) then
     f = 1.0 - (mu/r0)*c2
     g = dt - mu*c3
     fdot = -(mu/(fp*r0))*c1
     gdot = 1. - (mu/fp)*c2

     x = x0*f + vx0*g
     y = y0*f + vy0*g
     z = z0*f + vz0*g
     vx = x0*fdot + vx0*gdot
     vy = y0*fdot + vy0*gdot
     vz = z0*fdot + vz0*gdot

     x0 = x
     y0 = y
     z0 = z
     vx0 = vx
     vy0 = vy
     vz0 = vz
  endif

  return
end subroutine drift_dan   ! drift_dan

!********************************************************************#
!                  DRIFT_KEPMD
!********************************************************************#
!  Subroutine for solving kepler's equation in difference form for an
!  ellipse, given SMALL dm and SMALL eccentricity.  See DRIFT_DAN.F
!  for the criteria.
!  WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL
!  EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY
!  CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO.

!	Input:
!	    dm		==> increment in mean anomaly M (real*8 scalar)
!	    es,ec       ==> ecc. times sin and cos of E_0 (real*8 scalars)

!       Output:
!            x          ==> solution to Kepler's difference eqn (real*8 scalar)
!            s,c        ==> sin and cosine of x (real*8 scalars)


subroutine drift_kepmd(dm,es,ec,x,s,c)


  implicit none


  !...    Inputs
  real(double_precision) :: dm,es,ec

  !...	Outputs
  real(double_precision) :: x,s,c

  !...    Internals
  real(double_precision), parameter :: A0 = 39916800.d0
  real(double_precision), parameter :: A1 = 6652800.d0
  real(double_precision), parameter :: A2 = 332640.d0
  real(double_precision), parameter :: A3 = 7920.d0
  real(double_precision), parameter :: A4 = 110.d0

  real(double_precision) :: dx
  real(double_precision) :: fac1,fac2,q,y
  real(double_precision) :: f,fp,fpp,fppp


  !...    calc initial guess for root
  fac1 = 1.d0/(1.d0 - ec)
  q = fac1*dm
  fac2 = es*es*fac1 - ec/3.d0
  x = q*(1.d0 -0.5d0*fac1*q*(es -q*fac2))

  !...  excellent approx. to sin and cos of x for small x.
  y = x*x
  s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
  c = sqrt(1.d0 - s*s)

  !...    Compute better value for the root using quartic Newton method
  f = x - ec*s + es*(1.-c) - dm
  fp = 1. - ec*c + es*s
  fpp = ec*s + es*c
  fppp = ec*c - es*s
  dx = -f/fp
  dx = -f/(fp + 0.5*dx*fpp)
  dx = -f/(fp + 0.5*dx*fpp + 0.16666666666666666*dx*dx*fppp)
  x = x + dx

  !...  excellent approx. to sin and cos of x for small x.
  y = x*x
  s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
  c = sqrt(1.d0 - s*s)

  return
end subroutine drift_kepmd
!-----------------------------------------------------------------------------

!*************************************************************************
!                        DRIFT_KEPU.F
!*************************************************************************
! subroutine for solving kepler's equation using universal variables.

!             Input:
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 fp            ==>  f' from p170  
!                                       (real scalors)
!                 c1,c2,c3      ==>  c's from p171-172
!                                       (real scalors)
!                 iflg          ==>  =0 if converged; !=0 if not

! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 2/3/93

subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

  use mercury_constant

  implicit none


  !...  Inputs: 
  real(double_precision) :: dt,r0,mu,alpha,u

  !...  Outputs:
  real(double_precision) :: fp,c1,c2,c3
  integer :: iflg

  !...  Internals:
  real(double_precision) :: s,st,fo,fn

  !----
  !...  Executable code 

  call drift_kepu_guess(dt,r0,mu,alpha,u,s)

  st = s
  !..     store initial guess for possible use later in
  !..     laguerre's method, in case newton's method fails.

  call drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
  if(iflg.ne.0) then
     call drift_kepu_fchk(dt,r0,mu,alpha,u,st,fo)
     call drift_kepu_fchk(dt,r0,mu,alpha,u,s,fn)
     if(abs(fo).lt.abs(fn)) then
        s = st 
     endif
     call drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
  endif

  return
end subroutine drift_kepu    ! drift_kepu
!----------------------------------------------------------------------

!*************************************************************************
!                        DRIFT_KEPU_FCHK.F
!*************************************************************************
! Returns the value of the function f of which we are trying to find the root
! in universal variables.

!             Input:
!                 dt            ==>  time step (real scalar)
!                 r0            ==>  Distance between `Sun' and particle
!                                     (real scalar)
!                 mu            ==>  Reduced mass of system (real scalar)
!                 alpha         ==>  Twice the binding energy (real scalar)
!                 u             ==>  Vel. dot radial vector (real scalar)
!                 s             ==>  Approx. root of f 
!             Output:
!                 f             ==>  function value ( = 0 if O.K.) (integer)

! Author:  Martin Duncan  
! Date:    March 12/93
! Last revision: March 12/93

subroutine drift_kepu_fchk(dt,r0,mu,alpha,u,s,f)


  implicit none

  !...  Inputs: 
  real(double_precision) :: dt,r0,mu,alpha,u,s

  !...  Outputs:
  real(double_precision) :: f

  !...  Internals:
  real(double_precision) ::  x,c0,c1,c2,c3

  !----
  !...  Executable code 

  x=s*s*alpha
  call drift_kepu_stumpff(x,c0,c1,c2,c3)
  c1=c1*s
  c2 = c2*s*s
  c3 = c3*s*s*s
  f = r0*c1 + u*c2 + mu*c3 - dt

  return
end subroutine drift_kepu_fchk     !   drift_kepu_fchk
!-------------------------------------------------------------------

!*************************************************************************
!                        DRIFT_KEPU_GUESS.F
!*************************************************************************
! Initial guess for solving kepler's equation using universal variables.

!             Input:
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  initial guess for the value of 
!                                    universal variable

! Author:  Hal Levison & Martin Duncan 
! Date:    3/12/93
! Last revision: April 6/93
! Modified by JEC: 8/6/98

subroutine drift_kepu_guess(dt,r0,mu,alpha,u,s)

  use mercury_constant
  use utilities, only : mco_sine

  implicit none


  !...  Inputs: 
  real(double_precision) :: dt,r0,mu,alpha,u

  !...  Inputs and Outputs:
  real(double_precision) :: s

  !...  Internals:
  integer :: iflg
  real(double_precision) :: y,sy,cy,sigma,es
  real(double_precision) :: x,a
  real(double_precision) :: en,ec,e

  !----
  !...  Executable code 

  if (alpha.gt.0.0) then 
     !...       find initial guess for elliptic motion

     if( dt/r0 .le. 0.4)  then
        s = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0)
        return
     else
        a = mu/alpha
        en = sqrt(mu/(a*a*a))
        ec = 1.0 - r0/a
        es = u/(en*a*a)
        e = sqrt(ec*ec + es*es)
        y = en*dt - es
        
        call mco_sine (y,sy,cy)
        
        sigma = dsign(1.d0,(es*cy + ec*sy))
        x = y + sigma*.85*e
        s = x/sqrt(alpha)
     endif

  else
     !...       find initial guess for hyperbolic motion.
     call drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)
     if(iflg.ne.0) then
        s = dt/r0
     endif
  endif

  return
end subroutine drift_kepu_guess     !   drift_kepu_guess
!-------------------------------------------------------------------

!*************************************************************************
!                        DRIFT_KEPU_LAG.F
!*************************************************************************
! subroutine for solving kepler's equation in universal variables.
! using LAGUERRE'S METHOD

!             Input:
!                 s             ==>  inital value of universal variable
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  final value of universal variable
!                 fp            ==>  f' from p170  
!                                       (real scalors)
!                 c1,c2,c3      ==>  c's from p171-172
!                                       (real scalors)
!                 iflgn          ==>  =0 if converged; !=0 if not

! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 4/21/93

subroutine drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

  use mercury_constant

  implicit none


  !...  Inputs: 
  real(double_precision) :: s,dt,r0,mu,alpha,u

  !...  Outputs:
  real(double_precision) :: fp,c1,c2,c3
  integer :: iflg

  !...  Internals:
  integer :: nc,ncmax
  real(double_precision) :: ln
  real(double_precision) :: x,fpp,ds,c0,f
  real(double_precision) :: fdt

  integer, parameter :: NTMP=NLAG2+1


  !----
  !...  Executable code 

  !...    To get close approch needed to take lots of iterations if alpha<0
  if(alpha.lt.0.0) then
     ncmax = NLAG2
  else
     ncmax = NLAG2
  endif

  ln = 5.0
  !...    start laguere's method
  do nc =0,ncmax
     x = s*s*alpha
     call drift_kepu_stumpff(x,c0,c1,c2,c3)
     c1 = c1*s 
     c2 = c2*s*s 
     c3 = c3*s*s*s
     f = r0*c1 + u*c2 + mu*c3 - dt
     fp = r0*c0 + u*c1 + mu*c2
     fpp = (-40.0*alpha + mu)*c1 + u*c0
     ds = - ln*f/(fp + dsign(1.d0,fp)*sqrt(abs((ln - 1.0)*(ln - 1.0)*fp*fp - (ln - 1.0)*ln*f*fpp)))
     s = s + ds

     fdt = f/dt

     !..        quartic convergence
     if( fdt*fdt.lt.DANBYB*DANBYB) then 
        iflg = 0
        return
     endif
     !...      Laguerre's method succeeded
  enddo

  iflg = 2

  return

end subroutine drift_kepu_lag    !    drift_kepu_leg
!-----------------------------------------------------------------------

!*************************************************************************
!                        DRIFT_KEPU_NEW.F
!*************************************************************************
! subroutine for solving kepler's equation in universal variables.
! using NEWTON'S METHOD

!             Input:
!                 s             ==>  inital value of universal variable
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  final value of universal variable
!                 fp            ==>  f' from p170  
!                                       (real scalors)
!                 c1,c2,c3      ==>  c's from p171-172
!                                       (real scalors)
!                 iflgn          ==>  =0 if converged; !=0 if not

! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 4/21/93
! Modified by JEC: 31/3/98

subroutine drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflgn)

  use mercury_constant

  implicit none


  !...  Inputs: 
  real(double_precision) :: s,dt,r0,mu,alpha,u

  !...  Outputs:
  real(double_precision) :: fp,c1,c2,c3
  integer :: iflgn

  !...  Internals:
  integer :: nc
  real(double_precision) :: x,c0,ds,s2
  real(double_precision) :: f,fpp,fppp,fdt

  !----
  !...  Executable code 

  do nc=0,6
     s2 = s * s
     x = s2*alpha
     call drift_kepu_stumpff(x,c0,c1,c2,c3)
     c1 = c1*s 
     c2 = c2*s2 
     c3 = c3*s*s2
     f = r0*c1 + u*c2 + mu*c3 - dt
     fp = r0*c0 + u*c1 + mu*c2
     fpp = (mu - r0*alpha)*c1 + u*c0
     fppp = (mu - r0*alpha)*c0 - u*alpha*c1
     ds = - f/fp
     ds = - f/(fp + .5d0*ds*fpp)
     ds = -f/(fp + .5d0*ds*fpp + ds*ds*fppp*.1666666666666667d0)
     s = s + ds
     fdt = f/dt

     !..      quartic convergence
     if( fdt*fdt.lt.DANBYB*DANBYB) then 
        iflgn = 0
        return
     endif
     !...     newton's method succeeded

  enddo

  !..     newton's method failed
  iflgn = 1
  return

end subroutine drift_kepu_new  ! drift_kepu_new
!----------------------------------------------------------------------

!*************************************************************************
!                        DRIFT_KEPU_P3SOLVE.F
!*************************************************************************
! Returns the real root of cubic often found in solving kepler
! problem in universal variables.

!             Input:
!                 dt            ==>  time step (real scalar)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalar)
!                 mu            ==>  Reduced mass of system (real scalar)
!                 alpha         ==>  Twice the binding energy (real scalar)
!                 u             ==>  Vel. dot radial vector (real scalar)
!             Output:
!                 s             ==>  solution of cubic eqn for the  
!                                    universal variable
!                 iflg          ==>  success flag ( = 0 if O.K.) (integer)

! Author:  Martin Duncan  
! Date:    March 12/93
! Last revision: March 12/93

subroutine drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)

  implicit none

  !...  Inputs: 
  real(double_precision) :: dt,r0,mu,alpha,u

  !...  Outputs:
  integer :: iflg
  real(double_precision) :: s

  !...  Internals:
  real(double_precision) :: denom,a0,a1,a2,q,r,sq2,sq,p1,p2

  !----
  !...  Executable code 

  denom = (mu - alpha*r0)/6.d0
  a2 = 0.5*u/denom
  a1 = r0/denom
  a0 =-dt/denom

  q = (a1 - a2*a2/3.d0)/3.d0
  r = (a1*a2 -3.d0*a0)/6.d0 - (a2**3)/27.d0
  sq2 = q**3 + r**2

  if( sq2 .ge. 0.d0) then
     sq = sqrt(sq2)

     if ((r+sq) .le. 0.d0) then
        p1 =  -(-(r + sq))**(1.d0/3.d0)
     else
        p1 = (r + sq)**(1.d0/3.d0)
     endif
     if ((r-sq) .le. 0.d0) then
        p2 =  -(-(r - sq))**(1.d0/3.d0)
     else
        p2 = (r - sq)**(1.d0/3.d0)
     endif

     iflg = 0
     s = p1 + p2 - a2/3.d0

  else
     iflg = 1
     s = 0
  endif

  return
end subroutine drift_kepu_p3solve     !   drift_kepu_p3solve
!-------------------------------------------------------------------

!*************************************************************************
!                        DRIFT_KEPU_STUMPFF.F
!*************************************************************************
! subroutine for the calculation of stumpff functions
! see Danby p.172  equations 6.9.15

!             Input:
!                 x             ==>  argument
!             Output:
!                 c0,c1,c2,c3   ==>  c's from p171-172
!                                       (real scalors)
! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 2/3/93
! Modified by JEC: 31/3/98

subroutine drift_kepu_stumpff(x,c0,c1,c2,c3)

  use mercury_constant

  implicit none


  !...  Inputs: 
  real(double_precision) :: x

  !...  Outputs:
  real(double_precision) :: c0,c1,c2,c3

  !...  Internals:
  integer :: n,i
  real(double_precision) :: xm,x2,x3,x4,x5,x6

  !----
  !...  Executable code 

  n = 0
  xm = 0.1
  do while(abs(x).ge.xm)
     n = n + 1
     x = x * .25d0
  enddo
  
  x2 = x  * x
  x3 = x  * x2
  x4 = x2 * x2
  x5 = x2 * x3
  x6 = x3 * x3
  
  c2 = 1.147074559772972d-11*x6 - 2.087675698786810d-9*x5 + 2.755731922398589d-7*x4  - 2.480158730158730d-5*x3&
       + 1.388888888888889d-3*x2  - 4.166666666666667d-2*x + .5d0
  
  c3 = 7.647163731819816d-13*x6 - 1.605904383682161d-10*x5 + 2.505210838544172d-8*x4  - 2.755731922398589d-6*x3&
       + 1.984126984126984d-4*x2  - 8.333333333333333d-3*x + 1.666666666666667d-1
  
  c1 = 1. - x*c3
  c0 = 1. - x*c2
  
  if(n.ne.0) then
     do i=n,1,-1
        c3 = (c2 + c0*c3)*.25d0
        c2 = c1*c1*.5d0
        c1 = c0*c1
        c0 = 2.*c0*c0 - 1.
        x = x * 4.
     enddo
  endif

  return
end subroutine drift_kepu_stumpff     !   drift_kepu_stumpff


end module drift
