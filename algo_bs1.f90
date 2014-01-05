!******************************************************************************
! MODULE: algo_bs1
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that gather various functions about the BS
!! algorithm.\n\n
!! The Bulirsch-Stoer algorithms are described in W.H.Press et al. (1992)
!! ``Numerical Recipes in Fortran'', pub. Cambridge.
!
!******************************************************************************

module algo_bs1

  use types_numeriques

  private
  
  real(double_precision), parameter :: SHRINK=.55d0
  real(double_precision), parameter :: GROW=1.3d0
  
  public :: mdt_bs1
  
  contains
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MDT_BS1.FOR    (ErikSoft   2 March 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Integrates NBOD bodies (of which NBIG are Big) for one timestep H0
! using the Bulirsch-Stoer method. The accelerations are calculated using the 
! subroutine FORCE. The accuracy of the step is approximately determined 
! by the tolerance parameter TOL.

! N.B. Input/output must be in coordinates with respect to the central body.
! ===

!------------------------------------------------------------------------------

subroutine mdt_bs1 (time,h0,hdid,tol,jcen,nbod,nbig,mass,x0,v0,s,rphys,rcrit,ngf,stat,dtflag,ngflag,nce,ice,jce,force)
  use physical_constant
  use mercury_constant
  use mercury_globals

  implicit none
  
  ! Input/Output
  integer, intent(in) :: nbod
  integer, intent(in) :: nbig
  integer, intent(in) :: stat(nbod)
  integer, intent(in) :: dtflag
  integer, intent(in) :: ngflag
  integer, intent(in) :: nce
  integer, intent(in) :: ice(nce)
  integer, intent(in) :: jce(nce)
  real(double_precision), intent(in) :: time
  real(double_precision), intent(in) :: tol
  real(double_precision), intent(in) :: jcen(3)
  real(double_precision), intent(in) :: mass(nbod)
  real(double_precision), intent(in) :: s(3,nbod)
  real(double_precision), intent(in) :: ngf(4,nbod)
  real(double_precision), intent(in) :: rphys(nbod)
  real(double_precision), intent(in) :: rcrit(nbod)
  
  real(double_precision), intent(out) :: hdid
  
  real(double_precision), intent(inout) :: h0
  real(double_precision), intent(inout) :: x0(3,nbod)
  real(double_precision), intent(inout) :: v0(3,nbod)
  
  external force
  
  ! Local
  integer :: j, j1, k, n
  real(double_precision) :: tmp0,tmp1,tmp2,errmax,tol2,h,hx2,h2(8)
  real(double_precision) :: x(3,nb_bodies_initial),v(3,nb_bodies_initial),xend(3,nb_bodies_initial),vend(3,nb_bodies_initial)
  real(double_precision) :: a(3,nb_bodies_initial),a0(3,nb_bodies_initial),d(6,nb_bodies_initial,8)
  real(double_precision) :: xscal(nb_bodies_initial),vscal(nb_bodies_initial)
  
  !------------------------------------------------------------------------------
  
  tol2 = tol * tol
  
  ! Calculate arrays used to scale the relative error (R^2 for position and
  ! V^2 for velocity).
  do k = 2, nbod
     tmp1 = x0(1,k)*x0(1,k) + x0(2,k)*x0(2,k) + x0(3,k)*x0(3,k)
     tmp2 = v0(1,k)*v0(1,k) + v0(2,k)*v0(2,k) + v0(3,k)*v0(3,k)
     xscal(k) = 1.d0 / tmp1
     vscal(k) = 1.d0 / tmp2
  end do
  
  ! Calculate accelerations at the start of the step
  call force (time,jcen,nbod,nbig,mass,x0,v0,s,rcrit,a0,stat,ngf,ngflag,nce,ice,jce)
  
  do
     
     ! For each value of N, do a modified-midpoint integration with 2N substeps
     do n = 1, 8
        h = h0 / (2.d0 * float(n))
        h2(n) = .25d0 / (n*n)
        hx2 = h * 2.d0
        
        do k = 2, nbod
           x(1,k) = x0(1,k) + h*v0(1,k)
           x(2,k) = x0(2,k) + h*v0(2,k)
           x(3,k) = x0(3,k) + h*v0(3,k)
           v(1,k) = v0(1,k) + h*a0(1,k)
           v(2,k) = v0(2,k) + h*a0(2,k)
           v(3,k) = v0(3,k) + h*a0(3,k)
        end do
        call force (time,jcen,nbod,nbig,mass,x,v,s,rcrit,a,stat,ngf,ngflag,nce,ice,jce)
        do k = 2, nbod
           xend(1,k) = x0(1,k) + hx2*v(1,k)
           xend(2,k) = x0(2,k) + hx2*v(2,k)
           xend(3,k) = x0(3,k) + hx2*v(3,k)
           vend(1,k) = v0(1,k) + hx2*a(1,k)
           vend(2,k) = v0(2,k) + hx2*a(2,k)
           vend(3,k) = v0(3,k) + hx2*a(3,k)
        end do
        
        do j = 2, n
           call force (time,jcen,nbod,nbig,mass,xend,vend,s,rcrit,a,stat,ngf,ngflag,nce,ice,jce)
           do k = 2, nbod
              x(1,k) = x(1,k) + hx2*vend(1,k)
              x(2,k) = x(2,k) + hx2*vend(2,k)
              x(3,k) = x(3,k) + hx2*vend(3,k)
              v(1,k) = v(1,k) + hx2*a(1,k)
              v(2,k) = v(2,k) + hx2*a(2,k)
              v(3,k) = v(3,k) + hx2*a(3,k)
           end do
           call force (time,jcen,nbod,nbig,mass,x,v,s,rcrit,a,stat,ngf,ngflag,nce,ice,jce)
           do k = 2, nbod
              xend(1,k) = xend(1,k) + hx2*v(1,k)
              xend(2,k) = xend(2,k) + hx2*v(2,k)
              xend(3,k) = xend(3,k) + hx2*v(3,k)
              vend(1,k) = vend(1,k) + hx2*a(1,k)
              vend(2,k) = vend(2,k) + hx2*a(2,k)
              vend(3,k) = vend(3,k) + hx2*a(3,k)
           end do
        end do
        
        call force (time,jcen,nbod,nbig,mass,xend,vend,s,rcrit,a,stat,ngf,ngflag,nce,ice,jce)
        
        do k = 2, nbod
           d(1,k,n) = .5d0*(xend(1,k) + x(1,k) + h*vend(1,k))
           d(2,k,n) = .5d0*(xend(2,k) + x(2,k) + h*vend(2,k))
           d(3,k,n) = .5d0*(xend(3,k) + x(3,k) + h*vend(3,k))
           d(4,k,n) = .5d0*(vend(1,k) + v(1,k) + h*a(1,k))
           d(5,k,n) = .5d0*(vend(2,k) + v(2,k) + h*a(2,k))
           d(6,k,n) = .5d0*(vend(3,k) + v(3,k) + h*a(3,k))
        end do
        
        ! Update the D array, used for polynomial extrapolation
        do j = n - 1, 1, -1
           j1 = j + 1
           tmp0 = 1.d0 / (h2(j) - h2(n))
           tmp1 = tmp0 * h2(j1)
           tmp2 = tmp0 * h2(n)
           do k = 2, nbod
              d(1,k,j) = tmp1 * d(1,k,j1)  -  tmp2 * d(1,k,j)
              d(2,k,j) = tmp1 * d(2,k,j1)  -  tmp2 * d(2,k,j)
              d(3,k,j) = tmp1 * d(3,k,j1)  -  tmp2 * d(3,k,j)
              d(4,k,j) = tmp1 * d(4,k,j1)  -  tmp2 * d(4,k,j)
              d(5,k,j) = tmp1 * d(5,k,j1)  -  tmp2 * d(5,k,j)
              d(6,k,j) = tmp1 * d(6,k,j1)  -  tmp2 * d(6,k,j)
           end do
        end do
        
        ! After several integrations, test the relative error on extrapolated values
        if (n.gt.3) then
           errmax = 0.d0
           
           ! Maximum relative position and velocity errors (last D term added)
           do k = 2, nbod
              tmp1 = max( d(1,k,1)*d(1,k,1), d(2,k,1)*d(2,k,1),        d(3,k,1)*d(3,k,1) )
              tmp2 = max( d(4,k,1)*d(4,k,1), d(5,k,1)*d(5,k,1),        d(6,k,1)*d(6,k,1) )
              errmax = max(errmax, tmp1*xscal(k), tmp2*vscal(k))
           end do
           
           ! If error is smaller than TOL, update position and velocity arrays, and exit
           if (errmax.le.tol2) then
              do k = 2, nbod
                 x0(1,k) = d(1,k,1)
                 x0(2,k) = d(2,k,1)
                 x0(3,k) = d(3,k,1)
                 v0(1,k) = d(4,k,1)
                 v0(2,k) = d(5,k,1)
                 v0(3,k) = d(6,k,1)
              end do
              
              do j = 2, n
                 do k = 2, nbod
                    x0(1,k) = x0(1,k) + d(1,k,j)
                    x0(2,k) = x0(2,k) + d(2,k,j)
                    x0(3,k) = x0(3,k) + d(3,k,j)
                    v0(1,k) = v0(1,k) + d(4,k,j)
                    v0(2,k) = v0(2,k) + d(5,k,j)
                    v0(3,k) = v0(3,k) + d(6,k,j)
                 end do
              end do
              
              ! Save the actual stepsize used
              hdid = h0
              
              ! Recommend a new stepsize for the next call to this subroutine
              if (n.eq.8) h0 = h0 * SHRINK
              if (n.lt.7) h0 = h0 * GROW
              return
           end if
        end if
        
     end do
     
     ! If errors were too large, redo the step with half the previous step size.
     h0 = h0 * .5d0
  end do
  
  !------------------------------------------------------------------------------
  
end subroutine mdt_bs1
  
end module algo_bs1
