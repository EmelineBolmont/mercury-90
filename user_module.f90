module user_module

!*************************************************************
!** Modules that contains user defined modules. 
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
! ===
!------------------------------------------------------------------------------

subroutine mfo_user (time,jcen,nbod,nbig,m,x,v,a)
  
  use physical_constant
  use mercury_constant  

  implicit none

  
  ! Input/Output
  integer, intent(in) :: nbod, nbig
  real(double_precision),intent(in) :: time,jcen(3),m(nbod),x(3,nbod),v(3,nbod)
  real(double_precision),intent(out) :: a(3,nbod)
  
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
end subroutine mfo_user



end module user_module
