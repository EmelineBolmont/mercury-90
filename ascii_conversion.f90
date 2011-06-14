module ascii_conversion

!*************************************************************
!** Modules that allow conversion between reals and ascii strings
!** in order to store those strings in files and thus compress data
!** Version 1.0 - june 2011
!*************************************************************
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_C2FL.FOR    (ErikSoft  1 July 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Converts a CHARACTER*8 ASCII string into a REAL*8 variable.
!
! N.B. X will lie in the range -1.e112 < X < 1.e112
! ===
!
!------------------------------------------------------------------------------
!
function mio_c2fl (c)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  real(double_precision) :: mio_c2fl
  character*8 c
  !
  ! Local
  integer :: ex
  real(double_precision) :: x
  !
  !------------------------------------------------------------------------------
  !
  x = mio_c2re (c(1:8), 0.d0, 1.d0, 7)
  x = x * 2.d0 - 1.d0
  ex = ichar(c(8:8)) - 32 - 112
  mio_c2fl = x * (10.d0**dble(ex))
  !
  !------------------------------------------------------------------------------
  !
  return
end function mio_c2fl
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_C2RE.FOR    (ErikSoft  1 July 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts an ASCII string into a REAL*8 variable X, where XMIN <= X < XMAX,
! using the new format compression:
!
! X is assumed to be made up of NCHAR base-224 digits, each one represented
! by a character in the ASCII string. Each digit is given by the ASCII
! number of the character minus 32.
! The first 32 ASCII characters (CTRL characters) are avoided, because they
! cause problems when using some operating systems.
!
!------------------------------------------------------------------------------
!
function mio_c2re (c,xmin,xmax,nchar)
  !
  use types_numeriques

  implicit none

  !
  ! Input/output
  integer :: nchar
  real(double_precision) :: xmin,xmax,mio_c2re
  character*8 c
  !
  ! Local
  integer :: j
  real(double_precision) :: y
  !
  !------------------------------------------------------------------------------
  !
  y = 0
  do j = nchar, 1, -1
     y = (y + dble(ichar(c(j:j)) - 32)) / 224.d0
  end do
  !
  mio_c2re = xmin  +  y * (xmax - xmin)
  !
  !------------------------------------------------------------------------------
  !
  return
end function mio_c2re

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_FL2C.FOR    (ErikSoft  1 July 1998)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts a (floating point) REAL*8 variable X, into a CHARACTER*8 ASCII 
! string, using the new format compression:
!
! X is first converted to base 224, and then each base 224 digit is converted 
! to an ASCII character, such that 0 -> character 32, 1 -> character 33...
! and 223 -> character 255.
! The first 7 characters in the string are used to store the mantissa, and the
! eighth character is used for the exponent.
!
! ASCII characters 0 - 31 (CTRL characters) are not used, because they
! cause problems when using some operating systems.
!
! N.B. X must lie in the range -1.e112 < X < 1.e112
! ===
!
!------------------------------------------------------------------------------
!
function mio_fl2c (x)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  real(double_precision) :: x
  character*8 mio_fl2c
  !
  ! Local
  integer :: ex
  real(double_precision) :: ax,y

  !
  !------------------------------------------------------------------------------
  !
  if (x.eq.0) then
     y = .5d0
  else
     ax = abs(x)
     ex = int(log10(ax))
     if (ax.ge.1) ex = ex + 1
     y = ax*(10.d0**(-ex))
     if (y.eq.1) then
        y = y * .1d0
        ex = ex + 1
     end if
     y = sign(y,x) *.5d0 + .5d0
  end if
  !
  mio_fl2c(1:8) = mio_re2c (y, 0.d0, 1.d0)
  ex = ex + 112
  if (ex.gt.223) ex = 223
  if (ex.lt.0) ex = 0
  mio_fl2c(8:8) = char(ex+32)
  !
  !------------------------------------------------------------------------------
  !
  return
end function mio_fl2c

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_RE2C.FOR    (ErikSoft  27 June 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts a REAL*8 variable X, where XMIN <= X < XMAX, into an ASCII string
! of 8 characters, using the new format compression: 
!
! X is first converted to base 224, and then each base 224 digit is converted 
! to an ASCII character, such that 0 -> character 32, 1 -> character 33...
! and 223 -> character 255.
!
! ASCII characters 0 - 31 (CTRL characters) are not used, because they
! cause problems when using some operating systems.
!
!------------------------------------------------------------------------------
!
function mio_re2c (x,xmin,xmax)
  !
  use types_numeriques

  implicit none

  !
  ! Input/output
  real(double_precision) :: x,xmin,xmax
  character*8 mio_re2c
  !
  ! Local
  integer :: j
  real(double_precision) :: y,z
  !
  !------------------------------------------------------------------------------
  !
  mio_re2c(1:8) = '        '
  y = (x - xmin) / (xmax - xmin)
  !
  if (y.ge.1) then
     do j = 1, 8
        mio_re2c(j:j) = char(255)
     end do
  else if (y.gt.0) then
     z = y
     do j = 1, 8
        z = mod(z, 1.d0) * 224.d0
        mio_re2c(j:j) = char(int(z) + 32)
     end do
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end function mio_re2c


end module ascii_conversion
