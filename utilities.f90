module utilities

!*************************************************************
!** Modules that gather various functions about string manipulation
!** and things that are perfectly separated from mercury particuliar
!** behaviour
!**
!** Version 1.0 - june 2011
!*************************************************************
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MXX_SORT.FOR    (ErikSoft 24 May 1997)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Sorts an array X, of size N, using Shell's method. Also returns an array
! INDEX that gives the original index of every item in the sorted array X.
!
! N.B. The maximum array size is 29523.
! ===
!
!------------------------------------------------------------------------------
!
subroutine mxx_sort (n,x,index)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: n,index(n)
  real(double_precision) :: x(n)
  !
  ! Local
  integer :: i,j,k,l,m,inc,incarr(9),iy
  real(double_precision) :: y
  data incarr/1,4,13,40,121,364,1093,3280,9841/
  !
  !------------------------------------------------------------------------------
  !
  do i = 1, n
     index(i) = i
  end do
  !
  m = 0
10 m = m + 1
  if (incarr(m).lt.n) goto 10
  m = m - 1
  !
  do i = m, 1, -1
     inc = incarr(i)
     do j = 1, inc
        do k = inc, n - j, inc
           y = x(j+k)
           iy = index(j+k)
           do l = j + k - inc, j, -inc
              if (x(l).le.y) goto 20
              x(l+inc) = x(l)
              index(l+inc) = index(l)
           end do
20         x(l+inc) = y
           index(l+inc) = iy
        end do
     end do
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mxx_sort

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    MCE_BOX.FOR    (ErikSoft   30 September 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Given initial and final coordinates and velocities, the routine returns
! the X and Y coordinates of a box bounding the motion in between the
! end points.
!
! If the X or Y velocity changes sign, the routine performs a quadratic
! interpolation to estimate the corresponding extreme value of X or Y.
!
!------------------------------------------------------------------------------
!
subroutine mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod
  real(double_precision) :: h,x0(3,nbod), x1(3,nbod), v0(3,nbod),v1(3,nbod)
  real(double_precision) ::   xmin(nbod), xmax(nbod), ymin(nbod),ymax(nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: temp
  !
  !------------------------------------------------------------------------------
  !
  do j = 2, nbod
     xmin(j) = min (x0(1,j), x1(1,j))
     xmax(j) = max (x0(1,j), x1(1,j))
     ymin(j) = min (x0(2,j), x1(2,j))
     ymax(j) = max (x0(2,j), x1(2,j))
     !
     ! If velocity changes sign, do an interpolation
     if ((v0(1,j).lt.0.and.v1(1,j).gt.0).or.(v0(1,j).gt.0.and.v1(1,j).lt.0)) then
        temp = (v0(1,j)*x1(1,j) - v1(1,j)*x0(1,j)       - .5d0*h*v0(1,j)*v1(1,j)) / (v0(1,j) - v1(1,j))
        xmin(j) = min (xmin(j),temp)
        xmax(j) = max (xmax(j),temp)
     end if
     !
     if ((v0(2,j).lt.0.and.v1(2,j).gt.0).or.(v0(2,j).gt.0.and.v1(2,j).lt.0)) then
        temp = (v0(2,j)*x1(2,j) - v1(2,j)*x0(2,j)       - .5d0*h*v0(2,j)*v1(2,j)) / (v0(2,j) - v1(2,j))
        ymin(j) = min (ymin(j),temp)
        ymax(j) = max (ymax(j),temp)
     end if
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mce_box
!

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCE_MIN.FOR    (ErikSoft  1 December 1998)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates minimum value of a quantity D, within an interval H, given initial
! and final values D0, D1, and their derivatives D0T, D1T, using third-order
! (i.e. cubic) interpolation.
!
! Also calculates the value of the independent variable T at which D is a
! minimum, with respect to the epoch of D1.
!
! N.B. The routine assumes that only one minimum is present in the interval H.
! ===
!------------------------------------------------------------------------------
!
subroutine mce_min (d0,d1,d0t,d1t,h,d2min,tmin)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  real(double_precision) :: d0,d1,d0t,d1t,h,d2min,tmin
  !
  ! Local
  real(double_precision) :: a,b,c,temp,tau
  !
  !------------------------------------------------------------------------------
  !
  if (d0t*h.gt.0.or.d1t*h.lt.0) then
     if (d0.le.d1) then
        d2min = d0
        tmin = -h
     else
        d2min = d1
        tmin = 0.d0
     end if
  else
     temp = 6.d0*(d0 - d1)
     a = temp + 3.d0*h*(d0t + d1t)
     b = temp + 2.d0*h*(d0t + 2.d0*d1t)
     c = h * d1t
     !
     temp =-.5d0*(b + sign (sqrt(max(b*b - 4.d0*a*c,0.d0)), b) )
     if (temp.eq.0) then
        tau = 0.d0
     else
        tau = c / temp
     end if
     !
     ! Make sure TAU falls in the interval -1 < TAU < 0
     tau = min(tau, 0.d0)
     tau = max(tau, -1.d0)
     !
     ! Calculate TMIN and D2MIN
     tmin = tau * h
     temp = 1.d0 + tau
     d2min = tau*tau*((3.d0+2.d0*tau)*d0 + temp*h*d0t)    + temp*temp*((1.d0-2.d0*tau)*d1 + tau*h*d1t)
     !
     ! Make sure D2MIN is not negative
     d2min = max(d2min, 0.d0)
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mce_min

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_JD2Y.FOR    (ErikSoft  7 July 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts from Julian day number to Julian/Gregorian Calendar dates, assuming
! the dates are those used by the English calendar.
!
! Algorithm taken from `Practical Astronomy with your calculator' (1988)
! by Peter Duffett-Smith, 3rd edition, C.U.P.
!
! Algorithm for negative Julian day numbers (Julian calendar assumed) by
! J. E. Chambers.
!
! N.B. The output date is with respect to the Julian Calendar on or before
! ===  4th October 1582, and with respect to the Gregorian Calendar on or 
!      after 15th October 1582.
!
!
!------------------------------------------------------------------------------
!
subroutine mio_jd2y (jd0,year,month,day)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: year,month
  real(double_precision) :: jd0,day
  !
  ! Local
  integer :: i,a,b,c,d,e,g
  real(double_precision) :: jd,f,temp,x,y,z
  !
  !------------------------------------------------------------------------------
  !
  if (jd0.le.0) goto 50
  !
  jd = jd0 + 0.5d0
  i = sign( dint(dabs(jd)), jd )
  f = jd - 1.d0*i
  !
  ! If on or after 15th October 1582
  if (i.gt.2299160) then
     temp = (1.d0*i - 1867216.25d0) / 36524.25d0
     a = sign( dint(dabs(temp)), temp )
     temp = .25d0 * a
     b = i + 1 + a - sign( dint(dabs(temp)), temp )
  else
     b = i
  end if
  !
  c = b + 1524
  temp = (1.d0*c - 122.1d0) / 365.25d0
  d = sign( dint(dabs(temp)), temp )
  temp = 365.25d0 * d
  e = sign( dint(dabs(temp)), temp )
  temp = (c-e) / 30.6001d0
  g = sign( dint(dabs(temp)), temp )
  !
  temp = 30.6001d0 * g
  day = 1.d0*(c-e) + f - 1.d0*sign( dint(dabs(temp)), temp )
  !
  if (g.le.13) month = g - 1
  if (g.gt.13) month = g - 13
  !
  if (month.gt.2) year = d - 4716
  if (month.le.2) year = d - 4715
  !
  if (day.gt.32) then
     day = day - 32
     month = month + 1
  end if
  !
  if (month.gt.12) then
     month = month - 12
     year = year + 1
  end if
  return
  !
50 continue
  !
  ! Algorithm for negative Julian day numbers (Duffett-Smith doesn't work)
  x = jd0 - 2232101.5
  f = x - dint(x)
  if (f.lt.0) f = f + 1.d0
  y = dint(mod(x,1461.d0) + 1461.d0)
  z = dint(mod(y,365.25d0))
  month = int((z + 0.5d0) / 30.61d0)
  day = dint(z + 1.5d0 - 30.61d0*dble(month)) + f
  month = mod(month + 2, 12) + 1
  !
  year = 1399 + int (x / 365.25d0)
  if (x.lt.0) year = year - 1
  if (month.lt.3) year = year + 1
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mio_jd2y
!

!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_SPL.FOR    (ErikSoft  14 November 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Given a character string STRING, of length LEN bytes, the routine finds 
! the beginnings and ends of NSUB substrings present in the original, and 
! delimited by spaces. The positions of the extremes of each substring are 
! returned in the array DELIMIT.
! Substrings are those which are separated by spaces or the = symbol.
!
!------------------------------------------------------------------------------
!
subroutine mio_spl (length,string,nsub,delimit)

  implicit none

  !
  ! Input/Output
  integer, intent(in) :: length
  integer, intent(out) :: nsub!,delimit(2,100)
  integer, dimension(:,:), intent(out) :: delimit
  ! TODO virer les warnings à propos de 'Actual argument contains too few elements for dummy argument
  character*1 string(length)
  ! TODO character(len=length) :: string make an error on the outputs
  !
  ! Local
  integer :: j,k
  character*1 c
  
  ! If the first dimension of the array 'delimit' is not 2, return an error
  if (size(delimit,1) /= 2) then
    write(*,*) "mio_spl: The first dimension of 'delimit' must be of size 2"
    stop
  end if
  !
  !------------------------------------------------------------------------------
  !
  nsub = 0
  j = 0
  c = ' '
  delimit(1,1) = -1
  !
  ! Find the start of string
10 j = j + 1
  if (j.gt.length) goto 99
  c = string(j)
  if (c.eq.' '.or.c.eq.'=') goto 10
  !
  ! Find the end of string
  k = j
20 k = k + 1
  if (k.gt.length) goto 30
  c = string(k)
  if (c.ne.' '.and.c.ne.'=') goto 20
  !
  ! Store details for this string
30 nsub = nsub + 1
  delimit(1,nsub) = j
  delimit(2,nsub) = k - 1
  !
  if (k.lt.length) then
     j = k
     goto 10
  end if
  !
99 continue
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mio_spl

end module utilities
