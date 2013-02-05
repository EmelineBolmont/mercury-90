module contour

!*************************************************************
!** Modules that contains a subroutine to print a contour of a 3D plot for a 'z' value
!** HOW TO USE
!** EXAMPLE :
!** use get_contour
!** call get_contour(z,x,y,'contour.dat')
!** where 'z' is a matrix with first dimension the same length of the 'x' array and the second dimension of 'z' is the same dimension than 'y'
!** i.e z(1:10,1:5), x(1:10), y(1:5) will work
!**
!** Version 1.0 - sept 2011
!*************************************************************
  use types_numeriques

  implicit none
  
  private
  
  public :: get_contour, test_get_contour ! test_get_contour must be used only to test function, you should not use it out of the debug phase

contains

subroutine get_contour(matrix, x, y, filename, lvl)

implicit none

real(double_precision), dimension(:,:), intent(in) :: matrix
real(double_precision), dimension(:), intent(in) :: x, y
real(double_precision), intent(in) :: lvl
logical, dimension(:,:), allocatable :: not_tested ! in order to count only once each point
character(len=*), intent(in) :: filename

! Local
integer :: nb_x, nb_y
integer :: i,j
logical :: change_sign_L0, change_sign_L1
!-----

nb_x = size(x)
nb_y = size(y)

if (size(matrix,1).ne.nb_x) then
  stop 'Error: The first dimension of "matrix" must be the same size than the "x" array'
end if

if (size(matrix,2).ne.nb_y) then
  stop 'Error: The second dimension of "matrix" must be the same size than the "y" array'
end if

allocate(not_tested(nb_x, nb_y))
not_tested(1:nb_x,1:nb_y) = .True.

open(10, file=filename, status='replace')
close(10)
! First, we test right and top side separately because we only test one type of line each.

! We test the top side of the matrix
j=nb_y
do i=1,nb_x-1
  if (not_tested(i,j)) then
    change_sign_L0 = ((matrix(i,j) - lvl) * (matrix(i+1,j) - lvl)).lt.0.
    if (change_sign_L0) then
      call get_contour_line(i_max=nb_x, j_max=nb_y ,i1_start=i, j1_start=j,i2_start=i+1, j2_start=j, matrix=matrix, &
                              x=x, y=y, not_tested=not_tested, lvl=lvl, filename=filename)
    end if
    not_tested(i,j) = .False.
  end if
end do

! We test the right side of the matrix
i=nb_x
do j=1,nb_y-1
  if (not_tested(i,j)) then
    change_sign_L1 = ((matrix(i,j) - lvl) * (matrix(i,j+1) - lvl)).lt.0.
    if (change_sign_L1) then
      call get_contour_line(i_max=nb_x, j_max=nb_y ,i1_start=i, j1_start=j,i2_start=i, j2_start=j+1, matrix=matrix, &
                            x=x, y=y, not_tested=not_tested, lvl=lvl, filename=filename)
    end if
    not_tested(i,j) = .False.
  end if
end do

! We test the rest of the matrix.
do j=1,nb_y-1
  do i=1,nb_x-1
    if (not_tested(i,j)) then
      ! We search for a starting point of a contour line. Once we have found it, we call a routine to draw the entire contour 
      ! line, but since this is not included in the loop on i,j, we put in a logical array if we have drawn a point or not, in 
      ! order not to test twice the same point. if ((matrix(i-1,j).lt.lvl).or.(matrix(i+1,j).lt.lvl)) then
      change_sign_L0 = ((matrix(i,j) - lvl) * (matrix(i+1,j) - lvl)).lt.0.
      change_sign_L1 = ((matrix(i,j) - lvl) * (matrix(i,j+1) - lvl)).lt.0.
      if (change_sign_L0) then
        call get_contour_line(i_max=nb_x, j_max=nb_y ,i1_start=i, j1_start=j,i2_start=i+1, j2_start=j, matrix=matrix, &
                              x=x, y=y, not_tested=not_tested, lvl=lvl, filename=filename)
      else if (change_sign_L1) then
        call get_contour_line(i_max=nb_x, j_max=nb_y ,i1_start=i, j1_start=j,i2_start=i, j2_start=j+1, matrix=matrix, &
                              x=x, y=y, not_tested=not_tested, lvl=lvl, filename=filename)
      end if
      not_tested(i,j) = .False.
    end if
  end do
end do
!-----
end subroutine get_contour

subroutine get_contour_line(i_max, j_max, i1_start,j1_start,i2_start,j2_start, matrix, x, y, not_tested, lvl, filename)

implicit none
integer :: i_max, j_max
real(double_precision), dimension(i_max,j_max), intent(in) :: matrix
real(double_precision), dimension(i_max), intent(in) :: x
real(double_precision), dimension(j_max), intent(in) :: y
real(double_precision), intent(in) :: lvl
integer, intent(in) :: i1_start,j1_start,i2_start,j2_start
logical, dimension(i_max,j_max), intent(inout) :: not_tested
character(len=*), intent(in) :: filename


! Local 
real(double_precision) :: x_start, y_start ! the starting point of the contour line. Usefull to check if w've come at the end of the contour line (if (i,j) equal theses starting values)
integer :: i1, j1, i2, j2, i3, j3, i4, j4 ! Points that define a box for the check. see the drawing for details
integer :: i_temp, j_temp ! temporary variables to exchange two positions
! a vector to tell the direction to seek the isoline. This direction will of 
! course change during the search and allow us to test only 3 monolines.
integer, dimension(2) :: direction, tmp_dir 
logical :: end_of_line ! set to True if we must end the seek of the line because we have reached the end
real(double_precision) :: x_cont, y_cont ! the current (x,y) coordinates of the last contour point found
logical :: change_sign_L1, change_sign_L2 ! (boolean to tell if values of edges have differents signs around the lvl value)

! we store initial position
i1 = i1_start
j1 = j1_start
i2 = i2_start
j2 = j2_start

end_of_line = .False.


! We search for the starting point of the isoline
if (i1.eq.i2) then
  x_start = x(i1)
  y_start = y(j1) + (y(j2) - y(j1)) / (matrix(i2,j2) - matrix(i1,j1)) * (lvl - matrix(i1,j1))
else
  y_start = y(j1)
  x_start = x(i1) + (x(i2) - x(i1)) / (matrix(i2,j2) - matrix(i1,j1)) * (lvl - matrix(i1,j1))
end if

! We search for the direction at the starting line
if (i1.eq.i2) then
  direction(2) = 0
  if ((j2.gt.j1).and.(i1.ne.1)) then
    direction(1) = -1
  else if ((j2.gt.j1).and.(i1.eq.1)) then
  ! the line is on the left side of the matrix. We revert points 1 and 2 and change direction
    direction(1) = 1
    i_temp = i1
    j_temp = j1
    i1 = i2
    j1 = j2
    i2 = i_temp
    j2 = j_temp
  else if ((j1.gt.j2).and.(i1.ne.i_max)) then
    direction(1) = 1
  else
  ! the line is on the right side of the matrix. We revert points 1 and 2 and change direction
    direction(1) = -1
    i_temp = i1
    j_temp = j1
    i1 = i2
    j1 = j2
    i2 = i_temp
    j2 = j_temp
  end if
else 
  direction(1) = 0
  if ((i1.gt.i2).and.(j1.ne.1)) then
    direction(2) = -1
  else if ((i1.gt.i2).and.(j1.eq.1)) then
  ! the line is on the bottom side of the matrix. We revert points 1 and 2 and change direction
    direction(2) = 1
    i_temp = i1
    j_temp = j1
    i1 = i2
    j1 = j2
    i2 = i_temp
    j2 = j_temp
  else if ((i1.lt.i2).and.(j1.ne.j_max)) then
    direction(2) = 1
  else
  ! the line is on the top side of the matrix. We revert points 1 and 2 and change direction
    direction(2) = -1
    i_temp = i1
    j_temp = j1
    i1 = i2
    j1 = j2
    i2 = i_temp
    j2 = j_temp
  end if
end if

!                      3/ else : this line
!                         hereafter :
!                         LINE 3
!            (i3,j3) +-------------------+ (i4,j4)
!                    |                   |
!                    |                   |
!                    |                   |
!     1/ first line  |    Direction of   | 2/ second line
!        to check    |      research     |    to check
!        hereafter : |         .         |    hereafter :
!         LINE 1     |        / \        |    LINE 2
!                    |         |         |
!                    |         |         |
!            (i1,j1) +-------------------+ (i2,j2)
!                     Reference LINE 0 we
!                     know has a change in sign

open(10, file=filename, access='append')
write(10,*) x_start, y_start, lvl

do while (.not.end_of_line)
  i3 = i1 + direction(1)
  j3 = j1 + direction(2)
  
  i4 = i2 + direction(1)
  j4 = j2 + direction(2)
  
  ! We mark the 4 points as tested to avoid outputing the same contour several times
  not_tested(i1,j1) = .False.
  not_tested(i2,j2) = .False.
  not_tested(i3,j3) = .False.
  not_tested(i4,j4) = .False.

  change_sign_L1 = ((matrix(i1,j1) - lvl) * (matrix(i3,j3) - lvl)).lt.0.
  change_sign_L2 = ((matrix(i2,j2) - lvl) * (matrix(i4,j4) - lvl)).lt.0.
  
  ! We check first if the change in sign is on the LINE1
  if (change_sign_L1) then
    ! We test to know if its an x=cte or y=cte line
    if (direction(1).eq.0) then
      x_cont = x(i1)
      y_cont = y(j1) + (y(j3) - y(j1)) / (matrix(i3,j3) - matrix(i1,j1)) * (lvl - matrix(i1,j1))
    else
      x_cont = x(i1) + (x(i3) - x(i1)) / (matrix(i3,j3) - matrix(i1,j1)) * (lvl - matrix(i1,j1))
      y_cont = y(j1)
    end if
    
    ! We determine the new values for i1,j1,i2,j2
    ! NOTE that i1 and j1 doesn't change, but the direction does.
    i2 = i3
    j2 = j3
    ! We make a rotation of 90 degrees
    tmp_dir(1:2) = direction(1:2)
    direction(1) = -tmp_dir(2)
    direction(2) =  tmp_dir(1)
  ! We then check if the change in sign is on the LINE 2
  else if (change_sign_L2) then
    ! We test to know if its an x=cte or y=cte line
    if (direction(1).eq.0) then
      x_cont = x(i2)
      y_cont = y(j2) + (y(j4) - y(j2)) / (matrix(i4,j4) - matrix(i2,j2)) * (lvl - matrix(i2,j2))
    else
      x_cont = x(i2) + (x(i4) - x(i2)) / (matrix(i4,j4) - matrix(i2,j2)) * (lvl - matrix(i2,j2))
      y_cont = y(j2)
    end if
    
    ! We determine the new values for i1,j1,i2,j2
    ! NOTE that i2 and j2 doesn't change, but the direction does.
    i1 = i4
    j1 = j4
    ! We make a rotation of -90 degrees
    tmp_dir(1:2) = direction(1:2)
    direction(1) =  tmp_dir(2)
    direction(2) = -tmp_dir(1)
  ! So the change must be in the LINE 3
  else
    if (direction(1).eq.0) then
      x_cont = x(i4) + (x(i3) - x(i4)) / (matrix(i3,j3) - matrix(i4,j4)) * (lvl - matrix(i4,j4))
      y_cont = y(j3)
    else
      x_cont = x(i3)
      y_cont = y(j4) + (y(j3) - y(j4)) / (matrix(i3,j3) - matrix(i4,j4)) * (lvl - matrix(i4,j4))
    end if
    ! We determine the new values for i1,j1,i2,j2
    ! Note that the direction doesn't change here, we advance in straight line.
    i1 = i3
    j1 = j3
    i2 = i4
    j2 = j4
  end if
  ! We output in the file the current point of the contour
  write(10,*) x_cont, y_cont, lvl
  
  ! We test if we are at the end of the contour, i.e we are at the starting point again, or at an edge of the matrix.
  if ((i1.eq.i1_start).and.(i2.eq.i2_start).and.(j1.eq.j1_start).and.(j2.eq.j2_start)) then
    end_of_line = .True.
    ! In the case of the closed line, we redo the print of the first point to have a close line
    write(10,*) x_start, y_start, lvl
  else if ((i1.eq.i2).and.(i1.eq.1)) then
    ! Left border of the matrix
    end_of_line = .True.
  else if ((j1.eq.j2).and.(j1.eq.1)) then
    ! Low border of the matrix
    end_of_line = .True.
  else if ((i1.eq.i2).and.(i1.eq.i_max)) then
    ! right border of the matrix
    end_of_line = .True.
  else if ((j1.eq.j2).and.(j1.eq.j_max)) then
    ! Top border of the matrix
    end_of_line = .True.
  end if
end do
! At the end of each single contour line, we output a blank line to separate them in the output file (usefull for GNUPLOT in particular)
write(10,*) ''
close(10)

end subroutine get_contour_line

subroutine test_get_contour()
  implicit none
  integer, parameter :: dimx = 10
  integer, parameter :: dimy = 10
  real(double_precision), dimension(dimx, dimy) :: graph
  real(double_precision), dimension(dimx) :: x
  real(double_precision), dimension(dimy) :: y
  real(double_precision), parameter :: circle_radius = 10.
  integer :: i, j
  
  x = (/(i*1, i=1,dimx)/)
  y = (/(i*1., i=1,dimy)/)
  
  graph(:,:) = -1.
  
  graph(1:2,1:2)   = 1.  ! left bottom corner
  graph(1:2,5:6)   = 1.  ! left line
  graph(1:2,9:10)  = 1. ! left top corner
  graph(5:6,9:10)  = 1. ! top line
  graph(9:10,9:10) = 1. ! right top corner
  graph(9:10,5:6)  = 1. ! right line
  graph(9:10,1:2)  = 1. ! right bottom corner
  graph(5:6,1:2)   = 1.  ! bottom line
  graph(5:6,5:6)   = 1.  ! middle
  
  ! Contains the real map taht we try to contour in 'contour.dat'
  open(10, file='test_contour.dat')
  do i=1, dimx
    do j=1, dimy
      
      write(10,*) x(i), y(j), graph(i,j)
    end do
    write(10,*) ''
  end do
  close(10)
  
  ! Gnuplot script to launch and see if the contour match the map
  open(10, file='test_contour.gnuplot')
  write(10,*) 'set pm3d map'
  write(10,*) 'set pm3d explicit'
  write(10,*) 'set palette rgbformulae 22,13,-31'
  write(10,*) 'splot "test_contour.dat" with pm3d notitle, "contour.dat" with line linetype -1 notitle'
  write(10,*) 'pause -1 # wait until a carriage return is hit'
  
  close(10)
  
  call get_contour(matrix=graph, x=x, y=y, filename='contour.dat', lvl=0.d0)
end subroutine test_get_contour

end module contour