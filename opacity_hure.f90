module opacity_hure
  use types_numeriques
  use physical_constant
  
  implicit none
  
  private
  
  public :: get_opacity_hure, init_opacity_hure
  
  real(double_precision) :: dlog_r ! Space between two 'r' values in the opacity table (this width must be constant)
  real(double_precision) :: dlog_t ! Space between two temperature values in the opacity table (this width must be constant)
  
  integer :: nb_r_points ! number of density points
  integer :: nb_t_points ! number of temperature points

  real(double_precision), dimension(:), allocatable :: hure_r ! Array for R, which is a mix parameter between T and rho
  real(double_precision), dimension(:), allocatable :: hure_t ! temperature points
  real(double_precision), dimension(:,:), allocatable :: hure_k_table ! absorption coefficients (first dimension is T, the second is R
  
contains

  subroutine init_opacity_hure()
    ! Initialisation of the opacity table by reading the data file.
    implicit none
    !     ..local..
    integer :: q                 ! counter for density points (q<=mx)
    integer :: r                 ! counter for temperature points (r<=my)
    !     ..body..
    open(unit=10,file='opacity_table_hure.dat',status='old') 
    
    !     ..reads size of the grid..
    read(10,'(///,2(i6.6,1x))') nb_r_points,nb_t_points
    
    if (.not.allocated(hure_k_table)) then
      allocate(hure_r(nb_r_points))
      allocate(hure_t(nb_t_points))
      allocate(hure_k_table(nb_t_points, nb_r_points))
    end if
    
    read(10,'(4x, 24(1x,f6.3))') hure_r(1:nb_r_points)
    
    do r=1,nb_t_points
       !     ..reads temperature scale and absorption coefficients..
       read(10,'(f5.3,sp,24(1x,f6.3),ss)') hure_t(r),hure_k_table(r, 1:nb_r_points)
    enddo
    close(10)
    !     ..scale on density..
    
    ! We initialise the two width, for temperature and r, between each sample value
    ! /!\ these widths must be constant over the opacity table
    dlog_t = hure_t(2) - hure_t(1)
    dlog_r = hure_r(2) - hure_r(1)
    
    return
  end subroutine init_opacity_hure

  subroutine linear_interpolation(x1, y1, x2, y2, x, y)
  
  implicit none
  
  real(double_precision), intent(in) :: x1, y1, x2, y2 ! the two points used to interpolate
  real(double_precision), intent(in) :: x ! the x coordinate of the point where we want the y coordinate
  real(double_precision), intent(out) :: y
  
  y = (y2 - y1) * (x - x1) / (x2 - x1) + y1
  
  return 
  end subroutine linear_interpolation

  function get_opacity_hure(temperature, num_bulk_density)
  ! subroutine that return the opacity of the disk at the location of the planet given various parameters
  ! must have the structure defined in disk_parameters.f90 to allow pointers to be used
  
  implicit none
  real(double_precision), intent(in) :: temperature ! temperature of the disk [K]
  real(double_precision), intent(in) :: num_bulk_density ! bulk density of the gas disk [MSUN/AU^3] (in numerical units)
  
  real(double_precision) :: get_opacity_hure ! the opacity in numerical units (AU**2/MSUN)
  
  ! parameters
  real(double_precision), parameter :: num_to_phys_bulk_density = MSUN / AU**3
  real(double_precision), parameter :: phys_to_num_opacity = MSUN / AU**2
  
  ! locals
  real(double_precision) :: bulk_density
  real(double_precision) :: log_t, log_r
  real(double_precision) :: t1, t2, r1, r2
  integer :: id_t1, id_t2, id_r1, id_r2
  real(double_precision) :: p11, p12, p21, p22
  real(double_precision) :: p1, p2 ! temporary points for the interpolation done on temperature axis.
  
  !-----------------------------
  ! we convert the bulk_density from numerical units(AU, MS, DAY) to physical units (CGS)
  bulk_density = num_to_phys_bulk_density * num_bulk_density
  
  log_t = log10(temperature)
  log_r = log10(bulk_density) + 18.d0 - 3.d0 * log_t
  
  id_t1 = 1 + int((log_t - hure_t(1)) / dlog_t)
  id_r1 = 1 + int((log_r - hure_r(1)) / dlog_r)
  
  ! if exterior or equal to the boundaries, we consider the boundaries, without interpolation. 
  ! actually there will be an interpolation, but between equal values.
  if (log_t.lt.hure_t(1)) then
    id_t1 = 1
    id_t2 = id_t1
  else if (log_t.ge.hure_t(nb_t_points)) then
    id_t1 = nb_t_points
    id_t2 = id_t1
  else 
    id_t2 = id_t1 +1
  end if
  
  if (log_r.lt.hure_r(1)) then
    id_r1 = 1
    id_r2 = id_r1
  else if (log_r.ge.hure_t(nb_r_points)) then
    id_r1 = nb_r_points
    id_r2 = id_r1
  else
    id_r2 = id_r1 + 1
  end if
  
  t1 = hure_t(id_t1)
  t2 = hure_t(id_t2)
  r1 = hure_r(id_r1)
  r2 = hure_r(id_r2)
  
  
  !    p12                     p22
  ! r2 +-----------------------+
  !    |                       |
  !    |                       |
  !    |                       |
  !    |                       |
  !    |                       |
  !    |                       |
  ! r1 +-----------------------+
  !    p11                     p21
  !    t1                      t2
  
  p11 = hure_k_table(id_t1, id_r1)
  p12 = hure_k_table(id_t1, id_r2)
  p21 = hure_k_table(id_t2, id_r1)
  p22 = hure_k_table(id_t2, id_r2)
  
  if (id_t1.ne.id_t2) then
    call linear_interpolation(x1=t1, y1=p11, x2=t2, y2=p21, x=log_t, y=p1) ! for r=r1
    call linear_interpolation(x1=t1, y1=p12, x2=t2, y2=p22, x=log_t, y=p2) ! for r=r2
  else
    p1 = p11
    p2 = p22
  end if
  
  if (id_r1.ne.id_r2) then
    call linear_interpolation(x1=r1, y1=p1, x2=r2, y2=p2, x=log_r, y=get_opacity_hure)
  else
    get_opacity_hure = p1
  end if
  
!~   write(*,*) '----------------------------------'
!~   write(*,*) 'temperature=', temperature, "; bulk_density=", bulk_density
!~   write(*,*) 'log_T =', log_t, '; log_r =', log_r
!~   write(*,*) id_t1, id_t2, id_r1, id_r2
!~   write(*,*) p11, p12, p21, p22
!~   write(*,*) p1, p2
!~   write(*,*) get_opacity_hure
!~   write(*,*) '----------------------------------'
  
  ! we do not want log(opacity) but simply the opacity
  get_opacity_hure = 10.0d0**(get_opacity_hure)
  
!~   write(*,*) get_opacity_hure
  
  ! we change the opacity from physical units to numerical units
  get_opacity_hure = phys_to_num_opacity * get_opacity_hure
  
 
  end function get_opacity_hure

end module opacity_hure
