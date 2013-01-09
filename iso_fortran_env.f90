module iso_fortran_env

  ! Nonintrinsic version for Lahey/Fujitsu Fortran for Linux. 
  ! See Subclause 13.8.2 of the Fortran 2003 standard. 

  implicit NONE 
  public 

  integer, parameter :: Character_Storage_Size = 8 
  integer, parameter :: Error_Unit = 0 
  integer, parameter :: File_Storage_Size = 8 
  integer, parameter :: Input_Unit = 5 
  integer, parameter :: IOSTAT_END = -1 
  integer, parameter :: IOSTAT_EOR = -2 
  integer, parameter :: Numeric_Storage_Size = 32 
  integer, parameter :: Output_Unit = 6 

end module iso_fortran_env