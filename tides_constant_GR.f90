module tides_constant_GR
  use types_numeriques

  implicit none

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  !      USER_TIDES.INC  -- for use with mercury6_tides.f
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! Author: Sean Raymond, Emeline Bolmont
  ! Date: 07/08/09
  !
  ! Parameters that govern tidal evolution (from Golreich & Soter
  !         1966 and Jackson et al 2008): 
  !
  !   NTID = # of tidally evolving planets, in order in BIG file 
  !                     (must be 1 or greater)
  !   SIGMAST = Tidal dissipation factor (Msun-1.AU-2.day-1)
  !   RSun = Stellar radius IN AU
  !   K2P = Love number of degree 2 of a rocky planet
  !   k2pdeltap = k2_p*delta_p for a rocky planet in day
  !   REARTH = Radius of Earth in AU
  
  !
  !  M2EARTH = conversion from solar -> Earth masses
  !  BIGG = Gravitational constant in AU^3.Msun-1.day-2
  !  CMAU = conversion from cm to AU
  !  AUCM =     "        "  AU to cm
  !  MEARTH = Earth mass in solar mass
  !  C2 = speed of light (in AU/day)
  !
  integer, parameter :: ntid=2
  integer, parameter :: tidflag=1
  integer, parameter :: ceflag=0
  integer, parameter :: brown_dwarf=1
  integer, parameter :: M_dwarf=0
  integer, parameter :: Sun_like_star=0
  ! For an utilization of the code with no changing host body
  integer, parameter :: Rscst=0 !=1 : Rs = cst, rg2s = cst
  
  
  ! Integration stuff
  real(double_precision), parameter :: t_init = 1.d6*365.25!1.0d6*365.25
  ! If crash, write last line of spin.dat here : 
  integer, parameter :: crash=0
  real(double_precision), parameter :: t_crash = 0.0d0*365.25!1.0d6*365.25
  real(double_precision), parameter, dimension(3) :: rot_crash  = (/0.0,0.0,0.0/)
  real(double_precision), parameter, dimension(3) :: rot_crashp1 = (/0.0,0.0,0.0/)
  real(double_precision), parameter, dimension(3) :: rot_crashp2 = (/0.0,0.0,0.0/)
  real(double_precision), parameter, dimension(3) :: rot_crashp3 = (/0.0,0.0,0.0/)


  ! Planet dissipation, and caracteristics
  
  ! If pseudo_rot eq 0 : initial period as given by Pp0 (in hr)
  ! If pseudo_rot eq toto : initial period = toto*pseudo_synchronization period 
  real(double_precision), parameter, dimension(ntid) :: pseudo_rot = (/1,1/)
  real(double_precision), parameter, dimension(ntid) :: Pp0 = (/24.d0, 24.d0/) 	
  real(double_precision), parameter, dimension(ntid) :: dissplan = (/1.d0,10.d0/)
  ! Love number of degree 2: Jup: 0.38d0, E: 0.305d0
  real(double_precision), parameter, dimension(ntid) :: k2p = (/0.305d0,0.305d0/)
  ! Earth, k2pdeltap = 213s = 2.465278d-3 day
  ! Jup,  k2pdeltap = 7d-4s = 8.101852d-9 day
  real(double_precision), parameter, dimension(ntid) :: k2pdeltap = (/2.465278d-3,2.465278d-3/) 
  ! Radius of gyration: Jup: 2.54d-1, E:3.308d-1
  real(double_precision), parameter, dimension(ntid) :: rg2p = (/3.308d-1,3.308d-1/)
  ! Planets obliquities in rad
  real(double_precision), parameter, dimension(ntid) :: oblp = (/0.0d0,0.0d0/)				
  real(double_precision), parameter :: Rjup = 10.9d0 !rearth	
  ! Indicate if Planet is of known parameters.
  ! 0: Earth-like, 1: Jupiter, 2: Saturn
  ! 3,4: 55 Cnc e and b
  ! 5: radius_p (Jupiter radius)
  integer, parameter, dimension(ntid) :: jupiter = (/0,0/)
  real(double_precision), parameter, dimension(ntid) :: radius_p = (/0.954d0,0.954d0/)
  
  
  ! Star dissipation, and caracteristics in CGS
  real(double_precision), parameter :: dissstar = 1.0d0!1.0d0!1.d2
  
  
  ! Dissipation factors of allowed host body
  
  ! For R=cst, choose sigmast:
  real(double_precision), parameter :: sigma_what = 2.006*3.845764d4 !-60+64
  real(double_precision), parameter :: rg2_what = 2.0d-1
  real(double_precision), parameter :: k2st_what = 0.307d0 
  ! For R=cst, or dM or Suns
  real(double_precision), parameter :: Period_st   = 8.0d0    !day
  real(double_precision), parameter :: radius_star = 0.943 !Rsun
  
  ! No Need to chang stuff from here
  ! Radius of gyration for dM 
  real(double_precision), parameter :: rg2_dM = 2.0d-1
  real(double_precision), parameter :: k2st_dM = 0.307d0 
  ! Radius of gyration for Suns
  real(double_precision), parameter :: rg2_Sun = 5.9d-2
  real(double_precision), parameter :: k2st_Sun = 0.03d0 
  
  ! Sigma for BD, dM, Suns
  ! BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
  real(double_precision), parameter :: sigma_BD = 2.006*3.845764d4 !-60+64
  real(double_precision), parameter :: sigma_dM = 2.006*3.845764d4 !-60+64
  ! Sun-like-star: sigmast = 4.992d-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
  real(double_precision), parameter :: sigma_Sun = 4.992*3.845764d-2 !-66+64
  ! If planet not terrestrial, dissipation factor Gas Giant
  real(double_precision), parameter :: sigma_gg = 2.006*3.845764d4
  
  
  ! Some stuff, constants mainly
  real(double_precision), parameter :: rsun = 4.67920694d-3
  real(double_precision), parameter :: rearth = 4.25874677d-5
  real(double_precision), parameter :: m2earth = (1.9891d6/5.9794)
  real(double_precision), parameter :: minau = 6.68458d-12
  real(double_precision), parameter :: C2 = 1.731444830225d2


end module tides_constant_GR
