module tides_constant_GR
  use types_numeriques

  implicit none
  !
  ! Author: Emeline Bolmont
  ! Date: 04/04/13
  !
  ! Number of tidally evolving planets
  integer, parameter :: ntid=2
  ! Nature of host body
  integer, parameter :: brown_dwarf=1
  integer, parameter :: M_dwarf=0
  integer, parameter :: Sun_like_star=0
  ! For an utilization of the code with no changing host body
  integer, parameter :: Rscst=0 !=1 : Rs = cst, rg2s = cst
  real(double_precision), parameter :: Rjup = 10.9d0 !rearth	
  
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
  ! Planets obliquities in rad
  real(double_precision), parameter, dimension(ntid) :: oblp = (/0.0d0,0.0d0/)				
  
  ! Indicate if Planet is of known parameters.
  ! 0: Earth-like, 1: Terrestrial (no mass-radius relationship), 2: Gas giant
  ! 3: others (like Neptune)
  integer, parameter, dimension(ntid) :: jupiter = (/0,0/)
  ! If jupiter ne 0, then indicate radius in Rearth, for ex: 1 or 0.954d0*Rjup
  real(double_precision), parameter, dimension(ntid) :: radius_p = (/0.d0,0.d0/)
  ! Radius of gyration, love number and k2delta_t for other planets (jupiter=3)
  real(double_precision), parameter :: rg2p_what = 3.308d-1
  real(double_precision), parameter :: k2p_what = 0.305d0
  real(double_precision), parameter :: k2pdeltap_what = 2.465278d-3
  
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
 
  
  !*********************************************************************
  !*********************************************************************
  ! No Need to chang stuff from here
  ! Radius of gyration and love number for dM 
  real(double_precision), parameter :: rg2_dM = 2.0d-1
  real(double_precision), parameter :: k2st_dM = 0.307d0 
  ! Radius of gyration and love number for Suns
  real(double_precision), parameter :: rg2_Sun = 5.9d-2
  real(double_precision), parameter :: k2st_Sun = 0.03d0 
  ! Radius of gyration, love number and k2delta_t for terrestrial planets
  real(double_precision), parameter :: rg2p_terr = 3.308d-1
  real(double_precision), parameter :: k2p_terr = 0.305d0
  real(double_precision), parameter :: k2pdeltap_terr = 2.465278d-3
  ! Radius of gyration, love number and k2delta_t for gas giants
  real(double_precision), parameter :: rg2p_gg = 2.54d-1
  real(double_precision), parameter :: k2p_gg = 0.38d0
  real(double_precision), parameter :: k2pdeltap_gg = 8.101852d-9
 
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
  ! meter in AU
  real(double_precision), parameter :: minau = 6.68458d-12
  ! Speed of light
  real(double_precision), parameter :: C2 = 1.731444830225d2


end module tides_constant_GR
