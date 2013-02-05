module turbulence

!*************************************************************
!** Modules that contains routines useful for turbulences
!**
!** Version 1.0 - mai 2012
!*************************************************************

  use types_numeriques
  use physical_constant

  implicit none
  
  integer :: wavenumber_min = 1
  integer :: wavenumber_max = 96
  real(double_precision) :: lifetime_prefactor = 0.1
  
  ! We define a new type for the properties of the planet
  type TurbulenceMode
    ! Properties of the mode
    real(double_precision) :: r ! the radial position of the mode [AU]
    real(double_precision) :: phi ! the azimuthal position of the mode [radians]
    integer :: wavenumber ! The mode number associated to the mode
    real(double_precision) :: t_init ! The time at which the mode was initialised [days]
    real(double_precision) :: lifetime ! the duration of the mode since its initialisation [days]
    real(double_precision) :: chi ! I don't know exactly what it is. Seems like a prefactor [dimensionless]
    real(double_precision) :: radial_extent ! The radial extent of the mode [AU]

  end type TurbulenceMode

  contains


subroutine init_mode(time, mode)
	implicit none
	real(double_precision), intent(in) :: time
	type(TurbulenceMode), intent(out) :: mode
	
	! locals
	real(double_precision) :: random_number_01, random_number_02, random_number_03, random_number_normal
	
	!------------------------------------------------------------------------------
	call random_number(random_number_01)! par defaut, c'est entre 0 et 1
	call random_number(random_number_02)
	call random_number(random_number_03)
	
	random_number_normal = dble(wavenumber_min) * exp(log(dble(wavenumber_max) / dble(wavenumber_min)) * random_number_01) ! pour le mode m, on passe d'une distrib uniforme a une  distrib normale, modemax vaut 96 et modemin vaut 1
	mode%wavenumber = dble(int(random_number_normal))! on force le fait qeu ce soit un entier

	mode%r = rmax + (rmax - rmin) * random_number_02 ! c'est le rk du papier, c'est la position radiale de l'origine de la turbulence

	mode%phi = 2.d0 * pi * random_number_03 ! pareil pour phik
	call normal(mode%chi) ! ca genere xi_k, distribution normale de moyenne nulle et d'ecart type 1
	
	mode%radial_extent = pi * mode%r / (4.0d0 * mode%wavenumber) ! c'est l'etendu du mode de la turbulence
	
	mode%t_init = time ! origine de la generation du mode 
	! The factor lifetime_prefactor is here to have a better agreement with mhd simulations. lifemode value has 0.1
	mode%lifetime = lifetime_prefactor * 2.d0 * pi * mode%r**(1.5d0) / (mode%wavenumber * hsr)!c'est le temps de vie du mode
	
end subroutine init_mode

subroutine potentielturb
  implicit none
  integer i,j,k,m
  double precision rand,randl,rand2,rand3
  double precision tp,lturb
  double precision mode_center_r(k),chi(k),devturb(k),deltat(k),mode_init_time(k),mturb(k),mode_center_phi(k)
  double precision phim,phip
  character*2 lab


  ! initialisation
  if ((time.eq.0.0d0).AND.(initdisk.eq.0)) then
     do k=1,nb_modes
        mode_init_time(k)=0.0d0
        deltat(k)=0.0d0
     enddo
  endif

  ! le nombre de mode nmod=50
  do k=1,nb_modes
    if (abs(time-mode_init_time(k)).ge.deltat(k)) then!si le mode doit être détruit, on en créé un autre pour remplacer
			call random_number(rand)! par defaut, c'est entre 0 et 1
			call random_number(rand2)
			call random_number(rand3)
			randl=dble(modemin)*exp(log(dble(modemax)/dble(modemin))*rand) ! pour le mode m, on passe d'une distrib uniforme a une  distrib normale, modemax vaut 96 et modemin vaut 1
			mturb(k)=dble(int(randl))! on force le fait qeu ce soit un entier

			mode_center_r(k) = rmax + (rmax - rmin) * rand2 ! c'est le rk du papier, c'est la position radiale de l'origine de la turbulence

			mode_center_phi(k) = 2.d0 * pi * rand3 ! pareil pour phik
			call normal(chi(k)) ! ca genere xi_k, distribution normale de moyenne nulle et d'ecart type 1
			devturb(k)=pi*mode_center_r(k)/(4.0d0*mturb(k)) ! c'est l'etendu du mode de la turbulence
			mode_init_time(k)=time!origine de la generation du mode 
			! The factor lifemode is here to have a better agreement with mhd simulations. lifemode value has 0.1
			deltat(k)=lifemode * 2.d0 * pi * mode_center_r(k)**(1.5d0) / (mturb(k) * hsr)!c'est le temps de vie du mode
		endif
	enddo




potturb=0.0d0


! la on calcule le potentiel turbulent, on fait la smme de tous les
! modes en tenant compte de leur evolution teporemlle
do k=1,nb_modes
	if (mturb(k).gt.modecut) then
		 lturb=0.0d0
	else
		 tp = time - mode_init_time(k)!c'est le temps relatif au tempsd 'origine du mode
		 lturb = forc * chi(k) / rc(i)*!mettre r2omega^2 ici, vu que je l'ai
		 TODO
		 !forc, c'est le terme gamma, prefacte'ur commun a tous les modes
		 &    exp(-(rc(i)-mode_center_r(k))**2/devturb(k)**2)
		 &   *cos(mturb(k)*tc(j)-mode_center_phi(k)-mode_center_r(k)**(-1.5d0)*tp
		 &    +omega*tp )
		 &   *sin(pi*tp/deltat(k))!expression 14 
		 potturb = potturb + lturb
		 !         write(11,*) mode_init_time(k),mode_center_r(k),mode_center_phi(k),chi(k),
		 !     &   devturb(k),mturb(k),deltat(k),modemin,modemax,lifemode
		 !         call flush(11)
	endif
enddo
!         stop
! TODO il faut que je mette ls expressions pour l'acceleration
! turbulente, aec deux termes, en cosinus et sinus. Parc que je derive
! directement dans le code pour avoir l'acceleration qui derive du
! potentiel turbulence. Til faut que je fasse aussi la conversion de
! cylindrique a cartesien, j'ai le formules sur le bloc.

end subroutine potentielturb



FUNCTION logdis(idum)
INTEGER idum
REAL logdis
REAL ran2,v1

v1=1.0+ran2(idum)
logdis=(v1*log(v1)-v1+1.0)/(2.0*log(2.0)-1.0)
return
END FUNCTION logdis



SUBROUTINE normal(x)
IMPLICIT NONE
double precision v1,v2,v11,v22
double precision r,x

r=1.5d0
do while ((r.gt.1.d0).OR.(r.eq.0.d0))
  call random_number(v1)
  call random_number(v2)
  v11 = 2.d0 * v1 - 1.d0
  v22 = 2.d0 * v2 - 1.d0
  r=v11**2 + v22**2
enddo
x=v11*dsqrt(-2.d0*dlog(r)/r)
END SUBROUTINE normal

FUNCTION gasdev(idum,mu)
INTEGER idum
REAL gasdev,mu
INTEGER iset
REAL fac,gset,rsq,v1,v2,ran2
SAVE iset,gset
DATA iset/0/
if (idum.lt.0) iset=0
if (iset.eq.0) then
1 v1=2.*ran2(idum)-1.
  v2=2.*ran2(idum)-1.
  rsq=v1**2+v2**2
  if (rsq.ge.1..or.rsq.eq.0.)goto 1
  fac=sqrt(-2.*log(rsq)/rsq)
  gset=v1*fac
  gasdev=v2*fac
  iset=1
else
  gasdev=gset
  iset=0
endif
gasdev=mu+gasdev
return
end FUNCTION gasdev


FUNCTION ran2(idum)
INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
REAL ran2,AM,eps,RNMX
PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
*IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
*IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER idum2,j,k,iv(NTAB),iy
SAVE iv,iy,idum2
DATA idum2/123456789/, iv/NTAB*0/, iy/0/
if (idum.le.0) then
  idum=max(-idum,1)
  idum2=idum
  do 11 j=NTAB+8,1,-1
     k=idum/IQ1
     idum=IA1*(idum-k*IQ1)-k*IR1
     if (idum.lt.0) idum=idum+IM1
     if (j.le.NTAB) iv(j)=idum
11   continue
     iy=iv(1)
  endif
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum.lt.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if (iy.lt.1) iy=iy+IMM1
  ran2=min(AM*iy,RNMX)
  return
END function ran2

end module turbulence