module bessel

!*************************************************************
!** Modules that allow to compute bessel functions. Routines
!** are from numerical recipes
!** Version 1.1 - december 2011
!*************************************************************
  use types_numeriques
  use physical_constant
  
  implicit none
  
  private
  
  public :: bessik
  
  contains
  
	SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
	IMPLICIT NONE
	REAL(double_precision), INTENT(IN) :: x,xnu
	REAL(double_precision), INTENT(OUT) :: ri,rk,rip,rkp
	INTEGER, PARAMETER :: MAXIT=10000
	REAL(double_precision), PARAMETER :: XMIN=2.0
	REAL(double_precision), PARAMETER :: EPS=1.0d-10,FPMIN=1.0d-30
	INTEGER :: i,l,nl
	REAL(double_precision) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,&
		gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,&
		ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,rktemp,&
		s,sum,sum1,x2,xi,xi2,xmu,xmu2
	if (x<=0.0) then
		write (*,*) "in subroutine 'bessik' of module 'bessel', 'x' must be >0.0" 
		stop
	end if
	if (xnu<0.0) then
		write (*,*) "in subroutine 'bessik' of module 'bessel', 'xnu' must be >=0.0" 
		stop
	end if
!~ 	call assert(x > 0.0, xnu >= 0.0, 'bessik args')
	nl=int(xnu+0.5d0)
	xmu=xnu-nl
	xmu2=xmu*xmu
	xi=1.0d0/x
	xi2=2.0d0*xi
	h=xnu*xi
	if (h < FPMIN) h=FPMIN
	b=xi2*xnu
	d=0.0
	c=h
	do i=1,MAXIT
		b=b+xi2
		d=1.0d0/(b+d)
		c=b+1.0d0/c
		del=c*d
		h=del*h
		if (abs(del-1.0d0) < EPS) exit
	end do
	if (i > MAXIT) write(*,*) 'x too large in bessik; try asymptotic expansion'
	ril=FPMIN
	ripl=h*ril
	ril1=ril
	rip1=ripl
	fact=xnu*xi
	do l=nl,1,-1
		ritemp=fact*ril+ripl
		fact=fact-xi
		ripl=fact*ritemp+ril
		ril=ritemp
	end do
	f=ripl/ril
	if (x < XMIN) then
		x2=0.5d0*x
		pimu=PI*xmu
		if (abs(pimu) < EPS) then
			fact=1.0
		else
			fact=pimu/sin(pimu)
		end if
		d=-log(x2)
		e=xmu*d
		if (abs(e) < EPS) then
			fact2=1.0
		else
			fact2=sinh(e)/e
		end if
		call beschb_s(xmu,gam1,gam2,gampl,gammi)
		ff=fact*(gam1*cosh(e)+gam2*fact2*d)
		sum=ff
		e=exp(e)
		p=0.5d0*e/gampl
		q=0.5d0/(e*gammi)
		c=1.0
		d=x2*x2
		sum1=p
		do i=1,MAXIT
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*ff
			sum=sum+del
			del1=c*(p-i*ff)
			sum1=sum1+del1
			if (abs(del) < abs(sum)*EPS) exit
		end do
		if (i > MAXIT) write(*,*) 'bessk series failed to converge'
		rkmu=sum
		rk1=sum1*xi2
	else
		b=2.0d0*(1.0d0+x)
		d=1.0d0/b
		delh=d
		h=delh
		q1=0.0
		q2=1.0
		a1=0.25d0-xmu2
		c=a1
		q=c
		a=-a1
		s=1.0d0+q*delh
		do i=2,MAXIT
			a=a-2*(i-1)
			c=-a*c/i
			qnew=(q1-b*q2)/a
			q1=q2
			q2=qnew
			q=q+c*qnew
			b=b+2.0d0
			d=1.0d0/(b+a*d)
			delh=(b*d-1.0d0)*delh
			h=h+delh
			dels=q*delh
			s=s+dels
			if (abs(dels/s) < EPS) exit
		end do
		if (i > MAXIT) write(*,*) 'bessik: failure to converge in cf2'
		h=a1*h
		rkmu=sqrt(PI/(2.0d0*x))*exp(-x)/s
		rk1=rkmu*(xmu+x+0.5d0-h)*xi
	end if
	rkmup=xmu*xi*rkmu-rk1
	rimu=xi/(f*rkmu-rkmup)
	ri=(rimu*ril1)/ril
	rip=(rimu*rip1)/ril
	do i=1,nl
		rktemp=(xmu+i)*xi2*rk1+rkmu
		rkmu=rk1
		rk1=rktemp
	end do
	rk=rkmu
	rkp=xnu*xi*rkmu-rk1
	END SUBROUTINE bessik

	SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
	IMPLICIT NONE
	REAL(double_precision), INTENT(IN) :: x
	REAL(double_precision), INTENT(OUT) :: gam1,gam2,gampl,gammi
	INTEGER, PARAMETER :: NUSE1=5,NUSE2=5
	REAL(simple_precision) :: xx
	REAL(simple_precision), DIMENSION(7) :: c1=(/-1.142022680371168,&
		6.5165112670737e-3,3.087090173086e-4,-3.4706269649e-6,&
		6.9437664e-9,3.67795e-11,-1.356e-13/)
	REAL(simple_precision), DIMENSION(8) :: c2=(/1.843740587300905,&
		-7.68528408447867e-2,1.2719271366546e-3,&
		-4.9717367042e-6, -3.31261198e-8,2.423096e-10,&
		-1.702e-13,-1.49e-15/)
	xx=8.0d0*x*x-1.0d0
	gam1=chebev_s(-1.0,1.0,c1(1:NUSE1),xx)
	gam2=chebev_s(-1.0,1.0,c2(1:NUSE2),xx)
	gampl=gam2-x*gam1
	gammi=gam2+x*gam1
	END SUBROUTINE beschb_s

	FUNCTION chebev_s(a,b,c,x)
	IMPLICIT NONE
	REAL(simple_precision), INTENT(IN) :: a,b,x
	REAL(simple_precision), DIMENSION(:), INTENT(IN) :: c
	REAL(simple_precision) :: chebev_s
	INTEGER :: j,m
	REAL(simple_precision) :: d,dd,sv,y,y2
	if ((x-a)*(x-b) > 0.0) write(*,*) 'x not in range in chebev_s'
	m=size(c)
	d=0.0
	dd=0.0
	y=(2.0*x-a-b)/(b-a)
	y2=2.0*y
	do j=m,2,-1
		sv=d
		d=y2*d-dd+c(j)
		dd=sv
	end do
	chebev_s=y*d-dd+0.5*c(1)
	END FUNCTION chebev_s
end module bessel
