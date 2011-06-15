module mercury_outputs

!*************************************************************
!** Modules that write files, wheter they are outputs files, or errors
!**
!** Version 1.0 - june 2011
!*************************************************************
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_CE.FOR    (ErikSoft   1 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes details of close encounter minima to an output file, and decides how
! to continue the integration depending upon the close-encounter option
! chosen by the user. Close encounter details are stored until either 100
! have been accumulated, or a data dump is done, at which point the stored
! encounter details are also output.
!
! For each encounter, the routine outputs the time and distance of closest
! approach, the identities of the objects involved, and the output
! variables of the objects at this time. The output variables are:
! expressed as
!  r = the radial distance
!  theta = polar angle
!  phi = azimuthal angle
!  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and
!                             kinetic energies. (Note that 0 < fv < 1).
!  vtheta = polar angle of velocity vector
!  vphi = azimuthal angle of the velocity vector
!
!------------------------------------------------------------------------------
!
subroutine mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,nclo,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,&
     outfile,nstored,ceflush)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use ascii_conversion
  use orbital_elements
  use utilities

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,opt(8),stat(nbod),lmem(NMESS),stopflag
  integer :: nclo,iclo(nclo),jclo(nclo),nstored,ceflush
  real(double_precision) :: time,tstart,rcen,rmax,m(nbod),tclo(nclo),dclo(nclo)
  real(double_precision) :: ixvclo(6,nclo),jxvclo(6,nclo)
  character*80 outfile(3),mem(NMESS)
  character*8 id(nbod)
  !
  ! Local
  integer :: k,year,month
  real(double_precision) :: tmp0,t1,rfac,fr,fv,theta,phi,vtheta,vphi
  character*80 c(200)
  character*38 fstop
  character*6 tstring
  integer :: error
  !
  !------------------------------------------------------------------------------
  !
  save c
  !
  ! Scaling factor (maximum possible range) for distances
  rfac = log10 (rmax / rcen)
  !
  ! Store details of each new close-encounter minimum
  do k = 1, nclo
     nstored = nstored + 1
     c(nstored)(1:8)   = mio_fl2c(tclo(k))
     c(nstored)(9:16)  = mio_re2c(dble(iclo(k)-1),0.d0,11239423.99d0)
     c(nstored)(12:19) = mio_re2c(dble(jclo(k)-1),0.d0,11239423.99d0)
     c(nstored)(15:22) = mio_fl2c(dclo(k))
     !
     call mco_x2ov (rcen,rmax,m(1),0.d0,ixvclo(1,k),ixvclo(2,k),ixvclo(3,k),ixvclo(4,k),ixvclo(5,k),ixvclo(6,k),fr,theta,phi,fv,&
          vtheta,vphi)
     c(nstored)(23:30) = mio_re2c (fr    , 0.d0, rfac)
     c(nstored)(27:34) = mio_re2c (theta , 0.d0, PI)
     c(nstored)(31:38) = mio_re2c (phi   , 0.d0, TWOPI)
     c(nstored)(35:42) = mio_re2c (fv    , 0.d0, 1.d0)
     c(nstored)(39:46) = mio_re2c (vtheta, 0.d0, PI)
     c(nstored)(43:50) = mio_re2c (vphi  , 0.d0, TWOPI)
     !
     call mco_x2ov (rcen,rmax,m(1),0.d0,jxvclo(1,k),jxvclo(2,k),jxvclo(3,k),jxvclo(4,k),jxvclo(5,k),jxvclo(6,k),fr,theta,phi,fv,&
          vtheta,vphi)
     c(nstored)(47:54) = mio_re2c (fr    , 0.d0, rfac)
     c(nstored)(51:58) = mio_re2c (theta , 0.d0, PI)
     c(nstored)(55:62) = mio_re2c (phi   , 0.d0, TWOPI)
     c(nstored)(59:66) = mio_re2c (fv    , 0.d0, 1.d0)
     c(nstored)(63:74) = mio_re2c (vtheta, 0.d0, PI)
     c(nstored)(67:78) = mio_re2c (vphi  , 0.d0, TWOPI)
  end do
  !
  ! If required, output the stored close encounter details
  if (nstored.ge.100.or.ceflush.eq.0) then
     open (22, file=outfile(2), status='old', access='append',iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile(2))
        stop
     end if
     do k = 1, nstored
        write (22,'(a1,a2,a70)') char(12),'6b',c(k)(1:70)
     end do
     close (22)
     nstored = 0
  end if
  !
  ! If new encounter minima have occurred, decide whether to stop integration
  stopflag = 0
  if (opt(1).eq.1.and.nclo.gt.0) then
     open (23, file=outfile(3), status='old', access='append',iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile(3))
        stop
     end if
     ! If time style is Gregorian date then...
     tmp0 = tclo(1)
     if (opt(3).eq.1) then
        fstop = '(5a,/,9x,a,i10,1x,i2,1x,f4.1)'
        call mio_jd2y (tmp0,year,month,t1)
        write (23,fstop) mem(121)(1:lmem(121)),mem(126)(1:lmem(126)),id(iclo(1)),',',id(jclo(1)),mem(71)(1:lmem(71)),year,month,t1
        ! Otherwise...
     else
        if (opt(3).eq.3) then
           tstring = mem(2)
           fstop = '(5a,/,9x,a,f14.3,a)'
           t1 = (tmp0 - tstart) / 365.25d0
        else
           tstring = mem(1)
           fstop = '(5a,/,9x,a,f14.1,a)'
           if (opt(3).eq.0) t1 = tmp0
           if (opt(3).eq.2) t1 = tmp0 - tstart
        end if
        write (23,fstop) mem(121)(1:lmem(121)),mem(126)(1:lmem(126)),id(iclo(1)),',',id(jclo(1)),mem(71)(1:lmem(71)),t1,tstring
     end if
     stopflag = 1
     close(23)
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mio_ce
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_DUMP.FOR    (ErikSoft   21 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes masses, coordinates, velocities etc. of all objects, and integration
! parameters, to dump files. Also updates a restart file containing other
! variables used internally by MERCURY.
!
!------------------------------------------------------------------------------
!
subroutine mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,x,v,s,rho,rceh,stat,id,&
     ngf,epoch,opt,opflag,dumpfile,mem,lmem)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: algor,nbod,nbig,stat(nbod),opt(8),opflag,ndump,nfun
  integer :: lmem(NMESS)
  real(double_precision) :: time,tstart,tstop,dtout,h0,tol,rmax,en(3),am(3)
  real(double_precision) :: jcen(3),rcen,cefac,m(nbod),x(3,nbod),v(3,nbod)
  real(double_precision) :: s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod),epoch(nbod)
  character*80 dumpfile(4),mem(NMESS)
  character*8 id(nbod)
  !
  ! Local
  integer :: idp,i,j,k,len,j1,j2
  real(double_precision) :: rhocgs,k_2,rcen_2,rcen_4,rcen_6,x0(3,NMAX),v0(3,NMAX)
  character*150 c
  integer :: error
  !
  !------------------------------------------------------------------------------
  !
  rhocgs = AU * AU * AU * K2 / MSUN
  k_2 = 1.d0 / K2
  rcen_2 = 1.d0 / (rcen * rcen)
  rcen_4 = rcen_2 * rcen_2
  rcen_6 = rcen_4 * rcen_2
  !
  ! If using close-binary star, convert to user coordinates
  !      if (algor.eq.11) call mco_h2ub (time,jcen,nbod,nbig,h0,m,x,v,
  !     %   x0,v0)
  !
  ! Dump to temporary files (idp=1) and real dump files (idp=2)
  do idp = 1, 2
     !
     ! Dump data for the Big (i=1) and Small (i=2) bodies
     do i = 1, 2
        if (idp.eq.1) then
           if (i.eq.1) c(1:12) = 'big.tmp     '
           if (i.eq.2) c(1:12) = 'small.tmp   '
           open (31, file=c(1:12), status='unknown', iostat=error)
           if (error /= 0) then
              write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(c(1:12))
              stop
           end if
        else
           open (31, file=dumpfile(i), status='old', iostat=error)
           if (error /= 0) then
              write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(dumpfile(i))
              stop
           end if
        end if
        !
        ! Write header lines, data style (and epoch for Big bodies)
        write (31,'(a)') mem(151+i)(1:lmem(151+i))
        if (i.eq.1) then
           j1 = 2
           j2 = nbig
        else
           j1 = nbig + 1
           j2 = nbod
        end if
        write (31,'(a)') mem(154)(1:lmem(154))
        write (31,'(a)') mem(155)(1:lmem(155))
        write (31,*) mem(156)(1:lmem(156)),'Cartesian'
        if (i.eq.1) write (31,*) mem(157)(1:lmem(157)),time
        write (31,'(a)') mem(155)(1:lmem(155))
        !
        ! For each body...
        do j = j1, j2
           len = 37
           c(1:8) = id(j)
           write (c(9:37),'(1p,a3,e11.5,a3,e11.5)') ' r=',rceh(j),' d=',rho(j)/rhocgs
           if (m(j).gt.0) then
              write (c(len+1:len+25),'(a3,e22.15)') ' m=',m(j)*k_2
              len = len + 25
           end if
           do k = 1, 3
              if (ngf(k,j).ne.0) then
                 write (c(len+1:len+16),'(a2,i1,a1,e12.5)') ' a',k,'=',ngf(k,j)
                 len = len + 16
              end if
           end do
           if (ngf(4,j).ne.0) then
              write (c(len+1:len+15),'(a3,e12.5)') ' b=',ngf(4,j)
              len = len + 15
           end if
           write (31,'(a)') c(1:len)
           if (algor.eq.11) then
              write (31,312) x0(1,j), x0(2,j), x0(3,j)
              write (31,312) v0(1,j), v0(2,j), v0(3,j)
           else
              write (31,312) x(1,j), x(2,j), x(3,j)
              write (31,312) v(1,j), v(2,j), v(3,j)
           end if
           write (31,312) s(1,j)*k_2, s(2,j)*k_2, s(3,j)*k_2
        enddo
        close (31)
     end do
     !
     ! Dump the integration parameters
40   if (idp.eq.1) open (33,file='param.tmp',status='unknown',err=40)
45   if (idp.eq.2) open (33, file=dumpfile(3), status='old', err=45)
     !
     ! Important parameters
     write (33,'(a)') mem(151)(1:lmem(151))
     write (33,'(a)') mem(154)(1:lmem(154))
     write (33,'(a)') mem(155)(1:lmem(155))
     write (33,'(a)') mem(158)(1:lmem(158))
     write (33,'(a)') mem(155)(1:lmem(155))
     if (algor.eq.1) then
        write (33,*) mem(159)(1:lmem(159)),'MVS'
     else if (algor.eq.2) then
        write (33,*) mem(159)(1:lmem(159)),'BS'
     else if (algor.eq.3) then
        write (33,*) mem(159)(1:lmem(159)),'BS2'
     else if (algor.eq.4) then
        write (33,*) mem(159)(1:lmem(159)),'RADAU'
     else if (algor.eq.10) then
        write (33,*) mem(159)(1:lmem(159)),'HYBRID'
     else if (algor.eq.11) then
        write (33,*) mem(159)(1:lmem(159)),'CLOSE'
     else if (algor.eq.12) then
        write (33,*) mem(159)(1:lmem(159)),'WIDE'
     else
        write (33,*) mem(159)(1:lmem(159)),'0'
     end if
     write (33,*) mem(160)(1:lmem(160)),tstart
     write (33,*) mem(161)(1:lmem(161)),tstop
     write (33,*) mem(162)(1:lmem(162)),dtout
     write (33,*) mem(163)(1:lmem(163)),h0
     write (33,*) mem(164)(1:lmem(164)),tol
     !
     ! Integration options
     write (33,'(a)') mem(155)(1:lmem(155))
     write (33,'(a)') mem(165)(1:lmem(165))
     write (33,'(a)') mem(155)(1:lmem(155))
     if (opt(1).eq.0) then
        write (33,'(2a)') mem(166)(1:lmem(166)),mem(5)(1:lmem(5))
     else
        write (33,'(2a)') mem(166)(1:lmem(166)),mem(6)(1:lmem(6))
     end if
     if (opt(2).eq.0) then
        write (33,'(2a)') mem(167)(1:lmem(167)),mem(5)(1:lmem(5))
        write (33,'(2a)') mem(168)(1:lmem(168)),mem(5)(1:lmem(5))
     else if (opt(2).eq.2) then
        write (33,'(2a)') mem(167)(1:lmem(167)),mem(6)(1:lmem(6))
        write (33,'(2a)') mem(168)(1:lmem(168)),mem(6)(1:lmem(6))
     else
        write (33,'(2a)') mem(167)(1:lmem(167)),mem(6)(1:lmem(6))
        write (33,'(2a)') mem(168)(1:lmem(168)),mem(5)(1:lmem(5))
     end if
     if (opt(3).eq.0.or.opt(3).eq.2) then
        write (33,'(2a)') mem(169)(1:lmem(169)),mem(1)(1:lmem(1))
     else
        write (33,'(2a)') mem(169)(1:lmem(169)),mem(2)(1:lmem(2))
     end if
     if (opt(3).eq.2.or.opt(3).eq.3) then
        write (33,'(2a)') mem(170)(1:lmem(170)),mem(6)(1:lmem(6))
     else
        write (33,'(2a)') mem(170)(1:lmem(170)),mem(5)(1:lmem(5))
     end if
     if (opt(4).eq.1) then
        write (33,'(2a)') mem(171)(1:lmem(171)),mem(7)(1:lmem(7))
     else if (opt(4).eq.3) then
        write (33,'(2a)') mem(171)(1:lmem(171)),mem(9)(1:lmem(9))
     else
        write (33,'(2a)') mem(171)(1:lmem(171)),mem(8)(1:lmem(8))
     end if
     write (33,'(a)') mem(172)(1:lmem(172))
     if (opt(7).eq.1) then
        write (33,'(2a)') mem(173)(1:lmem(173)),mem(6)(1:lmem(6))
     else
        write (33,'(2a)') mem(173)(1:lmem(173)),mem(5)(1:lmem(5))
     end if
     if (opt(8).eq.1) then
        write (33,'(2a)') mem(174)(1:lmem(174)),mem(6)(1:lmem(6))
     else
        write (33,'(2a)') mem(174)(1:lmem(174)),mem(5)(1:lmem(5))
     end if
     !
     ! Infrequently-changed parameters
     write (33,'(a)') mem(155)(1:lmem(155))
     write (33,'(a)') mem(175)(1:lmem(175))
     write (33,'(a)') mem(155)(1:lmem(155))
     write (33,*) mem(176)(1:lmem(176)),rmax
     write (33,*) mem(177)(1:lmem(177)),rcen
     write (33,*) mem(178)(1:lmem(178)),m(1) * k_2
     write (33,*) mem(179)(1:lmem(179)),jcen(1) * rcen_2
     write (33,*) mem(180)(1:lmem(180)),jcen(2) * rcen_4
     write (33,*) mem(181)(1:lmem(181)),jcen(3) * rcen_6
     write (33,*) mem(182)(1:lmem(182))
     write (33,*) mem(183)(1:lmem(183))
     write (33,*) mem(184)(1:lmem(184)),cefac
     write (33,*) mem(185)(1:lmem(185)),ndump
     write (33,*) mem(186)(1:lmem(186)),nfun
     close (33)
     !
     ! Create new version of the restart file
60   if (idp.eq.1) open (35, file='restart.tmp', status='unknown',err=60)
65   if (idp.eq.2) open (35, file=dumpfile(4), status='old', err=65)
     write (35,'(1x,i2)') opflag
     write (35,*) en(1) * k_2
     write (35,*) am(1) * k_2
     write (35,*) en(3) * k_2
     write (35,*) am(3) * k_2
     write (35,*) s(1,1) * k_2
     write (35,*) s(2,1) * k_2
     write (35,*) s(3,1) * k_2
     close (35)
  end do
  !
  !------------------------------------------------------------------------------
  !
311 format (1x,a8,1x,a1,1p,e22.15,2(1x,e11.5))
312 format (1p,3(1x,e22.15),1x,i8)
313 format (1p,1x,e22.15,0p,2x,a)
314 format (1x,a8,1x,a1,1p,e22.15,4(1x,e12.5),1x,e22.15,2(1x,e11.5))
  return
end subroutine mio_dump
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_ERR.FOR    (ErikSoft  6 December 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes out an error message and terminates Mercury.
!
!------------------------------------------------------------------------------
!
subroutine mio_err (unit,s1,ls1,s2,ls2,s3,ls3,s4,ls4)

  implicit none

  !
  ! Input/Output
  integer :: unit,ls1,ls2,ls3,ls4
  character*80 s1,s2
  character(len=*) :: s3,s4
  !
  !------------------------------------------------------------------------------
  !
  write (*,'(/,2a)') ' ERROR: Programme terminated. See information',' file for details.'
  !
  write (unit,'(/,3a,/,2a)') s1(1:ls1),s2(1:ls2),s3(1:ls3),' ',s4(1:ls4)
  stop
  !
  !------------------------------------------------------------------------------
  !
end subroutine mio_err

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_LOG.FOR    (ErikSoft   25 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes a progress report to the log file (or the screen if you are running
! Mercury interactively).
!
!------------------------------------------------------------------------------
!
subroutine mio_log (time,tstart,en,am,opt,mem,lmem)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use utilities

  implicit none

  !
  ! Input/Output
  integer :: lmem(NMESS), opt(8)
  real(double_precision) :: time, tstart, en(3), am(3)
  character*80 mem(NMESS)
  !
  ! Local
  integer :: year, month
  real(double_precision) :: tmp0, tmp1, t1
  character*38 flog
  character*6 tstring
  !
  !------------------------------------------------------------------------------
  !
  if (opt(3).eq.0.or.opt(3).eq.2) then
     tstring = mem(1)
     flog = '(1x,a,f14.1,a,2(a,1p1e12.5))'
  else if (opt(3).eq.1) then
     flog = '(1x,a,i10,1x,i2,1x,f4.1,2(a,1p1e12.5))'
  else
     tstring = mem(2)
     flog = '(1x,a,f14.3,a,2(a,1p1e12.5))'
  end if
  !
  tmp0 = 0.d0
  tmp1 = 0.d0
  if (en(1).ne.0) tmp0 = (en(2) + en(3) - en(1)) / abs(en(1))
  if (am(1).ne.0) tmp1 = (am(2) + am(3) - am(1)) / abs(am(1))
  !
  if (opt(3).eq.1) then
     call mio_jd2y (time,year,month,t1)
     write (*,flog) mem(64)(1:lmem(64)), year, month, t1,mem(65)(1:lmem(65)), tmp0,mem(66)(1:lmem(66)), tmp1
  else
     if (opt(3).eq.0) t1 = time
     if (opt(3).eq.2) t1 = time - tstart
     if (opt(3).eq.3) t1 = (time - tstart) / 365.25d0
     write (*,flog) mem(63)(1:lmem(63)), t1, tstring,mem(65)(1:lmem(65)), tmp0, mem(66)(1:lmem(66)), tmp1
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mio_log
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_OUT.FOR    (ErikSoft   13 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes output variables for each object to an output file. Each variable
! is scaled between the minimum and maximum possible values and then
! written in a compressed format using ASCII characters.
! The output variables are:
!  r = the radial distance
!  theta = polar angle
!  phi = azimuthal angle
!  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and
!                             kinetic energies. (Note that 0 < fv < 1).
!  vtheta = polar angle of velocity vector
!  vphi = azimuthal angle of the velocity vector
!
! If this is the first output (OPFLAG = -1), or the first output since the 
! number of the objects or their masses have changed (OPFLAG = 1), then 
! the names, masses and spin components of all the objects are also output.
!
! N.B. Each object's distance must lie between RCEN < R < RMAX
! ===  
!
!------------------------------------------------------------------------------
!
subroutine mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho,stat,id,opt,opflag,algor,outfile)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use ascii_conversion
  use orbital_elements

  implicit none

  !
  ! Input/Output
  integer :: nbod, nbig, stat(nbod), opt(8), opflag, algor
  real(double_precision) :: time,jcen(3),rcen,rmax,m(nbod),xh(3,nbod),vh(3,nbod)
  real(double_precision) :: s(3,nbod),rho(nbod)
  character*80 outfile
  character*8 id(nbod)
  !
  ! Local
  integer :: k, len, nchar
  real(double_precision) :: rhocgs,k_2,rfac,rcen_2,fr,fv,theta,phi,vtheta,vphi
  character*80 header,c(NMAX)
  character*5 fout
  integer :: error
  !
  !------------------------------------------------------------------------------
  !
  rhocgs = AU * AU * AU * K2 / MSUN
  k_2 = 1.d0 / K2
  rcen_2 = 1.d0 / (rcen * rcen)
  !
  ! Scaling factor (maximum possible range) for distances
  rfac = log10 (rmax / rcen)
  !
  ! Create the format list, FOUT, used when outputting the orbital elements
  if (opt(4).eq.1) nchar = 2
  if (opt(4).eq.2) nchar = 4
  if (opt(4).eq.3) nchar = 7
  len = 3  +  6 * nchar
  fout(1:5) = '(a  )'
  if (len.lt.10) write (fout(3:3),'(i1)') len
  if (len.ge.10) write (fout(3:4),'(i2)') len
  !
  ! Open the orbital-elements output file
  open (21, file=outfile, status='old', access='append', iostat=error)
  if (error /= 0) then
     write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile)
     stop
  end if
  !
  !------------------------------------------------------------------------------
  !
  !  SPECIAL  OUTPUT  PROCEDURE
  !
  ! If this is a new integration or a complete output is required (e.g. because
  ! the number of objects has changed), then output object details & parameters.
  if (opflag.eq.-1.or.opflag.eq.1) then
     !
     ! Compose a header line with time, number of objects and relevant parameters
     header(1:8)   = mio_fl2c (time)
     header(9:16)  = mio_re2c (dble(nbig - 1),   0.d0, 11239423.99d0)
     header(12:19) = mio_re2c (dble(nbod - nbig),0.d0, 11239423.99d0)
     header(15:22) = mio_fl2c (m(1) * k_2)
     header(23:30) = mio_fl2c (jcen(1) * rcen_2)
     header(31:38) = mio_fl2c (jcen(2) * rcen_2 * rcen_2)
     header(39:46) = mio_fl2c (jcen(3) * rcen_2 * rcen_2 * rcen_2)
     header(47:54) = mio_fl2c (rcen)
     header(55:62) = mio_fl2c (rmax)
     !
     ! For each object, compress its index number, name, mass, spin components
     ! and density (some of these need to be converted to normal units).
     do k = 2, nbod
        c(k)(1:8) = mio_re2c (dble(k - 1), 0.d0, 11239423.99d0)
        c(k)(4:11) = id(k)
        c(k)(12:19) = mio_fl2c (m(k) * k_2)
        c(k)(20:27) = mio_fl2c (s(1,k) * k_2)
        c(k)(28:35) = mio_fl2c (s(2,k) * k_2)
        c(k)(36:43) = mio_fl2c (s(3,k) * k_2)
        c(k)(44:51) = mio_fl2c (rho(k) / rhocgs)
     end do
     !
     ! Write compressed output to file
     write (21,'(a1,a2,i2,a62,i1)') char(12),'6a',algor,header(1:62),opt(4)
     do k = 2, nbod
        write (21,'(a51)') c(k)(1:51)
     end do
  end if
  !
  !------------------------------------------------------------------------------
  !
  !  NORMAL  OUTPUT  PROCEDURE
  !
  ! Compose a header line containing the time and number of objects
  header(1:8)   = mio_fl2c (time)
  header(9:16)  = mio_re2c (dble(nbig - 1),    0.d0, 11239423.99d0)
  header(12:19) = mio_re2c (dble(nbod - nbig), 0.d0, 11239423.99d0)
  !
  ! Calculate output variables for each body and convert to compressed format
  do k = 2, nbod
     call mco_x2ov (rcen,rmax,m(1),m(k),xh(1,k),xh(2,k),xh(3,k),vh(1,k),vh(2,k),vh(3,k),fr,theta,phi,fv,vtheta,vphi)
     !
     ! Object's index number and output variables
     c(k)(1:8) = mio_re2c (dble(k - 1), 0.d0, 11239423.99d0)
     c(k)(4:11)                 = mio_re2c (fr,     0.d0, rfac)
     c(k)(4+  nchar:11+  nchar) = mio_re2c (theta,  0.d0, PI)
     c(k)(4+2*nchar:11+2*nchar) = mio_re2c (phi,    0.d0, TWOPI)
     c(k)(4+3*nchar:11+3*nchar) = mio_re2c (fv,     0.d0, 1.d0)
     c(k)(4+4*nchar:11+4*nchar) = mio_re2c (vtheta, 0.d0, PI)
     c(k)(4+5*nchar:11+5*nchar) = mio_re2c (vphi,   0.d0, TWOPI)
  end do
  !
  ! Write compressed output to file
  write (21,'(a1,a2,a14)') char(12),'6b',header(1:14)
  do k = 2, nbod
     write (21,fout) c(k)(1:len)
  end do
  !
  close (21)
  opflag = 0
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mio_out

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_AEI.FOR    (ErikSoft   31 January 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Creates a filename and opens a file to store aei information for an object.
! The filename is based on the name of the object.
!
!------------------------------------------------------------------------------
!
subroutine mio_clo (id,unitnum,header,lenhead,mem,lmem)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use utilities, only : mio_spl

  implicit none

  !
  ! Input/Output
  integer :: unitnum,lenhead,lmem(NMESS)
  character(len=4) :: extn = ".clo"
  character*8 id
  character*250 header
  character*80 mem(NMESS)
  !
  ! Local
  integer :: j,k,itmp,nsub,lim(2,4)
  logical test
  character*1 bad(5)
  character*250 filename
  !
  !------------------------------------------------------------------------------
  !
  data bad/ '*', '/', '.', ':', '&'/
  !
  ! Create a filename based on the object's name
  call mio_spl (8,id,nsub,lim)
  itmp = min(7,lim(2,1)-lim(1,1))
  filename(1:itmp+1) = id(1:itmp+1)
  filename(itmp+2:itmp+5) = extn
  do j = itmp + 6, 250
     filename(j:j) = ' '
  end do
  !
  ! Check for inappropriate characters in the filename
  do j = 1, itmp + 1
     do k = 1, 5
        if (filename(j:j).eq.bad(k)) filename(j:j) = '_'
     end do
  end do
  !
  ! If the file exists already, give a warning and don't overwrite it
  inquire (file=filename, exist=test)
  if (test) then
     write (*,'(/,3a)') mem(121)(1:lmem(121)),mem(87)(1:lmem(87)),filename(1:80)
     unitnum = -1
  else
     open (unitnum, file=filename, status='new')
     write (unitnum, '(/,30x,a8,//,a)') id,header(1:lenhead)
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mio_clo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_AEI.FOR    (ErikSoft   31 January 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Creates a filename and opens a file to store aei information for an object.
! The filename is based on the name of the object.
!
!------------------------------------------------------------------------------
!
subroutine mio_aei (id,unitnum,header,lenhead,mem,lmem)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use utilities, only : mio_spl

  implicit none

  !
  ! Input/Output
  integer :: unitnum,lenhead,lmem(NMESS)
  character(len=4) :: extn = ".aei"
  character*8 id
  character*250 header
  character*80 mem(NMESS)
  !
  ! Local
  integer :: j,k,itmp,nsub,lim(2,4)
  logical test
  character*1 bad(5)
  character*250 filename
  !
  !------------------------------------------------------------------------------
  !
  data bad/ '*', '/', '.', ':', '&'/
  !
  ! Create a filename based on the object's name
  call mio_spl (8,id,nsub,lim)
  itmp = min(7,lim(2,1)-lim(1,1))
  filename(1:itmp+1) = id(1:itmp+1)
  filename(itmp+2:itmp+5) = extn
  do j = itmp + 6, 250
     filename(j:j) = ' '
  end do
  !
  ! Check for inappropriate characters in the filename
  do j = 1, itmp + 1
     do k = 1, 5
        if (filename(j:j).eq.bad(k)) filename(j:j) = '_'
     end do
  end do
  !
  ! If the file exists already, give a warning and don't overwrite it
  inquire (file=filename, exist=test)
  if (test) then
     write (*,'(/,3a)') mem(121)(1:lmem(121)),mem(87)(1:lmem(87)),filename(1:80)
     unitnum = -1
  else
     open (unitnum, file=filename, status='new')
     write (unitnum, '(/,30x,a8,//,a)') id,header(1:lenhead)
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mio_aei

end module mercury_outputs
