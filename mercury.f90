!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MERCURY.F90    (ErikSoft   3 May 2002)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers (and Christophe Cossou)
!
! Mercury is a general-purpose N-body integration package for problems in
! celestial mechanics.
!
!------------------------------------------------------------------------------
! This package contains some subroutines taken from the Swift integration 
! package by H.F.Levison and M.J.Duncan (1994) Icarus, vol 108, pp18.
! Routines taken from Swift have names beginning with `drift' or `orbel'.
!
! The standard symplectic (MVS) algorithm is described in J.Widsom and
! M.Holman (1991) Astronomical Journal, vol 102, pp1528.
!
! The hybrid symplectic algorithm is described in J.E.Chambers (1999)
! Monthly Notices of the RAS, vol 304, pp793.
!
! RADAU is described in E.Everhart (1985) in ``The Dynamics of Comets:
! Their Origin and Evolution'' p185-202, eds. A.Carusi & G.B.Valsecchi,
! pub. Reidel.
!
! The Bulirsch-Stoer algorithms are described in W.H.Press et al. (1992)
! ``Numerical Recipes in Fortran'', pub. Cambridge.
!------------------------------------------------------------------------------
!
! Variables:
! ---------
!  M      = mass (in solar masses)
!  XH     = coordinates (x,y,z) with respect to the central body (AU)
!  VH     = velocities (vx,vy,vz) with respect to the central body (AU/day)
!  S      = spin angular momentum (solar masses AU^2/day)
!  RHO    = physical density (g/cm^3)
!  RCEH   = close-encounter limit (Hill radii)
!  STAT   = status (0 => alive, <>0 => to be removed)
!  ID     = name of the object (8 characters)
!  CE     = close encounter status
!  NGF    = (1-3) cometary non-gravitational (jet) force parameters
!   "     =  (4)  beta parameter for radiation pressure and P-R drag
!  EPOCH  = epoch of orbit (days)
!  NBOD  = current number of bodies (INCLUDING the central object)
!  NBIG  =    "       "    " big bodies (ones that perturb everything else)
!  TIME  = current epoch (days)
!  TOUT  = time of next output evaluation
!  TDUMP = time of next data dump
!  TFUN  = time of next periodic effect (e.g. next check for ejections)
!  H     = current integration timestep (days)
!  EN(1) = initial energy of the system
!  " (2) = current    "    "  "    "
!  " (3) = energy change due to collisions, ejections etc.
!  AM(1,2,3) = as above but for angular momentum
!
! Integration Parameters :
! ----------------------
!  ALGOR = 1  ->  Mixed-variable symplectic
!          2  ->  Bulirsch-Stoer integrator
!          3  ->         "           "      (conservative systems only)
!          4  ->  RA15 `radau' integrator
!          10 ->  Hybrid MVS/BS (democratic-heliocentric coords)
!          11 ->  Close-binary hybrid (close-binary coords)
!          12 ->  Wide-binary hybrid (wide-binary coords)
!
! TSTART = epoch of first required output (days)
! TSTOP  =   "      final required output ( "  )
! DTOUT  = data output interval           ( "  )
! DTDUMP = data-dump interval             ( "  )
! DTFUN  = interval for other periodic effects (e.g. check for ejections)
!  H0    = initial integration timestep (days)
!  TOL   = Integrator tolerance parameter (approx. error per timestep)
!  RMAX  = heliocentric distance at which objects are considered ejected (AU)
!  RCEN  = radius of central body (AU)
!  JCEN(1,2,3) = J2,J4,J6 for central body (units of RCEN^i for Ji)
!
! Options:
!  OPT(1) = close-encounter option (0=stop after an encounter, 1=continue)
!  OPT(2) = collision option (0=no collisions, 1=merge, 2=merge+fragment)
!  OPT(3) = time style (0=days 1=Greg.date 2/3=days/years w/respect to start)
!  OPT(4) = o/p precision (1,2,3 = 4,9,15 significant figures)
!  OPT(5) = < Not used at present >
!  OPT(6) = < Not used at present >
!  OPT(7) = apply post-Newtonian correction? (0=no, 1=yes)
!  OPT(8) = apply user-defined force routine mfo_user? (0=no, 1=yes)
!
! File variables :
! --------------
!  OUTFILE  (1) = osculating coordinates/velocities and masses
!     "     (2) = close encounter details
!     "     (3) = information file
!  DUMPFILE (1) = Big-body data
!     "     (2) = Small-body data
!     "     (3) = integration parameters
!     "     (4) = restart file
!
! Flags :
! -----
!  NGFLAG = do any bodies experience non-grav. forces?
!                            ( 0 = no non-grav forces)
!                              1 = cometary jets only
!                              2 = radiation pressure/P-R drag only
!                              3 = both
!  OPFLAG = integration mode (-2 = synchronising epochs)
!                             -1 = integrating towards start epoch
!                              0 = main integration, normal output
!                              1 = main integration, full output
!
!------------------------------------------------------------------------------
!

program mercury

  use physical_constant
  use mercury_constant
  use types_numeriques
  use system_properties
  use mercury_outputs
  use utilities

  implicit none

  !
  integer :: j,algor,nbod,nbig,opt(8),stat(NMAX),lmem(NMESS)
  integer :: opflag,ngflag,ndump,nfun
  integer :: error
  real(double_precision) :: m(NMAX),xh(3,NMAX),vh(3,NMAX),s(3,NMAX),rho(NMAX)
  real(double_precision) :: rceh(NMAX),epoch(NMAX),ngf(4,NMAX),rmax,rcen,jcen(3)
  real(double_precision) :: cefac,time,tstart,tstop,dtout,h0,tol,en(3),am(3)
  character*8 id(NMAX)
  character*80 outfile(3), dumpfile(4), mem(NMESS)
  external mdt_mvs, mdt_bs1, mdt_bs2, mdt_ra15, mdt_hy
  external mco_dh2h,mco_h2dh
  external mco_b2h,mco_h2b,mco_h2mvs,mco_mvs2h,mco_iden
  !
  data opt/0,1,1,2,0,1,0,0/
  !
  !------------------------------------------------------------------------------
  !
  ! Get initial conditions and integration parameters
  call mio_in (time,tstart,tstop,dtout,algor,h0,tol,rmax,rcen,jcen,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,id,&
       epoch,ngf,opt,opflag,ngflag,outfile,dumpfile,lmem,mem)
  !
  ! If this is a new integration, integrate all the objects to a common epoch.
  if (opflag.eq.-2) then
     open (23,file=outfile(3),status='old',access='append',iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile(3))
        stop
     end if
     write (23,'(/,a)') mem(55)(1:lmem(55))
     write (*,'(a)') mem(55)(1:lmem(55))
     call mxx_sync (time,tstart,h0,tol,jcen,nbod,nbig,m,xh,vh,s,rho,rceh,stat,id,epoch,ngf,opt,ngflag)
     write (23,'(/,a,/)') mem(56)(1:lmem(56))
     write (*,'(a)') mem(56)(1:lmem(56))
     opflag = -1
     close (23)
  end if
  !
  ! Main integration
  if (algor.eq.1) call mal_hcon (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,&
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,mdt_mvs,mco_h2mvs,mco_mvs2h)
  !
  if (algor.eq.9) call mal_hcon (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,&
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,mdt_mvs,mco_iden,mco_iden)
  !
  if (algor.eq.2) call mal_hvar (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,&
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,mdt_bs1)
  !
  if (algor.eq.3) call mal_hvar (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,&
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,mdt_bs2)
  !
  if (algor.eq.4) call mal_hvar (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,&
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,mdt_ra15)
  !
  if (algor.eq.10) call mal_hcon (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,&
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,mdt_hy,mco_h2dh,mco_dh2h)
  !
  ! Do a final data dump
  do j = 2, nbod
     epoch(j) = time
  end do
  call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,&
       id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
  !
  ! Calculate and record the overall change in energy and ang. momentum
  open  (23, file=outfile(3), status='old', access='append',iostat=error)
  if (error /= 0) then
     write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile(3))
     stop
  end if
  write (23,'(/,a)') mem(57)(1:lmem(57))
  call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
  !
  write (23,231) mem(58)(1:lmem(58)),abs((en(2) + en(3) - en(1)) / en(1))
  write (23,232) mem(59)(1:lmem(59)), abs((am(2) + am(3) - am(1)) / am(1))
  write (23,231) mem(60)(1:lmem(60)), abs(en(3) / en(1))
  write (23,232) mem(61)(1:lmem(61)), abs(am(3) / am(1))
  close (23)
  write (*,'(a)') mem(57)(1:lmem(57))
  !
  !------------------------------------------------------------------------------
  !
231 format (/,a,1p1e12.5)
232 format (a,1p1e12.5)
  stop
  
  contains
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_IN.FOR    (ErikSoft   4 May 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Reads names, masses, coordinates and velocities of all the bodies,
! and integration parameters for the MERCURY integrator package. 
! If DUMPFILE(4) exists, the routine assumes this is a continuation of
! an old integration, and reads all the data from the dump files instead
! of the input files.
!
! N.B. All coordinates are with respect to the central body!!
! ===
!
!------------------------------------------------------------------------------
!
subroutine mio_in (time,tstart,tstop,dtout,algor,h0,tol,rmax,rcen,jcen,en,am,cefac,ndump,nfun,nbod,nbig,m,x,v,s,rho,rceh,stat,&
     id,epoch,ngf,opt,opflag,ngflag,outfile,dumpfile,lmem,mem)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use orbital_elements

  implicit none

  !
  ! Input/Output
  integer :: algor,nbod,nbig,stat(NMAX),opt(8),opflag,ngflag
  integer :: lmem(NMESS),ndump,nfun
  real(double_precision) :: time,tstart,tstop,dtout,h0,tol,rmax,rcen,jcen(3)
  real(double_precision) :: en(3),am(3),m(NMAX),x(3,NMAX),v(3,NMAX),s(3,NMAX)
  real(double_precision) :: rho(NMAX),rceh(NMAX),epoch(NMAX),ngf(4,NMAX),cefac
  character*80 outfile(3),dumpfile(4), mem(NMESS)
  character*8 id(NMAX)
  !
  ! Local
  integer :: j,k,itmp,jtmp,informat,lim(2,10),nsub,year,month,lineno
  real(double_precision) :: q,a,e,i,p,n,l,temp,tmp2,tmp3,rhocgs,t1,tmp4,tmp5,tmp6
  !      real(double_precision) :: v0(3,NMAX),x0(3,NMAX)
  logical test,oldflag,flag1,flag2
  character*1 c1
  character*3 c3,alg(60)
  character*80 infile(3),filename,c80
  character*150 string
  integer :: error
  !
  !------------------------------------------------------------------------------
  !
  data alg/'MVS','Mvs','mvs','mvs','mvs','BS ','Bs ','bs ','Bul', 'bul','BS2','Bs2','bs2','Bu2','bu2',&
       'RAD','Rad','rad','RA ', 'ra ','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx', 'xxx',&
       'xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx', 'xxx','TES','Tes','tes','Tst','tst',&
       'HYB','Hyb','hyb','HY ', 'hy ','CLO','Clo','clo','CB ','cb ','WID','Wid','wid','WB ', 'wb '/
  !
  rhocgs = AU * AU * AU * K2 / MSUN
  do j = 1, 80
     filename(j:j) = ' '
  end do
  do j = 1, 3
     infile(j)   = filename
     outfile(j)  = filename
     dumpfile(j) = filename
  end do
  dumpfile(4) = filename
  !
  ! Read in output messages
  inquire (file='message.in', exist=test)
  if (.not.test) then
     write (*,'(/,2a)') ' ERROR: This file is needed to start',' the integration:  message.in'
     stop
  end if
  open (16, file='message.in', status='old')
  do
     read (16,'(i3,1x,i2,1x,a80)', iostat=error) j,lmem(j),mem(j)
     if (error /= 0) exit
  end do
  close (16)
  !
  ! Read in filenames and check for duplicate filenames
  inquire (file='files.in', exist=test)
  if (.not.test) call mio_err (6,mem(81),lmem(81),mem(88),lmem(88),' ',1,'files.in',8)
  open (15, file='files.in', status='old')
  !
  ! Input files
  do j = 1, 3
     read (15,'(a150)') string
     call mio_spl (150,string,nsub,lim)
     infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
     do k = 1, j - 1
        if (infile(j).eq.infile(k)) call mio_err (6,mem(81),lmem(81),mem(89),lmem(89),infile(j),80,mem(86),lmem(86))
     end do
  end do
  !
  ! Output files
  do j = 1, 3
     read (15,'(a150)') string
     call mio_spl (150,string,nsub,lim)
     outfile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
     do k = 1, j - 1
        if (outfile(j).eq.outfile(k)) call mio_err (6,mem(81),lmem(81),mem(89),lmem(89),outfile(j),80,mem(86),lmem(86))
     end do
     do k = 1, 3
        if (outfile(j).eq.infile(k)) call mio_err (6,mem(81),lmem(81),mem(89),lmem(89),outfile(j),80,mem(86),lmem(86))
     end do
  end do
  !
  ! Dump files
  do j = 1, 4
     read (15,'(a150)') string
     call mio_spl (150,string,nsub,lim)
     dumpfile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
     do k = 1, j - 1
        if (dumpfile(j).eq.dumpfile(k)) call mio_err (6,mem(81),lmem(81),mem(89),lmem(89),dumpfile(j),80,mem(86),lmem(86))
     end do
     do k = 1, 3
        if (dumpfile(j).eq.infile(k)) call mio_err (6,mem(81),lmem(81),mem(89),lmem(89),dumpfile(j),80,mem(86),lmem(86))
     end do
     do k = 1, 3
        if (dumpfile(j).eq.outfile(k)) call mio_err (6,mem(81),lmem(81),mem(89),lmem(89),dumpfile(j),80,mem(86),lmem(86))
     end do
  end do
  close (15)
  !
  ! Find out if this is an old integration (i.e. does the restart file exist)
  inquire (file=dumpfile(4), exist=oldflag)
  !
  ! Check if information file exists, and append a continuation message
  if (oldflag) then
     inquire (file=outfile(3), exist=test)
     if (.not.test) call mio_err (6,mem(81),lmem(81),mem(88),lmem(88),' ',1,outfile(3),80)
     open(23,file=outfile(3),status='old',access='append',iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile(3))
        stop
     end if
  else
     !
     ! If new integration, check information file doesn't exist, and then create it
     inquire (file=outfile(3), exist=test)
     if (test) call mio_err (6,mem(81),lmem(81),mem(87),lmem(87),' ',1,outfile(3),80)
     open(23, file = outfile(3), status = 'new', iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile(3))
        stop
     end if
  end if
  !
  !------------------------------------------------------------------------------
  !
  !  READ  IN  INTEGRATION  PARAMETERS
  !
  ! Check if the file containing integration parameters exists, and open it
  filename = infile(3)
  if (oldflag) filename = dumpfile(3)
  inquire (file=filename, exist=test)
  if (.not.test) call mio_err (23,mem(81),lmem(81),mem(88),lmem(88),' ',1,filename,80)
  open(13, file=filename, status='old', iostat=error)
  if (error /= 0) then
     write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(filename)
     stop
  end if
  !
  ! Read integration parameters
  lineno = 0
  do j = 1, 26
     do
        lineno = lineno + 1
        read (13,'(a150)') string
        if (string(1:1).ne.')') exit
     end do

     call mio_spl (150,string,nsub,lim)
     c80(1:3) = '   '
     c80 = string(lim(1,nsub):lim(2,nsub))
     if (j.eq.1) then
        algor = 0
        do k = 1, 60
           if (c80(1:3).eq.alg(k)) algor = (k + 4) / 5
        end do
        if (algor.eq.0) call mio_err (23,mem(81),lmem(81),mem(98),lmem(98),c80(lim(1,nsub):lim(2,nsub)),lim(2,nsub)-lim(1,nsub)+1,&
             mem(85),lmem(85))
     end if
     if (j.eq.2) read (c80,*,err=661) tstart
     if (j.eq.3) read (c80,*,err=661) tstop
     if (j.eq.4) read (c80,*,err=661) dtout
     if (j.eq.5) read (c80,*,err=661) h0
     if (j.eq.6) read (c80,*,err=661) tol
     c1 = c80(1:1)
     if ((j.eq.7).and.((c1.eq.'y').or.(c1.eq.'Y'))) opt(1) = 1
     if ((j.eq.8).and.((c1.eq.'n').or.(c1.eq.'N'))) opt(2) = 0
     if ((j.eq.9).and.((c1.eq.'y').or.(c1.eq.'Y'))) opt(2) = 2
     if ((j.eq.10).and.((c1.eq.'d').or.(c1.eq.'D'))) opt(3) = 0
     if ((j.eq.11).and.((c1.eq.'y').or.(c1.eq.'Y'))) opt(3) = opt(3) + 2
     if (j.eq.12) then
        if((c1.eq.'l').or.(c1.eq.'L')) then
           opt(4) = 1
        else if (j.eq.12.and.(c1.eq.'m'.or.c1.eq.'M')) then
           opt(4) = 2
        else if (j.eq.12.and.(c1.eq.'h'.or.c1.eq.'H')) then
           opt(4) = 3
        else
           goto 661
        end if
     end if
     if (j.eq.15.and.(c1.eq.'y'.or.c1.eq.'Y')) opt(8) = 1
     if (j.eq.16) read (c80,*,err=661) rmax
     if (j.eq.17) read (c80,*,err=661) rcen
     if (j.eq.18) read (c80,*,err=661) m(1)
     if (j.eq.19) read (c80,*,err=661) jcen(1)
     if (j.eq.20) read (c80,*,err=661) jcen(2)
     if (j.eq.21) read (c80,*,err=661) jcen(3)
     if (j.eq.24) read (c80,*,err=661) cefac
     if (j.eq.25) read (c80,*,err=661) ndump
     if (j.eq.26) read (c80,*,err=661) nfun
  end do
  h0 = abs(h0)
  tol = abs(tol)
  rmax = abs(rmax)
  rcen = abs(rcen)
  cefac = abs(cefac)
  close (13)
  !
  ! Change quantities for central object to suitable units
  m(1) = abs(m(1)) * K2
  jcen(1) = jcen(1) * rcen * rcen
  jcen(2) = jcen(2) * rcen * rcen * rcen * rcen
  jcen(3) = jcen(3) * rcen * rcen * rcen * rcen * rcen * rcen
  s(1,1) = 0.d0
  s(2,1) = 0.d0
  s(3,1) = 0.d0
  !
  ! Make sure that RCEN isn't too small, since it is used to scale the output
  ! (Minimum value corresponds to a central body with density 100g/cm^3).
  temp = 1.1235d-3 * m(1) ** .333333333333333d0
  if (rcen.lt.temp) then
     rcen = temp
     write (13,'(/,2a)') mem(121)(1:lmem(121)),mem(131)(1:lmem(131))
  end if
  !
  !------------------------------------------------------------------------------
  !
  !  READ  IN  DATA  FOR  BIG  AND  SMALL  BODIES
  !
  nbod = 1
  do j = 1, 2
     if (j.eq.2) nbig = nbod
     !
     ! Check if the file containing data for Big bodies exists, and open it
     filename = infile(j)
     if (oldflag) filename = dumpfile(j)
     inquire (file=filename, exist=test)
     if (.not.test) call mio_err (23,mem(81),lmem(81),mem(88),lmem(88),' ',1,filename,80)
     open (11, file=filename, status='old', iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(filename)
        stop
     end if
     !
     ! Read data style
     
     do
      read (11,'(a150)') string
      if (string(1:1).ne.')') exit
    end do
     call mio_spl (150,string,nsub,lim)
     c3 = string(lim(1,nsub):(lim(1,nsub)+2))
     if (c3.eq.'Car'.or.c3.eq.'car'.or.c3.eq.'CAR') then
        informat = 1
     else if (c3.eq.'Ast'.or.c3.eq.'ast'.or.c3.eq.'AST') then
        informat = 2
     else if (c3.eq.'Com'.or.c3.eq.'com'.or.c3.eq.'COM') then
        informat = 3
     else
        call mio_err (23,mem(81),lmem(81),mem(91),lmem(91),' ',1,mem(82+j),lmem(82+j))
     end if
     !
     ! Read epoch of Big bodies
     if (j.eq.1) then
      do
        read (11,'(a150)') string
        if (string(1:1).ne.')') exit
      end do
        call mio_spl (150,string,nsub,lim)
        read (string(lim(1,nsub):lim(2,nsub)),*,err=667) time
     end if
     !
     ! Read information for each object
     do
      read (11,'(a)',iostat=error) string
      if (error /= 0) exit
      
      if (string(1:1).eq.')') cycle
    
     call mio_spl (150,string,nsub,lim)
     
     if (lim(1,1).eq.-1) exit
     !
     ! Determine the name of the object
     nbod = nbod + 1
     if (nbod.gt.NMAX) call mio_err (23,mem(81),lmem(81),mem(90),lmem(90),' ',1,mem(82),lmem(82))
     !
     if ((lim(2,1)-lim(1,1)).gt.7) then
        write (23,'(/,3a)') mem(121)(1:lmem(121)),mem(122)(1:lmem(122)),string( lim(1,1):lim(2,1) )
     end if
     id(nbod) = string( lim(1,1):min(7+lim(1,1),lim(2,1)) )
     ! Check if another object has the same name
     do k = 1, nbod - 1
        if (id(k).eq.id(nbod)) call mio_err (23,mem(81),lmem(81),mem(103),lmem(103),id(nbod),8,' ',1)
     end do
     !
     ! Default values of mass, close-encounter limit, density etc.
     m(nbod) = 0.d0
     rceh(nbod) = 1.d0
     rho(nbod) = rhocgs
     epoch(nbod) = time
     do k = 1, 4
        ngf(k,nbod) = 0.d0
     end do
     !
     ! Read values of mass, close-encounter limit, density etc.
     do k = 3, nsub, 2
        c80 = string(lim(1,k-1):lim(2,k-1))
        read (string(lim(1,k):lim(2,k)),*,err=666) temp
        if (c80(1:1).eq.'m'.or.c80(1:1).eq.'M') then
           m(nbod) = temp * K2
        else if (c80(1:1).eq.'r'.or.c80(1:1).eq.'R') then
           rceh(nbod) = temp
        else if (c80(1:1).eq.'d'.or.c80(1:1).eq.'D') then
           rho(nbod) = temp * rhocgs
        else if (m(nbod).lt.0.or.rceh(nbod).lt.0.or.rho(nbod).lt.0) then
           call mio_err (23,mem(81),lmem(81),mem(97),lmem(97),id(nbod),8,mem(82+j),lmem(82+j))
        else if (c80(1:2).eq.'ep'.or.c80(1:2).eq.'EP'.or.c80(1:2).eq.'Ep') then
           epoch (nbod) = temp
        else if (c80(1:2).eq.'a1'.or.c80(1:2).eq.'A1') then
           ngf (1,nbod) = temp
        else if (c80(1:2).eq.'a2'.or.c80(1:2).eq.'A2') then
           ngf (2,nbod) = temp
        else if (c80(1:2).eq.'a3'.or.c80(1:2).eq.'A3') then
           ngf (3,nbod) = temp
        else if (c80(1:1).eq.'b'.or.c80(1:1).eq.'B') then
           ngf (4,nbod) = temp
        else
           goto 666
        end if
     end do
     !
     ! If required, read Cartesian coordinates, velocities and spins of the bodies
     jtmp = 100
     do
       read (11,'(a150)',end=666) string
       if (string(1:1).ne.')') exit
     end do 
     backspace 11
     if (informat.eq.1) then
        read (11,*,err=666) x(1,nbod),x(2,nbod),x(3,nbod),v(1,nbod),v(2,nbod),v(3,nbod),s(1,nbod),s(2,nbod),s(3,nbod)
     else
        read (11,*,err=666) a,e,i,p,n,l,s(1,nbod),s(2,nbod),s(3,nbod)
        i = i * DEG2RAD
        p = (p + n) * DEG2RAD
        n = n * DEG2RAD
        temp = m(nbod)  +  m(1)
        !
        ! Alternatively, read Cometary or asteroidal elements
        if (informat.eq.3) then
           q = a
           a = q / (1.d0 - e)
           l = mod (sqrt(temp/(abs(a*a*a))) * (epoch(nbod) - l), TWOPI)
        else
           q = a * (1.d0 - e)
           l = l * DEG2RAD
        end if
        if ((algor.eq.11).and.(nbod.ne.2)) temp = temp + m(2)
        call mco_el2x (temp,q,e,i,p,n,l,x(1,nbod),x(2,nbod),x(3,nbod),v(1,nbod),v(2,nbod),v(3,nbod))
     end if
     !
     s(1,nbod) = s(1,nbod) * K2
     s(2,nbod) = s(2,nbod) * K2
     s(3,nbod) = s(3,nbod) * K2
     !
     
     end do
      close (11)
  end do
  !
  ! Set non-gravitational-forces flag, NGFLAG
  ngflag = 0
  do j = 2, nbod
     if (ngf(1,j).ne.0.or.ngf(2,j).ne.0.or.ngf(3,j).ne.0) then
        if (ngflag.eq.0) ngflag = 1
        if (ngflag.eq.2) ngflag = 3
     else if (ngf(4,j).ne.0) then
        if (ngflag.eq.0) ngflag = 2
        if (ngflag.eq.1) ngflag = 3
     end if
  end do
  !
  !------------------------------------------------------------------------------
  !
  !  IF  CONTINUING  AN  OLD  INTEGRATION
  !
  if (oldflag) then
     if (opt(3).eq.1) then
        call mio_jd2y (time,year,month,t1)
        write (23,'(/,a,i10,i2,f8.5,/)') mem(62)(1:lmem(62)),year,month,t1
     else if (opt(3).eq.3) then
        t1 = (time - tstart) / 365.25d0
        write (23,'(/,a,f18.7,a,/)') mem(62)(1:lmem(62)),t1,mem(2)(1:lmem(2))
     else
        if (opt(3).eq.0) t1 = time
        if (opt(3).eq.2) t1 = time - tstart
        write (23,'(/,a,f18.5,a,/)') mem(62)(1:lmem(62)),t1,mem(1)(1:lmem(1))
     end if
     !
     ! Read in energy and angular momentum variables, and convert to internal units
     open (35, file=dumpfile(4), status='old', iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(dumpfile(4))
        stop
     end if
     read (35,*) opflag
     read (35,*) en(1),am(1),en(3),am(3)
     en(1) = en(1) * K2
     en(3) = en(3) * K2
     am(1) = am(1) * K2
     am(3) = am(3) * K2
     read (35,*) s(1,1),s(2,1),s(3,1)
     s(1,1) = s(1,1) * K2
     s(2,1) = s(2,1) * K2
     s(3,1) = s(3,1) * K2
     close (35)
     if (opflag.eq.0) opflag = 1
     !
     !------------------------------------------------------------------------------
     !
     !  IF  STARTING  A  NEW  INTEGRATION
     !
  else
     opflag = -2
     !
     ! Write integration parameters to information file
     write (23,'(/,a)') mem(11)(1:lmem(11))
     write (23,'(a)') mem(12)(1:lmem(12))
     j = algor + 13
     write (23,'(/,2a)') mem(13)(1:lmem(13)),mem(j)(1:lmem(j))
     if (tstart.ge.1.d11.or.tstart.le.-1.d10) then
        write (23,'(/,a,1p,e19.13,a)') mem(26)(1:lmem(26)),tstart,mem(1)(1:lmem(1))
     else
        write (23,'(/,a,f19.7,a)') mem(26)(1:lmem(26)),tstart,mem(1)(1:lmem(1))
     end if
     if (tstop.ge.1.d11.or.tstop.le.-1.d10) then
        write (23,'(a,1p,e19.13)') mem(27)(1:lmem(27)),tstop
     else
        write (23,'(a,f19.7)') mem(27)(1:lmem(27)),tstop
     end if
     write (23,'(a,f15.3)') mem(28)(1:lmem(28)),dtout
     if (opt(4).eq.1) write (23,'(2a)') mem(40)(1:lmem(40)),mem(7)(1:lmem(7))
     if (opt(4).eq.2) write (23,'(2a)') mem(40)(1:lmem(40)),mem(8)(1:lmem(8))
     if (opt(4).eq.3) write (23,'(2a)') mem(40)(1:lmem(40)),mem(9)(1:lmem(9))
     !
     write (23,'(/,a,f10.3,a)') mem(30)(1:lmem(30)),h0,mem(1)(1:lmem(1))
     write (23,'(a,1p1e10.4)') mem(31)(1:lmem(31)),tol
     write (23,'(a,1p1e10.4,a)') mem(32)(1:lmem(32)),m(1)/K2,mem(3)(1:lmem(3))
     write (23,'(a,1p1e11.4)') mem(33)(1:lmem(33)),jcen(1)/rcen**2
     write (23,'(a,1p1e11.4)') mem(34)(1:lmem(34)),jcen(2)/rcen**4
     write (23,'(a,1p1e11.4)') mem(35)(1:lmem(35)),jcen(3)/rcen**6
     write (23,'(a,1p1e10.4,a)') mem(36)(1:lmem(36)),rmax,mem (4)(1:lmem(4))
     write (23,'(a,1p1e10.4,a)') mem(37)(1:lmem(37)),rcen,mem (4)(1:lmem(4))
     !
     itmp = 5
     if (opt(2).eq.1.or.opt(2).eq.2) itmp = 6
     write (23,'(/,2a)') mem(41)(1:lmem(41)),mem(itmp)(1:lmem(itmp))
     itmp = 5
     if (opt(2).eq.2) itmp = 6
     write (23,'(2a)') mem(42)(1:lmem(42)),mem(itmp)(1:lmem(itmp))
     itmp = 5
     if (opt(7).eq.1) itmp = 6
     write (23,'(2a)') mem(45)(1:lmem(45)),mem(itmp)(1:lmem(itmp))
     itmp = 5
     if (opt(8).eq.1) itmp = 6
     write (23,'(2a)') mem(46)(1:lmem(46)),mem(itmp)(1:lmem(itmp))
     !
     ! Check that element and close-encounter files don't exist, and create them
     do j = 1, 2
        inquire (file=outfile(j), exist=test)
        if (test) call mio_err (23,mem(81),lmem(81),mem(87),lmem(87),' ',1,outfile(j),80)
        open  (20+j, file=outfile(j), status='new', iostat=error)
        if (error /= 0) then
           write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(outfile(j))
           stop
        end if
        close (20+j)
     end do
     !
     ! Check that dump files don't exist, and then create them
     do j = 1, 4
        inquire (file=dumpfile(j), exist=test)
        if (test) call mio_err (23,mem(81),lmem(81),mem(87),lmem(87),' ',1,dumpfile(j),80)
        open  (30+j, file=dumpfile(j), status='new', iostat=error)
        if (error /= 0) then
           write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(dumpfile(j))
           stop
        end if
        close (30+j)
     end do
     !
     ! Write number of Big bodies and Small bodies to information file
     write (23,'(/,a,i4)') mem(38)(1:lmem(38)), nbig - 1
     write (23,'(a,i4)') mem(39)(1:lmem(39)), nbod - nbig
     !
     ! Calculate initial energy and angular momentum and write to information file
     s(1,1) = 0.d0
     s(2,1) = 0.d0
     s(3,1) = 0.d0
     call mxx_en (jcen,nbod,nbig,m,x,v,s,en(1),am(1))
     write (23,'(//,a)') mem(51)(1:lmem(51))
     write (23,'(a)')    mem(52)(1:lmem(52))
     write (23,'(/,a,1p1e12.5,a)') mem(53)(1:lmem(53)),en(1)/K2,mem(72)(1:lmem(72))
     write (23,'(a,1p1e12.5,a)')   mem(54)(1:lmem(54)),am(1)/K2,mem(73)(1:lmem(73))
     !
     ! Initialize lost energy and angular momentum
     en(3) = 0.d0
     am(3) = 0.d0
     !
     ! Write warning messages if necessary
     if (tstop.lt.tstart) write (23,'(/,2a)') mem(121)(1:lmem(121)),mem(123)(1:lmem(123))
     if (nbig.le.0) write (23,'(/,2a)') mem(121)(1:lmem(121)),mem(124)(1:lmem(124))
     if (nbig.eq.nbod) write (23,'(/,2a)') mem(121)(1:lmem(121)),mem(125)(1:lmem(125))
  end if
  !
  !------------------------------------------------------------------------------
  !
  !  CHECK  FOR  ATTEMPTS  TO  DO  INCOMPATIBLE  THINGS
  !
  ! If using close-binary algorithm, set radius of central body to be no less
  ! than the periastron of binary star.
  if (algor.eq.11) then
     temp = m(1) + m(2)
     call mco_x2el (temp,x(1,2),x(2,2),x(3,2),v(1,2),v(2,2),v(3,2),a,tmp2,tmp3,tmp4,tmp5,tmp6)
     rcen = max (rcen, a)
  end if
  !
  ! Check if non-grav forces are being used with an incompatible algorithm
  if ((ngflag.ne.0).and.((algor.eq.3).or.(algor.eq.11).or.(algor.eq.12))) then
     call mio_err (23,mem(81),lmem(81),mem(92),lmem(92),' ',1,mem(85),lmem(85))
  endif
  !
  ! Check if user-defined force routine is being used with wrong algorithm
  if ((opt(8).eq.1).and.((algor.eq.11).or.(algor.eq.12))) call mio_err(23,mem(81),lmem(81),mem(93),lmem(93),' ',1,mem(85),lmem(85))
  !
  ! Check whether MVS is being used to integrate massive Small bodies,
  ! or whether massive Small bodies have different epochs than Big bodies.
  flag1 = .false.
  flag2 = .false.
  do j = nbig + 1, nbod
     if (m(j).ne.0) then
        if (algor.eq.1) call mio_err (23,mem(81),lmem(81),mem(94),lmem(94),' ',1,mem(85),lmem(85))
        flag1 = .true.
     end if
     if (epoch(j).ne.time) flag2 = .true.
  end do
  if (flag1.and.flag2) call mio_err (23,mem(81),lmem(81),mem(95),  lmem(95),' ',1,mem(84),lmem(84))
  !
  ! Check if central oblateness is being used with close-binary algorithm
  if (algor.eq.11.and.(jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0)) then
     call mio_err (23,mem(81),lmem(81),mem(102),lmem(102),' ',1,mem(85),lmem(85))
  endif
  !
  ! Check whether RCEN > RMAX or RMAX/RCEN is very large
  if (rcen.gt.rmax) call mio_err (23,mem(81),lmem(81),mem(105),lmem(105),' ',1,mem(85),lmem(85))
  if (rmax/rcen.ge.1.d12) write (23,'(/,2a,/a)') mem(121)(1:lmem(121)),mem(106)(1:lmem(106)),mem(85)(1:lmem(85))
  close (23)
  return
  !
  ! Error reading from the input file containing integration parameters
661 write (c3,'(i3)') lineno
  call mio_err (23,mem(81),lmem(81),mem(99),lmem(99),c3,3,mem(85),lmem(85))
  !
  ! Error reading from the input file for Big or Small bodies
666 call mio_err (23,mem(81),lmem(81),mem(100),lmem(100),id(nbod),8,mem(82+j),lmem(82+j))
  !
  ! Error reading epoch of Big bodies
667 call mio_err (23,mem(81),lmem(81),mem(101),lmem(101),' ',1,mem(83),lmem(83))
  !
  !------------------------------------------------------------------------------
  !
end subroutine mio_in

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MAL_HVAR.FOR    (ErikSoft   4 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Does an integration using a variable-timestep integration algorithm. The
! particular integrator routine is ONESTEP and the algorithm must use
! coordinates with respect to the central body.
!
! N.B. This routine is also called by the synchronisation routine mxx_sync,
! ===  in which case OPFLAG = -2. Beware when making changes involving OPFLAG.
!
!------------------------------------------------------------------------------
!
subroutine mal_hvar (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,&
     id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,onestep)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use dynamic

  implicit none

  !
  ! Input/Output
  integer :: algor,nbod,nbig,stat(nbod),opt(8),opflag,ngflag,ndump,nfun
  integer :: lmem(NMESS)
  real(double_precision) :: time,tstart,tstop,dtout,h0,tol,jcen(3),rcen,rmax
  real(double_precision) :: en(3),am(3),cefac,m(nbod),xh(3,nbod),vh(3,nbod)
  real(double_precision) :: s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod)
  character*8 id(nbod)
  character*80 outfile(3),dumpfile(4),mem(NMESS)
  !
  ! Local
  integer :: i,j,k,n,itmp,nhit,ihit(CMAX),jhit(CMAX),chit(CMAX)
  integer :: dtflag,ejflag,nowflag,stopflag,nstored,ce(NMAX)
  integer :: nclo,iclo(CMAX),jclo(CMAX),nce,ice(NMAX),jce(NMAX)
  real(double_precision) :: tmp0,h,hdid,tout,tdump,tfun,tlog,tsmall,dtdump,dtfun
  real(double_precision) :: thit(CMAX),dhit(CMAX),thit1,x0(3,NMAX),v0(3,NMAX)
  real(double_precision) :: rce(NMAX),rphys(NMAX),rcrit(NMAX),a(NMAX)
  real(double_precision) :: dclo(CMAX),tclo(CMAX),epoch(NMAX)
  real(double_precision) :: ixvclo(6,CMAX),jxvclo(6,CMAX)
  external mfo_all,onestep
  !
  !------------------------------------------------------------------------------
  !
  ! Initialize variables. DTFLAG = 0 implies first ever call to ONESTEP
  dtout  = abs(dtout)
  dtdump = abs(h0) * ndump
  dtfun  = abs(h0) * nfun
  dtflag = 0
  nstored = 0
  tsmall = h0 * 1.d-8
  h = h0
  do j = 2, nbod
     ce(j) = 0.d0
  end do
  !
  ! Calculate close-encounter limits and physical radii for massive bodies
  call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
  !
  ! Set up time of next output, times of previous dump, log and periodic effect
  if (opflag.eq.-1) then
     tout = tstart
  else
     n = int (abs (time - tstart) / dtout) + 1
     tout = tstart  +  dtout * sign (dble(n), tstop - tstart)
     if ((tstop - tstart)*(tout - tstop).gt.0) tout = tstop
  end if
  tdump = time
  tfun  = time
  tlog  = time
  !
  !------------------------------------------------------------------------------
  !
  !  MAIN  LOOP  STARTS  HERE
  !
  do
     !
     ! Is it time for output ?
     if (abs(tout-time).lt.abs(tsmall).and.opflag.ge.-1) then
        !
        ! Beware: the integration may change direction at this point!!!!
        if (opflag.eq.-1) dtflag = 0
        !
        ! Output data for all bodies
        call mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho,stat,id,opt,opflag,algor,outfile(1))
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,outfile,&
             nstored,0)
        tmp0 = tstop - tout
        tout = tout + sign( min( abs(tmp0), abs(dtout) ), tmp0 )
        !
        ! Update the data dump files
        do j = 2, nbod
           epoch(j) = time
        end do
        call mio_dump (time,tstart,tstop,dtout,algor,h,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,&
             id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
     end if
     !
     ! If integration has finished return to the main part of programme
     if (abs(tstop-time).le.abs(tsmall).and.opflag.ne.-1) return
     !
     ! Set the timestep
     if (opflag.eq.-1) tmp0 = tstart - time
     if (opflag.eq.-2) tmp0 = tstop  - time
     if (opflag.ge.0)  tmp0 = tout   - time
     h = sign ( max( min( abs(tmp0), abs(h) ), tsmall), tmp0 )
     !
     ! Save the current coordinates and velocities
     call mco_iden (time,jcen,nbod,nbig,h,m,xh,vh,x0,v0,ngf,ngflag,opt)
     !
     ! Advance one timestep
     call onestep (time,h,hdid,tol,jcen,nbod,nbig,m,xh,vh,s,rphys,rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce,mfo_all)
     time = time + hdid
     !
     ! Check if close encounters or collisions occurred
     nclo = 0
     call mce_stat (time,h,rcen,nbod,nbig,m,x0,v0,xh,vh,rce,rphys,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit,chit,dhit,&
          thit,thit1,nowflag,stat,outfile(3),mem,lmem)
     !
     !------------------------------------------------------------------------------
     !
     !  CLOSE  ENCOUNTERS
     !
     ! If encounter minima occurred, output details and decide whether to stop
     if (nclo.gt.0.and.opflag.ge.-1) then
        itmp = 1
        if (nhit.ne.0) itmp = 0
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,nclo,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,&
             mem,lmem,outfile,nstored,itmp)
        if (stopflag.eq.1) return
     end if
     !
     !------------------------------------------------------------------------------
     !
     !  COLLISIONS
     !
     ! If a collision occurred, output details and resolve the collision
     if (nhit.gt.0.and.opt(2).ne.0) then
        do k = 1, nhit
           if (chit(k).eq.1) then
              i = ihit(k)
              j = jhit(k)
              call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
           end if
        end do
        !
        ! Remove lost objects, reset flags and recompute Hill and physical radii
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat,id,mem,lmem,outfile(3),itmp)
        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
     end if
     !
     !------------------------------------------------------------------------------
     !
     !  COLLISIONS  WITH  CENTRAL  BODY
     !
     ! Check for collisions
     call mce_cent (time,hdid,rcen,jcen,2,nbod,nbig,m,x0,v0,xh,vh,nhit,jhit,thit,dhit,algor,ngf,ngflag)
     !
     ! Resolve the collisions
     if (nhit.gt.0) then
        do k = 1, nhit
           i = 1
           j = jhit(k)
           call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
        end do
        !
        ! Remove lost objects, reset flags and recompute Hill and physical radii
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat,id,mem,lmem,outfile(3),itmp)
        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
     end if
     !
     !------------------------------------------------------------------------------
     !
     !  DATA  DUMP  AND  PROGRESS  REPORT
     !
     ! Do the data dump
     if (abs(time-tdump).ge.abs(dtdump).and.opflag.ge.-1) then
        do j = 2, nbod
           epoch(j) = time
        end do
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,outfile,&
             nstored,0)
        call mio_dump (time,tstart,tstop,dtout,algor,h,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,&
             id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
     end if
     !
     ! Write a progress report to the log file
     if (abs(time-tlog).ge.abs(dtdump).and.opflag.ge.0) then
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
        call mio_log (time,tstart,en,am,opt,mem,lmem)
        tlog = time
     end if
     !
     !------------------------------------------------------------------------------
     !
     !  CHECK  FOR  EJECTIONS  AND  DO  OTHER  PERIODIC  EFFECTS
     !
     if (abs(time-tfun).ge.abs(dtfun).and.opflag.ge.-1) then
        !
        ! Recompute close encounter limits, to allow for changes in Hill radii
        call mce_hill (nbod,m,xh,vh,rce,a)
        do j = 2, nbod
           rce(j) = rce(j) * rceh(j)
        end do
        !
        ! Check for ejections
        call mxx_ejec (time,tstart,rmax,en,am,jcen,2,nbod,nbig,m,xh,vh,s,stat,id,opt,ejflag,outfile(3),mem,lmem)
        !
        ! Remove lost objects, reset flags and recompute Hill and physical radii
        if (ejflag.ne.0) then
           call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat,id,mem,lmem,outfile(3),itmp)
           dtflag = 1
           if (opflag.ge.0) opflag = 1
           call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
        end if
        tfun = time
     end if
     !
     ! Go on to the next time step
  end do
  !
  !------------------------------------------------------------------------------
  !
end subroutine mal_hvar
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MAL_HCON.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Does an integration using an integrator with a constant stepsize H.
! Input and output to this routine use coordinates XH, and velocities VH,
! with respect to the central body, but the integration algorithm uses
! its own internal coordinates X, and velocities V.
!
! The programme uses the transformation routines COORD and BCOORD to change
! to and from the internal coordinates, respectively.
!
!------------------------------------------------------------------------------
!
subroutine mal_hcon (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,&
     id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,onestep,coord,bcoord)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use dynamic

  implicit none

  !
  ! Input/Output
  integer :: algor,nbod,nbig,stat(nbod),opt(8),opflag,ngflag
  integer :: lmem(NMESS),ndump,nfun
  real(double_precision) :: time,tstart,tstop,dtout,h0,tol,jcen(3),rcen,rmax
  real(double_precision) :: en(3),am(3),cefac,m(nbod),xh(3,nbod),vh(3,nbod)
  real(double_precision) :: s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod)
  character*8 id(nbod)
  character*80 outfile(3),dumpfile(4),mem(NMESS)
  !
  ! Local
  integer :: i,j,k,n,itmp,nclo,nhit,jhit(CMAX),iclo(CMAX),jclo(CMAX)
  integer :: dtflag,ejflag,stopflag,colflag,nstored
  real(double_precision) :: x(3,NMAX),v(3,NMAX),xh0(3,NMAX),vh0(3,NMAX)
  real(double_precision) :: rce(NMAX),rphys(NMAX),rcrit(NMAX),epoch(NMAX)
  real(double_precision) :: hby2,tout,tmp0,tdump,tfun,tlog,dtdump,dtfun
  real(double_precision) :: dclo(CMAX),tclo(CMAX),dhit(CMAX),thit(CMAX)
  real(double_precision) :: ixvclo(6,CMAX),jxvclo(6,CMAX),a(NMAX)
  external onestep,coord,bcoord
  !
  !------------------------------------------------------------------------------
  !
  ! Initialize variables. DTFLAG = 0/2: first call ever/normal
  dtout  = abs(dtout)
  dtdump = abs(h0) * ndump
  dtfun  = abs(h0) * nfun
  dtflag = 0
  nstored = 0
  hby2 = 0.500001d0 * abs(h0)
  !
  ! Calculate close-encounter limits and physical radii
  call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
  !
  ! Set up time of next output, times of previous dump, log and periodic effect
  if (opflag.eq.-1) then
     tout = tstart
  else
     n = int (abs (time-tstart) / dtout) + 1
     tout = tstart  +  dtout * sign (dble(n), tstop - tstart)
     if ((tstop-tstart)*(tout-tstop).gt.0) tout = tstop
  end if
  tdump = time
  tfun  = time
  tlog  = time
  !
  ! Convert to internal coordinates and velocities
  call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  !------------------------------------------------------------------------------
  !
  !  MAIN  LOOP  STARTS  HERE
  !
  do
     !
     ! Is it time for output ?
     if (abs(tout-time).le.hby2.and.opflag.ge.-1) then
        !
        ! Beware: the integration may change direction at this point!!!!
        if (opflag.eq.-1.and.dtflag.ne.0) dtflag = 1
        !
        ! Convert to heliocentric coordinates and output data for all bodies
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho,stat,id,opt,opflag,algor,outfile(1))
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,outfile,&
             nstored,0)
        tmp0 = tstop - tout
        tout = tout + sign( min( abs(tmp0), abs(dtout) ), tmp0 )
        !
        ! Update the data dump files
        do j = 2, nbod
           epoch(j) = time
        end do
        call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,&
             stat,id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
     end if
     !
     ! If integration has finished, convert to heliocentric coords and return
     if (abs(tstop-time).le.hby2.and.opflag.ge.0) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        return
     end if
     !
     ! Make sure the integration is heading in the right direction
     ! The timestep will be redo if there are collisions with central body. Else, we exit the loop in the last 'if' statement.
     do
        tmp0 = tstop - time
        if (opflag.eq.-1) tmp0 = tstart - time
        h0 = sign (h0, tmp0)
        !
        ! Save the current heliocentric coordinates and velocities
        if (algor.eq.1) then
           call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,xh0,vh0,ngf,ngflag,opt)
        else
           call bcoord(time,jcen,nbod,nbig,h0,m,x,v,xh0,vh0,ngf,ngflag,opt)
        end if
        call onestep (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag,&
             ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile,mem,lmem)
        time = time + h0
        !
        !------------------------------------------------------------------------------
        !
        !  CLOSE  ENCOUNTERS
        !
        ! If encounter minima occurred, output details and decide whether to stop
        if ((nclo.gt.0).and.(opflag.ge.-1)) then
           itmp = 1
           if (colflag.ne.0) itmp = 0
           call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,nclo,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,&
                outfile,nstored,itmp)
           if (stopflag.eq.1) return
        end if
        !
        !------------------------------------------------------------------------------
        !
        !  COLLISIONS
        !
        ! If collisions occurred, output details and remove lost objects
        if (colflag.ne.0) then
           !
           ! Reindex the surviving objects
           call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
           call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat,id,mem,lmem,outfile(3),itmp)
           !
           ! Reset flags, and calculate new Hill radii and physical radii
           dtflag = 1
           if (opflag.ge.0) opflag = 1
           call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
           call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
        end if
        !
        !------------------------------------------------------------------------------
        !
        !  COLLISIONS  WITH  CENTRAL  BODY
        !
        ! Check for collisions with the central body
        if (algor.eq.1) then
           call mco_iden(time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        else
           call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        end if
        itmp = 2
        if ((algor.eq.11).or.(algor.eq.12)) itmp = 3
        call mce_cent (time,h0,rcen,jcen,itmp,nbod,nbig,m,xh0,vh0,xh,vh,nhit,jhit,thit,dhit,algor,ngf,ngflag)
        !
        ! If something hit the central body, restore the coords prior to this step
        if (nhit.eq.0) then
           exit
        else 
           ! Redo that integration time step
           call mco_iden (time,jcen,nbod,nbig,h0,m,xh0,vh0,xh,vh,ngf,ngflag,opt)
           time = time - h0
           !
           ! Merge the object(s) with the central body
           do k = 1, nhit
              i = 1
              j = jhit(k)
              call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
           end do
           !
           ! Remove lost objects, reset flags and recompute Hill and physical radii
           call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat,id,mem,lmem,outfile(3),itmp)
           if (opflag.ge.0) opflag = 1
           dtflag = 1
           call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
           if (algor.eq.1) then
              call mco_iden (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
           else
              call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
           end if
           !
        end if
     end do
     !
     !------------------------------------------------------------------------------
     !
     !  DATA  DUMP  AND  PROGRESS  REPORT
     !
     ! Convert to heliocentric coords and do the data dump
     if (abs(time-tdump).ge.abs(dtdump).and.opflag.ge.-1) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        do j = 2, nbod
           epoch(j) = time
        end do
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,outfile,&
             nstored,0)
        call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,&
             stat,id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
     end if
     !
     ! Convert to heliocentric coords and write a progress report to the log file
     if (abs(time-tlog).ge.abs(dtdump).and.opflag.ge.0) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
        call mio_log (time,tstart,en,am,opt,mem,lmem)
        tlog = time
     end if
     !
     !------------------------------------------------------------------------------
     !
     !  CHECK  FOR  EJECTIONS  AND  DO  OTHER  PERIODIC  EFFECTS
     !
     if (abs(time-tfun).ge.abs(dtfun).and.opflag.ge.-1) then
        if (algor.eq.1) then
           call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        else
           call bcoord(time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        end if
        !
        ! Recompute close encounter limits, to allow for changes in Hill radii
        call mce_hill (nbod,m,xh,vh,rce,a)
        do j = 2, nbod
           rce(j) = rce(j) * rceh(j)
        end do
        !
        ! Check for ejections
        itmp = 2
        if ((algor.eq.11).or.(algor.eq.12)) itmp = 3
        call mxx_ejec (time,tstart,rmax,en,am,jcen,itmp,nbod,nbig,m,xh,vh,s,stat,id,opt,ejflag,outfile(3),mem,lmem)
        !
        ! Remove ejected objects, reset flags, calculate new Hill and physical radii
        if (ejflag.ne.0) then
           call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat,id,mem,lmem,outfile(3),itmp)
           if (opflag.ge.0) opflag = 1
           dtflag = 1
           call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
           if (algor.eq.1) then
              call mco_iden (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
           else
              call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
           end if
        end if
        tfun = time
     end if
     !
     ! Go on to the next time step
  end do
  !
  !------------------------------------------------------------------------------
  !
end subroutine mal_hcon

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MXX_SYNC.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Synchronizes the epochs of NBIG Big bodies (having a common epoch) and
! NBOD-NBIG Small bodies (possibly having differing epochs), for an 
! integration using MERCURY.
! The Small bodies are picked up in order starting with the one with epoch
! furthest from the time, TSTART, at which the main integration will begin
! producing output.
!
! N.B. The synchronization integrations use Everhart's RA15 routine.
! ---
!
!------------------------------------------------------------------------------
!
subroutine mxx_sync (time,tstart,h0,tol,jcen,nbod,nbig,m,x,v,s,rho,rceh,stat,id,epoch,ngf,opt,ngflag)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8),stat(nbod)
  real(double_precision) :: time,tstart,h0,tol,jcen(3),m(nbod),x(3,nbod),v(3,nbod)
  real(double_precision) :: s(3,nbod),rceh(nbod),rho(nbod),epoch(nbod),ngf(4,nbod)
  character*8 id(nbod)
  !
  ! Local
  integer :: j,k,l,nsml,nsofar,indx(NMAX),itemp,jtemp(NMAX)
  integer :: raflag,nce,ice(NMAX),jce(NMAX)
  real(double_precision) :: temp,epsml(NMAX),rtemp(NMAX)
  real(double_precision) :: h,hdid,tsmall,rphys(NMAX),rcrit(NMAX)
  character*8 ctemp(NMAX)
  external mfo_all
  !
  !------------------------------------------------------------------------------
  !
  ! Reorder Small bodies by epoch so that ep(1) is furthest from TSTART
  nsml = nbod - nbig
  do j = nbig + 1, nbod
     epsml(j-nbig) = epoch(j)
  end do
  call mxx_sort (nsml,epsml,indx)
  !
  if (abs(epsml(1)-tstart).lt.abs(epsml(nsml)-tstart)) then
     k = nsml + 1
     do j = 1, nsml / 2
        l = k - j
        temp = epsml(j)
        epsml(j) = epsml (l)
        epsml(l) = temp
        itemp = indx(j)
        indx(j) = indx (l)
        indx(l) = itemp
     end do
  end if
  !
  do j = nbig + 1, nbod
     epoch(j) = epsml(j-nbig)
  end do
  !
  ! Reorder the other arrays associated with each Small body
  do k = 1, 3
     do j = 1, nsml
        rtemp(j) = x(k,j+nbig)
     end do
     do j = 1, nsml
        x(k,j+nbig) = rtemp(indx(j))
     end do
     do j = 1, nsml
        rtemp(j) = v(k,j+nbig)
     end do
     do j = 1, nsml
        v(k,j+nbig) = rtemp(indx(j))
     end do
     do j = 1, nsml
        rtemp(j) = s(k,j+nbig)
     end do
     do j = 1, nsml
        s(k,j+nbig) = rtemp(indx(j))
     end do
  end do
  !
  do j = 1, nsml
     rtemp(j) = m(j+nbig)
  end do
  do j = 1, nsml
     m(j+nbig) = rtemp(indx(j))
  end do
  do j = 1, nsml
     rtemp(j) = rceh(j+nbig)
  end do
  do j = 1, nsml
     rceh(j+nbig) = rtemp(indx(j))
  end do
  do j = 1, nsml
     rtemp(j) = rho(j+nbig)
  end do
  do j = 1, nsml
     rho(j+nbig) = rtemp(indx(j))
  end do
  !
  do j = 1, nsml
     ctemp(j) = id(j+nbig)
     jtemp(j) = stat(j+nbig)
  end do
  do j = 1, nsml
     id(j+nbig) = ctemp(indx(j))
     stat(j+nbig) = jtemp(indx(j))
  end do
  !
  ! Integrate Small bodies up to the same epoch
  h = h0
  tsmall = h0 * 1.d-12
  raflag = 0
  !
  do j = nbig + 1, nbod
     nsofar = j - 1
     do while (abs(time-epoch(j)).gt.tsmall)
        temp = epoch(j) - time
        h = sign(max(min(abs(temp),abs(h)),tsmall),temp)
        call mdt_ra15 (time,h,hdid,tol,jcen,nsofar,nbig,m,x,v,s,rphys,rcrit,ngf,stat,raflag,ngflag,opt,nce,ice,jce,mfo_all)
        time = time + hdid
     end do
     raflag = 1
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mxx_sync

end program mercury
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_USER.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Applies an arbitrary force, defined by the user.
!
! If using with the symplectic algorithm MAL_MVS, the force should be
! small compared with the force from the central object.
! If using with the conservative Bulirsch-Stoer algorithm MAL_BS2, the
! force should not be a function of the velocities.
!
! N.B. All coordinates and velocities must be with respect to central body
! ===
!------------------------------------------------------------------------------
!
subroutine mfo_user (time,jcen,nbod,nbig,m,x,v,a)
  !
  use physical_constant
  use mercury_constant  
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod, nbig
  real(double_precision) :: time,jcen(3),m(nbod),x(3,nbod),v(3,nbod),a(3,nbod)
  !
  ! Local
  integer :: j
  !
  !------------------------------------------------------------------------------
  !
  do j = 1, nbod
     a(1,j) = 0.d0
     a(2,j) = 0.d0
     a(3,j) = 0.d0
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_user
!

!

!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_ACSH.FOR    (ErikSoft  2 March 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates inverse hyperbolic cosine of an angle X (in radians).
!
!------------------------------------------------------------------------------
!
function mco_acsh (x)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  real(double_precision) :: x,mco_acsh
  !
  !------------------------------------------------------------------------------
  !
  if (x.ge.1.d0) then
     mco_acsh = log (x + sqrt(x*x - 1.d0))
  else
     mco_acsh = 0.d0
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end function mco_acsh
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_B2H.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts barycentric coordinates to coordinates with respect to the central
! body.
!
!------------------------------------------------------------------------------
!
subroutine mco_b2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag,opt)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
  real(double_precision) :: vh(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j
  !
  !------------------------------------------------------------------------------
  !
  do j = 2, nbod
     xh(1,j) = x(1,j) - x(1,1)
     xh(2,j) = x(2,j) - x(2,1)
     xh(3,j) = x(3,j) - x(3,1)
     vh(1,j) = v(1,j) - v(1,1)
     vh(2,j) = v(2,j) - v(2,1)
     vh(3,j) = v(3,j) - v(3,1)
  enddo
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_b2h
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_DH2H.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts democratic heliocentric coordinates to coordinates with respect
! to the central body.
!
!------------------------------------------------------------------------------
!
subroutine mco_dh2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag,opt)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
  real(double_precision) :: vh(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: mvsum(3),temp
  !
  !------------------------------------------------------------------------------
  !
  mvsum(1) = 0.d0
  mvsum(2) = 0.d0
  mvsum(3) = 0.d0
  !
  do j = 2, nbod
     xh(1,j) = x(1,j)
     xh(2,j) = x(2,j)
     xh(3,j) = x(3,j)
     mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
     mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
     mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
  end do
  !
  temp = 1.d0 / m(1)
  mvsum(1) = temp * mvsum(1)
  mvsum(2) = temp * mvsum(2)
  mvsum(3) = temp * mvsum(3)
  !
  do j = 2, nbod
     vh(1,j) = v(1,j) + mvsum(1)
     vh(2,j) = v(2,j) + mvsum(2)
     vh(3,j) = v(3,j) + mvsum(3)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_dh2h
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_IDEN.FOR    (ErikSoft   2 November 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Makes a new copy of a set of coordinates.
!
!------------------------------------------------------------------------------
!
subroutine mco_iden (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
  real(double_precision) :: vh(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j
  !
  !------------------------------------------------------------------------------
  !
  do j = 1, nbod
     x(1,j) = xh(1,j)
     x(2,j) = xh(2,j)
     x(3,j) = xh(3,j)
     v(1,j) = vh(1,j)
     v(2,j) = vh(2,j)
     v(3,j) = vh(3,j)
  enddo
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_iden
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_MVS2H.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Applies a symplectic corrector, which converts coordinates for a second-
! order mixed-variable symplectic integrator to ones with respect to the
! central body.
!
!------------------------------------------------------------------------------
!
subroutine mco_mvs2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag,opt)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use drift
  
  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
  real(double_precision) :: vh(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j,k,iflag,stat(NMAX)
  real(double_precision) :: minside,msofar,gm(NMAX),a(3,NMAX),xj(3,NMAX),vj(3,NMAX)
  real(double_precision) :: ha(2),hb(2),rt10,angf(3,NMAX),ausr(3,NMAX)
  !
  !------------------------------------------------------------------------------
  !
  rt10 = sqrt(10.d0)
  ha(1) =  h * rt10 * 3.d0 / 10.d0
  hb(1) = -h * rt10 / 72.d0
  ha(2) =  h * rt10 / 5.d0
  hb(2) =  h * rt10 / 24.d0
  do j = 2, nbod
     angf(1,j) = 0.d0
     angf(2,j) = 0.d0
     angf(3,j) = 0.d0
     ausr(1,j) = 0.d0
     ausr(2,j) = 0.d0
     ausr(3,j) = 0.d0
  end do
  call mco_iden (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag,opt)
  !
  ! Calculate effective central masses for Kepler drifts
  minside = m(1)
  do j = 2, nbig
     msofar = minside + m(j)
     gm(j) = m(1) * msofar / minside
     minside = msofar
  end do
  !
  ! Two step corrector
  do k = 1, 2
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j(time,jcen,nbig,nbig,h,m,xh,vh,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),ha(k),iflag)
     end do
     !
     ! Advance Interaction Hamiltonian
     call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,xh,xj,a,stat)
     !
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,xh,vh,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,xh,vh,angf,ngf)
     !
     do j = 2, nbod
        vh(1,j) = vh(1,j)  -  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        vh(2,j) = vh(2,j)  -  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        vh(3,j) = vh(3,j)  -  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j(time,jcen,nbig,nbig,h,m,xh,vh,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),-2.d0*ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),-2.d0*ha(k),iflag)
     end do
     !
     ! Advance Interaction Hamiltonian
     call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,xh,xj,a,stat)
     !
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,xh,vh,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,xh,vh,angf,ngf)
     !
     do j = 2, nbod
        vh(1,j) = vh(1,j)  +  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        vh(2,j) = vh(2,j)  +  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        vh(3,j) = vh(3,j)  +  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j(time,jcen,nbig,nbig,h,m,xh,vh,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j),vh(3,j),ha(k),iflag)
     end do
     call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag,opt)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_mvs2h
!

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2B.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts coordinates with respect to the central body to barycentric
! coordinates.
!
!------------------------------------------------------------------------------
!
subroutine mco_h2b (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
  real(double_precision) :: v(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: mtot,temp
  !
  !------------------------------------------------------------------------------
  !
  mtot = 0.d0
  x(1,1) = 0.d0
  x(2,1) = 0.d0
  x(3,1) = 0.d0
  v(1,1) = 0.d0
  v(2,1) = 0.d0
  v(3,1) = 0.d0
  !
  ! Calculate coordinates and velocities of the central body
  do j = 2, nbod
     mtot = mtot  +  m(j)
     x(1,1) = x(1,1)  +  m(j) * xh(1,j)
     x(2,1) = x(2,1)  +  m(j) * xh(2,j)
     x(3,1) = x(3,1)  +  m(j) * xh(3,j)
     v(1,1) = v(1,1)  +  m(j) * vh(1,j)
     v(2,1) = v(2,1)  +  m(j) * vh(2,j)
     v(3,1) = v(3,1)  +  m(j) * vh(3,j)
  enddo
  !
  temp = -1.d0 / (mtot + m(1))
  x(1,1) = temp * x(1,1)
  x(2,1) = temp * x(2,1)
  x(3,1) = temp * x(3,1)
  v(1,1) = temp * v(1,1)
  v(2,1) = temp * v(2,1)
  v(3,1) = temp * v(3,1)
  !
  ! Calculate the barycentric coordinates and velocities
  do j = 2, nbod
     x(1,j) = xh(1,j) + x(1,1)
     x(2,j) = xh(2,j) + x(2,1)
     x(3,j) = xh(3,j) + x(3,1)
     v(1,j) = vh(1,j) + v(1,1)
     v(2,j) = vh(2,j) + v(2,1)
     v(3,j) = vh(3,j) + v(3,1)
  enddo
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_h2b
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2MVS.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Applies an inverse symplectic corrector, which converts coordinates with
! respect to the central body to integrator coordinates for a second-order
! mixed-variable symplectic integrator.
!
!------------------------------------------------------------------------------
!
subroutine mco_h2mvs (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use drift
  
  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
  real(double_precision) :: v(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j,k,iflag,stat(NMAX)
  real(double_precision) :: minside,msofar,gm(NMAX),a(3,NMAX),xj(3,NMAX),vj(3,NMAX)
  real(double_precision) :: ha(2),hb(2),rt10,angf(3,NMAX),ausr(3,NMAX)
  !
  !------------------------------------------------------------------------------
  !
  rt10 = sqrt(10.d0)
  ha(1) = -h * rt10 / 5.d0
  hb(1) = -h * rt10 / 24.d0
  ha(2) = -h * rt10 * 3.d0 / 10.d0
  hb(2) =  h * rt10 / 72.d0
  do j = 2, nbod
     angf(1,j) = 0.d0
     angf(2,j) = 0.d0
     angf(3,j) = 0.d0
     ausr(1,j) = 0.d0
     ausr(2,j) = 0.d0
     ausr(3,j) = 0.d0
  end do
  call mco_iden (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  ! Calculate effective central masses for Kepler drifts
  minside = m(1)
  do j = 2, nbig
     msofar = minside + m(j)
     gm(j) = m(1) * msofar / minside
     minside = msofar
  end do
  !
  do k = 1, 2
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j (time,jcen,nbig,nbig,h,m,x,v,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),ha(k),iflag)
     end do
     !
     ! Advance Interaction Hamiltonian
     call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
     !
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
     !
     do j = 2, nbod
        v(1,j) = v(1,j)  +  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  +  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  +  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j (time,jcen,nbig,nbig,h,m,x,v,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),-2.d0*ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),-2.d0*ha(k),iflag)
     end do
     !
     ! Advance Interaction Hamiltonian
     call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
     !
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
     !
     do j = 2, nbod
        v(1,j) = v(1,j)  -  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  -  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  -  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
     end do
     !
     ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
     call mco_h2j (time,jcen,nbig,nbig,h,m,x,v,xj,vj,ngf,ngflag,opt)
     do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),ha(k),iflag)
     end do
     do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),ha(k),iflag)
     end do
     call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag,opt)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_h2mvs
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2DH.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Convert coordinates with respect to the central body to democratic
! heliocentric coordinates.
!
!------------------------------------------------------------------------------
!
subroutine mco_h2dh (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
  real(double_precision) :: v(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: mtot,temp,mvsum(3)
  !
  !------------------------------------------------------------------------------
  !
  mtot = 0.d0
  mvsum(1) = 0.d0
  mvsum(2) = 0.d0
  mvsum(3) = 0.d0
  !
  do j = 2, nbod
     x(1,j) = xh(1,j)
     x(2,j) = xh(2,j)
     x(3,j) = xh(3,j)
     mtot = mtot + m(j)
     mvsum(1) = mvsum(1)  +  m(j) * vh(1,j)
     mvsum(2) = mvsum(2)  +  m(j) * vh(2,j)
     mvsum(3) = mvsum(3)  +  m(j) * vh(3,j)
  end do
  !
  temp = 1.d0 / (m(1) + mtot)
  !
  mvsum(1) = temp * mvsum(1)
  mvsum(2) = temp * mvsum(2)
  mvsum(3) = temp * mvsum(3)
  !
  do j = 2, nbod
     v(1,j) = vh(1,j) - mvsum(1)
     v(2,j) = vh(2,j) - mvsum(2)
     v(3,j) = vh(3,j) - mvsum(3)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_h2dh
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2J.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts coordinates with respect to the central body to Jacobi coordinates.
!
! N.B. The coordinates respect to the central body for the small bodies
! ===  are assumed to be equal to their Jacobi coordinates.
!
!------------------------------------------------------------------------------
!
subroutine mco_h2j (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,opt)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbig),xh(3,nbig),vh(3,nbig),x(3,nbig)
  real(double_precision) :: v(3,nbig),ngf(4,nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: mtot, mx, my, mz, mu, mv, mw, temp
  !
  !------------------------------------------------------------------------------c
  mtot = m(2)
  x(1,2) = xh(1,2)
  x(2,2) = xh(2,2)
  x(3,2) = xh(3,2)
  v(1,2) = vh(1,2)
  v(2,2) = vh(2,2)
  v(3,2) = vh(3,2)
  mx = m(2) * xh(1,2)
  my = m(2) * xh(2,2)
  mz = m(2) * xh(3,2)
  mu = m(2) * vh(1,2)
  mv = m(2) * vh(2,2)
  mw = m(2) * vh(3,2)
  !
  do j = 3, nbig - 1
     temp = 1.d0 / (mtot + m(1))
     mtot = mtot + m(j)
     x(1,j) = xh(1,j)  -  temp * mx
     x(2,j) = xh(2,j)  -  temp * my
     x(3,j) = xh(3,j)  -  temp * mz
     v(1,j) = vh(1,j)  -  temp * mu
     v(2,j) = vh(2,j)  -  temp * mv
     v(3,j) = vh(3,j)  -  temp * mw
     mx = mx  +  m(j) * xh(1,j)
     my = my  +  m(j) * xh(2,j)
     mz = mz  +  m(j) * xh(3,j)
     mu = mu  +  m(j) * vh(1,j)
     mv = mv  +  m(j) * vh(2,j)
     mw = mw  +  m(j) * vh(3,j)
  enddo
  !
  if (nbig.gt.2) then
     temp = 1.d0 / (mtot + m(1))
     x(1,nbig) = xh(1,nbig)  -  temp * mx
     x(2,nbig) = xh(2,nbig)  -  temp * my
     x(3,nbig) = xh(3,nbig)  -  temp * mz
     v(1,nbig) = vh(1,nbig)  -  temp * mu
     v(2,nbig) = vh(2,nbig)  -  temp * mv
     v(3,nbig) = vh(3,nbig)  -  temp * mw
  end if
  !
  do j = nbig + 1, nbod
     x(1,j) = xh(1,j)
     x(2,j) = xh(2,j)
     x(3,j) = xh(3,j)
     v(1,j) = vh(1,j)
     v(2,j) = vh(2,j)
     v(3,j) = vh(3,j)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_h2j
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_J2H.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts Jacobi coordinates to coordinates with respect to the central
! body.
!
! N.B. The Jacobi coordinates of the small bodies are assumed to be equal
! ===  to their coordinates with respect to the central body.
!
!------------------------------------------------------------------------------
!
subroutine mco_j2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag,opt)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,opt(8)
  real(double_precision) :: time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
  real(double_precision) :: vh(3,nbod),ngf(4,nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: mtot, mx, my, mz, mu, mv, mw, temp
  !
  !------------------------------------------------------------------------------
  !
  xh(1,2) = x(1,2)
  xh(2,2) = x(2,2)
  xh(3,2) = x(3,2)
  vh(1,2) = v(1,2)
  vh(2,2) = v(2,2)
  vh(3,2) = v(3,2)
  mtot = m(2)
  temp = m(2) / (mtot + m(1))
  mx = temp * x(1,2)
  my = temp * x(2,2)
  mz = temp * x(3,2)
  mu = temp * v(1,2)
  mv = temp * v(2,2)
  mw = temp * v(3,2)
  !
  do j = 3, nbig - 1
     xh(1,j) = x(1,j) + mx
     xh(2,j) = x(2,j) + my
     xh(3,j) = x(3,j) + mz
     vh(1,j) = v(1,j) + mu
     vh(2,j) = v(2,j) + mv
     vh(3,j) = v(3,j) + mw
     mtot = mtot + m(j)
     temp = m(j) / (mtot + m(1))
     mx = mx  +  temp * x(1,j)
     my = my  +  temp * x(2,j)
     mz = mz  +  temp * x(3,j)
     mu = mu  +  temp * v(1,j)
     mv = mv  +  temp * v(2,j)
     mw = mw  +  temp * v(3,j)
  enddo
  !
  if (nbig.gt.2) then
     xh(1,nbig) = x(1,nbig) + mx
     xh(2,nbig) = x(2,nbig) + my
     xh(3,nbig) = x(3,nbig) + mz
     vh(1,nbig) = v(1,nbig) + mu
     vh(2,nbig) = v(2,nbig) + mv
     vh(3,nbig) = v(3,nbig) + mw
  end if
  !
  do j = nbig + 1, nbod
     xh(1,j) = x(1,j)
     xh(2,j) = x(2,j)
     xh(3,j) = x(3,j)
     vh(1,j) = v(1,j)
     vh(2,j) = v(2,j)
     vh(3,j) = v(3,j)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_j2h
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_SINE.FOR    (ErikSoft  17 April 1997)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates sin and cos of an angle X (in radians).
!
!------------------------------------------------------------------------------
!
subroutine mco_sine (x,sx,cx)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  real(double_precision) :: x,sx,cx
  !
  ! Local
  real(double_precision) :: pi,twopi
  !
  !------------------------------------------------------------------------------
  !
  pi = 3.141592653589793d0
  twopi = 2.d0 * pi
  !
  if (x.gt.0) then
     x = mod(x,twopi)
  else
     x = mod(x,twopi) + twopi
  end if
  !
  cx = cos(x)
  !
  if (x.gt.pi) then
     sx = -sqrt(1.d0 - cx*cx)
  else
     sx =  sqrt(1.d0 - cx*cx)
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_sine
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_SINH.FOR    (ErikSoft  12 June 1998)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculates sinh and cosh of an angle X (in radians)
!
!------------------------------------------------------------------------------
!
subroutine mco_sinh (x,sx,cx)
  !
  use types_numeriques

  implicit none

  !
  ! Input/Output
  real(double_precision) :: x,sx,cx
  !
  !------------------------------------------------------------------------------
  !
  sx = sinh(x)
  cx = sqrt (1.d0 + sx*sx)
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mco_sinh
!

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_BS1.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H0
! using the Bulirsch-Stoer method. The accelerations are calculated using the 
! subroutine FORCE. The accuracy of the step is approximately determined 
! by the tolerance parameter TOL.
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
subroutine mdt_bs1 (time,h0,hdid,tol,jcen,nbod,nbig,mass,x0,v0,s,rphys,rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce,force)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  real(double_precision) :: SHRINK,GROW
  parameter (SHRINK=.55d0,GROW=1.3d0)
  !
  ! Input/Output
  integer :: nbod, nbig, opt(8), stat(nbod), dtflag, ngflag
  integer :: nce, ice(nce), jce(nce)
  real(double_precision) :: time,h0,hdid,tol,jcen(3),mass(nbod),x0(3,nbod),v0(3,nbod)
  real(double_precision) :: s(3,nbod),ngf(4,nbod),rphys(nbod),rcrit(nbod)
  external force
  !
  ! Local
  integer :: j, j1, k, n
  real(double_precision) :: tmp0,tmp1,tmp2,errmax,tol2,h,hx2,h2(8)
  real(double_precision) :: x(3,NMAX),v(3,NMAX),xend(3,NMAX),vend(3,NMAX)
  real(double_precision) :: a(3,NMAX),a0(3,NMAX),d(6,NMAX,8),xscal(NMAX),vscal(NMAX)
  !
  !------------------------------------------------------------------------------
  !
  tol2 = tol * tol
  !
  ! Calculate arrays used to scale the relative error (R^2 for position and
  ! V^2 for velocity).
  do k = 2, nbod
     tmp1 = x0(1,k)*x0(1,k) + x0(2,k)*x0(2,k) + x0(3,k)*x0(3,k)
     tmp2 = v0(1,k)*v0(1,k) + v0(2,k)*v0(2,k) + v0(3,k)*v0(3,k)
     xscal(k) = 1.d0 / tmp1
     vscal(k) = 1.d0 / tmp2
  end do
  !
  ! Calculate accelerations at the start of the step
  call force (time,jcen,nbod,nbig,mass,x0,v0,s,rcrit,a0,stat,ngf,ngflag,opt,nce,ice,jce)
  !
  do
     !
     ! For each value of N, do a modified-midpoint integration with 2N substeps
     do n = 1, 8
        h = h0 / (2.d0 * float(n))
        h2(n) = .25d0 / (n*n)
        hx2 = h * 2.d0
        !
        do k = 2, nbod
           x(1,k) = x0(1,k) + h*v0(1,k)
           x(2,k) = x0(2,k) + h*v0(2,k)
           x(3,k) = x0(3,k) + h*v0(3,k)
           v(1,k) = v0(1,k) + h*a0(1,k)
           v(2,k) = v0(2,k) + h*a0(2,k)
           v(3,k) = v0(3,k) + h*a0(3,k)
        end do
        call force (time,jcen,nbod,nbig,mass,x,v,s,rcrit,a,stat,ngf,ngflag,opt,nce,ice,jce)
        do k = 2, nbod
           xend(1,k) = x0(1,k) + hx2*v(1,k)
           xend(2,k) = x0(2,k) + hx2*v(2,k)
           xend(3,k) = x0(3,k) + hx2*v(3,k)
           vend(1,k) = v0(1,k) + hx2*a(1,k)
           vend(2,k) = v0(2,k) + hx2*a(2,k)
           vend(3,k) = v0(3,k) + hx2*a(3,k)
        end do
        !
        do j = 2, n
           call force (time,jcen,nbod,nbig,mass,xend,vend,s,rcrit,a,stat,ngf,ngflag,opt,nce,ice,jce)
           do k = 2, nbod
              x(1,k) = x(1,k) + hx2*vend(1,k)
              x(2,k) = x(2,k) + hx2*vend(2,k)
              x(3,k) = x(3,k) + hx2*vend(3,k)
              v(1,k) = v(1,k) + hx2*a(1,k)
              v(2,k) = v(2,k) + hx2*a(2,k)
              v(3,k) = v(3,k) + hx2*a(3,k)
           end do
           call force (time,jcen,nbod,nbig,mass,x,v,s,rcrit,a,stat,ngf,ngflag,opt,nce,ice,jce)
           do k = 2, nbod
              xend(1,k) = xend(1,k) + hx2*v(1,k)
              xend(2,k) = xend(2,k) + hx2*v(2,k)
              xend(3,k) = xend(3,k) + hx2*v(3,k)
              vend(1,k) = vend(1,k) + hx2*a(1,k)
              vend(2,k) = vend(2,k) + hx2*a(2,k)
              vend(3,k) = vend(3,k) + hx2*a(3,k)
           end do
        end do
        !
        call force (time,jcen,nbod,nbig,mass,xend,vend,s,rcrit,a,stat,ngf,ngflag,opt,nce,ice,jce)
        !
        do k = 2, nbod
           d(1,k,n) = .5d0*(xend(1,k) + x(1,k) + h*vend(1,k))
           d(2,k,n) = .5d0*(xend(2,k) + x(2,k) + h*vend(2,k))
           d(3,k,n) = .5d0*(xend(3,k) + x(3,k) + h*vend(3,k))
           d(4,k,n) = .5d0*(vend(1,k) + v(1,k) + h*a(1,k))
           d(5,k,n) = .5d0*(vend(2,k) + v(2,k) + h*a(2,k))
           d(6,k,n) = .5d0*(vend(3,k) + v(3,k) + h*a(3,k))
        end do
        !
        ! Update the D array, used for polynomial extrapolation
        do j = n - 1, 1, -1
           j1 = j + 1
           tmp0 = 1.d0 / (h2(j) - h2(n))
           tmp1 = tmp0 * h2(j1)
           tmp2 = tmp0 * h2(n)
           do k = 2, nbod
              d(1,k,j) = tmp1 * d(1,k,j1)  -  tmp2 * d(1,k,j)
              d(2,k,j) = tmp1 * d(2,k,j1)  -  tmp2 * d(2,k,j)
              d(3,k,j) = tmp1 * d(3,k,j1)  -  tmp2 * d(3,k,j)
              d(4,k,j) = tmp1 * d(4,k,j1)  -  tmp2 * d(4,k,j)
              d(5,k,j) = tmp1 * d(5,k,j1)  -  tmp2 * d(5,k,j)
              d(6,k,j) = tmp1 * d(6,k,j1)  -  tmp2 * d(6,k,j)
           end do
        end do
        !
        ! After several integrations, test the relative error on extrapolated values
        if (n.gt.3) then
           errmax = 0.d0
           !
           ! Maximum relative position and velocity errors (last D term added)
           do k = 2, nbod
              tmp1 = max( d(1,k,1)*d(1,k,1), d(2,k,1)*d(2,k,1),        d(3,k,1)*d(3,k,1) )
              tmp2 = max( d(4,k,1)*d(4,k,1), d(5,k,1)*d(5,k,1),        d(6,k,1)*d(6,k,1) )
              errmax = max(errmax, tmp1*xscal(k), tmp2*vscal(k))
           end do
           !
           ! If error is smaller than TOL, update position and velocity arrays, and exit
           if (errmax.le.tol2) then
              do k = 2, nbod
                 x0(1,k) = d(1,k,1)
                 x0(2,k) = d(2,k,1)
                 x0(3,k) = d(3,k,1)
                 v0(1,k) = d(4,k,1)
                 v0(2,k) = d(5,k,1)
                 v0(3,k) = d(6,k,1)
              end do
              !
              do j = 2, n
                 do k = 2, nbod
                    x0(1,k) = x0(1,k) + d(1,k,j)
                    x0(2,k) = x0(2,k) + d(2,k,j)
                    x0(3,k) = x0(3,k) + d(3,k,j)
                    v0(1,k) = v0(1,k) + d(4,k,j)
                    v0(2,k) = v0(2,k) + d(5,k,j)
                    v0(3,k) = v0(3,k) + d(6,k,j)
                 end do
              end do
              !
              ! Save the actual stepsize used
              hdid = h0
              !
              ! Recommend a new stepsize for the next call to this subroutine
              if (n.eq.8) h0 = h0 * SHRINK
              if (n.lt.7) h0 = h0 * GROW
              return
           end if
        end if
        !
     end do
     !
     ! If errors were too large, redo the step with half the previous step size.
     h0 = h0 * .5d0
  end do
  !
  !------------------------------------------------------------------------------
  !
end subroutine mdt_bs1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_BS2.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H0
! using the Bulirsch-Stoer method. The accelerations are calculated using the 
! subroutine FORCE. The accuracy of the step is approximately determined 
! by the tolerance parameter TOL.
!
! N.B. This version only works for conservative systems (i.e. force is a
! ===  function of position only) !!!! Hence, non-gravitational forces
!      and post-Newtonian corrections cannot be used.
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
subroutine mdt_bs2 (time,h0,hdid,tol,jcen,nbod,nbig,mass,x0,v0,s,rphys,rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce,force)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  real(double_precision) :: SHRINK,GROW
  parameter (SHRINK=.55d0,GROW=1.3d0)
  !
  ! Input/Output
  integer :: nbod, nbig, opt(8), stat(nbod), dtflag, ngflag
  real(double_precision) :: time,h0,hdid,tol,jcen(3),mass(nbod),x0(3,nbod),v0(3,nbod)
  real(double_precision) :: s(3,nbod),ngf(4,nbod),rphys(nbod),rcrit(nbod)
  integer :: nce,ice(nce),jce(nce)
  external force
  !
  ! Local
  integer :: j,j1,k,n
  real(double_precision) :: tmp0,tmp1,tmp2,errmax,tol2,h,h2(12),hby2,h2by2
  real(double_precision) :: xend(3,NMAX),b(3,NMAX),c(3,NMAX)
  real(double_precision) :: a(3,NMAX),a0(3,NMAX),d(6,NMAX,12),xscal(NMAX),vscal(NMAX)
  !
  !------------------------------------------------------------------------------
  !
  tol2 = tol * tol
  !
  ! Calculate arrays used to scale the relative error (R^2 for position and
  ! V^2 for velocity).
  do k = 2, nbod
     tmp1 = x0(1,k)*x0(1,k) + x0(2,k)*x0(2,k) + x0(3,k)*x0(3,k)
     tmp2 = v0(1,k)*v0(1,k) + v0(2,k)*v0(2,k) + v0(3,k)*v0(3,k)
     xscal(k) = 1.d0 / tmp1
     vscal(k) = 1.d0 / tmp2
  end do
  !
  ! Calculate accelerations at the start of the step
  call force (time,jcen,nbod,nbig,mass,x0,v0,s,rcrit,a0,stat,ngf,ngflag,opt,nce,ice,jce)
  !
  do
     !
     ! For each value of N, do a modified-midpoint integration with N substeps
     do n = 1, 12
        h = h0 / (dble(n))
        hby2  = .5d0 * h
        h2(n) = h * h
        h2by2 = .5d0 * h2(n)
        !
        do k = 2, nbod
           b(1,k) = .5d0*a0(1,k)
           b(2,k) = .5d0*a0(2,k)
           b(3,k) = .5d0*a0(3,k)
           c(1,k) = 0.d0
           c(2,k) = 0.d0
           c(3,k) = 0.d0
           xend(1,k) = h2by2 * a0(1,k)  +  h * v0(1,k)  +  x0(1,k)
           xend(2,k) = h2by2 * a0(2,k)  +  h * v0(2,k)  +  x0(2,k)
           xend(3,k) = h2by2 * a0(3,k)  +  h * v0(3,k)  +  x0(3,k)
        end do
        !
        do j = 2, n
           call force (time,jcen,nbod,nbig,mass,xend,v0,s,rcrit,a,stat,ngf,ngflag,opt,nce,ice,jce)
           tmp0 = h * dble(j)
           do k = 2, nbod
              b(1,k) = b(1,k) + a(1,k)
              b(2,k) = b(2,k) + a(2,k)
              b(3,k) = b(3,k) + a(3,k)
              c(1,k) = c(1,k) + b(1,k)
              c(2,k) = c(2,k) + b(2,k)
              c(3,k) = c(3,k) + b(3,k)
              xend(1,k) = h2(n)*c(1,k) + h2by2*a0(1,k) + tmp0*v0(1,k)      + x0(1,k)
              xend(2,k) = h2(n)*c(2,k) + h2by2*a0(2,k) + tmp0*v0(2,k)      + x0(2,k)
              xend(3,k) = h2(n)*c(3,k) + h2by2*a0(3,k) + tmp0*v0(3,k)      + x0(3,k)
           end do
        end do
        !
        call force (time,jcen,nbod,nbig,mass,xend,v0,s,rcrit,a,stat,ngf,ngflag,opt,nce,ice,jce)
        !
        do k = 2, nbod
           d(1,k,n) = xend(1,k)
           d(2,k,n) = xend(2,k)
           d(3,k,n) = xend(3,k)
           d(4,k,n) = h*b(1,k) + hby2*a(1,k) + v0(1,k)
           d(5,k,n) = h*b(2,k) + hby2*a(2,k) + v0(2,k)
           d(6,k,n) = h*b(3,k) + hby2*a(3,k) + v0(3,k)
        end do
        !
        ! Update the D array, used for polynomial extrapolation
        do j = n - 1, 1, -1
           j1 = j + 1
           tmp0 = 1.d0 / (h2(j) - h2(n))
           tmp1 = tmp0 * h2(j1)
           tmp2 = tmp0 * h2(n)
           do k = 2, nbod
              d(1,k,j) = tmp1 * d(1,k,j1)  -  tmp2 * d(1,k,j)
              d(2,k,j) = tmp1 * d(2,k,j1)  -  tmp2 * d(2,k,j)
              d(3,k,j) = tmp1 * d(3,k,j1)  -  tmp2 * d(3,k,j)
              d(4,k,j) = tmp1 * d(4,k,j1)  -  tmp2 * d(4,k,j)
              d(5,k,j) = tmp1 * d(5,k,j1)  -  tmp2 * d(5,k,j)
              d(6,k,j) = tmp1 * d(6,k,j1)  -  tmp2 * d(6,k,j)
           end do
        end do
        !
        ! After several integrations, test the relative error on extrapolated values
        if (n.gt.3) then
           errmax = 0.d0
           !
           ! Maximum relative position and velocity errors (last D term added)
           do k = 2, nbod
              tmp1 = max( d(1,k,1)*d(1,k,1), d(2,k,1)*d(2,k,1),        d(3,k,1)*d(3,k,1) )
              tmp2 = max( d(4,k,1)*d(4,k,1), d(5,k,1)*d(2,k,1),        d(6,k,1)*d(6,k,1) )
              errmax = max( errmax, tmp1*xscal(k), tmp2*vscal(k) )
           end do
           !
           ! If error is smaller than TOL, update position and velocity arrays and exit
           if (errmax.le.tol2) then
              do k = 2, nbod
                 x0(1,k) = d(1,k,1)
                 x0(2,k) = d(2,k,1)
                 x0(3,k) = d(3,k,1)
                 v0(1,k) = d(4,k,1)
                 v0(2,k) = d(5,k,1)
                 v0(3,k) = d(6,k,1)
              end do
              !
              do j = 2, n
                 do k = 2, nbod
                    x0(1,k) = x0(1,k) + d(1,k,j)
                    x0(2,k) = x0(2,k) + d(2,k,j)
                    x0(3,k) = x0(3,k) + d(3,k,j)
                    v0(1,k) = v0(1,k) + d(4,k,j)
                    v0(2,k) = v0(2,k) + d(5,k,j)
                    v0(3,k) = v0(3,k) + d(6,k,j)
                 end do
              end do
              !
              ! Save the actual stepsize used
              hdid = h0
              !
              ! Recommend a new stepsize for the next call to this subroutine
              if (n.ge.8) h0 = h0 * SHRINK
              if (n.lt.7) h0 = h0 * GROW
              return
           end if
        end if
        !
     end do
     !
     ! If errors were too large, redo the step with half the previous step size.
     h0 = h0 * .5d0
  end do
  !
  !------------------------------------------------------------------------------
  !
end subroutine mdt_bs2
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_HY.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H
! using a second-order hybrid-symplectic integrator algorithm
!
! DTFLAG = 0 implies first ever call to this subroutine, 
!        = 1 implies first call since number/masses of objects changed.
!        = 2 normal call
!
! N.B. Input/output must be in democratic heliocentric coordinates.
! ===
!
!------------------------------------------------------------------------------
!
subroutine mdt_hy (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag,ngflag,&
     opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile,mem,lmem)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use drift
  use dynamic
  
  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,stat(nbod),algor,opt(8),dtflag,ngflag,opflag
  integer :: colflag,lmem(NMESS),nclo,iclo(CMAX),jclo(CMAX)
  real(double_precision) :: time,tstart,h0,tol,rmax,en(3),am(3),jcen(3),rcen
  real(double_precision) :: m(nbod),x(3,nbod),v(3,nbod),s(3,nbod),rphys(nbod)
  real(double_precision) :: rce(nbod),rcrit(nbod),ngf(4,nbod),tclo(CMAX),dclo(CMAX)
  real(double_precision) :: ixvclo(6,CMAX),jxvclo(6,CMAX)
  character*80 outfile(3),mem(NMESS)
  character*8 id(nbod)
  !
  ! Local
  integer :: j,nce,ice(NMAX),jce(NMAX),ce(NMAX),iflag
  real(double_precision) :: a(3,NMAX),hby2,hrec,x0(3,NMAX),v0(3,NMAX),mvsum(3),temp
  real(double_precision) :: angf(3,NMAX),ausr(3,NMAX)
  external mfo_hkce
  !
  !------------------------------------------------------------------------------
  !
  save a, hrec, angf, ausr
  hby2 = h0 * .5d0
  nclo = 0
  colflag = 0
  !
  ! If accelerations from previous call are not valid, calculate them now
  if (dtflag.ne.2) then
     if (dtflag.eq.0) hrec = h0
     call mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
     dtflag = 2
     do j = 2, nbod
        angf(1,j) = 0.d0
        angf(2,j) = 0.d0
        angf(3,j) = 0.d0
        ausr(1,j) = 0.d0
        ausr(2,j) = 0.d0
        ausr(3,j) = 0.d0
     end do
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
  end if
  !
  ! Advance interaction Hamiltonian for H/2
  do j = 2, nbod
     v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
     v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
     v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
  end do
  !
  ! Advance solar Hamiltonian for H/2
  mvsum(1) = 0.d0
  mvsum(2) = 0.d0
  mvsum(3) = 0.d0
  do j = 2, nbod
     mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
     mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
     mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
  end do
  !
  temp = hby2 / m(1)
  mvsum(1) = temp * mvsum(1)
  mvsum(2) = temp * mvsum(2)
  mvsum(3) = temp * mvsum(3)
  do j = 2, nbod
     x(1,j) = x(1,j)  +  mvsum(1)
     x(2,j) = x(2,j)  +  mvsum(2)
     x(3,j) = x(3,j)  +  mvsum(3)
  end do
  !
  ! Save the current coordinates and velocities
  call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,x0,v0,ngf,ngflag,opt)
  !
  ! Advance H_K for H
  do j = 2, nbod
     iflag = 0
     call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),h0,iflag)
  end do
  !
  ! Check whether any object separations were < R_CRIT whilst advancing H_K
  call mce_snif (h0,2,nbod,nbig,x0,v0,x,v,rcrit,ce,nce,ice,jce)
  !
  ! If objects had close encounters, advance H_K using Bulirsch-Stoer instead
  if (nce.gt.0) then
     do j = 2, nbod
        if (ce(j).ne.0) then
           x(1,j) = x0(1,j)
           x(2,j) = x0(2,j)
           x(3,j) = x0(3,j)
           v(1,j) = v0(1,j)
           v(2,j) = v0(2,j)
           v(3,j) = v0(3,j)
        end if
     end do
     call mdt_hkce (time,tstart,h0,hrec,tol,rmax,en(3),jcen,rcen,nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,ngflag,&
          colflag,ce,nce,ice,jce,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile,mem,lmem,mfo_hkce)
  end if
  !
  ! Advance solar Hamiltonian for H/2
  mvsum(1) = 0.d0
  mvsum(2) = 0.d0
  mvsum(3) = 0.d0
  do j = 2, nbod
     mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
     mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
     mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
  end do
  !
  temp = hby2 / m(1)
  mvsum(1) = temp * mvsum(1)
  mvsum(2) = temp * mvsum(2)
  mvsum(3) = temp * mvsum(3)
  do j = 2, nbod
     x(1,j) = x(1,j)  +  mvsum(1)
     x(2,j) = x(2,j)  +  mvsum(2)
     x(3,j) = x(3,j)  +  mvsum(3)
  end do
  !
  ! Advance interaction Hamiltonian for H/2
  call mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
  if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
  if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
  !
  do j = 2, nbod
     v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
     v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
     v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mdt_hy
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_HKCE.FOR    (ErikSoft   1 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H under
! the Hamiltonian H_K, including close-encounter terms.
!
!------------------------------------------------------------------------------
!
subroutine mdt_hkce (time,tstart,h0,hrec,tol,rmax,elost,jcen,rcen,nbod,nbig,m,x,v,s,rphy,rcrit,rce,stat,id,ngf,algor,opt,ngflag,&
     colflag,ce,nce,ice,jce,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile,mem,lmem,force)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use dynamic

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,nce,ice(nce),jce(nce),stat(nbod),ngflag,ce(nbod)
  integer :: algor,opt(8),colflag,lmem(NMESS),nclo,iclo(CMAX)
  integer :: jclo(CMAX)
  real(double_precision) :: time,tstart,h0,hrec,tol,rmax,elost,jcen(3),rcen
  real(double_precision) :: m(nbod),x(3,nbod),v(3,nbod),s(3,nbod)
  real(double_precision) :: rce(nbod),rphy(nbod),rcrit(nbod),ngf(4,nbod)
  real(double_precision) :: tclo(CMAX),dclo(CMAX),ixvclo(6,CMAX),jxvclo(6,CMAX)
  character*80 outfile(3),mem(NMESS)
  character*8 id(nbod)
  external force
  !
  ! Local
  integer :: iback(NMAX),index(NMAX),ibs(NMAX),jbs(NMAX),nclo_old
  integer :: i,j,k,nbs,nbsbig,statbs(NMAX)
  integer :: nhit,ihit(CMAX),jhit(CMAX),chit(CMAX),nowflag,dtflag
  real(double_precision) :: tlocal,hlocal,hdid,tmp0
  real(double_precision) :: mbs(NMAX),xbs(3,NMAX),vbs(3,NMAX),sbs(3,NMAX)
  real(double_precision) :: rcritbs(NMAX),rcebs(NMAX),rphybs(NMAX)
  real(double_precision) :: ngfbs(4,NMAX),x0(3,NMAX),v0(3,NMAX)
  real(double_precision) :: thit(CMAX),dhit(CMAX),thit1,temp
  character*8 idbs(NMAX)
  !
  !------------------------------------------------------------------------------
  !
  ! N.B. Don't set nclo to zero!!
  nbs = 1
  nbsbig = 0
  mbs(1) = m(1)
  if (algor.eq.11) mbs(1) = m(1) + m(2)
  sbs(1,1) = s(1,1)
  sbs(2,1) = s(2,1)
  sbs(3,1) = s(3,1)
  !
  ! Put data for close-encounter bodies into local arrays for use with BS routine
  do j = 2, nbod
     if (ce(j).ne.0) then
        nbs = nbs + 1
        if (j.le.nbig) nbsbig = nbs
        mbs(nbs)   = m(j)
        xbs(1,nbs) = x(1,j)
        xbs(2,nbs) = x(2,j)
        xbs(3,nbs) = x(3,j)
        vbs(1,nbs) = v(1,j)
        vbs(2,nbs) = v(2,j)
        vbs(3,nbs) = v(3,j)
        sbs(1,nbs) = s(1,j)
        sbs(2,nbs) = s(2,j)
        sbs(3,nbs) = s(3,j)
        rcebs(nbs) = rce(j)
        rphybs(nbs) = rphy(j)
        statbs(nbs) = stat(j)
        rcritbs(nbs) = rcrit(j)
        idbs(nbs) = id(j)
        index(nbs) = j
        iback(j) = nbs
     end if
  end do
  !
  do k = 1, nce
     ibs(k) = iback(ice(k))
     jbs(k) = iback(jce(k))
  end do
  !
  tlocal = 0.d0
  hlocal = sign(hrec,h0)
  !
  ! Begin the Bulirsch-Stoer integration
50 continue
  tmp0 = abs(h0) - abs(tlocal)
  hrec = hlocal
  if (abs(hlocal).gt.tmp0) hlocal = sign (tmp0, h0)
  !
  ! Save old coordinates and integrate
  call mco_iden (time,jcen,nbs,0,h0,mbs,xbs,vbs,x0,v0,ngf,ngflag,opt)
  call mdt_bs2 (time,hlocal,hdid,tol,jcen,nbs,nbsbig,mbs,xbs,vbs,sbs,rphybs,rcritbs,ngfbs,statbs,dtflag,ngflag,opt,nce,ibs,jbs,&
       force)
  tlocal = tlocal + hdid
  !
  ! Check for close-encounter minima
  nclo_old = nclo
  temp = time + tlocal
  call mce_stat (temp,hdid,rcen,nbs,nbsbig,mbs,x0,v0,xbs,vbs,rcebs,rphybs,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit,&
       chit,dhit,thit,thit1,nowflag,statbs,outfile(3),mem,lmem)
  !
  ! If collisions occurred, resolve the collision and return a flag
  if ((nhit.gt.0).and.(opt(2).ne.0)) then
     do k = 1, nhit
        if (chit(k).eq.1) then
           i = ihit(k)
           j = jhit(k)
           call mce_coll (thit(k),tstart,elost,jcen,i,j,nbs,nbsbig,mbs,xbs,vbs,sbs,rphybs,statbs,idbs,opt,mem,lmem,outfile(3))
           colflag = colflag + 1
        end if
     end do
  end if
  !
  ! If necessary, continue integrating objects undergoing close encounters
  if ((tlocal - h0)*h0.lt.0) goto 50
  !
  ! Return data for the close-encounter objects to global arrays
  do k = 2, nbs
     j = index(k)
     m(j)   = mbs(k)
     x(1,j) = xbs(1,k)
     x(2,j) = xbs(2,k)
     x(3,j) = xbs(3,k)
     v(1,j) = vbs(1,k)
     v(2,j) = vbs(2,k)
     v(3,j) = vbs(3,k)
     s(1,j) = sbs(1,k)
     s(2,j) = sbs(2,k)
     s(3,j) = sbs(3,k)
     stat(j) = statbs(k)
  end do
  do k = 1, nclo
     iclo(k) = index(iclo(k))
     jclo(k) = index(jclo(k))
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mdt_hkce
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_MVS.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H
! using a second-order mixed-variable symplectic integrator.
!
! DTFLAG = 0 implies first ever call to this subroutine, 
!        = 1 implies first call since number/masses of objects changed.
!        = 2 normal call
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
subroutine mdt_mvs (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag,ngflag,&
     opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile,mem,lmem)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques
  use drift
  use dynamic
  
  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,stat(nbod),algor,opt(8),dtflag,ngflag,opflag
  integer :: colflag,lmem(NMESS),nclo,iclo(CMAX),jclo(CMAX)
  real(double_precision) :: time,tstart,h0,tol,rmax,en(3),am(3),jcen(3),rcen
  real(double_precision) :: m(nbod),x(3,nbod),v(3,nbod),s(3,nbod),rphys(nbod)
  real(double_precision) :: rce(nbod),rcrit(nbod),ngf(4,nbod),tclo(CMAX),dclo(CMAX)
  real(double_precision) :: ixvclo(6,CMAX),jxvclo(6,CMAX)
  character*80 outfile(3),mem(NMESS)
  character*8 id(nbod)
  !
  ! Local
  integer :: j,iflag,nhit,ihit(CMAX),jhit(CMAX),chit(CMAX),nowflag
  real(double_precision) :: xj(3,NMAX),vj(3,NMAX),a(3,NMAX),gm(NMAX),hby2,thit1,temp
  real(double_precision) :: msofar,minside,x0(3,NMAX),v0(3,NMAX),dhit(CMAX),thit(CMAX)
  real(double_precision) :: angf(3,NMAX),ausr(3,NMAX)
  !
  !------------------------------------------------------------------------------
  !
  save a, xj, gm, angf, ausr
  hby2 = .5d0 * h0
  nclo = 0
  !
  ! If accelerations from previous call are not valid, calculate them now,
  ! and also the Jacobi coordinates XJ, and effective central masses GM.
  if (dtflag.ne.2) then
     dtflag = 2
     call mco_h2j (time,jcen,nbig,nbig,h0,m,x,v,xj,vj,ngf,ngflag,opt)
     call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
     !
     minside = m(1)
     do j = 2, nbig
        msofar = minside + m(j)
        gm(j) = m(1) * msofar / minside
        minside = msofar
        angf(1,j) = 0.d0
        angf(2,j) = 0.d0
        angf(3,j) = 0.d0
        ausr(1,j) = 0.d0
        ausr(2,j) = 0.d0
        ausr(3,j) = 0.d0
     end do
     ! If required, apply non-gravitational and user-defined forces
     if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
     if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
  end if
  !
  ! Advance interaction Hamiltonian for H/2
  do j = 2, nbod
     v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
     v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
     v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
  end do
  !
  ! Save current coordinates and velocities
  call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,x0,v0,ngf,ngflag,opt)
  !
  ! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
  call mco_h2j (time,jcen,nbig,nbig,h0,m,x,v,xj,vj,ngf,ngflag,opt)
  do j = 2, nbig
     iflag = 0
     call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j),vj(2,j),vj(3,j),h0,iflag)
  end do
  do j = nbig + 1, nbod
     iflag = 0
     call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),h0,iflag)
  end do
  call mco_j2h (time,jcen,nbig,nbig,h0,m,xj,vj,x,v,ngf,ngflag,opt)
  !
  ! Check for close-encounter minima during drift step
  temp = time + h0
  call mce_stat (temp,h0,rcen,nbod,nbig,m,x0,v0,x,v,rce,rphys,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit,chit,dhit,thit,&
       thit1,nowflag,stat,outfile(3),mem,lmem)
  !
  ! Advance interaction Hamiltonian for H/2
  call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
  if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
  if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
  !
  do j = 2, nbod
     v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
     v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
     v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mdt_mvs
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_RA15.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H0 using
! Everhart's RA15 integrator algorithm. The accelerations are calculated
! using the subroutine FORCE. The accuracy of the step is approximately 
! determined by the tolerance parameter TOL.
!
! Based on RADAU by E. Everhart, Physics Department, University of Denver.
! Comments giving equation numbers refer to Everhart (1985) ``An Efficient
! Integrator that Uses Gauss-Radau Spacings'', in The Dynamics of Comets:
! Their Origin and Evolution, p185-202, eds. A. Carusi & G. B. Valsecchi,
! pub Reidel. (A listing of the original subroutine is also given in this 
! paper.)
!
! DTFLAG = 0 implies first ever call to this subroutine, 
!        = 1 implies first call since number/masses of objects changed.
!        = 2 normal call
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
subroutine mdt_ra15 (time,t,tdid,tol,jcen,nbod,nbig,mass,x1,v1,spin,rphys,rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce,force)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,dtflag,ngflag,opt(8),stat(nbod)
  integer :: nce,ice(nce),jce(nce)
  real(double_precision) :: time,t,tdid,tol,jcen(3),mass(nbod)
  real(double_precision) :: x1(3*nbod),v1(3*nbod),spin(3*nbod)
  real(double_precision) :: ngf(4,nbod),rphys(nbod),rcrit(nbod)
  external force
  !
  ! Local
  integer :: nv,niter,j,k,n
  real(double_precision) :: x(3*NMAX),v(3*NMAX),a(3*NMAX),a1(3*NMAX)
  real(double_precision) :: g(7,3*NMAX),b(7,3*NMAX),e(7,3*NMAX)
  real(double_precision) :: h(8),xc(8),vc(7),c(21),d(21),r(28),s(9)
  real(double_precision) :: q,q2,q3,q4,q5,q6,q7,temp,gk
  !
  !------------------------------------------------------------------------------
  !
  save h,xc,vc,c,d,r,b,e
  !
  ! Gauss-Radau spacings for substeps within a sequence, for the 15th order 
  ! integrator. The sum of the H values should be 3.733333333333333
  !
  data h/ 0.d0,.0562625605369221d0,.1802406917368924d0,.3526247171131696d0,.5471536263305554d0,.7342101772154105d0,&
       .8853209468390958d0,.9775206135612875d0/
  !
  ! Constant coefficients used in series expansions for X and V
  !  XC: 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72
  !  VC: 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8
  data xc/.5d0,.1666666666666667d0,.08333333333333333d0,.05d0,.03333333333333333d0,.02380952380952381d0,&
       .01785714285714286d0,.01388888888888889d0/
  data vc/.5d0,.3333333333333333d0,.25d0,.2d0,.1666666666666667d0,.1428571428571429d0,.125d0/
  !
  ! If this is first call to the subroutine, set values of the constant arrays
  ! (R = R21, R31, R32, R41, R42, R43 in Everhart's paper.)
  if (dtflag.eq.0) then
     n = 0
     do j = 2, 8
        do k = 1, j - 1
           n = n + 1
           r(n) = 1.d0 / (h(j) - h(k))
        end do
     end do
     !
     ! Constants to convert between B and G arrays (C = C21, C31, C32, C41, C42...)
     c(1) = - h(2)
     d(1) =   h(2)
     n = 1
     do j = 3, 7
        n = n + 1
        c(n) = -h(j) * c(n-j+2)
        d(n) =  h(2) * d(n-j+2)
        do k = 3, j - 1
           n = n + 1
           c(n) = c(n-j+1)  -  h(j) * c(n-j+2)
           d(n) = d(n-j+1)  +  h(k) * d(n-j+2)
        end do
        n = n + 1
        c(n) = c(n-j+1) - h(j)
        d(n) = d(n-j+1) + h(j)
     end do
     !
     dtflag = 1
  end if
  !
  nv = 3 * nbod
100 continue
  !
  ! If this is first call to subroutine since number/masses of objects changed
  ! do 6 iterations and initialize B, E arrays, otherwise do 2 iterations.
  if (dtflag.eq.1) then
     niter = 6
     do j = 4, nv
        do k = 1, 7
           b (k,j) = 0.d0
           e (k,j) = 0.d0
        end do
     end do
  else
     niter = 2
  end if
  !
  ! Calculate forces at the start of the sequence
  call force (time,jcen,nbod,nbig,mass,x1,v1,spin,rcrit,a1,stat,ngf,ngflag,opt,nce,ice,jce)
  !
  ! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
  do k = 4, nv
     g(1,k) = b(7,k)*d(16) + b(6,k)*d(11) + b(5,k)*d(7) + b(4,k)*d(4) + b(3,k)*d(2) + b(2,k)*d(1) + b(1,k)
     g(2,k) = b(7,k)*d(17) + b(6,k)*d(12) + b(5,k)*d(8) + b(4,k)*d(5) + b(3,k)*d(3) + b(2,k)
     g(3,k) = b(7,k)*d(18) + b(6,k)*d(13) + b(5,k)*d(9) + b(4,k)*d(6) + b(3,k)
     g(4,k) = b(7,k)*d(19) + b(6,k)*d(14) + b(5,k)*d(10) + b(4,k)
     g(5,k) = b(7,k)*d(20) + b(6,k)*d(15) + b(5,k)
     g(6,k) = b(7,k)*d(21) + b(6,k)
     g(7,k) = b(7,k)
  end do
  !
  !------------------------------------------------------------------------------
  !
  !  MAIN  LOOP  STARTS  HERE
  !
  ! For each iteration (six for first call to subroutine, two otherwise)...
  do n = 1, niter
     !
     ! For each substep within a sequence...
     do j = 2, 8
        !
        ! Calculate position predictors using Eqn. 9 of Everhart
        s(1) = t * h(j)
        s(2) = s(1) * s(1) * .5d0
        s(3) = s(2) * h(j) * .3333333333333333d0
        s(4) = s(3) * h(j) * .5d0
        s(5) = s(4) * h(j) * .6d0
        s(6) = s(5) * h(j) * .6666666666666667d0
        s(7) = s(6) * h(j) * .7142857142857143d0
        s(8) = s(7) * h(j) * .75d0
        s(9) = s(8) * h(j) * .7777777777777778d0
        !
        do k = 4, nv
           x(k) = s(9)*b(7,k) + s(8)*b(6,k) + s(7)*b(5,k) + s(6)*b(4,k) + s(5)*b(3,k) + s(4)*b(2,k) &
                + s(3)*b(1,k) + s(2)*a1(k)  + s(1)*v1(k) + x1(k)
        end do
        !
        ! If necessary, calculate velocity predictors too, from Eqn. 10 of Everhart
        if (ngflag.ne.0) then
           s(1) = t * h(j)
           s(2) = s(1) * h(j) * .5d0
           s(3) = s(2) * h(j) * .6666666666666667d0
           s(4) = s(3) * h(j) * .75d0
           s(5) = s(4) * h(j) * .8d0
           s(6) = s(5) * h(j) * .8333333333333333d0
           s(7) = s(6) * h(j) * .8571428571428571d0
           s(8) = s(7) * h(j) * .875d0
           !
           do k = 4, nv
              v(k) = s(8)*b(7,k) + s(7)*b(6,k) + s(6)*b(5,k) + s(5)*b(4,k) + s(4)*b(3,k)&
                   + s(3)*b(2,k) + s(2)*b(1,k) + s(1)*a1(k)  + v1(k)
           end do
        end if
        !
        ! Calculate forces at the current substep
        call force (time,jcen,nbod,nbig,mass,x,v,spin,rcrit,a,stat,ngf,ngflag,opt,nce,ice,jce)
        !
        ! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
        select case (j)
        case (2)
           do k = 4, nv
              temp = g(1,k)
              g(1,k) = (a(k) - a1(k)) * r(1)
              b(1,k) = b(1,k) + g(1,k) - temp
           end do
        case (3)
           do k = 4, nv
              temp = g(2,k)
              gk = a(k) - a1(k)
              g(2,k) = (gk*r(2) - g(1,k))*r(3)
              temp = g(2,k) - temp
              b(1,k) = b(1,k)  +  temp * c(1)
              b(2,k) = b(2,k)  +  temp
           end do
        case (4)
           do k = 4, nv
              temp = g(3,k)
              gk = a(k) - a1(k)
              g(3,k) = ((gk*r(4) - g(1,k))*r(5) - g(2,k))*r(6)
              temp = g(3,k) - temp
              b(1,k) = b(1,k)  +  temp * c(2)
              b(2,k) = b(2,k)  +  temp * c(3)
              b(3,k) = b(3,k)  +  temp
           end do
        case (5)
           do k = 4, nv
              temp = g(4,k)
              gk = a(k) - a1(k)
              g(4,k) = (((gk*r(7) - g(1,k))*r(8) - g(2,k))*r(9)- g(3,k))*r(10)
              temp = g(4,k) - temp
              b(1,k) = b(1,k)  +  temp * c(4)
              b(2,k) = b(2,k)  +  temp * c(5)
              b(3,k) = b(3,k)  +  temp * c(6)
              b(4,k) = b(4,k)  +  temp
           end do
        case (6)
           do k = 4, nv
              temp = g(5,k)
              gk = a(k) - a1(k)
              g(5,k) =  ((((gk*r(11) - g(1,k))*r(12) - g(2,k))*r(13)- g(3,k))*r(14) - g(4,k))*r(15)
              temp = g(5,k) - temp
              b(1,k) = b(1,k)  +  temp * c(7)
              b(2,k) = b(2,k)  +  temp * c(8)
              b(3,k) = b(3,k)  +  temp * c(9)
              b(4,k) = b(4,k)  +  temp * c(10)
              b(5,k) = b(5,k)  +  temp
           end do
        case (7)
           do k = 4, nv
              temp = g(6,k)
              gk = a(k) - a1(k)
              g(6,k) = (((((gk*r(16) - g(1,k))*r(17) - g(2,k))*r(18)- g(3,k))*r(19) - g(4,k))*r(20) - g(5,k))*r(21)
              temp = g(6,k) - temp
              b(1,k) = b(1,k)  +  temp * c(11)
              b(2,k) = b(2,k)  +  temp * c(12)
              b(3,k) = b(3,k)  +  temp * c(13)
              b(4,k) = b(4,k)  +  temp * c(14)
              b(5,k) = b(5,k)  +  temp * c(15)
              b(6,k) = b(6,k)  +  temp
           end do
        case (8)
           do k = 4, nv
              temp = g(7,k)
              gk = a(k) - a1(k)
              g(7,k) = ((((((gk*r(22) - g(1,k))*r(23) - g(2,k))*r(24) - g(3,k))*r(25) &
                   - g(4,k))*r(26) - g(5,k))*r(27) - g(6,k))*r(28)
              temp = g(7,k) - temp
              b(1,k) = b(1,k)  +  temp * c(16)
              b(2,k) = b(2,k)  +  temp * c(17)
              b(3,k) = b(3,k)  +  temp * c(18)
              b(4,k) = b(4,k)  +  temp * c(19)
              b(5,k) = b(5,k)  +  temp * c(20)
              b(6,k) = b(6,k)  +  temp * c(21)
              b(7,k) = b(7,k)  +  temp
           end do
        end select
     end do
  end do
  !
  !------------------------------------------------------------------------------
  !
  !  END  OF  MAIN  LOOP
  !
  ! Estimate suitable sequence size for the next call to subroutine (Eqs. 15, 16)
  temp = 0.d0
  do k = 4, nv
     temp = max( temp, abs( b(7,k) ) )
  end do
  temp = temp / (72.d0 * abs(t)**7)
  tdid = t
  if (temp.eq.0) then
     t = tdid * 1.4d0
  else
     t = sign( (tol/temp)**(1.d0/9.d0), tdid )
  end if
  !
  ! If sequence size for the first subroutine call is too big, go back and redo
  ! the sequence using a smaller size.
  if ((dtflag.eq.1).and.(abs(t/tdid).lt.1)) then
     t = t * .8d0
     goto 100
  end if
  !
  ! If new sequence size is much bigger than the current one, reduce it
  if (abs(t/tdid).gt.1.4d0) t = tdid * 1.4d0
  !
  ! Find new position and velocity values at end of the sequence (Eqs. 11, 12)
  temp = tdid * tdid
  do k = 4 , nv
     x1(k) = (xc(8)*b(7,k) + xc(7)*b(6,k) + xc(6)*b(5,k) + xc(5)*b(4,k) + xc(4)*b(3,k) + xc(3)*b(2,k)&
          + xc(2)*b(1,k) + xc(1)*a1(k))*temp + v1(k)*tdid + x1(k)
     !
     v1(k) = (vc(7)*b(7,k) + vc(6)*b(6,k) + vc(5)*b(5,k) + vc(4)*b(4,k) + vc(3)*b(3,k) + vc(2)*b(2,k)&
          + vc(1)*b(1,k) + a1(k))*tdid + v1(k)
  end do
  !
  ! Predict new B values to use at the start of the next sequence. The predicted
  ! values from the last call are saved as E. The correction, BD, between the
  ! actual and predicted values of B is applied in advance as a correction.
  q = t / tdid
  q2 = q  * q
  q3 = q  * q2
  q4 = q2 * q2
  q5 = q2 * q3
  q6 = q3 * q3
  q7 = q3 * q4
  !
  do k = 4, nv
     s(1) = b(1,k) - e(1,k)
     s(2) = b(2,k) - e(2,k)
     s(3) = b(3,k) - e(3,k)
     s(4) = b(4,k) - e(4,k)
     s(5) = b(5,k) - e(5,k)
     s(6) = b(6,k) - e(6,k)
     s(7) = b(7,k) - e(7,k)
     !
     ! Estimate B values for the next sequence (Eqs. 13 of Everhart).
     e(1,k) = q* (b(7,k)* 7.d0 + b(6,k)* 6.d0 + b(5,k)* 5.d0 + b(4,k)* 4.d0 + b(3,k)* 3.d0 + b(2,k)*2.d0 + b(1,k))
     e(2,k) = q2*(b(7,k)*21.d0 + b(6,k)*15.d0 + b(5,k)*10.d0 + b(4,k)* 6.d0 + b(3,k)* 3.d0 + b(2,k))
     e(3,k) = q3*(b(7,k)*35.d0 + b(6,k)*20.d0 + b(5,k)*10.d0 + b(4,k)*4.d0 + b(3,k))
     e(4,k) = q4*(b(7,k)*35.d0 + b(6,k)*15.d0 + b(5,k)*5.d0 + b(4,k))
     e(5,k) = q5*(b(7,k)*21.d0 + b(6,k)*6.d0 + b(5,k))
     e(6,k) = q6*(b(7,k)*7.d0 + b(6,k))
     e(7,k) = q7* b(7,k)
     !
     b(1,k) = e(1,k) + s(1)
     b(2,k) = e(2,k) + s(2)
     b(3,k) = e(3,k) + s(3)
     b(4,k) = e(4,k) + s(4)
     b(5,k) = e(5,k) + s(5)
     b(6,k) = e(6,k) + s(6)
     b(7,k) = e(7,k) + s(7)
  end do
  dtflag = 2
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mdt_ra15
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_ALL.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations on a set of NBOD bodies (of which NBIG are Big)
! due to Newtonian gravitational perturbations, post-Newtonian
! corrections (if required), cometary non-gravitational forces (if required)
! and user-defined forces (if required).
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
subroutine mfo_all (time,jcen,nbod,nbig,m,x,v,s,rcrit,a,stat,ngf,ngflag,opt,nce,ice,jce)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,ngflag,stat(nbod),opt(8),nce,ice(nce),jce(nce)
  real(double_precision) :: time,jcen(3),m(nbod),x(3,nbod),v(3,nbod),s(3,nbod)
  real(double_precision) :: a(3,nbod),ngf(4,nbod),rcrit(nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: acor(3,NMAX),acen(3)
  !
  !------------------------------------------------------------------------------
  !
  ! Newtonian gravitational forces
  call mfo_grav (nbod,nbig,m,x,v,a,stat)
  !
  ! Correct for oblateness of the central body
  if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
     call mfo_obl (jcen,nbod,m,x,acor,acen)
     do j = 2, nbod
        a(1,j) = a(1,j) + (acor(1,j) - acen(1))
        a(2,j) = a(2,j) + (acor(2,j) - acen(2))
        a(3,j) = a(3,j) + (acor(3,j) - acen(3))
     end do
  end if
  !
  ! Include non-gravitational (cometary jet) accelerations if necessary
  if (ngflag.eq.1.or.ngflag.eq.3) then
     call mfo_ngf (nbod,x,v,acor,ngf)
     do j = 2, nbod
        a(1,j) = a(1,j) + acor(1,j)
        a(2,j) = a(2,j) + acor(2,j)
        a(3,j) = a(3,j) + acor(3,j)
     end do
  end if
  !
  ! Include radiation pressure/Poynting-Robertson drag if necessary
  if (ngflag.eq.2.or.ngflag.eq.3) then
     call mfo_pr (nbod,nbig,m,x,v,acor,ngf)
     do j = 2, nbod
        a(1,j) = a(1,j) + acor(1,j)
        a(2,j) = a(2,j) + acor(2,j)
        a(3,j) = a(3,j) + acor(3,j)
     end do
  end if
  !
  ! Include post-Newtonian corrections if required
  if (opt(7).eq.1) then
     call mfo_pn (nbod,nbig,m,x,v,acor)
     do j = 2, nbod
        a(1,j) = a(1,j) + acor(1,j)
        a(2,j) = a(2,j) + acor(2,j)
        a(3,j) = a(3,j) + acor(3,j)
     end do
  end if
  !
  ! Include user-defined accelerations if required
  if (opt(8).eq.1) then
     call mfo_user (time,jcen,nbod,nbig,m,x,v,acor)
     do j = 2, nbod
        a(1,j) = a(1,j) + acor(1,j)
        a(2,j) = a(2,j) + acor(2,j)
        a(3,j) = a(3,j) + acor(3,j)
     end do
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_all
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_GRAV.FOR    (ErikSoft   3 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations on a set of NBOD bodies (NBIG of which are Big)
! due to gravitational perturbations by all the other bodies, except that
! Small bodies do not interact with one another.
!
! The positions and velocities are stored in arrays X, V with the format
! (x,y,z) and (vx,vy,vz) for each object in succession. The accelerations 
! are stored in the array A (ax,ay,az).
!
! N.B. All coordinates and velocities must be with respect to central body!!!!
! ===
!------------------------------------------------------------------------------
!
subroutine mfo_grav (nbod,nbig,m,x,v,a,stat)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod, nbig, stat(nbod)
  real(double_precision) :: m(nbod), x(3,nbod), v(3,nbod), a(3,nbod)
  !
  ! Local
  integer :: i, j
  real(double_precision) :: sx, sy, sz, dx, dy, dz, tmp1, tmp2, s_1, s2, s_3, r3(NMAX)
  !
  !------------------------------------------------------------------------------
  !
  sx = 0.d0
  sy = 0.d0
  sz = 0.d0
  do i = 2, nbod
     a(1,i) = 0.d0
     a(2,i) = 0.d0
     a(3,i) = 0.d0
     s2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
     s_1  = 1.d0 / sqrt(s2)
     r3(i) = s_1 * s_1 * s_1
  end do
  !
  do i = 2, nbod
     tmp1 = m(i) * r3(i)
     sx = sx  -  tmp1 * x(1,i)
     sy = sy  -  tmp1 * x(2,i)
     sz = sz  -  tmp1 * x(3,i)
  end do
  !
  ! Direct terms
  do i = 2, nbig
     do j = i + 1, nbod
        dx = x(1,j) - x(1,i)
        dy = x(2,j) - x(2,i)
        dz = x(3,j) - x(3,i)
        s2 = dx*dx + dy*dy + dz*dz
        s_1 = 1.d0 / sqrt(s2)
        s_3 = s_1 * s_1 * s_1
        tmp1 = s_3 * m(i)
        tmp2 = s_3 * m(j)
        a(1,j) = a(1,j)  -  tmp1 * dx
        a(2,j) = a(2,j)  -  tmp1 * dy
        a(3,j) = a(3,j)  -  tmp1 * dz
        a(1,i) = a(1,i)  +  tmp2 * dx
        a(2,i) = a(2,i)  +  tmp2 * dy
        a(3,i) = a(3,i)  +  tmp2 * dz
     end do
  end do
  !
  ! Indirect terms (add these on last to reduce roundoff error)
  do i = 2, nbod
     tmp1 = m(1) * r3(i)
     a(1,i) = a(1,i)  +  sx  -  tmp1 * x(1,i)
     a(2,i) = a(2,i)  +  sy  -  tmp1 * x(2,i)
     a(3,i) = a(3,i)  +  sz  -  tmp1 * x(3,i)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_grav
!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_DRCT.FOR    (ErikSoft   27 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates direct accelerations between bodies in the interaction part
! of the Hamiltonian of a symplectic integrator that partitions close
! encounter terms (e.g. hybrid symplectic algorithms or SyMBA).
! The routine calculates accelerations between all pairs of bodies with
! indices I >= I0.
!
!------------------------------------------------------------------------------
!
subroutine mfo_drct (i0,nbod,nbig,m,x,rcrit,a,stat)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: i0, nbod, nbig, stat(nbod)
  real(double_precision) :: m(nbod), x(3,nbod), a(3,nbod), rcrit(nbod)
  !
  ! Local
  integer :: i,j
  real(double_precision) :: dx,dy,dz,s,s_1,s2,s_3,rc,rc2,q,q2,q3,q4,q5,tmp2,faci,facj
  !
  !------------------------------------------------------------------------------
  !
  if (i0.le.0) i0 = 2
  !
  do i = i0, nbig
     do j = i + 1, nbod
        dx = x(1,j) - x(1,i)
        dy = x(2,j) - x(2,i)
        dz = x(3,j) - x(3,i)
        s2 = dx * dx  +  dy * dy  +  dz * dz
        rc = max(rcrit(i), rcrit(j))
        rc2 = rc * rc
        !
        if (s2.ge.rc2) then
           s_1 = 1.d0 / sqrt(s2)
           tmp2 = s_1 * s_1 * s_1
        else if (s2.le.0.01*rc2) then
           tmp2 = 0.d0
        else
           s_1 = 1.d0 / sqrt(s2)
           s   = 1.d0 / s_1
           s_3 = s_1 * s_1 * s_1
           q = (s - 0.1d0*rc) / (0.9d0 * rc)
           q2 = q  * q
           q3 = q  * q2
           q4 = q2 * q2
           q5 = q2 * q3
           tmp2 = (10.d0*q3 - 15.d0*q4 + 6.d0*q5) * s_3
        end if
        !
        faci = tmp2 * m(i)
        facj = tmp2 * m(j)
        a(1,j) = a(1,j)  -  faci * dx
        a(2,j) = a(2,j)  -  faci * dy
        a(3,j) = a(3,j)  -  faci * dz
        a(1,i) = a(1,i)  +  facj * dx
        a(2,i) = a(2,i)  +  facj * dy
        a(3,i) = a(3,i)  +  facj * dz
     end do
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_drct
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_HY.FOR    (ErikSoft   2 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations due to the Interaction part of the Hamiltonian 
! of a hybrid symplectic integrator for a set of NBOD bodies (NBIG of which 
! are Big), where Small bodies do not interact with one another.
!
!------------------------------------------------------------------------------
!
subroutine mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod, nbig, stat(nbod)
  real(double_precision) :: jcen(3), m(nbod), x(3,nbod), a(3,nbod), rcrit(nbod)
  !
  ! Local
  integer :: k
  real(double_precision) :: aobl(3,NMAX),acen(3)
  !
  !------------------------------------------------------------------------------
  !
  ! Initialize accelerations to zero
  do k = 1, nbod
     a(1,k) = 0.d0
     a(2,k) = 0.d0
     a(3,k) = 0.d0
  end do
  !
  ! Calculate direct terms
  call mfo_drct (2,nbod,nbig,m,x,rcrit,a,stat)
  !
  ! Add accelerations due to oblateness of the central body
  if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
     call mfo_obl (jcen,nbod,m,x,aobl,acen)
     do k = 2, nbod
        a(1,k) = a(1,k) + aobl(1,k) - acen(1)
        a(2,k) = a(2,k) + aobl(2,k) - acen(2)
        a(3,k) = a(3,k) + aobl(3,k) - acen(3)
     end do
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_hy
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_HKCE.FOR    (ErikSoft   27 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations due to the Keplerian part of the Hamiltonian 
! of a hybrid symplectic integrator, when close encounters are taking place,
! for a set of NBOD bodies (NBIG of which are Big). Note that Small bodies
! do not interact with one another.
!
!------------------------------------------------------------------------------
!
subroutine mfo_hkce (time,jcen,nbod,nbig,m,x,v,spin,rcrit,a,stat,ngf,ngflag,opt,nce,ice,jce)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod,nbig,stat(nbod),ngflag,opt(8),nce,ice(nce),jce(nce)
  real(double_precision) :: time,jcen(3),rcrit(nbod),ngf(4,nbod),m(nbod)
  real(double_precision) :: x(3,nbod),v(3,nbod),a(3,nbod),spin(3,nbod)
  !
  ! Local
  integer :: i, j, k
  real(double_precision) :: tmp2,dx,dy,dz,s,s_1,s2,s_3,faci,facj,rc,rc2,q,q2,q3,q4,q5
  !
  !------------------------------------------------------------------------------
  !
  ! Initialize accelerations
  do j = 1, nbod
     a(1,j) = 0.d0
     a(2,j) = 0.d0
     a(3,j) = 0.d0
  end do
  !
  ! Direct terms
  do k = 1, nce
     i = ice(k)
     j = jce(k)
     dx = x(1,j) - x(1,i)
     dy = x(2,j) - x(2,i)
     dz = x(3,j) - x(3,i)
     s2 = dx * dx  +  dy * dy  +  dz * dz
     rc = max (rcrit(i), rcrit(j))
     rc2 = rc * rc
     !
     if (s2.lt.rc2) then
        s_1 = 1.d0 / sqrt(s2)
        s_3 = s_1 * s_1 * s_1
        if (s2.le.0.01*rc2) then
           tmp2 = s_3
        else
           s = 1.d0 / s_1
           q = (s - 0.1d0*rc) / (0.9d0 * rc)
           q2 = q * q
           q3 = q * q2
           q4 = q2 * q2
           q5 = q2 * q3
           tmp2 = (1.d0 - 10.d0*q3 + 15.d0*q4 - 6.d0*q5) * s_3
        end if
        !
        faci = tmp2 * m(i)
        facj = tmp2 * m(j)
        a(1,j) = a(1,j)  -  faci * dx
        a(2,j) = a(2,j)  -  faci * dy
        a(3,j) = a(3,j)  -  faci * dz
        a(1,i) = a(1,i)  +  facj * dx
        a(2,i) = a(2,i)  +  facj * dy
        a(3,i) = a(3,i)  +  facj * dz
     end if
  end do
  !
  ! Solar terms
  do i = 2, nbod
     s2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
     s_1 = 1.d0 / sqrt(s2)
     tmp2 = m(1) * s_1 * s_1 * s_1
     a(1,i) = a(1,i)  -  tmp2 * x(1,i)
     a(2,i) = a(2,i)  -  tmp2 * x(2,i)
     a(3,i) = a(3,i)  -  tmp2 * x(3,i)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_hkce
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_MVS.FOR    (ErikSoft   2 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations on a set of NBOD bodies (of which NBIG are Big)
! due to gravitational perturbations by all the other bodies.
! This routine is designed for use with a mixed-variable symplectic
! integrator using Jacobi coordinates.
!
! Based upon routines from Levison and Duncan's SWIFT integrator.
!
!------------------------------------------------------------------------------
!
subroutine mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod, nbig, stat(nbod)
  real(double_precision) :: jcen(3), m(nbod), x(3,nbod), xj(3,nbod), a(3,nbod)
  !
  ! Local
  integer :: i,j,k,k1
  real(double_precision) :: fac0,fac1,fac12,fac2,minside,dx,dy,dz,s_1,s2,s_3,faci,facj
  real(double_precision) :: a0(3),a0tp(3),a1(3,NMAX),a2(3,NMAX),a3(3,NMAX),aobl(3,NMAX)
  real(double_precision) :: r,r2,r3,rj,rj2,rj3,q,q2,q3,q4,q5,q6,q7,acen(3)
  !
  !------------------------------------------------------------------------------
  !
  ! Initialize variables
  a0(1) = 0.d0
  a0(2) = 0.d0
  a0(3) = 0.d0
  a1(1,2) = 0.d0
  a1(2,2) = 0.d0
  a1(3,2) = 0.d0
  a2(1,2) = 0.d0
  a2(2,2) = 0.d0
  a2(3,2) = 0.d0
  minside = 0.d0
  !
  ! Calculate acceleration terms
  do k = 3, nbig
     k1 = k - 1
     minside = minside + m(k1)
     r2   = x(1,k)  * x(1,k)  +  x(2,k) * x(2,k)  +  x(3,k) * x(3,k)
     rj2  = xj(1,k) * xj(1,k) + xj(2,k) * xj(2,k) + xj(3,k) * xj(3,k)
     r  = 1.d0 / sqrt(r2)
     rj = 1.d0 / sqrt(rj2)
     r3  = r  * r  * r
     rj3 = rj * rj * rj
     !
     fac0 = m(k) * r3
     fac12 = m(1) * rj3
     fac2 = m(k) * fac12 / (minside + m(1))
     q = (r2 - rj2) * .5d0 / rj2
     q2 = q  * q
     q3 = q  * q2
     q4 = q2 * q2
     q5 = q2 * q3
     q6 = q3 * q3
     q7 = q3 * q4
     fac1 = 402.1875d0*q7 - 187.6875d0*q6 + 86.625d0*q5   - 39.375d0*q4 + 17.5d0*q3 - 7.5d0*q2 + 3.d0*q - 1.d0
     !
     ! Add to A0 term
     a0(1) = a0(1)  -  fac0 * x(1,k)
     a0(2) = a0(2)  -  fac0 * x(2,k)
     a0(3) = a0(3)  -  fac0 * x(3,k)
     !
     ! Calculate A1 for this body
     a1(1,k) = fac12 * (xj(1,k) + fac1*x(1,k))
     a1(2,k) = fac12 * (xj(2,k) + fac1*x(2,k))
     a1(3,k) = fac12 * (xj(3,k) + fac1*x(3,k))
     !
     ! Calculate A2 for this body
     a2(1,k) = a2(1,k1)  +  fac2 * xj(1,k)
     a2(2,k) = a2(2,k1)  +  fac2 * xj(2,k)
     a2(3,k) = a2(3,k1)  +  fac2 * xj(3,k)
  end do
  !
  r2   = x(1,2)  * x(1,2)  +  x(2,2) * x(2,2)  +  x(3,2) * x(3,2)
  r  = 1.d0 / sqrt(r2)
  r3  = r  * r  * r
  fac0 = m(2) * r3
  a0tp(1) = a0(1)  -  fac0 * x(1,2)
  a0tp(2) = a0(2)  -  fac0 * x(2,2)
  a0tp(3) = a0(3)  -  fac0 * x(3,2)
  !
  ! Calculate A3 (direct terms)
  do k = 2, nbod
     a3(1,k) = 0.d0
     a3(2,k) = 0.d0
     a3(3,k) = 0.d0
  end do
  do i = 2, nbig
     do j = i + 1, nbig
        dx = x(1,j) - x(1,i)
        dy = x(2,j) - x(2,i)
        dz = x(3,j) - x(3,i)
        s2 = dx*dx + dy*dy + dz*dz
        s_1 = 1.d0 / sqrt(s2)
        s_3 = s_1 * s_1 * s_1
        faci = m(i) * s_3
        facj = m(j) * s_3
        a3(1,j) = a3(1,j)  -  faci * dx
        a3(2,j) = a3(2,j)  -  faci * dy
        a3(3,j) = a3(3,j)  -  faci * dz
        a3(1,i) = a3(1,i)  +  facj * dx
        a3(2,i) = a3(2,i)  +  facj * dy
        a3(3,i) = a3(3,i)  +  facj * dz
     end do
     !
     do j = nbig + 1, nbod
        dx = x(1,j) - x(1,i)
        dy = x(2,j) - x(2,i)
        dz = x(3,j) - x(3,i)
        s2 = dx*dx + dy*dy + dz*dz
        s_1 = 1.d0 / sqrt(s2)
        s_3 = s_1 * s_1 * s_1
        faci = m(i) * s_3
        a3(1,j) = a3(1,j)  -  faci * dx
        a3(2,j) = a3(2,j)  -  faci * dy
        a3(3,j) = a3(3,j)  -  faci * dz
     end do
  end do
  !
  ! Big-body accelerations
  do k = 2, nbig
     a(1,k) = a0(1) + a1(1,k) + a2(1,k) + a3(1,k)
     a(2,k) = a0(2) + a1(2,k) + a2(2,k) + a3(2,k)
     a(3,k) = a0(3) + a1(3,k) + a2(3,k) + a3(3,k)
  end do
  !
  ! Small-body accelerations
  do k = nbig + 1, nbod
     a(1,k) = a0tp(1) + a3(1,k)
     a(2,k) = a0tp(2) + a3(2,k)
     a(3,k) = a0tp(3) + a3(3,k)
  end do
  !
  ! Correct for oblateness of the central body
  if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
     call mfo_obl (jcen,nbod,m,x,aobl,acen)
     do k = 2, nbod
        a(1,k) = a(1,k) + (aobl(1,k) - acen(1))
        a(2,k) = a(2,k) + (aobl(2,k) - acen(2))
        a(3,k) = a(3,k) + (aobl(3,k) - acen(3))
     end do
  end if
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_mvs
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_NGF.FOR    (ErikSoft  29 November 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations on a set of NBOD bodies due to cometary
! non-gravitational jet forces. The positions and velocities are stored in
! arrays X, V with the format (x,y,z) and (vx,vy,vz) for each object in
! succession. The accelerations are stored in the array A (ax,ay,az). The
! non-gravitational accelerations follow a force law described by Marsden
! et al. (1973) Astron. J. 211-225, with magnitude determined by the
! parameters NGF(1,2,3) for each object.
!
! N.B. All coordinates and velocities must be with respect to central body!!!!
! ===
!------------------------------------------------------------------------------
!
subroutine mfo_ngf (nbod,x,v,a,ngf)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod
  real(double_precision) :: x(3,nbod), v(3,nbod), a(3,nbod), ngf(4,nbod)
  !
  ! Local
  integer :: j
  real(double_precision) :: r2,r,rv,q,g,tx,ty,tz,nx,ny,nz,a1,a2,a3
  !
  !------------------------------------------------------------------------------
  !
  do j = 2, nbod
     r2 = x(1,j)*x(1,j) + x(2,j)*x(2,j) +x(3,j)*x(3,j)
     !
     ! Only calculate accelerations if body is close to the Sun (R < 9.36 AU), 
     ! or if the non-gravitational force parameters are exceptionally large.
     if (r2.lt.88.d0.or.abs(ngf(1,j)).gt.1d-7.or.abs(ngf(2,j)).gt.1d-7.or.abs(ngf(3,j)).gt.1d-7) then
        r = sqrt(r2)
        rv = x(1,j)*v(1,j) + x(2,j)*v(2,j) + x(3,j)*v(3,j)
        !
        ! Calculate Q = R / R0, where R0 = 2.808 AU
        q = r * .3561253561253561d0
        g = .111262d0 * q**(-2.15d0) * (1.d0+q**5.093d0)**(-4.6142d0)
        !
        ! Within-orbital-plane transverse vector components
        tx = r2*v(1,j) - rv*x(1,j)
        ty = r2*v(2,j) - rv*x(2,j)
        tz = r2*v(3,j) - rv*x(3,j)
        !
        ! Orbit-normal vector components
        nx = x(2,j)*v(3,j) - x(3,j)*v(2,j)
        ny = x(3,j)*v(1,j) - x(1,j)*v(3,j)
        nz = x(1,j)*v(2,j) - x(2,j)*v(1,j)
        !
        ! Multiplication factors
        a1 = ngf(1,j) * g / r
        a2 = ngf(2,j) * g / sqrt(tx*tx + ty*ty + tz*tz)
        a3 = ngf(3,j) * g / sqrt(nx*nx + ny*ny + nz*nz)
        !
        ! X,Y and Z components of non-gravitational acceleration
        a(1,j) = a1*x(1,j) + a2*tx + a3*nx
        a(2,j) = a1*x(2,j) + a2*ty + a3*ny
        a(3,j) = a1*x(3,j) + a2*tz + a3*nz
     else
        a(1,j) = 0.d0
        a(2,j) = 0.d0
        a(3,j) = 0.d0
     end if
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_ngf
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_OBL.FOR    (ErikSoft   2 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates barycentric accelerations of NBOD bodies due to oblateness of
! the central body. Also returns the corresponding barycentric acceleration
! of the central body.
!
! N.B. All coordinates must be with respect to the central body!!!!
! ===
!------------------------------------------------------------------------------
!
subroutine mfo_obl (jcen,nbod,m,x,a,acen)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod
  real(double_precision) :: jcen(3), m(nbod), x(3,nbod), a(3,nbod), acen(3)
  !
  ! Local
  integer :: i
  real(double_precision) :: jr2,jr4,jr6,r2,r_1,r_2,r_3,u2,u4,u6,tmp1,tmp2,tmp3,tmp4
  !
  !------------------------------------------------------------------------------
  !
  acen(1) = 0.d0
  acen(2) = 0.d0
  acen(3) = 0.d0
  !
  do i = 2, nbod
     !
     ! Calculate barycentric accelerations on the objects
     r2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
     r_1 = 1.d0 / sqrt(r2)
     r_2 = r_1 * r_1
     r_3 = r_2 * r_1
     jr2 = jcen(1) * r_2
     jr4 = jcen(2) * r_2 * r_2
     jr6 = jcen(3) * r_2 * r_2 * r_2
     u2 = x(3,i) * x(3,i) * r_2
     u4 = u2 * u2
     u6 = u4 * u2
     !
     tmp1 = m(1) * r_3
     tmp2 = jr2*(7.5d0*u2 - 1.5d0) + jr4*(39.375d0*u4 - 26.25d0*u2 + 1.875d0) + &
          jr6*(187.6875d0*u6 -216.5625d0*u4 +59.0625d0*u2 -2.1875d0)
     tmp3 = jr2*3.d0 + jr4*(17.5d0*u2 - 7.5d0) + jr6*(86.625d0*u4 - 78.75d0*u2 + 13.125d0)
     !
     a(1,i) = x(1,i) * tmp1 * tmp2
     a(2,i) = x(2,i) * tmp1 * tmp2
     a(3,i) = x(3,i) * tmp1 * (tmp2 - tmp3)
     !
     ! Calculate barycentric accelerations on the central body
     tmp4 = m(i) / m(1)
     acen(1) = acen(1)  -  tmp4 * a(1,i)
     acen(2) = acen(2)  -  tmp4 * a(2,i)
     acen(3) = acen(3)  -  tmp4 * a(3,i)
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_obl
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_PN.FOR    (ErikSoft   3 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! ****** To be completed at a later date ******
!
! Calculates post-Newtonian relativistic corrective accelerations for a set
! of NBOD bodies (NBIG of which are Big).
!
! This routine should not be called from the symplectic algorithm MAL_MVS
! or the conservative Bulirsch-Stoer algorithm MAL_BS2.
!
! N.B. All coordinates and velocities must be with respect to central body!!!!
! ===
!------------------------------------------------------------------------------
!
subroutine mfo_pn (nbod,nbig,m,x,v,a)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod, nbig
  real(double_precision) :: m(nbod), x(3,nbod), v(3,nbod), a(3,nbod)
  !
  ! Local
  integer :: j
  !
  !------------------------------------------------------------------------------
  !
  do j = 1, nbod
     a(1,j) = 0.d0
     a(2,j) = 0.d0
     a(3,j) = 0.d0
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_pn
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_PR.FOR    (ErikSoft   3 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! ****** To be completed at a later date ******
!
! Calculates radiation pressure and Poynting-Robertson drag for a set
! of NBOD bodies (NBIG of which are Big).
!
! This routine should not be called from the symplectic algorithm MAL_MVS
! or the conservative Bulirsch-Stoer algorithm MAL_BS2.
!
! N.B. All coordinates and velocities must be with respect to central body!!!!
! ===
!------------------------------------------------------------------------------
!
subroutine mfo_pr (nbod,nbig,m,x,v,a,ngf)
  !
  use physical_constant
  use mercury_constant
  use types_numeriques

  implicit none

  !
  ! Input/Output
  integer :: nbod, nbig
  real(double_precision) :: m(nbod), x(3,nbod), v(3,nbod), a(3,nbod), ngf(4,nbod)
  !
  ! Local
  integer :: j
  !
  !------------------------------------------------------------------------------
  !
  do j = 1, nbod
     a(1,j) = 0.d0
     a(2,j) = 0.d0
     a(3,j) = 0.d0
  end do
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine mfo_pr

