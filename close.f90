!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      CLOSE6.FOR    (ErikSoft   5 June 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Makes output files containing details of close encounters that occurred
! during an integration using Mercury6 or higher.
!
! The user specifies the names of the required objects in the file close.in
!
!------------------------------------------------------------------------------
!

program close
  use physical_constant
  use mercury_constant
  use types_numeriques
  use ascii_conversion
  use orbital_elements
  use mercury_outputs
  use utilities

  implicit none

  !
  integer :: itmp,i,j,k,l,iclo,jclo,precision,lenin
  integer :: nmaster,nopen,nwait,nbig,nsml,nsub,lim(2,100)
  integer :: year,month,timestyle,line_num,lenhead,lmem(NMESS)
  integer :: nchar,algor,allflag,firstflag,ninfile
  integer :: unit(NMAX),master_unit(NMAX)
  real(double_precision) :: time,t0,t1,rmax,rcen,rfac,dclo,mcen,jcen(3)
  real(double_precision) :: fr,theta,phi,fv,vtheta,vphi,gm
  real(double_precision) :: x1(3),x2(3),v1(3),v2(3),m(NMAX)
  real(double_precision) :: a1,a2,e1,e2,i1,i2,p1,p2,n1,n2,l1,l2,q1,q2
  logical test
  character*250 string,fout,header,infile(50)
  character*80 mem(NMESS),cc,c(NMAX)
  character*8 master_id(NMAX),id(NMAX)
  character*5 fin
  character*1 check,style,type,c1
  !
  !------------------------------------------------------------------------------
  !
  allflag = 0
  !
  ! Read in output messages
  inquire (file='message.in', exist=test)
  if (.not.test) then
     write (*,'(/,2a)') ' ERROR: This file is needed to continue: ',' message.in'
     stop
  end if
  open (14, file='message.in', status='old')
10 continue
  read (14,'(i3,1x,i2,1x,a80)',end=20) j,lmem(j),mem(j)
  goto 10
20 close (14)
  !
  ! Open file containing parameters for this programme
  inquire (file='close.in', exist=test)
  if (test) then
     open (10, file='close.in', status='old')
  else
     call mio_err (6,mem(81),lmem(81),mem(88),lmem(88),' ',1,'close.in',9)
  end if
  !
  ! Read number of input files
30 read (10,'(a250)') string
  if (string(1:1).eq.')') goto 30
  call mio_spl (250,string,nsub,lim)
  read (string(lim(1,nsub):lim(2,nsub)),*) ninfile
  !
  ! Make sure all the input files exist
  do j = 1, ninfile
40   read (10,'(a250)') string
     if (string(1:1).eq.')') goto 40
     call mio_spl (250,string,nsub,lim)
     infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
     inquire (file=infile(j), exist=test)
     if (.not.test) call mio_err (6,mem(81),lmem(81),mem(88),lmem(88),' ',1,infile(j),80)
  end do
  !
  ! Read parameters used by this programme
  timestyle = 1
  do j = 1, 2
50   read (10,'(a250)') string
     if (string(1:1).eq.')') goto 50
     call mio_spl (250,string,nsub,lim)
     c1 = string(lim(1,nsub):lim(2,nsub))
     if (j.eq.1.and.(c1.eq.'d'.or.c1.eq.'D')) timestyle = 0
     if (j.eq.2.and.(c1.eq.'y'.or.c1.eq.'Y')) timestyle = timestyle+2
  end do
  !
  ! Read in the names of the objects for which orbital elements are required
  nopen = 0
  nwait = 0
  nmaster = 0
  call get_clo_format (timestyle,fout,header,lenhead)
60 continue
  read (10,'(a250)',end=70) string
  call mio_spl (250,string,nsub,lim)
  if (string(1:1).eq.')'.or.lim(1,1).eq.-1) goto 60
  !
  ! Either open an aei file for this object or put it on the waiting list
  nmaster = nmaster + 1
  itmp = min(7,lim(2,1)-lim(1,1))
  master_id(nmaster)='        '
  master_id(nmaster)(1:itmp+1) = string(lim(1,1):lim(1,1)+itmp)
  if (nopen.lt.NFILES) then
     nopen = nopen + 1
     master_unit(nmaster) = 10 + nopen
     call mio_clo (master_id(nmaster),master_unit(nmaster),  header,lenhead,mem,lmem)
  else
     nwait = nwait + 1
     master_unit(nmaster) = -2
  end if
  goto 60
  !
70 continue
  ! If no objects are listed in CLOSE.IN assume that all objects are required
  if (nopen.eq.0) allflag = 1
  close (10)
  !
  !------------------------------------------------------------------------------
  !
  !  LOOP  OVER  EACH  INPUT  FILE  CONTAINING  INTEGRATION  DATA
  !
90 continue
  firstflag = 0
  do i = 1, ninfile
     line_num = 0
     open (10, file=infile(i), status='old')
     !
     ! Loop over each time slice
100  continue
     line_num = line_num + 1
     read (10,'(3a1)',end=900,err=666) check,style,type
     line_num = line_num - 1
     backspace 10
     !
     ! Check if this is an old style input file
     if (ichar(check).eq.12.and.(style.eq.'0'.or.style.eq.'1'.or.style.eq.'2'.or.style.eq.'3'.or.style.eq.'4')) then
        write (*,'(/,2a)') ' ERROR: This is an old style data file','        Try running m_close5.for instead.'
        stop
     end if
     if (ichar(check).ne.12) goto 666
     !
     !------------------------------------------------------------------------------
     !
     !  IF  SPECIAL  INPUT,  READ  TIME,  PARAMETERS,  NAMES,  MASSES  ETC.
     !
     if (type.eq.'a') then
        line_num = line_num + 1
        read (10,'(3x,i2,a62,i1)') algor,cc(1:62),precision
        !
        ! Decompress the time, number of objects, central mass and J components etc.
        time = mio_c2fl (cc(1:8))
        if (firstflag.eq.0) then
           t0 = time
           firstflag = 1
        end if
        nbig = int(.5d0 + mio_c2re(cc(9:16), 0.d0, 11239424.d0, 3))
        nsml = int(.5d0 + mio_c2re(cc(12:19),0.d0, 11239424.d0, 3))
        mcen = mio_c2fl (cc(15:22)) * K2
        jcen(1) = mio_c2fl (cc(23:30))
        jcen(2) = mio_c2fl (cc(31:38))
        jcen(3) = mio_c2fl (cc(39:46))
        rcen = mio_c2fl (cc(47:54))
        rmax = mio_c2fl (cc(55:62))
        rfac = log10 (rmax / rcen)
        !
        ! Read in strings containing compressed data for each object
        do j = 1, nbig + nsml
           line_num = line_num + 1
           read (10,'(a)',err=666) c(j)(1:51)
        end do
        !
        ! Create input format list
        if (precision.eq.1) nchar = 2
        if (precision.eq.2) nchar = 4
        if (precision.eq.3) nchar = 7
        lenin = 3  +  6 * nchar
        fin(1:5) = '(a00)'
        write (fin(3:4),'(i2)') lenin
        !
        ! For each object decompress its name, code number, mass, spin and density
        do j = 1, nbig + nsml
           k = int(.5d0 + mio_c2re(c(j)(1:8),0.d0,11239424.d0,3))
           id(k) = c(j)(4:11)
           m(k)  = mio_c2fl (c(j)(12:19)) * K2
           !
           ! Find the object on the master list
           unit(k) = 0
           do l = 1, nmaster
              if (id(k).eq.master_id(l)) unit(k) = master_unit(l)
           end do
           !
           ! If object is not on the master list, add it to the list now
           if (unit(k).eq.0) then
              nmaster = nmaster + 1
              master_id(nmaster) = id(k)
              !
              ! Either open an aei file for this object or put it on the waiting list
              if (allflag.eq.1) then
                 if (nopen.lt.NFILES) then
                    nopen = nopen + 1
                    master_unit(nmaster) = 10 + nopen
                    call mio_clo (master_id(nmaster),master_unit(nmaster),header,lenhead,mem,lmem)
                 else
                    nwait = nwait + 1
                    master_unit(nmaster) = -2
                 end if
              else
                 master_unit(nmaster) = -1
              end if
              unit(k) = master_unit(nmaster)
           end if
        end do
        !
        !------------------------------------------------------------------------------
        !
        !  IF  NORMAL  INPUT,  READ  COMPRESSED  DATA  ON  THE  CLOSE  ENCOUNTER
        !
     else if (type.eq.'b') then
        line_num = line_num + 1
        read (10,'(3x,a70)',err=666) cc(1:70)
        !
        ! Decompress time, distance and orbital variables for each object
        time = mio_c2fl (cc(1:8))
        iclo = int(.5d0 + mio_c2re(cc(9:16),  0.d0, 11239424.d0, 3))
        jclo = int(.5d0 + mio_c2re(cc(12:19), 0.d0, 11239424.d0, 3))
        if (iclo.gt.NMAX.or.jclo.gt.NMAX) then
           write (*,'(/,2a)') mem(81)(1:lmem(81)),mem(90)(1:lmem(90))
           stop
        end if
        dclo = mio_c2fl (cc(15:22))
        fr     = mio_c2re (cc(23:30), 0.d0, rfac,  4)
        theta  = mio_c2re (cc(27:34), 0.d0, PI,    4)
        phi    = mio_c2re (cc(31:38), 0.d0, TWOPI, 4)
        fv     = mio_c2re (cc(35:42), 0.d0, 1.d0,  4)
        vtheta = mio_c2re (cc(39:46), 0.d0, PI,    4)
        vphi   = mio_c2re (cc(43:50), 0.d0, TWOPI, 4)
        call mco_ov2x (rcen,rmax,mcen,m(iclo),fr,theta,phi,fv,vtheta,vphi,x1(1),x1(2),x1(3),v1(1),v1(2),v1(3))
        !
        fr     = mio_c2re (cc(47:54), 0.d0, rfac,  4)
        theta  = mio_c2re (cc(51:58), 0.d0, PI,    4)
        phi    = mio_c2re (cc(55:62), 0.d0, TWOPI, 4)
        fv     = mio_c2re (cc(59:66), 0.d0, 1.d0,  4)
        vtheta = mio_c2re (cc(63:70), 0.d0, PI,    4)
        vphi   = mio_c2re (cc(67:74), 0.d0, TWOPI, 4)
        call mco_ov2x (rcen,rmax,mcen,m(jclo),fr,theta,phi,fv,vtheta,vphi,x2(1),x2(2),x2(3),v2(1),v2(2),v2(3))
        !
        ! Convert to Keplerian elements
        gm = mcen + m(iclo)
        call mco_x2el (gm,x1(1),x1(2),x1(3),v1(1),v1(2),v1(3),q1,e1,i1,p1,n1,l1)
        a1 = q1 / (1.d0 - e1)
        gm = mcen + m(jclo)
        call mco_x2el (gm,x2(1),x2(2),x2(3),v2(1),v2(2),v2(3),q2,e2,i2,p2,n2,l2)
        a2 = q2 / (1.d0 - e2)
        i1 = i1 * RAD2DEG
        i2 = i2 * RAD2DEG
        !
        ! Convert time to desired format
        if (timestyle.eq.0) t1 = time
        if (timestyle.eq.1) call mio_jd2y (time,year,month,t1)
        if (timestyle.eq.2) t1 = time - t0
        if (timestyle.eq.3) t1 = (time - t0) / 365.25d0
        !
        ! Write encounter details to appropriate files
        if (timestyle.eq.1) then
           if (unit(iclo).ge.10) write (unit(iclo),fout) year,month,t1,id(jclo),dclo,a1,e1,i1,a2,e2,i2
           !
           if (unit(jclo).ge.10) write (unit(jclo),fout) year,month,t1,id(iclo),dclo,a2,e2,i2,a1,e1,i1
        else
           if (unit(iclo).ge.10) write (unit(iclo),fout) t1,id(jclo),dclo,a1,e1,i1,a2,e2,i2
           if (unit(jclo).ge.10) write (unit(jclo),fout) t1,id(iclo),dclo,a2,e2,i2,a1,e1,i1
        end if
        !
        !------------------------------------------------------------------------------
        !
        !  IF  TYPE  IS  NOT  'a'  OR  'b',  THE  INPUT  FILE  IS  CORRUPTED
        !
     else
        goto 666
     end if
     !
     ! Move on to the next time slice
     goto 100
     !
     ! If input file is corrupted, try to continue from next uncorrupted time slice
666  continue
     write (*,'(2a,/,a,i10)') mem(121)(1:lmem(121)),infile(i)(1:60),mem(104)(1:lmem(104)),line_num
     c1 = ' '
     do while (ichar(c1).ne.12)
        line_num = line_num + 1
        read (10,'(a1)',end=900) c1
     end do
     line_num = line_num - 1
     backspace 10
     !
     ! Move on to the next file containing close encounter data
900  continue
     close (10)
  end do
  !
  ! Close clo files
  do j = 1, nopen
     close (10+j)
  end do
  nopen = 0
  !
  ! If some objects remain on waiting list, read through input files again
  if (nwait.gt.0) then
     do j = 1, nmaster
        if (master_unit(j).ge.10) master_unit(j) = -1
        if (master_unit(j).eq.-2.and.nopen.lt.NFILES) then
           nopen = nopen + 1
           nwait = nwait - 1
           master_unit(j) = 10 + nopen
           call mio_clo (master_id(j),master_unit(j),header,lenhead,mem,lmem)
        end if
     end do
     goto 90
  end if
  !
end program close

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      get_clo_format.FOR    (ErikSoft  30 November 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! This routine gives the header of the output .clo file regarding the timestyle used.
! The len of the header, and the output format is also returned.
!
!------------------------------------------------------------------------------
!
subroutine get_clo_format (timestyle,fout,header,lenhead)

  implicit none

  !
  ! Input
  integer,intent(in) :: timestyle
  
  ! Output
  integer,intent(out) :: lenhead ! length of the header
  character(len=250), intent(out) :: fout,& ! The output format for the .clo file
                                     header ! The header of the .clo file
  !
  !------------------------------------------------------------------------------
  !
  select case (timestyle)
    case (0,2)
      header(1:19) = '    Time (days)    '
      header(20:58) = '  Object   dmin (AU)     a1       e1    '
      header(59:90) = '   i1       a2       e2       i2'
      lenhead = 90
      fout = '(1x,f18.5,1x,a8,1x,f10.8,2(1x,f9.4,1x,f8.6,1x,f7.3))'
    case (1)
      header(1:23) = '     Year/Month/Day    '
      header(24:62) = '  Object   dmin (AU)     a1       e1    '
      header(63:94) = '   i1       a2       e2       i2'
      lenhead = 94
      fout(1:37) = '(1x,i10,1x,i2,1x,f8.5,1x,a8,1x,f10.8,'
      fout(38:64) = '2(1x,f9.4,1x,f8.6,1x,f7.3))'
    case default
      header(1:19) = '    Time (years)   '
      header(20:58) = '  Object   dmin (AU)     a1       e1    '
      header(59:90) = '   i1       a2       e2       i2'
      fout = '(1x,f18.7,1x,a8,1x,f10.8,2(1x,f9.4,1x,f8.6,1x,f7.3))'
      lenhead = 90
  end select
  !
  !------------------------------------------------------------------------------
  !
  return
end subroutine get_clo_format