!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      ELEMENT6.FOR    (ErikSoft   5 June 2001)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers

! Makes output files containing Keplerian orbital elements from data created
! by Mercury6 and higher.

! The user specifies the names of the required objects in the file elements.in
! See subroutine get_aei_format for the identities of each element in the EL array
! e.g. el(1)=a, el(2)=e etc.

!------------------------------------------------------------------------------


program element
  use types_numeriques
  use physical_constant
  use mercury_constant
  use ascii_conversion
  use orbital_elements
  use utilities
  use mercury_outputs
  use mercury_globals

  use system_properties
  use algo_mvs, only : mco_h2j


  implicit none


  integer :: itmp,i,j,k,l,precision,lenin
  integer :: nmaster,nopen,nwait,nbig,nsml,nbod,nsub,lim(2,100)
  integer :: year,month,timestyle,line_num,lenhead
  integer :: nchar,centre,allflag,firstflag,ninfile,nel,iel(22)
  integer :: nbod1,nbig1
  real(double_precision) :: time,teval,t0,t1,tprevious,rcen,rfac,rhocgs,temp
  real(double_precision) :: mcen,jcen(3), s(3)
  real(double_precision) :: fr,theta,phi,fv,vtheta,vphi,gm
  logical test
  character(len=250) :: string,fout,header,infile(50)
  character(len=80) :: cc
  character(len=5) :: fin
  character(len=1) :: check,style,type,c1
  character(len=2) :: c2
  
  ! fake variables to be able to use mco_iden
  integer :: ngflag !< do any bodies experience non-grav. forces?
!!\n                            ( 0 = no non-grav forces)
!!\n                              1 = cometary jets only
!!\n                              2 = radiation pressure/P-R drag only
!!\n                              3 = both
  real(double_precision), dimension(1,1) :: ngf
  
  integer :: error ! to store error when we allocate
  integer, dimension(:), allocatable :: unit, master_unit, code, iback ! (Number of bodies)
  character(len=80), dimension(:), allocatable :: c ! (Number of bodies)
  character(len=8), dimension(:), allocatable :: master_id, id ! (Number of bodies)
  real(double_precision), dimension(:), allocatable :: m, a, ns, is ! (Number of bodies)
  real(double_precision), dimension(:,:), allocatable :: el ! (22,Number of bodies)
  real(double_precision), dimension(:,:), allocatable :: x, v, xh, vh ! (3,Number of bodies)

  !------------------------------------------------------------------------------
  call getNumberOfBodies(nb_big_bodies=nbig, nb_bodies=nb_bodies_initial)
  
    ! We allocate and initialize arrays now that we know their sizes
  allocate(unit(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "unit" array'
  end if
  unit(1:nb_bodies_initial) = 0
  
  allocate(master_unit(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "master_unit" array'
  end if
  master_unit(1:nb_bodies_initial) = 0
  
  allocate(code(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "code" array'
  end if
  code(1:nb_bodies_initial) = 0
  
  allocate(iback(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "iback" array'
  end if
  iback(1:nb_bodies_initial) = 0
  
  allocate(c(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "c" array'
  end if
  c(1:nb_bodies_initial) = ''
  
  allocate(master_id(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "master_id" array'
  end if
  master_id(1:nb_bodies_initial) = ''
  
  allocate(id(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "id" array'
  end if
  id(1:nb_bodies_initial) = ''
  
  allocate(m(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "m" array'
  end if
  m(1:nb_bodies_initial) = 0.d0
  
  allocate(a(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "a" array'
  end if
  a(1:nb_bodies_initial) = 0.d0
  
  allocate(ns(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "ns" array'
  end if
  ns(1:nb_bodies_initial) = 0.d0
  
  allocate(is(nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "is" array'
  end if
  is(1:nb_bodies_initial) = 0.d0
  
  allocate(el(22,nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "el" array'
  end if
  el(1:22,1:nb_bodies_initial) = 0.d0
  
  allocate(x(3,nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "x" array'
  end if
  x(1:3,1:nb_bodies_initial) = 0.d0
  
  allocate(v(3,nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "v" array'
  end if
  v(1:3,1:nb_bodies_initial) = 0.d0
  
  allocate(vh(3,nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "vh" array'
  end if
  vh(1:3,1:nb_bodies_initial) = 0.d0
  
  allocate(xh(3,nb_bodies_initial), stat=error)
  if (error.ne.0) then
    write(*,*) 'Error: failed to allocate "xh" array'
  end if
  xh(1:3,1:nb_bodies_initial) = 0.d0
  
  
  !------------------------------------------------------------------------------
  allflag = 0
  tprevious = 0.d0
  rhocgs = AU_cgs * AU_cgs * AU_cgs * K2 / MSUN_cgs

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

  ! Open file containing parameters for this programme
  inquire (file='element.in', exist=test)
  if (test) then
     open (10, file='element.in', status='old')
  else
     call mio_err (6,mem(81),lmem(81),mem(88),lmem(88),' ',1,'element.in',9)
  end if

  ! Read number of input files
30 read (10,'(a250)') string
  if (string(1:1).eq.')') goto 30
  call mio_spl (250,string,nsub,lim)
  read (string(lim(1,nsub):lim(2,nsub)),*) ninfile

  ! Make sure all the input files exist
  do j = 1, ninfile
40   read (10,'(a250)') string
     if (string(1:1).eq.')') goto 40
     call mio_spl (250,string,nsub,lim)
     infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
     inquire (file=infile(j), exist=test)
     if (.not.test) call mio_err (6,mem(81),lmem(81),mem(88),lmem(88),' ',1,infile(j),80)
  end do

  ! What type elements does the user want?
  centre = 0
45 read (10,'(a250)') string
  if (string(1:1).eq.')') goto 45
  call mio_spl (250,string,nsub,lim)
  c2 = string(lim(1,nsub):(lim(1,nsub)+1))
  if (c2.eq.'ce'.or.c2.eq.'CE'.or.c2.eq.'Ce') then
     centre = 0
  else if (c2.eq.'ba'.or.c2.eq.'BA'.or.c2.eq.'Ba') then
     centre = 1
  else if (c2.eq.'ja'.or.c2.eq.'JA'.or.c2.eq.'Ja') then
     centre = 2
  else
     call mio_err (6,mem(81),lmem(81),mem(107),lmem(107),' ',1,'       Check element.in',23)
  end if

  ! Read parameters used by this programme
  timestyle = 1
  do j = 1, 4
50   read (10,'(a250)') string
     if (string(1:1).eq.')') goto 50
     call mio_spl (250,string,nsub,lim)
     c1 = string(lim(1,nsub):lim(2,nsub))
     if (j.eq.1) read (string(lim(1,nsub):lim(2,nsub)),*) teval
     teval = abs(teval) * .999d0
     if (j.eq.2.and.(c1.eq.'d'.or.c1.eq.'D')) timestyle = 0
     if (j.eq.3.and.(c1.eq.'y'.or.c1.eq.'Y')) timestyle = timestyle+2
     if (j.eq.4) call get_aei_format (string,timestyle,nel,iel,fout,header,lenhead)
  end do

  ! Read in the names of the objects for which orbital elements are required
  nopen = 0
  nwait = 0
  nmaster = 0
60 continue
  read (10,'(a250)',end=70) string
  call mio_spl (250,string,nsub,lim)
  if (string(1:1).eq.')'.or.lim(1,1).eq.-1) goto 60

  ! Either open an aei file for this object or put it on the waiting list
  nmaster = nmaster + 1
  itmp = min(7,lim(2,1)-lim(1,1))
  master_id(nmaster)='        '
  master_id(nmaster)(1:itmp+1) = string(lim(1,1):lim(1,1)+itmp)
  if (nopen.lt.NFILES) then
     nopen = nopen + 1
     master_unit(nmaster) = 10 + nopen
     call mio_aei (master_id(nmaster),master_unit(nmaster),  header,lenhead)
  else
     nwait = nwait + 1
     master_unit(nmaster) = -2
  end if
  goto 60

70 continue
  ! If no objects are listed in ELEMENT.IN assume that all objects are required
  if (nopen.eq.0) allflag = 1
  close (10)

  !------------------------------------------------------------------------------

  !  LOOP  OVER  EACH  INPUT  FILE  CONTAINING  INTEGRATION  DATA

90 continue
  firstflag = 0
  do i = 1, ninfile
     line_num = 0
     open (10, file=infile(i), status='old')

     ! Loop over each time slice
100  continue
     line_num = line_num + 1
     read (10,'(3a1)',end=900,err=666) check,style,type
     line_num = line_num - 1
     backspace 10

     ! Check if this is an old style input file
     if (ichar(check).eq.12.and.(style.eq.'0'.or.style.eq.'1'.or.style.eq.'2'.or.style.eq.'3'.or.style.eq.'4')) then
        write (*,'(/,2a)') ' ERROR: This is an old style data file','        Try running m_elem5.for instead.'
        stop
     end if
     if (ichar(check).ne.12) goto 666

     !------------------------------------------------------------------------------

     !  IF  SPECIAL  INPUT,  READ  TIME,  PARAMETERS,  NAMES,  MASSES  ETC.

     if (type.eq.'a') then
        line_num = line_num + 1
        read (10,'(3x,i2,a62,i1)') algor,cc(1:62),precision

        ! Decompress the time, number of objects, central mass and J components etc.
        time = mio_c2fl (cc(1:8))
        nbig = int(.5d0 + mio_c2re(cc(9:16), 0.d0, 11239424.d0, 3))
        nsml = int(.5d0 + mio_c2re(cc(12:19),0.d0, 11239424.d0, 3))
        mcen = mio_c2fl (cc(15:22))
        jcen(1) = mio_c2fl (cc(23:30))
        jcen(2) = mio_c2fl (cc(31:38))
        jcen(3) = mio_c2fl (cc(39:46))
        rcen = mio_c2fl (cc(47:54))
        rmax = mio_c2fl (cc(55:62))
        rfac = log10 (rmax / rcen)

        ! Read in strings containing compressed data for each object
        do j = 1, nbig + nsml
           line_num = line_num + 1
           read (10,'(a)',err=666) c(j)(1:51)
        end do

        ! Create input format list
        if (precision.eq.1) nchar = 2
        if (precision.eq.2) nchar = 4
        if (precision.eq.3) nchar = 7
        lenin = 3  +  6 * nchar
        fin(1:5) = '(a00)'
        write (fin(3:4),'(i2)') lenin

        ! For each object decompress its name, code number, mass, spin and density
        do j = 1, nbig + nsml
           k = int(.5d0 + mio_c2re(c(j)(1:8),0.d0,11239424.d0,3))
           id(k) = c(j)(4:11)
           el(18,k) = mio_c2fl (c(j)(12:19))
           s(1) = mio_c2fl (c(j)(20:27))
           s(2) = mio_c2fl (c(j)(28:35))
           s(3) = mio_c2fl (c(j)(36:43))
           el(21,k) = mio_c2fl (c(j)(44:51))

           ! Calculate spin rate and longitude & inclination of spin vector
           temp = sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
           if (temp.gt.0) then
              call mce_spin (1.d0,el(18,k)*K2,temp*K2,el(21,k)*rhocgs,el(20,k))
              temp = s(3) / temp
              if (abs(temp).lt.1.d0) then
                 is(k) = acos (temp)
                 ns(k) = atan2 (s(1), -s(2))
              else
                 if (temp.gt.0) is(k) = 0.d0
                 if (temp.lt.0) is(k) = PI
                 ns(k) = 0.d0
              end if
           else
              el(20,k) = 0.d0
              is(k) = 0.d0
              ns(k) = 0.d0
           end if

           ! Find the object on the master list
           unit(k) = 0
           do l = 1, nmaster
              if (id(k).eq.master_id(l)) unit(k) = master_unit(l)
           end do

           ! If object is not on the master list, add it to the list now
           if (unit(k).eq.0) then
              nmaster = nmaster + 1
              master_id(nmaster) = id(k)

              ! Either open an aei file for this object or put it on the waiting list
              if (allflag.eq.1) then
                 if (nopen.lt.NFILES) then
                    nopen = nopen + 1
                    master_unit(nmaster) = 10 + nopen
                    call mio_aei (master_id(nmaster),master_unit(nmaster),header,lenhead)
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

        !------------------------------------------------------------------------------

        !  IF  NORMAL  INPUT,  READ  COMPRESSED  ORBITAL  VARIABLES  FOR  ALL  OBJECTS

     else if (type.eq.'b') then
        line_num = line_num + 1
        read (10,'(3x,a14)',err=666) cc(1:14)

        ! Decompress the time and the number of objects
        time = mio_c2fl (cc(1:8))
        nbig = int(.5d0 + mio_c2re(cc(9:16),  0.d0, 11239424.d0, 3))
        nsml = int(.5d0 + mio_c2re(cc(12:19), 0.d0, 11239424.d0, 3))
        nbod = nbig + nsml
        if (firstflag.eq.0) t0 = time

        ! Read in strings containing compressed data for each object
        do j = 1, nbod
           line_num = line_num + 1
           read (10,fin,err=666) c(j)(1:lenin)
        end do

        ! Look for objects for which orbital elements are required
        m(1) = mcen * K2
        do j = 1, nbod
           code(j) = int(.5d0 + mio_c2re(c(j)(1:8), 0.d0,11239424.d0, 3))
           if (code(j).gt.nb_bodies_initial) then
              write (*,'(/,2a)') mem(81)(1:lmem(81)),  mem(90)(1:lmem(90))
              stop
           end if

           ! Decompress orbital variables for each object
           l = j + 1
           m(l) = el(18,code(j)) * K2
           fr     = mio_c2re (c(j)(4:11), 0.d0, rfac,  nchar)
           theta  = mio_c2re (c(j)(4+  nchar:11+  nchar), 0.d0, PI, nchar)
           phi    = mio_c2re (c(j)(4+2*nchar:11+2*nchar), 0.d0, TWOPI,  nchar)
           fv     = mio_c2re (c(j)(4+3*nchar:11+3*nchar), 0.d0, 1.d0, nchar)
           vtheta = mio_c2re (c(j)(4+4*nchar:11+4*nchar), 0.d0, PI, nchar)
           vphi   = mio_c2re (c(j)(4+5*nchar:11+5*nchar), 0.d0, TWOPI, nchar)
           call mco_ov2x (rcen,m(1),m(l),fr,theta,phi,fv,vtheta,vphi,x(1,l),x(2,l),x(3,l),v(1,l),v(2,l),v(3,l))
           el(16,code(j)) = sqrt(x(1,l)*x(1,l) + x(2,l)*x(2,l)    + x(3,l)*x(3,l))
        end do

        ! Convert to barycentric, Jacobi or close-binary coordinates if desired
        nbod1 = nbod + 1
        nbig1 = nbig + 1
!~         xh(:,:) = x(:,:)
!~         vh(:,:) = v(:,:)
        call mco_iden (time,jcen,nbod1,nbig1,temp,m,x,v,xh,vh,ngf,ngflag)
        if (centre.eq.1) call mco_h2b (jcen,nbod1,nbig1,temp,m,xh,vh,x,v)
        if (centre.eq.2) call mco_h2j (jcen,nbod1,nbig1,temp,m,xh,vh,x,v)
        if ((centre.eq.0).and.(algor.eq.11)) call mco_h2cb (jcen,nbod1,nbig1,temp,m,xh,vh,x,v)

        ! Put Cartesian coordinates into element arrays
        do j = 1, nbod
           k = code(j)
           l = j + 1
           el(10,k) = x(1,l)
           el(11,k) = x(2,l)
           el(12,k) = x(3,l)
           el(13,k) = v(1,l)
           el(14,k) = v(2,l)
           el(15,k) = v(3,l)

           ! Convert to Keplerian orbital elements
           gm = (mcen + el(18,k)) * K2
           call mco_x2el (gm,el(10,k),el(11,k),el(12,k),el(13,k),el(14,k),el(15,k),el(8,k),el(2,k),el(3,k),el(7,k),el(5,k),el(6,k))
           el(1,k) = el(8,k) / (1.d0 - el(2,k))
           el(9,k) = el(1,k) * (1.d0 + el(2,k))
           el(4,k) = mod(el(7,k) - el(5,k) + TWOPI, TWOPI)
           ! Calculate true anomaly
           if (el(2,k).eq.0) then
              el(17,k) = el(6,k)
           else
              temp = (el(8,k)*(1.d0 + el(2,k))/el(16,k) - 1.d0) /el(2,k)
              temp = sign (min(abs(temp), 1.d0), temp)
              el(17,k) = acos(temp)
              if (sin(el(6,k)).lt.0.d0) el(17,k) = TWOPI - el(17,k)
           end if
           ! Calculate obliquity
           el(19,k) = acos (cos(el(3,k))*cos(is(k))+ sin(el(3,k))*sin(is(k))*cos(ns(k) - el(5,k)))

           ! Convert angular elements from radians to degrees
           do l = 3, 7
              el(l,k) = mod(el(l,k) * RAD2DEG, 360.d0)
           end do
           el(17,k) = el(17,k) * RAD2DEG
           el(19,k) = el(19,k) * RAD2DEG
        end do

        ! Convert time to desired format
        if (timestyle.eq.0) t1 = time
        if (timestyle.eq.1) call mio_jd2y (time,year,month,t1)
        if (timestyle.eq.2) t1 = time - t0
        if (timestyle.eq.3) t1 = (time - t0) / 365.25d0

        ! If output is required at this epoch, write elements to appropriate files
        if ((firstflag.eq.0).or.(abs(time-tprevious).ge.teval)) then
           firstflag = 1
           tprevious = time

           ! Write required elements to the appropriate aei file
           do j = 1, nbod
              k = code(j)
              if (unit(k).ge.10) then
                 if (timestyle.eq.1) then
                    write (unit(k),fout) year,month,t1,(el(iel(l),k),l=1,nel)
                 else
                    write (unit(k),fout) t1,(el(iel(l),k),l=1,nel)
                 end if
              end if
           end do
        end if

        !------------------------------------------------------------------------------

        !  IF  TYPE  IS  NOT  'a'  OR  'b',  THE  INPUT  FILE  IS  CORRUPTED

     else
        goto 666
     end if

     ! Move on to the next time slice
     goto 100

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

     ! Move on to the next file containing integration data
900  continue
     close (10)
  end do

  ! Close aei files
  do j = 1, nopen
     close (10+j)
  end do
  nopen = 0

  ! If some objects remain on waiting list, read through input files again
  if (nwait.gt.0) then
     do j = 1, nmaster
        if (master_unit(j).ge.10) master_unit(j) = -1
        if ((master_unit(j).eq.-2).and.(nopen.lt.NFILES)) then
           nopen = nopen + 1
           nwait = nwait - 1
           master_unit(j) = 10 + nopen
           call mio_aei (master_id(j),master_unit(j),header,lenhead)
        end if
     end do
     goto 90
  end if

  !------------------------------------------------------------------------------

  !  CREATE  A  SUMMARY  OF  FINAL  MASSES  AND  ELEMENTS

  open (10, file='element.out', status='unknown')
  rewind 10

  if ((timestyle.eq.0).or.(timestyle.eq.2)) then
     write (10,'(/,a,f18.5,/)') ' Time (days): ',t1
  else if (timestyle.eq.1) then
     write (10,'(/,a,i10,1x,i2,1x,f8.5,/)') ' Date: ',year,month,t1
  else if (timestyle.eq.3) then
     write (10,'(/,a,f18.7,/)') ' Time (years): ',t1
  end if
  write (10,'(2a,/)') '              a        e       i      mass','    Rot/day  Obl'

  ! Sort surviving objects in order of increasing semi-major axis
  do j = 1, nbod
     k = code(j)
     a(j) = el(1,k)
  end do
  call mxx_sort (nbod,a,iback)

  ! Write values of a, e, i and m for surviving objects in an output file
  do j = 1, nbod
     k = code(iback(j))
     write (10,213) id(k),el(1,k),el(2,k),el(3,k),el(18,k) / EARTH_MASS,el(20,k),  el(19,k)
  end do

  !------------------------------------------------------------------------------

  ! Format statements
213 format (1x,a8,1x,f8.4,1x,f7.5,1x,f7.3,1p,e11.4,0p,1x,f6.3,1x,f6.2)

end program element

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author John E. Chambers
!> 
!
!> @date 31 January 2001
!
! DESCRIPTION: 
!> @brief Makes an output format list and file header for the orbital-element files
!! created by M_ELEM3.FOR\n
!! Also identifies which orbital elements will be output for each object.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_aei_format (string,timestyle,nel,iel,fout,header,lenhead)

  use physical_constant
  use mercury_constant
  use utilities, only : mio_spl

  implicit none


  ! Input
  integer, intent(in) :: timestyle
  character(len=250), intent(in) :: string
  
  ! Output
  integer, intent(out) :: nel
  integer, intent(out) :: iel(22)
  integer, intent(out) :: lenhead
  character(len=250), intent(out) :: header
  character(len=250), intent(out) :: fout

  ! Local
  integer :: i,j,pos,nsub,lim(2,20),formflag,lenfout,f1,f2,itmp
  character(len=1), dimension(22), parameter :: elcode = (/ 'a','e','i','g','n','l','p','q','b','x','y','z','u','v',&
       'w','r','f','m','o','s','d','c'/)
  character(len=4), dimension(22), parameter :: elhead = (/ '  a ','  e ','  i ','peri','node','  M ','long','  q ',&
       '  Q ','  x ','  y ','  z ',' vx ',' vy ',' vz ','  r ','  f ','mass','oblq','spin','dens','comp'/)

  !------------------------------------------------------------------------------

  ! Initialize header to a blank string
  do i = 1, 250
     header(i:i) = ' '
  end do

  ! Create part of the format list and header for the required time style
  select case (timestyle)
  case (0,2)
     fout(1:9) = '(1x,f18.5'
     lenfout = 9
     header(1:19) = '    Time (days)    '
     lenhead = 19
  case (1)
     fout(1:21) = '(1x,i10,1x,i2,1x,f8.5'
     lenfout = 21
     header(1:23) = '    Year/Month/Day     '
     lenhead = 23
  case (3)
     fout(1:9) = '(1x,f18.7'
     lenfout = 9
     header(1:19) = '    Time (years)   '
     lenhead = 19
  case default
    lenfout = 0
    lenhead = 0
  end select

  ! Identify the required elements
  call mio_spl (250,string,nsub,lim)
  do i = 1, nsub
     do j = 1, 22
        if (string(lim(1,i):lim(1,i)).eq.elcode(j)) iel(i) = j
     end do
  end do
  nel = nsub

  ! For each element, see whether normal or exponential notation is required
  do i = 1, nsub
     formflag = 0
     do j = lim(1,i)+1, lim(2,i)
        if (formflag.eq.0) pos = j
        if (string(j:j).eq.'.') formflag = 1
        if (string(j:j).eq.'e') formflag = 2
     end do

     ! Create the rest of the format list and header
     if (formflag.eq.1) then
        read (string(lim(1,i)+1:pos-1),*) f1
        read (string(pos+1:lim(2,i)),*) f2
        write (fout(lenfout+1:lenfout+10),'(a10)') ',1x,f  .  '
        write (fout(lenfout+6:lenfout+7),'(i2)') f1
        write (fout(lenfout+9:lenfout+10),'(i2)') f2
        lenfout = lenfout + 10
     else if (formflag.eq.2) then
        read (string(lim(1,i)+1:pos-1),*) f1
        write (fout(lenfout+1:lenfout+16),'(a16)') ',1x,1p,e  .  ,0p'
        write (fout(lenfout+9:lenfout+10),'(i2)') f1
        write (fout(lenfout+12:lenfout+13),'(i2)') f1 - 7
        lenfout = lenfout + 16
     end if
     itmp = (f1 - 4) / 2
     header(lenhead+itmp+2:lenhead+itmp+5) = elhead(iel(i))
     lenhead = lenhead + f1 + 1
  end do

  lenfout = lenfout + 1
  fout(lenfout:lenfout) = ')'

  !------------------------------------------------------------------------------

  return
end subroutine get_aei_format

