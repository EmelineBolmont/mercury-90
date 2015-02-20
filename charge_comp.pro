@functions
@constants

;! Mass star (kg):
Ms = 0.08*Msun
;! Number planets
nbp = 2
;! Number of planets tidally evolving
n_tid = 2

k2pDeltap = [2.465278d-3,2.465278d-3] ;day

;! If comparison with IDL simulations 1, if not 0
idl = 0
if idl eq 1 then begin
   nbp_idl =2
   nline_idl = dblarr(nbp_idl)
   filename_idl = strarr(nbp_idl)
   filename_idl(0) = 'datatides_18.0000_0.00000_12_2_3_3_0.0261799_0.00000_0.dat'
   if nbp_idl eq 2 then filename_idl(1) = 'datatides_100.000_50.0000_12_5_3_4_0.0174533_0.00000_0.dat'
endif

;************************************************************

nlineheader = 4                 ;! number of header lines in the data files

nlineb = fltarr(nbp)
for i=0,nbp-1 do begin
   filename = 'PLANET'+strtrim(i+1,2)+'.aei'
   nlineb(i) = file_lines(filename)-nlineheader
endfor

if idl eq 1 then begin
   for i=0,nbp_idl-1 do begin
      filename = filename_idl(i)
      print,filename
      nline_idl(i) = file_lines(filename)-1
   endfor
endif

if n_tid ge 1 then begin
   filenames = 'spins.out'
   print,filenames
   ; In Mercury spin is in day-1, but later on it is converted to s-1
   readcol,filenames,toto1,spinstx,spinsty,spinstz,Rst,rg2s,k2s,sigmas,format='F,F,F,F,F,F,F,F'

   spinpx = dblarr(n_tid,n_elements(toto1))
   spinpy = dblarr(n_tid,n_elements(toto1))
   spinpz = dblarr(n_tid,n_elements(toto1))
   horbx  = dblarr(n_tid,n_elements(toto1))
   horby  = dblarr(n_tid,n_elements(toto1))
   horbz  = dblarr(n_tid,n_elements(toto1))
   Rp     = dblarr(n_tid,n_elements(toto1)) 
   rg2p   = dblarr(n_tid,n_elements(toto1))
   dEdt   = dblarr(n_tid,n_elements(toto1))

   for i=0,n_tid-1 do begin 
      filenamep = 'spinp'+strtrim(i+1,2)+'.out'
      print,filenamep
      ; in Mercury spin is in day-1, later converted to s-1
      ; Rp is here in Rsun, rg2p does not have unit
      readcol,filenamep,toto1,spinp1x,spinp1y,spinp1z,Rp1,rg2p1,format='F,F,F,F,F,F'
      spinpx(i,*) = spinp1x & spinpy(i,*) = spinp1y & spinpz(i,*) = spinp1z
      Rp(i,*) = Rp1 & rg2p(i,*) = rg2p1

      filenameh = 'horb'+strtrim(i+1,2)+'.out'
      print,filenameh
      ; The unit here does not matter, we always normalize later
      readcol,filenameh,toto1,horb1x,horb1y,horb1z,format='F,F,F,F'
      horbx(i,*) = horb1x & horby(i,*) = horb1y & horbz(i,*) = horb1z  

      filenamee = 'dEdt'+strtrim(i+1,2)+'.out'
      print,filenamee
      ; Mercury gives dE/dt in Msun.AU^2.day^-3, we convert here in W
      readcol,filenamee,toto1,tmp,format='F,F'
      dEdt(i,*) = tmp*6.90125d37 ;conversation from Msun.AU^2.day^-3 to W
   endfor
endif 

; n line maximum
nmaxb = max(nlineb)
if idl eq 1 then nline = max(nline_idl)

;; mass = [0.01,0.012,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.072,0.075,0.08]

;; nlines = dblarr(n_elements(mass))
;; for j = 0,n_elements(mass)-1 do begin
;;    filename = 'mass_'+strtrim(1000*mass(j),2)+'.dat' 
;;    nlines(j) = file_lines(filename)
;; endfor
;; nmax = max(nlines)

;; toto   = dblarr(n_elements(mass),nmax)
;; radius = dblarr(n_elements(mass),nmax)
;; lumi   = dblarr(n_elements(mass),nmax)
;; HZinGJ = dblarr(n_elements(mass),nmax)
;; HZoutGJ= dblarr(n_elements(mass),nmax)
;; HZinb  = dblarr(n_elements(mass),nmax)
;; HZoutb = dblarr(n_elements(mass),nmax)

;; for j = 0,n_elements(mass)-1 do begin
;;    filename = 'mass_'+strtrim(1000*mass(j),2)+'.dat' 
;;    nt = file_lines(filename)
;;    read_array = dblarr(7,nt)
;;    openr,1,filename
;;    readf,1,read_array
;;    close,1
   
;;    toto[j,0:nt-1] = read_array(0,*)
;;    radius[j,0:nt-1] = read_array(1,*)
;;    lumi[j,0:nt-1] = read_array(2,*)
;;    HZinGJ[j,0:nt-1] = read_array(3,*)
;;    HZoutGJ[j,0:nt-1] = read_array(4,*)
;;    HZinb[j,0:nt-1] = read_array(5,*)
;;    HZoutb[j,0:nt-1] = read_array(6,*)
;; endfor

header = strarr(nlineheader)
headeri = strarr(1)

; Table for the big planets (PLANET) : 

tb     =  dblarr(nbp,nmaxb)
ab     =  dblarr(nbp,nmaxb)
eb     =  dblarr(nbp,nmaxb)
incb   =  dblarr(nbp,nmaxb)
perib  =  dblarr(nbp,nmaxb)
nodeb  =  dblarr(nbp,nmaxb)
manomb =  dblarr(nbp,nmaxb)
lambda =  dblarr(nbp,nmaxb)
mb     =  dblarr(nbp,nmaxb)
xb     =  dblarr(nbp,nmaxb)
yb     =  dblarr(nbp,nmaxb)
zb     =  dblarr(nbp,nmaxb)
ub     =  dblarr(nbp,nmaxb)
vb     =  dblarr(nbp,nmaxb)
wb     =  dblarr(nbp,nmaxb)

;; print,nmaxb,nlineb1,nlineb2
;**************************************
; Table filling for big planets (PLANET)
for i=0,nbp-1 do begin
   filename = 'PLANET'+strtrim(i+1,2)+'.aei'
   nlineb1 = file_lines(filename)-nlineheader
   read_array = dblarr(14,nlineb1)
   openr,1,filename
   readf,1,header
   readf,1,read_array
   close,1
   ; time in years
   tb[i,0:nlineb1-1] = read_array(0,*)
   ; semi-major axis in AU
   ab[i,0:nlineb1-1] = read_array(1,*)
   eb[i,0:nlineb1-1] = read_array(2,*)
   ; inclination in degrees
   incb[i,0:nlineb1-1] = read_array(3,*)
   ; argument of pericenter, longitude of ascending node and mean anomaly in degrees
   perib[i,0:nlineb1-1] = read_array(4,*)
   nodeb[i,0:nlineb1-1] = read_array(5,*)
   manomb[i,0:nlineb1-1] = read_array(6,*)
   ; mass in Msun
   mb[i,0:nlineb1-1] = read_array(7,*)
   ; Positions in AU
   xb[i,0:nlineb1-1] = read_array(8,*)
   yb[i,0:nlineb1-1] = read_array(9,*)
   zb[i,0:nlineb1-1] = read_array(10,*)
   ; Velocities in AU/day? 
   ub[i,0:nlineb1-1] = read_array(11,*)
   vb[i,0:nlineb1-1] = read_array(12,*)
   wb[i,0:nlineb1-1] = read_array(13,*)
endfor
;zozo

; Angular momentum calculation
; calculated with sma, ecc and inclination
horb_sec = dblarr(nbp,nmaxb)
; circular angular momentum (what it should be with ecc = 0, inc = 0)
horb_circ = dblarr(nbp,nmaxb)
; Sum over all planets
horb_s = dblarr(nmaxb)
horb_c = dblarr(nmaxb)

; Angular momentum deficit calculation 
AMD_sec = dblarr(nmaxb)

for j=0,nbp-1 do begin
   for i=0,nmaxb-1 do begin
      horb_sec(j,i)  = Ms*mb(j,i)*Msun/(Ms+mb(j,i)*Msun)*sqrt(G*(Ms+mb(j,i)*Msun)*ab(j,i)*AU)*sqrt(1-eb(j,i)^2)*cos(incb(j,i)*!pi/180.d0) ; kg.m^2.s-1
      horb_circ(j,i) = Ms*mb(j,i)*Msun/(Ms+mb(j,i)*Msun)*sqrt(G*(Ms+mb(j,i)*Msun)*ab(j,i)*AU) ; kg.m^2.s-1
      ; Sum on number of planets to have total orbital momentum
      horb_s(i) = horb_s(i) + horb_sec(j,i)
      horb_c(i) = horb_c(i) + horb_circ(j,i)
      ; Calculation of angular momentum deficit
      AMD_sec(i) = horb_c(i) - horb_s(i)
   endfor
endfor

if n_tid ge 1 then begin   


   ;IDL Data
   if idl eq 1 then begin
       ;! Table for the idl planet (DATATIDES) : 
       ti         =  dblarr(nbp_idl,nline)
       ai         =  dblarr(nbp_idl,nline)
       ei         =  dblarr(nbp_idl,nline)
       oblpi      =  dblarr(nbp_idl,nline)
       oblsi      =  dblarr(nbp_idl,nline)
       rotpi      =  dblarr(nbp_idl,nline)
       rotsi      =  dblarr(nbp_idl,nline)
       Rpi        =  dblarr(nbp_idl,nline)
       Rsi        =  dblarr(nbp_idl,nline)
       rg2si      =  dblarr(nbp_idl,nline)
       tidefluxi  =  dblarr(nbp_idl,nline)
       Ip         =  dblarr(nbp_idl)  

       for i = 0,nbp_idl-1 do begin
           headeri = strarr(1)
           nline = file_lines(filename_idl(i))-1
           read_array = dblarr(10,nline)
           openr,1,filename_idl(i)
           readf,1,headeri
           readf,1,read_array
           close,1

           ti[i,0:nline-1] = read_array(0,*);+toto1(0)
           ai[i,0:nline-1] = read_array(1,*)
           ei[i,0:nline-1] = read_array(2,*)
           oblpi[i,0:nline-1] = read_array(3,*)
           oblsi[i,0:nline-1] = read_array(4,*)
           rotpi[i,0:nline-1] = read_array(5,*)
           rotsi[i,0:nline-1] = read_array(6,*)
           Rpi[i,0:nline-1] = read_array(7,*)
           Rsi[i,0:nline-1] = read_array(8,*)
           rg2si[i,0:nline-1] = read_array(9,*)
        endfor

        for i = 0,nbp_idl-1 do begin
            for j = 0,nline-1 do begin
                tidefluxi(i,j) = enerdot(ai(i,j)*AU,ei(i,j),rotpi(i,j),oblpi(i,j)*!Pi/180.d0,G $
  		          	,mb(i,0)*Msun,Ms,Rpi(i,j)*Rearth,k2pDeltap(i)*day)/(4.d0*!Pi*(Rpi(i,j)*Rearth)^2)
            endfor
            Ip(i) = rg2p(i,0)*mb(i,0)*Msun*(Rp(i,0)*rsun)^2
        endfor 
    endif

    ;! Obliquities calculations
    tmp              = dblarr(nmaxb)
    ; Obliquities in degrees, precession angle in degrees
    oblpm            = dblarr(n_tid,nmaxb)
    oblsm            = dblarr(n_tid,nmaxb)
    precession_angle = dblarr(n_tid,nmaxb)
    ; That's where the spin is converted to s-1
    spinp            = dblarr(n_tid,nmaxb)

    for i = 0,n_tid-1 do begin
       for bou = 0,nmaxb-1 do begin

          tmp(bou)=(horbx(i,bou)*spinpx(i,bou) $
                    +horby(i,bou)*spinpy(i,bou) $
                    +horbz(i,bou)*spinpz(i,bou)) $
                   /(sqrt(horbx(i,bou)^2+horby(i,bou)^2+horbz(i,bou)^2) $
                     *sqrt(spinpx(i,bou)^2+spinpy(i,bou)^2+spinpz(i,bou)^2))
          if abs(tmp(bou)) le 1.d0 then $ 
             oblpm(i,bou) = acos(tmp(bou))*180.d0/!Pi
          if abs(tmp(bou)) gt 1.d0 then oblpm(i,bou) = 1.0d-6

          tmp(bou)=(horbx(i,bou)*spinstx(bou) $
                    +horby(i,bou)*spinsty(bou) $
                    +horbz(i,bou)*spinstz(bou)) $
                   /(sqrt(horbx(i,bou)^2+horby(i,bou)^2+horbz(i,bou)^2) $
                     *sqrt(spinstx(bou)^2+spinsty(bou)^2+spinstz(bou)^2))
          if abs(tmp(bou)) le 1.d0 then $ 
             oblsm(i,bou) = acos(tmp(bou))*180.d0/!Pi
          if abs(tmp(bou)) gt 1.d0 then oblsm(i,bou) = 1.0d-6

          tmp(bou)=1.d0/sin(oblpm(i,bou)) $
                    *(xb(i,bou)*spinpx(i,bou) $
                    +yb(i,bou)*spinpy(i,bou) $
                    +zb(i,bou)*spinpz(i,bou)) $
                    /(sqrt(xb(i,bou)^2+yb(i,bou)^2+zb(i,bou)^2) $
                     *sqrt(spinpx(i,bou)^2+spinpy(i,bou)^2+spinpz(i,bou)^2))
          if abs(tmp(bou)) le 1.d0 then $
             precession_angle(i,bou) = acos(tmp(bou))*180.d0/!Pi
          if abs(tmp(bou)) gt 1.d0 then precession_angle(i,bou) = 1.0d-6

          ; if the angle between horb and spinp is more than 90: retrograde rotation
          if ((horbx(i,bou)*spinpx(i,bou) $
                    +horby(i,bou)*spinpy(i,bou) $
                    +horbz(i,bou)*spinpz(i,bou)) ge 0) $
              spinp(i,bou) = sqrt(spinpx(i,bou)^2 $
                 +spinpy(i,bou)^2+spinpz(i,bou)^2)/day
          if ((horbx(i,bou)*spinpx(i,bou) $
                    +horby(i,bou)*spinpy(i,bou) $
                    +horbz(i,bou)*spinpz(i,bou)) lt 0) $
              spinp(i,bou) = -sqrt(spinpx(i,bou)^2 $
                 +spinpy(i,bou)^2+spinpz(i,bou)^2)/day
       endfor
    endfor

    ; Angular momentum calculation
    ; calculated with horb coming from mercury
    horb_vec = dblarr(nbp,nmaxb)
    ; calculated with sma, ecc and inclination
    horb_sec = dblarr(nbp,nmaxb)
    ; circular angular momentum (what it should be with ecc = 0, inc = 0)
    horb_circ = dblarr(nbp,nmaxb)
    ; Sum over all planets
    horb_v = dblarr(nmaxb)
    horb_s = dblarr(nmaxb)
    horb_c = dblarr(nmaxb)

    ; Angular momentum deficit calculation 
    AMD     = dblarr(nmaxb)
    AMD_sec = dblarr(nmaxb)

    spinst = dblarr(nmaxb)
    momspin = dblarr(n_tid,nmaxb)
    momspitot = dblarr(nmaxb)
    momstar = dblarr(nmaxb)

    ; Tidal flux
    tidalflux = dblarr(n_tid,nmaxb)
    inst_tidalflux = dblarr(n_tid,nmaxb)

    for j=0,nbp-1 do begin
       for i=0,nmaxb-1 do begin
           horb_vec(j,i)  = Ms*mb(j,i)*Msun/(Ms+mb(j,i)*Msun)*sqrt(horbx(j,i)^2+horby(j,i)^2+horbz(j,i)^2)*(AU^2/day) ; kg.m^2.s-1
           horb_sec(j,i)  = Ms*mb(j,i)*Msun/(Ms+mb(j,i)*Msun)*sqrt(G*(Ms+mb(j,i)*Msun)*ab(j,i)*AU)*sqrt(1-eb(j,i)^2)*cos(incb(j,i)*!pi/180.d0) ; kg.m^2.s-1
           horb_circ(j,i) = Ms*mb(j,i)*Msun/(Ms+mb(j,i)*Msun)*sqrt(G*(Ms+mb(j,i)*Msun)*ab(j,i)*AU) ; kg.m^2.s-1
           ; Sum on number of planets to have total orbital momentum
           horb_v(i) = horb_v(i) + horb_vec(j,i)
           horb_s(i) = horb_s(i) + horb_sec(j,i)
           horb_c(i) = horb_c(i) + horb_circ(j,i)
           ; Calculation of angular momentum deficit
           AMD(i)     = horb_c(i) - horb_v(i)
           AMD_sec(i) = horb_c(i) - horb_s(i)

           spinst(i) = sqrt(spinstx(i)^2+spinsty(i)^2+spinstz(i)^2)/day ; s-1
       endfor
    endfor
    for i=0,nmaxb-1 do begin
        momstar(i)=rg2s(i)*Ms*(Rst(i)*Rsun)^2*spinst(i) 
    endfor

    for j=0,n_tid-1 do begin
       for i=0,nmaxb-1 do begin
           momspin(j,i) = rg2p(j,i)*mb(j,i)*Msun*(Rp(j,i)*rsun)^2*spinp(j,i)

           ; Calculation of energydot and tidal flux, in W/m2
           tidalflux(j,i) = enerdot(ab(j,i)*AU,eb(j,i),spinp(j,i),oblpm(j,i)*!Pi/180.d0,G,mb(j,i)*Msun $
                 ,Ms,Rp(j,i)*rsun,k2pDeltap(j)*day)/(4*!Pi*(Rp(j,i)*rsun)^2)
           inst_tidalflux(j,i) = dEdt(j,i)/(4*!Pi*(Rp(j,i)*rsun)^2)
       endfor
    endfor

    for i=0,nmaxb-1 do begin
        for j=0,n_tid-1 do begin
            momspitot(i)=momspitot(i)+momspin(j,i)
        endfor
    endfor

    indicend = dblarr(2,nbp)
    for j = 0,nbp-1 do begin
        for i = 0,n_elements(tb(j,*))-1 do begin
            if tb(j,i) le 1.d10 then begin
                indicend(0,j) = i
            endif else begin
                break
            endelse
        endfor
    endfor

    for i = 0,n_elements(toto1(*))-1 do begin
        if toto1(i) le 1.d10 then begin
            indicend(1,*) = i
        endif else begin
            break
        endelse
    endfor 
endif

END
