;PRO charge


;! Number planets
nbp = 2
;! Number of planets tidally evolving
n_tid = 2
;! If comparison with IDL simulations 1, if not 0
idl = 1

nlineheader = 4                 ;! number of header lines in the data files

nlineb = fltarr(nbp)
for i=0,nbp-1 do begin
   filename = 'PLANET'+strtrim(i+1,2)+'.aei'
   nlineb(i) = file_lines(filename)-nlineheader
endfor

if idl eq 1 then begin
   filename4 = 'datatides_18.0000_0.00000_12_2_3_3_0.0261799_0.00000_0.dat'
   print,filename4
   nlineb4 = file_lines(filename4)-1
   filename5 = 'datatides_100.000_50.0000_12_5_3_4_0.0174533_0.00000_0.dat'
   print,filename5
   nlineb5 = file_lines(filename5)-1
endif

if n_tid ge 1 then begin
   filenames = 'spins.dat'
   print,filenames
   readcol,filenames,sss,toto1,spinstx,spinsty,spinstz,Rst,format='A,F,F,F,F,F'
   filenamep1 = 'spinp1.dat'
   print,filenamep1
   readcol,filenamep1,ppp,toto1,spinp1x,spinp1y,spinp1z,format='A,F,F,F,F'
   filenameh1 = 'horb1.dat'
   print,filenameh1
   readcol,filenameh1,hhh,toto1,horb1x,horb1y,horb1z,format='A,F,F,F,F'
endif
if n_tid ge 2 then begin
   filenamep2 = 'spinp2.dat'
   print,filenamep2
   readcol,filenamep2,ppp,toto1,spinp2x,spinp2y,spinp2z,format='A,F,F,F,F'
   filenameh2 = 'horb2.dat'
   print,filenameh2
   readcol,filenameh2,hhh,toto1,horb2x,horb2y,horb2z,format='A,F,F,F,F'
endif
if n_tid ge 3 then begin
   filenamep3 = 'spinp3.dat'
   print,filenamep3
   readcol,filenamep3,ppp,toto1,spinp3x,spinp3y,spinp3z,format='A,F,F,F,F'
   filenameh3 = 'horb3.dat'
   print,filenameh3
   readcol,filenameh3,hhh,toto1,horb3x,horb3y,horb3z,format='A,F,F,F,F'
endif


; n line maximum
nmaxb = max(nlineb)
if idl eq 1 then nline = max([nlineb4,nlineb5])

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

   tb[i,0:nlineb1-1] = read_array(0,*)
   ab[i,0:nlineb1-1] = read_array(1,*)
   eb[i,0:nlineb1-1] = read_array(2,*)
   incb[i,0:nlineb1-1] = read_array(3,*)
   perib[i,0:nlineb1-1] = read_array(4,*)
   nodeb[i,0:nlineb1-1] = read_array(5,*)
   manomb[i,0:nlineb1-1] = read_array(6,*)
   mb[i,0:nlineb1-1] = read_array(7,*)
   xb[i,0:nlineb1-1] = read_array(8,*)
   yb[i,0:nlineb1-1] = read_array(9,*)
   zb[i,0:nlineb1-1] = read_array(10,*)
   ub[i,0:nlineb1-1] = read_array(11,*)
   vb[i,0:nlineb1-1] = read_array(12,*)
   wb[i,0:nlineb1-1] = read_array(13,*)
endfor

if idl eq 1 then begin
   ;!**************************************
   ;! Table for the idl planet (DATATIDES) : 
   ti         =  dblarr(n_tid,nline)
   ai         =  dblarr(n_tid,nline)
   ei         =  dblarr(n_tid,nline)
   oblpi      =  dblarr(n_tid,nline)
   oblsi      =  dblarr(n_tid,nline)
   rotpi      =  dblarr(n_tid,nline)
   rotsi      =  dblarr(n_tid,nline)
   Rpi        =  dblarr(n_tid,nline)
   Rsi        =  dblarr(n_tid,nline)
   rg2si      =  dblarr(n_tid,nline)
   
   headeri = strarr(1)
   nlineb4 = file_lines(filename4)-1
   read_array = dblarr(10,nlineb4)
   openr,1,filename4
   readf,1,headeri
   readf,1,read_array
   close,1
   
   ti[0,0:nlineb4-1] = read_array(0,*);+toto1(0)
   ai[0,0:nlineb4-1] = read_array(1,*)
   ei[0,0:nlineb4-1] = read_array(2,*)
   oblpi[0,0:nlineb4-1] = read_array(3,*)
   oblsi[0,0:nlineb4-1] = read_array(4,*)
   rotpi[0,0:nlineb4-1] = read_array(5,*)
   rotsi[0,0:nlineb4-1] = read_array(6,*)
   Rpi[0,0:nlineb4-1] = read_array(7,*)
   Rsi[0,0:nlineb4-1] = read_array(8,*)
   rg2si[0,0:nlineb4-1] = read_array(9,*)
   
   nlineb5 = file_lines(filename5)-1
   read_array = dblarr(10,nlineb5)
   openr,1,filename5
   readf,1,headeri
   readf,1,read_array
   close,1
   
   ti[1,0:nlineb5-1] = read_array(0,*);+toto1(0)
   ai[1,0:nlineb5-1] = read_array(1,*)
   ei[1,0:nlineb5-1] = read_array(2,*)
   oblpi[1,0:nlineb5-1] = read_array(3,*)
   oblsi[1,0:nlineb5-1] = read_array(4,*)
   rotpi[1,0:nlineb5-1] = read_array(5,*)
   rotsi[1,0:nlineb5-1] = read_array(6,*)
   Rpi[1,0:nlineb5-1] = read_array(7,*)
   Rsi[1,0:nlineb5-1] = read_array(8,*)
   rg2si[1,0:nlineb5-1] = read_array(9,*)
endif

if n_tid ge 1 then begin
;! Obliquities calculations

tmp     = dblarr(n_elements(horb1x))
oblp1m  = dblarr(n_elements(horb1x))
obls1m  = dblarr(n_elements(horb1x))
if n_tid ge 2 then begin
   oblp2m  = dblarr(n_elements(horb1x))
   obls2m  = dblarr(n_elements(horb1x))
endif
if n_tid ge 3 then begin
   oblp3m  = dblarr(n_elements(horb1x))
   obls3m  = dblarr(n_elements(horb1x))
endif

for bou = 0,n_elements(horb1x)-1 do begin

   tmp(bou)=(horb1x(bou)*spinp1x(bou) $
             +horb1y(bou)*spinp1y(bou) $
             +horb1z(bou)*spinp1z(bou)) $
            /(sqrt(horb1x(bou)^2+horb1y(bou)^2+horb1z(bou)^2) $
              *sqrt(spinp1x(bou)^2+spinp1y(bou)^2+spinp1z(bou)^2))
   if abs(tmp(bou)) le 1.d0 then $ 
      oblp1m(bou) = acos(tmp(bou))*180.d0/!Pi
   if abs(tmp(bou)) gt 1.d0 then oblp1m(bou) = 1.0d-6
  
   tmp(bou)=(horb1x(bou)*spinstx(bou) $
             +horb1y(bou)*spinsty(bou) $
             +horb1z(bou)*spinstz(bou)) $
            /(sqrt(horb1x(bou)^2+horb1y(bou)^2+horb1z(bou)^2) $
              *sqrt(spinstx(bou)^2+spinsty(bou)^2+spinstz(bou)^2))
   if abs(tmp(bou)) le 1.d0 then $ 
      obls1m(bou) = acos(tmp(bou))*180.d0/!Pi
   if abs(tmp(bou)) gt 1.d0 then obls1m(bou) = 1.0d-6
   
   if n_tid ge 2 then begin
      tmp(bou)=(horb2x(bou)*spinp2x(bou) $
             +horb2y(bou)*spinp2y(bou) $
             +horb2z(bou)*spinp2z(bou)) $
            /(sqrt(horb2x(bou)^2+horb2y(bou)^2+horb2z(bou)^2) $
              *sqrt(spinp2x(bou)^2+spinp2y(bou)^2+spinp2z(bou)^2))
      if abs(tmp(bou)) le 1.d0 then $ 
         oblp2m(bou) = acos(tmp(bou))*180.d0/!Pi 
      if abs(tmp(bou)) gt 1.d0 then oblp2m(bou) = 1.0d-6
      
      tmp(bou)=(horb2x(bou)*spinstx(bou) $
                +horb2y(bou)*spinsty(bou) $
                +horb2z(bou)*spinstz(bou)) $
               /(sqrt(horb2x(bou)^2+horb2y(bou)^2+horb2z(bou)^2) $
                 *sqrt(spinstx(bou)^2+spinsty(bou)^2+spinstz(bou)^2))
      if abs(tmp(bou)) le 1.d0 then $ 
         obls2m(bou) = acos(tmp(bou))*180.d0/!Pi
      if abs(tmp(bou)) gt 1.d0 then obls2m(bou) = 1.0d-6
   endif
   
   if n_tid ge 3 then begin
      tmp(bou)=(horb3x(bou)*spinp3x(bou) $
             +horb3y(bou)*spinp3y(bou) $
             +horb3z(bou)*spinp3z(bou)) $
            /(sqrt(horb3x(bou)^2+horb3y(bou)^2+horb3z(bou)^2) $
              *sqrt(spinp3x(bou)^2+spinp3y(bou)^2+spinp3z(bou)^2))
      if abs(tmp(bou)) le 1.d0 then $ 
         oblp3m(bou) = acos(tmp(bou))*180.d0/!Pi 
      if abs(tmp(bou)) gt 1.d0 then oblp3m(bou) = 1.0d-6
      
      tmp(bou)=(horb3x(bou)*spinstx(bou) $
                +horb3y(bou)*spinsty(bou) $
                +horb3z(bou)*spinstz(bou)) $
               /(sqrt(horb3x(bou)^2+horb3y(bou)^2+horb3z(bou)^2) $
                 *sqrt(spinstx(bou)^2+spinsty(bou)^2+spinstz(bou)^2))
      if abs(tmp(bou)) le 1.d0 then $ 
         obls3m(bou) = acos(tmp(bou))*180.d0/!Pi
      if abs(tmp(bou)) gt 1.d0 then obls3m(bou) = 1.0d-6
   endif
endfor
endif
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
if n_tid ge 1 then begin
   for i = 0,n_elements(toto1(*))-1 do begin
      if toto1(i) le 1.d10 then begin
         indicend(1,0) = i
         indicend(1,1) = i
      endif else begin
         break
      endelse
   endfor 
endif

END
