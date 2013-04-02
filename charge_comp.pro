;PRO charge


; Indicate here the number of big planets you took into account: 
; Compare idlcode results with 
; mercury_tides and mercury_tides_GR

nbp = 2

nlineheader = 4                 ; number of header lines in the data files
nlineb1 = fltarr(1)
nlineb2 = fltarr(1)

filename1 = 'PLANET1.aei'
nlineb1 = file_lines(filename1)-nlineheader
filename2 = 'PLANET2.aei'
nlineb2 = file_lines(filename2)-nlineheader
;;! filename3 = 'PLANET1.aei'
;;! nlineb3 = file_lines(filename3)-nlineheader

filename4 = 'datatides_18.0000_0.00000_12_2_3_3_0.0261799_0.00000_0.dat'
print,filename4
nlineb4 = file_lines(filename4)-1
filename5 = 'datatides_100.000_50.0000_12_5_3_4_0.0174533_0.00000_0.dat'
print,filename5
nlineb5 = file_lines(filename5)-1

filenames = 'spins.dat'
print,filenames
readcol,filenames,sss,toto1,spinstx,spinsty,spinstz,Rst,format='A,F,F,F,F,F'
filenamep1 = 'spinp1.dat'
print,filenamep1
readcol,filenamep1,ppp,toto1,spinp1x,spinp1y,spinp1z,format='A,F,F,F,F'
filenamep2 = 'spinp2.dat'
print,filenamep2
readcol,filenamep2,ppp,toto1,spinp2x,spinp2y,spinp2z,format='A,F,F,F,F'
filenameh1 = 'horb1.dat'
print,filenameh1
readcol,filenameh1,hhh,toto1,horb1x,horb1y,horb1z,format='A,F,F,F,F'
filenameh2 = 'horb2.dat'
print,filenameh2
readcol,filenameh2,hhh,toto1,horb2x,horb2y,horb2z,format='A,F,F,F,F'

; n line maximum
nmaxb = max([nlineb1,nlineb2]);,nlineb3
nline = max([nlineb4,nlineb5])

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

nlineb1 = file_lines(filename1)-nlineheader
read_array = dblarr(14,nlineb1)
openr,1,filename1
readf,1,header
readf,1,read_array
close,1

tb[0,0:nlineb1-1] = read_array(0,*)
ab[0,0:nlineb1-1] = read_array(1,*)
eb[0,0:nlineb1-1] = read_array(2,*)
incb[0,0:nlineb1-1] = read_array(3,*)
perib[0,0:nlineb1-1] = read_array(4,*)
nodeb[0,0:nlineb1-1] = read_array(5,*)
manomb[0,0:nlineb1-1] = read_array(6,*)
mb[0,0:nlineb1-1] = read_array(7,*)
xb[0,0:nlineb1-1] = read_array(8,*)
yb[0,0:nlineb1-1] = read_array(9,*)
zb[0,0:nlineb1-1] = read_array(10,*)
ub[0,0:nlineb1-1] = read_array(11,*)
vb[0,0:nlineb1-1] = read_array(12,*)
wb[0,0:nlineb1-1] = read_array(13,*)

;;;;**************************************
nlineb2 = file_lines(filename2)-nlineheader
read_array = dblarr(14,nlineb2)
openr,1,filename2
readf,1,header
readf,1,read_array
close,1

tb[1,0:nlineb2-1] = read_array(0,*)
ab[1,0:nlineb2-1] = read_array(1,*)
eb[1,0:nlineb2-1] = read_array(2,*)
incb[1,0:nlineb2-1] = read_array(3,*)
perib[1,0:nlineb2-1] = read_array(4,*)
nodeb[1,0:nlineb2-1] = read_array(5,*)
manomb[1,0:nlineb2-1] = read_array(6,*)
mb[1,0:nlineb2-1] = read_array(7,*)
xb[1,0:nlineb2-1] = read_array(8,*)
yb[1,0:nlineb2-1] = read_array(9,*)
zb[1,0:nlineb2-1] = read_array(10,*)
ub[1,0:nlineb2-1] = read_array(11,*)
vb[1,0:nlineb2-1] = read_array(12,*)
wb[1,0:nlineb2-1] = read_array(13,*)



;; ;;;;**************************************
;; nlineb3 = file_lines(filename3)-nlineheader
;; read_array = dblarr(8,nlineb3)
;; openr,1,filename3
;; readf,1,header
;; readf,1,read_array
;; close,1

;; tb[2,0:nlineb3-1] = read_array(0,*)
;; ab[2,0:nlineb3-1] = read_array(1,*)
;; eb[2,0:nlineb3-1] = read_array(2,*)
;; incb[2,0:nlineb3-1] = read_array(3,*)
;; perib[2,0:nlineb3-1] = read_array(4,*)
;; nodeb[2,0:nlineb3-1] = read_array(5,*)
;; manomb[2,0:nlineb3-1] = read_array(6,*)
;; mb[2,0:nlineb3-1] = read_array(7,*)

;**************************************
; Table for the idl planet (DATATIDES) : 
ti         =  dblarr(2,nline)
ai         =  dblarr(2,nline)
ei         =  dblarr(2,nline)
oblpi      =  dblarr(2,nline)
oblsi      =  dblarr(2,nline)
rotpi      =  dblarr(2,nline)
rotsi      =  dblarr(2,nline)
Rpi        =  dblarr(2,nline)
Rsi        =  dblarr(2,nline)
rg2si      =  dblarr(2,nline)

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


;! Obliquities calculations

tmp     = dblarr(n_elements(horb1x))
oblp1m  = dblarr(n_elements(horb1x))
oblp2m  = dblarr(n_elements(horb1x))
obls1m  = dblarr(n_elements(horb1x))
obls2m  = dblarr(n_elements(horb1x))

for bou = 0,n_elements(horb1x)-1 do begin
   tmp(bou)=(horb1x(bou)*spinp1x(bou) $
             +horb1y(bou)*spinp1y(bou) $
             +horb1z(bou)*spinp1z(bou)) $
            /(sqrt(horb1x(bou)^2+horb1y(bou)^2+horb1z(bou)^2) $
              *sqrt(spinp1x(bou)^2+spinp1y(bou)^2+spinp1z(bou)^2))
   if tmp(bou) le 1 then $ 
      oblp1m(bou) = acos(tmp(bou))*180.d0/!Pi
   if tmp(bou) gt 1 then oblp1m(bou) = 0.0d0
   
   tmp(bou)=(horb2x(bou)*spinp2x(bou) $
             +horb2y(bou)*spinp2y(bou) $
             +horb2z(bou)*spinp2z(bou)) $
            /(sqrt(horb2x(bou)^2+horb2y(bou)^2+horb2z(bou)^2) $
              *sqrt(spinp2x(bou)^2+spinp2y(bou)^2+spinp2z(bou)^2))
   if tmp(bou) le 1 then $ 
      oblp2m(bou) = acos(tmp(bou))*180.d0/!Pi 
   if tmp(bou) gt 1 then oblp2m(bou) = 0.0d0
   
   tmp(bou)=(horb1x(bou)*spinstx(bou) $
             +horb1y(bou)*spinsty(bou) $
             +horb1z(bou)*spinstz(bou)) $
            /(sqrt(horb1x(bou)^2+horb1y(bou)^2+horb1z(bou)^2) $
              *sqrt(spinstx(bou)^2+spinsty(bou)^2+spinstz(bou)^2))
   if tmp(bou) le 1 then $ 
      obls1m(bou) = acos(tmp(bou))*180.d0/!Pi
   if tmp(bou) gt 1 then obls1m(bou) = 0.0d0
   
   tmp(bou)=(horb2x(bou)*spinstx(bou) $
             +horb2y(bou)*spinsty(bou) $
             +horb2z(bou)*spinstz(bou)) $
            /(sqrt(horb2x(bou)^2+horb2y(bou)^2+horb2z(bou)^2) $
              *sqrt(spinstx(bou)^2+spinsty(bou)^2+spinstz(bou)^2))
   if tmp(bou) le 1 then $ 
      obls2m(bou) = acos(tmp(bou))*180.d0/!Pi
   if tmp(bou) gt 1 then obls2m(bou) = 0.0d0
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
      indicend(1,0) = i
      indicend(1,1) = i
   endif else begin
      break
   endelse
endfor   
END
