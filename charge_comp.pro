rg2p      = [3.308d-1,2.54d-1]
Msun      =  1.98892d30               ;kg

;! Mass star:
Ms = 0.08*Msun
;! Number planets
nbp = 2
;! Number of planets tidally evolving
n_tid = 2
;! if terrestrial: jupiter =0, if gas giant: jupiter=1
jupiter = [0,0]

;! If comparison with IDL simulations 1, if not 0
idl = 1
nbp_idl =2
nline_idl = dblarr(nbp_idl)
filename_idl = strarr(nbp_idl)
filename_idl(0) = 'datatides_18.0000_0.00000_12_2_3_3_0.0261799_0.00000_0.dat'
filename_idl(1) = 'datatides_100.000_50.0000_12_5_3_4_0.0174533_0.00000_0.dat'

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

filenames = 'spins.dat'
print,filenames
readcol,filenames,sss,toto1,spinstx,spinsty,spinstz,Rst,format='A,F,F,F,F,F'
spinpx = dblarr(n_tid,n_elements(toto1))
spinpy = dblarr(n_tid,n_elements(toto1))
spinpz = dblarr(n_tid,n_elements(toto1))
horbx  = dblarr(n_tid,n_elements(toto1))
horby  = dblarr(n_tid,n_elements(toto1))
horbZ  = dblarr(n_tid,n_elements(toto1))

for i=0,n_tid-1 do begin 
   filenamep = 'spinp'+strtrim(i+1,2)+'.dat'
   print,filenamep
   readcol,filenamep,ppp,toto1,spinp1x,spinp1y,spinp1z,format='A,F,F,F,F'
   spinpx(i,*) = spinp1x & spinpy(i,*) = spinp1y & spinpz(i,*) = spinp1z
   
   filenameh = 'horb'+strtrim(i+1,2)+'.dat'
   print,filenameh
   readcol,filenameh,hhh,toto1,horb1x,horb1y,horb1z,format='A,F,F,F,F'
   horbx(i,*) = horb1x & horby(i,*) = horb1y & horbz(i,*) = horb1z  
endfor

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
 
endif

if n_tid ge 1 then begin
   ;! Obliquities calculations
   tmp     = dblarr(n_elements(horb1x))
   oblpm  = dblarr(n_tid,n_elements(horb1x))
   oblsm  = dblarr(n_tid,n_elements(horb1x))
   spinp  = dblarr(n_tid,n_elements(horb1x))
   
   for i = 0,n_tid-1 do begin
      for bou = 0,n_elements(horb1x)-1 do begin
      
         tmp(bou)=(horbx(i,bou)*spinpx(i,bou) $
                   +horby(i,bou)*spinpy(i,bou) $
                   +horbz(i,bou)*spinpz(i,bou)) $
                  /(sqrt(horbx(i,bou)^2+horby(i,bou)^2+horbz(i,bou)^2) $
                    *sqrt(spinpx(i,bou)^2+spinpy(i,bou)^2+spinpz(i,bou)^2))
         if abs(tmp(bou)) le 1.d0 then $ 
            oblpm(i,bou) = acos(tmp(bou))*180.d0/!Pi
         if abs(tmp(bou)) gt 1.d0 then oblpm(i,bou) = 1.0d-6
        
         tmp(bou)=(horbx(bou)*spinstx(bou) $
                   +horby(bou)*spinsty(bou) $
                   +horbz(bou)*spinstz(bou)) $
                  /(sqrt(horbx(bou)^2+horby(bou)^2+horbz(bou)^2) $
                    *sqrt(spinstx(bou)^2+spinsty(bou)^2+spinstz(bou)^2))
         if abs(tmp(bou)) le 1.d0 then $ 
            oblsm(i,bou) = acos(tmp(bou))*180.d0/!Pi
         if abs(tmp(bou)) gt 1.d0 then oblsm(i,bou) = 1.0d-6
         
         spinp(i,bou) = sqrt(spinpx(i,bou)^2 $
                +spinpy(i,bou)^2+spinpz(i,bou)^2)
        
      endfor
   endfor
endif

; Angular momentum calculation
horb_vec = dblarr(nbp,n_elements(horb1x))
horbx= dblarr(nbp,n_elements(horb1x))
horby= dblarr(nbp,n_elements(horb1x))
horbz= dblarr(nbp,n_elements(horb1x))
horb = dblarr(n_elements(horb1x))
spinst = dblarr(n_elements(horb1x))
momspin = dblarr(n_tid,n_elements(horb1x))
momspitot = dblarr(n_elements(horb1x))
momstar = dblarr(n_elements(horb1x))

Ip   = dblarr(n_tid)
if Ms le 0.3*Msun then Is = 0.254d0*Ms*(Rst*Rsun)^2
if Ms gt 0.3*Msun then Is = 5.9d-2*Ms*(Rst*Rsun)^2
Isi = rg2si*Ms*(Rsi*Rsun)^2

for j=0,nbp-1 do begin
   for i=0,n_elements(horb1x)-1 do begin
      horbx(j,i)=(yb(j,i)*wb(j,i)-zb(j,i)*vb(j,i))
      horby(j,i)=(zb(j,i)*ub(j,i)-xb(j,i)*wb(j,i))
      horbz(j,i)=(xb(j,i)*vb(j,i)-yb(j,i)*ub(j,i))
      
      ; With x,y,z,u,v,w
      horb_vec(j,i) = Ms*mb(j,i)*Msun/(Ms+mb(j,i)*Msun)*sqrt(horbx(j,i)^2+horby(j,i)^2+horbz(j,i)^2)*(AU^2/day) ; kg.m^2.s-1
      horb(i) = horb(i)+horb_vec(j,i)
      spinst(i) = sqrt(spinstx(i)^2+spinsty(i)^2+spinstz(i)^2)/day ; s-1
   endfor
endfor
for i=0,n_elements(horb1x)-1 do begin
   momstar(i)=Is(i)*spinst(i) 
endfor
 
for j=0,n_tid-1 do begin
   Ip(j) = 3.308d-1*mb(j,0)*Msun*(1.01034d0*Rearth)^2
   for i=0,n_elements(horb1x)-1 do begin
      momspin(j,i) = Ip(j)*spinp(j,i)/day
   endfor
endfor
for i=0,n_elements(horb1x)-1 do begin
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
if n_tid ge 1 then begin
   for i = 0,n_elements(toto1(*))-1 do begin
      if toto1(i) le 1.d10 then begin
         indicend(1,*) = i
      endif else begin
         break
      endelse
   endfor 
endif

END
