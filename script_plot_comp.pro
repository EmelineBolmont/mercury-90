;PRO script_plot
; ratio between orbital angular momentum and spin angular momentum for
; pseudo synchronisation
Function alpha,e
  return, (1+15/2.*e^2+45/8.*e^4+5/16.*e^6)* $
     1./(1+3*e^2+3/8.*e^4)*1./(1-e^2)^1.5
end
; pseudo synchronisation rotation without the a dependance:
Function pseudorot,e,G,Mp,Ms
  return, alpha(e)*sqrt(G*(Mp+Ms))
end
; Constants 
G         =  6.6742367d-11            ;m^3.kg^-1.s^-2
AU        =  1.49588d11               ;m
yr        =  365.25*24*3600           ;s
day       =  24*3600.                 ;s
hr        =  3600.                    ;s

Msun      =  1.98892d30               ;kg
Mjup      =  9.5511d-4 * Msun         ;kg
Mearth    =  3.d-6 * Msun             ;kg
Rjup      =  69173.d3                 ;m
Rearth    =  6371.0d3                 ;m
Rsun      =  6.96d8                   ;m

Mp = [1.,10.]*Mearth
Ms = 0.08*Msun

device,decomposed=0
loadct,13

; Now, we plot everything !!

;; ops,file='11267_evolution_ae.eps',form=1;,/landscape
;; ops,file='11267_evolution_obl_inc_spin.eps',form=1;,/landscape

;;!p.multi=[0,2,2]
!P.Font=1
!p.thick = 6
!x.thick = 6
!y.thick = 6

Tinf    = 1.d2;8d6;
Tsup    = 1.d6;1.1d7;
amin    = 5d-3;5.d-3;0.019;
amax    = 3d-1;100d-3;6.2d-2;1.0d-2;0.021
emin    = 1.0d-4;8	
emax    = 1.0d0;2
oblmin  = -5d-2
oblmax  = 2d-1
incmin  = 1d-1
incmax  = 5d0
Trotmin = 10;20
Trotmax = 600;200

indcolor = 200
incolor  = 50;50

ae = 1
if ae eq 1 then begin
;! semi-major axis with respect to t
multiplot,[1,2],ygap=0.01
plot,tb(0,*),ab(0,*) $
     ,/nodata $
     ,xrange=[Tinf,Tsup],yrange=[amin,amax] $
     ,charsize=1.8,charthick=3 $
     ;; ,title='SMA' $
     ;; ,xtitle='t (years)' $
     ,ytitle='a (AU)' $
     ;; ,ytickname=['1.6e-2','2e-2','2.4e-2','2.6e-2','3d-2' $
     ;;             ,'3.4e-2'] $
     ;; ,yticks=5 $
     ;; ,ytickv=[1.6d-2,2d-2,2.4d-2,2.6d-2,3d-2,3.4d-2] $;
     ,xGRIDSTYLE=1,xTICKLEN=0.5 $
     ,xstyle=1,ystyle=1 $
     ,/xlog,/ylog

;;! Mercury
for i = 0,nbp-1 do begin
   oplot,tb(i,0:indicend(0,i)),ab(i,0:indicend(0,i)),color=incolor+i*indcolor,thick=7
   oplot,tb(i,0:indicend(0,i)),ab(i,0:indicend(0,i))*(1-eb(i,0:indicend(0,i))),color=incolor+i*indcolor,thick=7,linestyle=2
   oplot,tb(i,0:indicend(0,i)),ab(i,0:indicend(0,i))*(1+eb(i,0:indicend(0,i))),color=incolor+i*indcolor,thick=7,linestyle=2
endfor

;; oplot,toto1(*),Rst(*)*Rsun/AU,thick=5,color=255
oplot,toto1(*), $
      (G*(Ms+Mp(0)))^(1/3.)*(sqrt(spinstz(*)^2+spinstx(*)^2+spinsty(*)^2)/86400.d0)^(-2./3.)/AU $
            ,thick=7;,color=255
            
print,'Porb =',2.d0*!Pi/sqrt(G*(Ms+Mp))*(ab(0,0)*AU)^(3.d0/2.d0)*1.0d0/(86400.d0),' days'
print,'min(e) =',min(eb(0,*)),min(eb(1,*))
print,'max(e) =',max(eb(0,*)),max(eb(1,*))

;;! IDL
for i = 0,nbp-1 do begin
   oplot,ti(i,*),ai(i,*) $
      ,thick=5,linestyle=2,color=incolor+i*indcolor
   oplot,ti(i,*),(G*(Ms+Mp(0)))^(1/3.) $
      *(rotsi(i,*))^(-2./3.)/AU $
      ,thick=5,linestyle=2;,color=255
endfor
;;! oplot,ti(nbp-1,*),Rsi(nbp-1,*),thick=5,color=255,linestyle=2    ;,thick=5,linestyle=2
;;! print,'Rsi=',Rsi(nbp-1,0)
;; xyouts,4.1d7,1.55d-2,'1 M!dearth!n,'+greek('sigma')+'!dp!n' $
;;        ,charsize=1.8,charthick=3,charthick=1.5,color=incolor
;; xyouts,4.1d7,2.3d-2,'10 M!dearth!n, 100 '+greek('sigma')+'!dp!n' $
;;        ,charsize=1.8,charthick=3,charthick=1.5,color=incolor+1*indcolor
;; xyouts,1d5,2.55d-2,'0.08 M!dsun!n, '+greek('sigma')+'!dBD!n' $
;;        ,charsize=1.8,charthick=3,charthick=1.5


;;! eccentricity with respect to t
multiplot
plot,tb(0,*),eb(0,*) $
     ,/nodata $
     ,xrange=[Tinf,Tsup],yrange=[emin,emax] $
     ,charsize=1.8,charthick=3 $
     ;; ,title='Eccentricity' $
     ,xtitle='Age of BD - t!d0!n (years)' $
     ,ytitle='e' $
     ,xGRIDSTYLE=1,xTICKLEN=0.5 $
     ,xstyle=1,ystyle=1 $
     ,/ylog,/xlog
for i = 0,nbp-1 do begin
   oplot,tb(i,0:indicend(0,i)),eb(i,0:indicend(0,i)),color=incolor+i*indcolor,thick=7,linestyle=0;,psym=2
   oplot,ti(i,*),ei(i,*) $
      ,linestyle=2,color=incolor+i*indcolor,thick=5
endfor
multiplot,/reset
endif else begin
;! Inclination
multiplot,[1,3],ygap=0.01
plot,tb(0,*),incb(0,*) $
     ,/nodata $
     ,xrange=[Tinf,Tsup],yrange=[oblmin,oblmax] $
     ,charsize=1.8,charthick=3 $
     ,xGRIDSTYLE=1,xTICKLEN=0.5 $
     ;; ,title='Spin' $
     ;; ,xtitle='t (years)' $
     ,ytitle='obl (deg)' $
     ,xstyle=1,ystyle=1 $
     ,/xlog;,/ylog

;! Mercury
oplot,toto1(0:indicend(1,0)),oblp1m(0:indicend(1,0)),color=incolor+0*indcolor,thick=7
;;! print,'oblp1m(0:7) =',oblp1m(0:7)
;;! print,'oblp2m(0:7) =',oblp2m(0:7)
oplot,toto1(0:indicend(1,0)),oblp2m(0:indicend(1,0)),color=incolor+1*indcolor,thick=7

;! idl
for i=0,nbp-1 do begin
 oplot,ti(i,*),oblpi(i,*) $
       ,linestyle=2,color=incolor+i*indcolor,thick=5
endfor

multiplot
plot,tb(0,*),incb(0,*) $
     ,/nodata $
     ,xrange=[Tinf,Tsup],yrange=[incmin,incmax] $
     ,charsize=1.8,charthick=3 $
     ,xGRIDSTYLE=1,xTICKLEN=0.5 $
     ;; ,title='Spin' $
     ;; ,xtitle='t (years)' $
     ,ytitle='inc (deg)' $
     ,xstyle=1,ystyle=1 $
     ,/xlog,/ylog

;! Mercury
oplot,toto1(0:indicend(1,0)),obls1m(0:indicend(1,0)),color=incolor+0*indcolor,thick=7
oplot,toto1(0:indicend(1,0)),obls2m(0:indicend(1,0)),color=incolor+1*indcolor,thick=7
;! idl
for i=0,nbp-1 do begin
 oplot,ti(i,*),oblsi(i,*) $
       ,linestyle=2,color=incolor+i*indcolor,thick=5
endfor
;; for i = 0,nbp-1 do begin
;;    oplot,tb(i,*),incb(i,*),color=255,thick=5,linestyle=5
;; endfor

;! Rotation of bodies
multiplot
plot,tb(0,*),ab(0,*) $
     ,/nodata $
     ,xrange=[Tinf,Tsup],yrange=[Trotmin,Trotmax] $
     ,charsize=1.8,charthick=3 $
     ;; ,title='Rotation Period' $
     ,xtitle='Age of BD - t!d0!n (years)' $
     ,ytitle='P (hr)' $
     ,xGRIDSTYLE=1,xTICKLEN=0.5 $
     ;; ,ytickname=['20','30','40','50','70','100' $
     ;;             ] $
     ;; ,yticks=5 $
     ;; ,ytickv=[20,30,40,50,70,100] $
     ,xstyle=1,ystyle=1 $
     ,/ylog,/xlog
 
;! Mercury
i=0
oplot,toto1(0:indicend(1,i)),24.*2*!Pi $
   /(sqrt(spinp1x(0:indicend(1,i))^2+spinp1y(0:indicend(1,i))^2+spinp1z(0:indicend(1,i))^2)) $
      ,thick=7,linestyle=0,color=incolor+0*indcolor
oplot,toto1(0:indicend(1,i)),24.*2*!Pi $
	/(sqrt(spinp2x(0:indicend(1,i))^2+spinp2y(0:indicend(1,i))^2+spinp2z(0:indicend(1,i))^2)) $
      ,thick=7,linestyle=0,color=incolor+1*indcolor      
oplot,toto1(0:indicend(1,i)),24.*2*!Pi $
   /(sqrt(spinstx(0:indicend(1,i))^2+spinsty(0:indicend(1,i))^2+spinstz(0:indicend(1,i))^2)) $
      ,thick=7,linestyle=0;,color=incolor+1*indcolor

;! idl
for i = 0,nbp-1 do begin
   oplot,ti(i,*),2*!Pi/(rotpi(i,*)*hr) $
         ,linestyle=2,color=incolor+i*indcolor,thick=5
   oplot,ti(i,*),2*!Pi/(rotsi(i,*)*hr) $
         ,linestyle=2,thick=5
   ;; oplot,ti(i,*),2*!Pi/(pseudorot(ei(i,*),G,Mp(i),Ms)*(ai(i,*)*AU)^(-3./2.)*hr) $
   ;;       ,color=255,thick=3,linestyle=5     
   ;; oplot,tb(i,*),2*!Pi/(pseudorot(eb(i,*),G,Mp(i),Ms)*(ab(i,*)*AU)^(-3./2.)*hr) $
   ;;       ,color=255,thick=3,linestyle=5
endfor
 
multiplot,/reset
endelse
END
