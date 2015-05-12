
Tinf    = 1.d2     &  Tsup    = 1.d7

amin    = 1d-3     &  amax    = 1d0

emin    = 1.0d-6   &  emax    = 1.0d0

flxmin  = 1d-6     &  flxmax  = 1d5

oblmin  = 1d-3     &  oblmax  = 4.5d1

incmin  = 1d-3     &  incmax  = 4.5d1

Trotmin = 10       &  Trotmax = 10000

pr_angl_min = 0    &  pr_angl_max = 180.

AMD_min = 1d-30    &  AMD_max = 1d0

cons_min = 1d-10   &  cons_max = 1d-2

indcolor = 200   &  incolor  = 255
idlcol   = 200   &  idlicol  = 255

; Obliquity, inclination, rotation period: ae = 0
; semi-major axis, ecc, tidal flux: ae = 1
; Precession angle : ae = 2
; Angular momentum deficit: ae = 3
; Conservation of angular momentum: ae = 4
ae = 1
ops_plot=0

if ops_plot eq 1 then begin
   if ae eq 1 then ops,file='XXX_aef.eps',form=1;,/landscape
   if ae eq 0 then ops,file='XXX_ois.eps',form=1;,/landscape
   if ae eq 3 then ops,file='XXX_AMD.eps',form=1,/landscape
   if ae eq 4 then ops,file='XXX_H.eps',form=1,/landscape
endif

device,decomposed=0
loadct,13

; Now, we plot everything !!

!P.Font=1
!p.thick = 6
!x.thick = 6
!y.thick = 6

if ae eq 1 then begin
    ;! semi-major axis with respect to t
    if n_tid ge 1 then multiplot,[1,3],ygap=0.01
    if n_tid eq 0 then multiplot,[1,2],ygap=0.01
    plot,tb(0,*),ab(0,*) $
         ,/nodata $
         ,xrange=[Tinf,Tsup],yrange=[amin,amax] $
         ,charsize=1.8,charthick=3 $
         ,ytitle='Semi-major axis (AU)' $
         ;; ,ytickname=['1.6e-2','2e-2','2.4e-2','2.6e-2','3d-2' $
         ;;             ,'3.4e-2'] $
         ;; ,yticks=5 $
         ;; ,ytickv=[1.6d-2,2d-2,2.4d-2,2.6d-2,3d-2,3.4d-2] $;
         ,xGRIDSTYLE=1,xTICKLEN=0.5 $
         ,xstyle=1,ystyle=1 $
         ,/xlog,/ylog

    ;;! Mercury
    for j = 0,nbp-1 do begin
       ;i = nbp-1-j
       i = j
       oplot,tb(i,*),ab(i,*),color=incolor-i*indcolor,thick=5
       oplot,tb(i,*),ab(i,*)*(1-eb(i,*)),color=incolor-i*indcolor,thick=2,linestyle=2
       oplot,tb(i,*),ab(i,*)*(1+eb(i,*)),color=incolor-i*indcolor,thick=2,linestyle=2
    endfor

    ;; oplot,toto1(*),Rst(*)*Rsun/AU,thick=5,color=255
    if n_tid ge 1 then $
       oplot,toto1(*), $
             (G*(Ms+mb(0,*)*Msun))^(1/3.)*(sqrt(spinstz(*)^2+spinstx(*)^2+spinsty(*)^2)/86400.d0)^(-2./3.)/AU $
                   ,thick=7;,color=255

    print,'Porb =',2.d0*!Pi/sqrt(G*(Ms+mb(0,0)*msun))*(ab(0,0)*AU)^(3.d0/2.d0)*1.0d0/(86400.d0),' days'
    print,'min(e) =',min(eb(0,*));!,min(eb(1,*))
    print,'max(e) =',max(eb(0,*));!,max(eb(1,*))

    ;;! IDL
    if idl eq 1 then begin
       for i = 0,nbp_idl-1 do begin
          oplot,ti(i,*),ai(i,*) $
             ,color=idlicol-i*idlcol,thick=4,linestyle=2
          oplot,ti(i,*),(G*(Ms+mb(i,*)*msun))^(1/3.) $
             *(rotsi(i,*))^(-2./3.)/AU $
             ,thick=4,linestyle=2;,color=255
       endfor
    endif

    ;;! eccentricity with respect to t
    multiplot
    if n_tid eq 0 then plot,tb(0,*),eb(0,*) $
         ,/nodata $
         ,xrange=[Tinf,Tsup],yrange=[emin,emax] $
         ,charsize=1.8,charthick=3 $
         ,xtitle='Time (years)' $
         ;;,xtitle='Age of BD - t!d0!n (years)' $
         ,ytitle='Eccentricity' $
         ,xGRIDSTYLE=1,xTICKLEN=0.5 $
         ,xstyle=1,ystyle=1 $
         ,/xlog,/ylog
    if n_tid ge 1 then plot,tb(0,*),eb(0,*) $
         ,/nodata $
         ,xrange=[Tinf,Tsup],yrange=[emin,emax] $
         ,charsize=1.8,charthick=3 $
         ,ytitle='Eccentricity' $
         ,xGRIDSTYLE=1,xTICKLEN=0.5 $
         ,xstyle=1,ystyle=1 $
         ,/xlog,/ylog

    for j = 0,nbp-1 do begin
       ;i = nbp-1-j
       i = j
       oplot,tb(i,*),eb(i,*),color=incolor-i*indcolor,thick=5,linestyle=0;,psym=2
    endfor

    if idl eq 1 then begin
       for i = 0,nbp_idl-1 do begin
          oplot,ti(i,*),ei(i,*) $
             ,color=idlicol-i*idlcol,thick=4,linestyle=2
       endfor
    endif

    if n_tid ge 1 then begin
        ;;! Tidal flux with respect to t
        multiplot
        plot,tb(0,*),tidalflux(0,*) $
             ,/nodata $
             ,xrange=[Tinf,Tsup],yrange=[flxmin,flxmax] $
             ,charsize=1.8,charthick=3 $
             ;,xtitle='Age of BD - t!d0!n (years)' $
             ,xtitle='Time (years)' $
             ,ytitle='Tidal flux (W/m!u2!n)' $
             ,xGRIDSTYLE=1,xTICKLEN=0.5 $
             ,xstyle=1,ystyle=1 $
             ,/xlog,/ylog
        for j = 0,nbp-1 do begin
           ;i = nbp-1-j
           i = j
           oplot,tb(i,*),inst_tidalflux(i,*),color=incolor-i*indcolor,thick=1,linestyle=2;,psym=2
           oplot,tb(i,*),tidalflux(i,*),color=incolor-i*indcolor,thick=5,linestyle=0;,psym=2
        endfor

        if idl eq 1 then begin
           for i = 0,nbp_idl-1 do begin
              oplot,ti(i,*),tidefluxi(i,*) $
                 ,color=idlicol-i*idlcol,thick=4,linestyle=2
           endfor
        endif

        ;; oplot,[Tinf,Tsup],[2.4,2.4],linestyle=2,thick=8,color=1
        ;; oplot,[Tinf,Tsup],[4.8,4.8],linestyle=2,thick=8,color=1
        ;; oplot,[Tinf,Tsup],[300,300],linestyle=3,thick=8,color=1
    endif
    multiplot,/reset
endif 

if ae eq 0 then begin
    ;! Inclination
    multiplot,[1,3],ygap=0.01
    plot,tb(0,*),incb(0,*) $
         ,/nodata $
         ,xrange=[Tinf,Tsup],yrange=[oblmin,oblmax] $
         ,charsize=1.8,charthick=3 $
         ,xGRIDSTYLE=1,xTICKLEN=0.5 $
         ,ytitle='Obliquity (deg)' $
         ,xstyle=1,ystyle=1 $
         ,/xlog,/ylog

    ;! Mercury
    for j=0,n_tid-1 do begin
       ;i = n_tid-1-j
       i = j
       oplot,toto1(*),oblpm(i,*),color=incolor-i*indcolor,thick=5
    endfor

    ;! idl
    if idl eq 1 then begin
       for i=0,nbp_idl-1 do begin
        oplot,ti(i,*),oblpi(i,*) $
              ,color=idlicol-i*idlcol,thick=4,linestyle=2
       endfor
    endif

    multiplot
    plot,tb(0,*),incb(0,*) $
         ,/nodata $
         ,xrange=[Tinf,Tsup],yrange=[incmin,incmax] $
         ,charsize=1.8,charthick=3 $
         ,xGRIDSTYLE=1,xTICKLEN=0.5 $
         ,ytitle='Inclination (deg)' $
         ,xstyle=1,ystyle=1 $
         ,/xlog,/ylog

    ;! Mercury
    for j=0,n_tid-1 do begin
       ;i = n_tid-1-j
       i = j
       oplot,toto1(*),oblsm(i,*),color=incolor-i*indcolor,thick=5
    endfor

    ;! idl
    if idl eq 1 then begin
       for i=0,nbp_idl-1 do begin
        oplot,ti(i,*),oblsi(i,*) $
              ,color=idlicol-i*idlcol,thick=4,linestyle=2
       endfor
    endif

    ;; for i = 0,nbp-1 do begin
    ;;    oplot,tb(i,*),incb(i,*),color=255,thick=5,linestyle=5
    ;; endfor


    ;! Rotation of bodies
    multiplot
    plot,tb(0,*),ab(0,*) $
         ,/nodata $
         ,xrange=[Tinf,Tsup],yrange=[Trotmin,Trotmax] $
         ,charsize=1.8,charthick=3 $
         ;,xtitle='Age of BD - t!d0!n (years)' $
         ,xtitle='Time (years)' $
         ,ytitle='Rotation period (hr)' $
         ,xGRIDSTYLE=1,xTICKLEN=0.5 $
         ;; ,ytickname=['20','30','40','50','70','100' $
         ;;             ] $
         ;; ,yticks=5 $
         ;; ,ytickv=[20,30,40,50,70,100] $
         ,xstyle=1,ystyle=1 $
         ,/xlog,/ylog

    ;! Mercury
    i=0
    for j=0,n_tid-1 do begin
       ;i = n_tid-1-j
       i = j
       oplot,toto1(*),2*!Pi/(spinp(i,*)*hr) $
          ,thick=5,linestyle=0,color=incolor-i*indcolor
    endfor

    oplot,toto1(*),2*!Pi/(spinst(*)*hr) $
          ,thick=4,linestyle=0;,color=incolor-1*indcolor

    for j=0,n_tid-1 do begin
       ;i = n_tid-1-j
       i = j
       oplot,tb(i,*),2.d0*!Pi $
          /(pseudorot(eb(i,*),G,mb(i,*)*msun,Ms)*(ab(i,*)*AU)^(-3./2.)*hr) $
    		,color=incolor-i*indcolor,thick=3,linestyle=5
    endfor

    ;! idl
    if idl eq 1 then begin
       for i = 0,nbp_idl-1 do begin
          oplot,ti(i,*),2*!Pi/(rotpi(i,*)*hr) $
                ,color=idlicol-i*idlcol,thick=4,linestyle=2
          oplot,ti(i,*),2*!Pi/(rotsi(i,*)*hr) $
                ,color=255,linestyle=2,thick=4
          ;; oplot,ti(i,*),2*!Pi/(pseudorot(ei(i,*),G,mb(i,*)*Msun,Ms)*(ai(i,*)*AU)^(-3./2.)*hr) $
          ;;       ,color=255,thick=3,linestyle=5     
          ;; oplot,tb(i,*),2*!Pi/(pseudorot(eb(i,*),G,mb(i,*)*Msun,Ms)*(ab(i,*)*AU)^(-3./2.)*hr) $
          ;;       ,color=255,thick=3,linestyle=5
       endfor

    endif
multiplot,/reset
endif

if ae eq 2 then begin
    plot,tb(0,*),ab(0,*) $
         ,/nodata $
         ,xrange=[Tinf,Tsup],yrange=[pr_angl_min,pr_angl_max] $
         ,charsize=1.8,charthick=3 $
         ;,xtitle='Age of BD - t!d0!n (years)' $
         ,xtitle='Time (years)' $
         ,ytitle='Precession angle (deg)' $
         ,xGRIDSTYLE=1,xTICKLEN=0.5 $
         ,xstyle=1,ystyle=1

    ;! Mercury
    i=0
    for j=0,n_tid-1 do begin
       ;i = n_tid-1-j
       i = j
       oplot,toto1(*),precession_angle(j,*) $
          ,thick=5,linestyle=0,color=incolor-i*indcolor
    endfor
endif

if ae eq 3 then begin
    plot,tb(0,*),ab(0,*) $
         ,/nodata $
         ,xrange=[Tinf,Tsup],yrange=[AMD_min,AMD_max] $
         ,charsize=2.5,charthick=3.5 $
         ,xtitle='Time (years)' $
         ,ytitle='Angular momentum deficit (Msun.AU!u2!n.yr!u-1!n)' $
         ,xGRIDSTYLE=1,xTICKLEN=0.5 $
         ,xstyle=1,ystyle=1 $
         ,/xlog,/ylog

    if n_tid ge 1 then oplot,toto1(*),AMD(*)*yr/(Msun*AU^2),linestyle=0,thick=5
    oplot,tb(0,*),AMD_sec(*)*yr/(Msun*AU^2),linestyle=5,thick=5
    if n_tid ge 1 then print,'Average of AMD',mean(AMD(*)*yr/(Msun*AU^2))
    print,'Average of AMD_sec',mean(AMD_sec(*)*yr/(Msun*AU^2))
endif

if ae eq 4 then begin
    plot,tb(0,*),ab(0,*) $
         ,/nodata $
         ,xrange=[Tinf,Tsup],yrange=[cons_min,cons_max] $
         ,charsize=2.5,charthick=3.5 $ ;,charsize=1.8,charthick=3
         ,xtitle='Time (years)' $
         ,ytitle=''+greek('Delta')+'L/L!d0!n' $
         ,xGRIDSTYLE=1,xTICKLEN=0.5 $
         ,xstyle=1,ystyle=1 $
         ,/xlog,/ylog

    oplot,toto1(*),abs((Lorb_tot_x(*)+momspitot_x(*)+momstar_x(*)-(Lorb_tot_x(0)+momspitot_x(0)+momstar_x(0))) $
        /(Lorb_tot_x(0)+momspitot_x(0)+momstar_x(0))),linestyle=2,thick=3,color=80
    oplot,toto1(*),abs((Lorb_tot_y(*)+momspitot_y(*)+momstar_y(*)-(Lorb_tot_y(0)+momspitot_y(0)+momstar_y(0))) $
        /(Lorb_tot_y(0)+momspitot_y(0)+momstar_y(0))),linestyle=2,thick=3,color=100
    oplot,toto1(*),abs((Lorb_tot_z(*)+momspitot_z(*)+momstar_z(*)-(Lorb_tot_z(0)+momspitot_z(0)+momstar_z(0))) $
        /(Lorb_tot_z(0)+momspitot_z(0)+momstar_z(0))),linestyle=2,thick=3,color=120
    
    oplot,toto1(*),abs((Lorb_tot(*)+momspitot(*)+momstar(*)-(Lorb_tot(0)+momspitot(0)+momstar(0))) $
        /(Lorb_tot(0)+momspitot(0)+momstar(0))),linestyle=0,thick=3

    oplot,toto1(*),abs((horb_v(*)+momspitot(*)+momstar(*)-(horb_v(0)+momspitot(0)+momstar(0))) $
          /(horb_v(0)+momspitot(0)+momstar(0))),linestyle=3,thick=5,color=incolor-0*indcolor
    oplot,toto1(*),abs((horb_s(*)+momspitot(*)+momstar(*)-(horb_s(0)+momspitot(0)+momstar(0))) $
          /(horb_s(0)+momspitot(0)+momstar(0))),linestyle=5,thick=5,color=incolor-0*indcolor

    if idl eq 1 then begin
       for i = 0,nbp_idl-1 do begin
          oplot,ti(i,*),abs((mb(i,*)*msun*Ms*sqrt(G*ai(i,*)*AU*(1-ei(i,*)^2)/(mb(i,*)*msun+Ms)) $
    	+Ip(i)*rotpi(i,*)+rg2si(i,*)*Ms*(Rsi(i,*)*AU)^2*rotsi(i,*)-(mb(i,0)*msun*Ms*sqrt(G*ai(i,0)*AU*(1-ei(i,0)^2)/(mb(i,0)*msun+Ms)) $
    	+Ip(i)*rotpi(i,0)+rg2si(i,0)*Ms*(Rsi(i,0)*AU)^2*rotsi(i,0)))/(mb(i,0)*msun*Ms*sqrt(G*ai(i,0)*AU*(1-ei(i,0)^2)/(mb(i,0)*msun+Ms)) $
    	+Ip(i)*rotpi(i,0)+rg2si(i,0)*Ms*(Rsi(i,0)*AU)^2*rotsi(i,0))) $
                ,color=idlicol-i*idlcol,thick=4,linestyle=2
       endfor
    endif

endif

END

