

; The 2 eccentricity dependant factors in the equation in a 
Function Na1, e
  return, (1.d0+31.d0/2.d0*e^2+255.d0/8.d0*e^4 $
           +185.d0/16.d0*e^6+25.d0/64.d0*e^8)/(1.d0-e^2)^(15.d0/2.d0)
end
Function Na2, e
  return, (1d0+15/2.d0*e^2+45/8.d0*e^4+5/16.d0*e^6)/(1.d0-e^2)^6
end
Function No2, e
  return, (1.d0+3.d0*e^2+3.d0/8.d0*e^4)/(1.d0-e^2)^5
end
; Mean orbital angular velocity without the a dependance         (m^3/2.s-1)
Function norb,G,Mp,Ms
  return, sqrt(G)*sqrt(Mp+Ms)
end

;Jeremy's Ki factor : 
Function Kplan,k2deltat_plan,G,Mp,Ms,Rp,a
  return, 3.d0/2.d0*k2deltat_plan*(G*Mp^2/Rp)*(Ms/Mp)^2 $
          *(Rp/a)^6*(norb(G,Mp,Ms)*a^(-1.5d0))^2.
end

;Energy due to tidal dissipation in bodies Jeremy's exact formula
Function enerdot,a,e,rotp,oblp,G,Mp,Ms,Rp,k2deltat_plan
  return, 2.d0* Kplan(k2deltat_plan,G,Mp,Ms,Rp,a) $
          * (Na1(e)-2.d0*Na2(e)*cos(oblp)*(rotp/(norb(G,Mp,Ms)*a^(-1.5d0))) $
             + (1.0d0+cos(oblp)^2)/2.d0*No2(e)*sqrt(1.d0-e^2)*(rotp/(norb(G,Mp,Ms)*a^(-1.5d0)))^2)
end

; ratio between orbital angular momentum and spin angular momentum for
; pseudo synchronisation
Function alpha,e
  return, (1.d0+15.d0/2.d0*e^2+45.d0/8.d0*e^4+5.d0/16.d0*e^6)* $
     1.d0/(1.d0+3.d0*e^2+3.d0/8.d0*e^4)*1.d0/(1.d0-e^2)^1.5d0
end

; pseudo synchronisation rotation without the a dependance:
Function pseudorot,e,G,Mp,Ms
  return, alpha(e)*sqrt(G*(Mp+Ms))
end
