

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

; Cross product
Function cross_product,x,y,z,u,v,w
    a = y*w-z*v
    b = z*u-x*w
    c = x*v-y*u
    cross = [a,b,c]
    return,cross
end

; Orbital angular momentum for 4 bodies: one star, 3 planets
Function Lorb4bodies,ms,mp1,mp2,mp3,x1,y1,z1,u1,v1,w1,x2,y2,z2,u2,v2,w2,x3,y3,z3,u3,v3,w3

    a1 = ms*mp1/(ms+mp1)
    a2 = (ms+mp1)*mp2/(ms+mp1+mp2)
    a3 = (ms+mp1+mp2)*mp3/(ms+mp1+mp2+mp3)

    r1    = [x1,y1,z1]
    v1    = [u1,v1,w1]

    r23  = [x2,y2,z2]-mp1/(ms+mp1)*[x1,y1,z1]
    v23  = [u2,v2,w2]-mp1/(ms+mp1)*[u1,v1,w1]

    r234 = [x3,y3,z3]-mp1/(ms+mp1)*[x1,y1,z1]-mp2/(ms+mp1+mp2)*r23
    v234 = [u3,v3,w3]-mp1/(ms+mp1)*[u1,v1,w1]-mp2/(ms+mp1+mp2)*v23
    
    r1xv1     = cross_product(x1,y1,z1,u1,v1,w1)
    r23xv23   = cross_product(r23(0),r23(1),r23(2),v23(0),v23(1),v23(2))
    r234xv234 = cross_product(r234(0),r234(1),r234(2),v234(0),v234(1),v234(2))

    Lorb = a1 * r1xv1 + a2 * r23xv23 + a3 * r234xv234
    return,Lorb
end

