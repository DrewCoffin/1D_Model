;-----------------------------------------------------------------------
FUNCTION ion_ion, h1, h2

h_prime = h1 * h2 / sqrt(h1^2 + h2^2)

return, h_prime
END
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------

FUNCTION ion_neutral, hi, hn, z0

a = (hi^2 + hn^2)/(hi^2 * hn^2)
b = -2 * z0/hn^2
c = z0^2/hn^2

return, sqrt(1./a) * exp((b^2 - 4.*a*c)/(4.*a))
END
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
PRO cm3_reactions, r_ind, r_dep, h, n, ftint
common info,runt,dt,zoff,rdist
rootpi = sqrt(!dpi)

; Ion and Neutral Reactants
s_sp  = ion_neutral(h.sp,  h.s, zoff)
s_s2p = ion_neutral(h.s2p, h.s, zoff)
s_s3p = ion_neutral(h.s3p, h.s, zoff)
s_op  = ion_neutral(h.op,  h.s, zoff)
s_o2p = ion_neutral(h.o2p, h.s, zoff)

o_sp  = ion_neutral(h.sp,  h.o, zoff)
o_s2p = ion_neutral(h.s2p, h.o, zoff)
o_s3p = ion_neutral(h.s3p, h.o, zoff)
o_op  = ion_neutral(h.op,  h.o, zoff)
o_o2p = ion_neutral(h.o2p, h.o, zoff)

; Ion-Ion Reactants
sp_sp   = ion_ion(h.sp,  h.sp)
sp_s2p  = ion_ion(h.sp,  h.s2p)
sp_s3p  = ion_ion(h.sp,  h.s3p)
sp_op   = ion_ion(h.sp,  h.op)
sp_o2p  = ion_ion(h.sp,  h.o2p)

s2p_s2p = ion_ion(h.s2p, h.s2p)
s2p_s3p = ion_ion(h.s2p, h.s3p)
s2p_op  = ion_ion(h.s2p, h.op)
s2p_o2p = ion_ion(h.s2p, h.o2p)

s3p_s3p = ion_ion(h.s3p, h.s3p)
s3p_op  = ion_ion(h.s3p, h.op)
s3p_o2p = ion_ion(h.s3p, h.o2p)

op_op   = ion_ion(h.op,  h.op)
op_o2p  = ion_ion(h.op,  h.o2p)

o2p_o2p = ion_ion(h.o2p, h.o2p)

; Electron-Neutral Reactants
e_s = (1. * n.sp  * s_sp  + $
       2. * n.s2p * s_s2p + $
       3. * n.s3p * s_s3p + $
       1. * n.op  * s_op  + $
       2. * n.o2p * s_o2p)

e_o = (1. * n.sp  * o_sp  + $
       2. * n.s2p * o_s2p + $
       3. * n.s3p * o_s3p + $
       1. * n.op  * o_op  + $
       2. * n.o2p * o_o2p)

; Electron-Ion Reactants
e_sp = (1. * n.sp  * sp_sp  + $
        2. * n.s2p * sp_s2p + $
        3. * n.s3p * sp_s3p + $
        1. * n.op  * sp_op  + $
        2. * n.o2p * sp_o2p)

e_s2p = (1. * n.sp  * sp_s2p  + $
         2. * n.s2p * s2p_s2p + $
         3. * n.s3p * s2p_s3p + $
         1. * n.op  * s2p_op  + $
         2. * n.o2p * s2p_o2p)

e_s3p = (1. * n.sp  * sp_s3p  + $
         2. * n.s2p * s2p_s3p + $
         3. * n.s3p * s3p_s3p + $
         1. * n.op  * s3p_op  + $
         2. * n.o2p * s3p_o2p)

e_op = (1. * n.sp  * sp_op  + $
        2. * n.s2p * s2p_op + $
        3. * n.s3p * s3p_op + $
        1. * n.op  * op_op  + $
        2. * n.o2p * op_o2p)

e_o2p = (1. * n.sp  * sp_o2p  + $
         2. * n.s2p * s2p_o2p + $
         3. * n.s3p * s3p_o2p + $
         1. * n.op  * op_o2p  + $
         2. * n.o2p * o2p_o2p)

; Charge Exchange reactions
ftint.cx_k0  = rootpi * r_ind.cx_k0  * n.sp  * n.s2p * sp_s2p
ftint.cx_k1  = rootpi * r_ind.cx_k1  * n.s   * n.sp  * s_sp
ftint.cx_k2  = rootpi * r_ind.cx_k2  * n.s   * n.s2p * s_s2p
ftint.cx_k3  = rootpi * r_ind.cx_k3  * n.s   * n.s2p * s_s2p
ftint.cx_k4  = rootpi * r_ind.cx_k4  * n.s   * n.s3p * s_s3p
ftint.cx_k5  = rootpi * r_ind.cx_k5  * n.o   * n.op  * o_op
ftint.cx_k6  = rootpi * r_ind.cx_k6  * n.o   * n.o2p * o_o2p
ftint.cx_k7  = rootpi * r_ind.cx_k7  * n.o   * n.o2p * o_o2p
ftint.cx_k8  = rootpi * r_ind.cx_k8  * n.o   * n.sp  * o_sp
ftint.cx_k9  = rootpi * r_ind.cx_k9  * n.s   * n.op  * s_op
ftint.cx_k10 = rootpi * r_ind.cx_k10 * n.s   * n.o2p * s_o2p
ftint.cx_k11 = rootpi * r_ind.cx_k11 * n.s   * n.o2p * s_o2p
ftint.cx_k12 = rootpi * r_ind.cx_k12 * n.o   * n.s2p * o_s2p
ftint.cx_k13 = rootpi * r_ind.cx_k13 * n.o2p * n.sp  * sp_o2p
ftint.cx_k14 = rootpi * r_ind.cx_k14 * n.o   * n.s3p * o_s3p
ftint.cx_k15 = rootpi * r_ind.cx_k15 * n.o2p * n.s2p * s2p_o2p
ftint.cx_k16 = rootpi * r_ind.cx_k16 * n.s3p * n.sp  * sp_s3p

; Electron Impact Ionization (Thermal)
ftint.is   = rootpi * r_dep.is   * n.s   * e_s
ftint.isp  = rootpi * r_dep.isp  * n.sp  * e_sp
ftint.is2p = rootpi * r_dep.is2p * n.s2p * e_s2p
ftint.is3p = rootpi * r_dep.is3p * n.s3p * e_s3p

ftint.io   = rootpi * r_dep.io   * n.o   * e_o
ftint.iop  = rootpi * r_dep.iop  * n.op  * e_op
ftint.io2p = rootpi * r_dep.io2p * n.o2p * e_o2p

; Electron Impact Ionization (Hot Population)
ftint.ish   = rootpi * r_ind.ish   * n.elh * n.s   * h.s
ftint.isph  = rootpi * r_ind.isph  * n.elh * n.sp  * h.sp
ftint.is2ph = rootpi * r_ind.is2ph * n.elh * n.s2p * h.s2p
ftint.is3ph = rootpi * r_ind.is3ph * n.elh * n.s3p * h.s3p

ftint.ioh   = rootpi * r_ind.ioh   * n.elh * n.o   * h.o
ftint.ioph  = rootpi * r_ind.ioph  * n.elh * n.op  * h.op
ftint.io2ph = rootpi * r_ind.io2ph * n.elh * n.o2p * h.o2p

; Recombination (thermal only)
ftint.rsp = rootpi  * r_dep.rsp  * n.sp  * e_sp
ftint.rs2p = rootpi * r_dep.rs2p * n.s2p * e_s2p
ftint.rs3p = rootpi * r_dep.rs3p * n.s3p * e_s3p

ftint.rop = rootpi  * r_dep.rop  * n.op  * e_op
ftint.ro2p = rootpi * r_dep.ro2p * n.o2p * e_o2p

END
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
pro get_scale_heights,h,T,n
;all h in km
;-----------------------------------------------------------------------
if check_math(mask=211,/noclear) gt 0 then stop
@cm3_model_common

; The plasma isn't rigidly corotating anymore. It seems that instead
; of the corotation angular velocity, we should use the actual angular
; velocity as a function of System III longitude.
v_corot = 2. * !pi * rdist * !rjkm / (9.925 * 3600.)
v = v_corot - lag_const + lag_amp * cos((l3-lag_phase)*!dtor)
omega = v/(rdist*!rjkm) 

mp = 1.67e-27
ms = 32.*mp 
mo = 16.*mp

; Scale height from equation on p. 11,049 of Bagenal 1994
h.sp  = sqrt(2.0*T.sp *1.6e-19*(1+1.*T.el/T.sp)/(3.0*ms*omega^2))/1e3
h.s2p = sqrt(2.0*T.s2p*1.6e-19*(1+2.*T.el/T.s2p)/(3.0*ms*omega^2))/1e3
h.s3p = sqrt(2.0*T.s3p*1.6e-19*(1+3.*T.el/T.s3p)/(3.0*ms*omega^2))/1e3
;h.s4p = sqrt(2.0*T.s4p*1.6e-19*(1+4.*T.el/T.s4p)/(3.0*ms*omega^2))/1e3
h.op  = sqrt(2.0*T.op *1.6e-19*(1+1.*T.el/T.op)/(3.0*mo*omega^2))/1e3
h.o2p = sqrt(2.0*T.o2p*1.6e-19*(1+2.*T.el/T.o2p)/(3.0*mo*omega^2))/1e3

if check_math(mask=211,/noclear) gt 0 then stop
return
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
PRO lat_distribution, n, h, lat_dist
;-----------------------------------------------------------------------
common info,runt,dt,zoff,rdist

; Calculate the densities as a function of z in increments of .2 Rj.
unity = replicate(1., n_elements(lat_dist[0].z))

nsp0  = unity # n.sp
ns2p0 = unity # n.s2p
ns3p0 = unity # n.s3p
;ns4p0 = unity # n.s4p
nop0  = unity # n.op
no2p0 = unity # n.o2p

hsp  = unity # h.sp
hs2p = unity # h.s2p
hs3p = unity # h.s3p
;hs4p = unity # h.s4p
hop  = unity # h.op
ho2p = unity # h.o2p

; Ions
lat_dist.sp  = nsp0  * exp(-lat_dist.z^2/hsp^2)
lat_dist.s2p = ns2p0 * exp(-lat_dist.z^2/hs2p^2)
lat_dist.s3p = ns3p0 * exp(-lat_dist.z^2/hs3p^2)
;lat_dist.s4p = ns4p0 * exp(-lat_dist.z^2/hs4p^2)
lat_dist.op  = nop0  * exp(-lat_dist.z^2/hop^2)
lat_dist.o2p = no2p0 * exp(-lat_dist.z^2/ho2p^2)

; Electrons
lat_dist.el = lat_dist.sp + 2 * lat_dist.s2p + 3 * lat_dist.s3p + 4 * lat_dist.s4p + $
              lat_dist.op + 2 * lat_dist.o2p

END

;-----------------------------------------------------------------------
pro get_dependent_rates,r_dep,r_ind,n,T,h
;-----------------------------------------------------------------------
@cm3_model_common

;Electron impact ionization (Schreier et al. 1998)
;r_dep.is = 7.4e-8*sqrt(T.el_el)*exp(-10.36/T.el_el)     ;check this!!!
;r_dep.isp = 1.274e-8*sqrt(T.el_sp)*exp(-23.1/T.el_sp)
;r_dep.is2p = 5.8e-9*sqrt(T.el_s2p)*exp(-34.88/T.el_s2p)
;r_dep.io = 1.015e-8*sqrt(T.el_el)*exp(-14.54/T.el_el)
;r_dep.iop = 3.78e-9*sqrt(T.el_op)*exp(-33.42/T.el_op)

;Electron impact ionization (G. S. Voronov, 1997, ADNDT, 65, 1)
r_dep.is = cfit(16,16,T.el)
r_dep.isp = cfit(16,15,T.el)
r_dep.is2p = cfit(16,14,T.el)
r_dep.is3p = cfit(16,13,T.el)
r_dep.io = cfit(8,8,T.el)
r_dep.iop = cfit(8,7,T.el)
r_dep.io2p = cfit(8,6,T.el)

; Total electron recombination rates from the work of Sultana Nahar
;     S II, S III, O III: (ApJS, 101, 423 [1995]; ApJS 106, 213 [1996])
;     O I, O II:          (ApJS, 120, 131 [1999])
; Recombination rates for S I and S IV come from Mazzotta et al 1998, which
; provides formulae for the dielectronic recombination rate. The
; radiative recombination rate comes from Dima Verner's rrfit.f code,
; which uses: Shull & Van Steenberg, (1982, ApJS, 48, 95) for Sulfur 
;             Pequignot et al. (1991, A&A, 251, 680) for Oxygen
recombination_rates,T.el,rsp,rs2p,rs3p,rop,ro2p
r_dep.rsp  = rsp 
r_dep.rs2p = rs2p
r_dep.rs3p = rs3p
r_dep.rop  = rop
r_dep.ro2p = ro2p

fs = 1.0/(1.0+r_ind.o_to_s)
fo = r_ind.o_to_s/(1.0+r_ind.o_to_s)

; If transport_type=0, then the transport rate is inversely
; proprtional to the neutral source to the trans_exp power. If the
; trans_type keyword is set, then the transport rate will be inversely
; proportional to the ion density to the trans_exp power.
if keyword_set(trans_type) then BEGIN
   iondens=n.sp+n.s2p+n.s3p+n.s4p+n.op+n.o2p
   tau = tau0 * (iondens0/iondens)^trans_exp 
ENDIF ELSE tau = tau0 * (net_source0/net_source)^trans_exp
r_dep.transport = 1/(tau*8.64e4)

; The neutral source rate at the Equator is given by:
;  n(0) = N_tot / (sqrt(pi) x H)
; i.e. if the neutrals are more spread out along the field line, less
; will be produced at the equator, and more at higher latitudes. 
; Since scale heights are in km, and densities in cm need to multiply
; by 1e5.
r_ind.S_production = fs*net_source/(sqrt(!pi)*(h[0].s*1e5))
r_ind.O_production = fo*net_source/(sqrt(!pi)*(h[0].o*1e5))

end

;-----------------------------------------------------------------------
pro get_independent_rates,r_ind,T,h
;-----------------------------------------------------------------------
@cm3_model_common
; The hot electron temperature is fixed, so this should be done
; outside of any loops.
r_ind.ish = cfit(16,16,T[0].elh,c)
r_ind.isph = cfit(16,15,T[0].elh,c)
r_ind.is2ph = cfit(16,14,T[0].elh,c)
r_ind.is3ph = cfit(16,13,T[0].elh,c)
r_ind.ioh = cfit(8,8,T[0].elh,c)
r_ind.ioph = cfit(8,7,T[0].elh,c)
r_ind.io2ph = cfit(8,6,T[0].elh,c)

;Charge exchange
r_ind.cx_k0 =  1.2e-8       ;S+ + S++ -> S++ + S+ from McGrath & Johnson figure
;r_ind.cx_k0 =  8.1e-9       ;S+ + S++ -> S++ + S+ from Smith and Strobel
r_ind.cx_k1 =  2.4e-8
r_ind.cx_k2 =  3.0e-10
r_ind.cx_k3 =  7.8e-9
r_ind.cx_k4 =  1.32e-8
r_ind.cx_k5 =  1.32e-8
r_ind.cx_k6 =  5.2e-10
r_ind.cx_k7 =  5.4e-9
;r_ind.cx_k8 =  9.0e-11   ;Schreier
r_ind.cx_k8 =  6.0e-11    ;McGrath
r_ind.cx_k9 =  3.12e-9
r_ind.cx_k10 = 2.34e-8
r_ind.cx_k11 = 1.62e-8
r_ind.cx_k12 = 2.28e-9   ;new
;r_ind.cx_k12 = 7.3e-9  ;old
r_ind.cx_k13 = 1.38e-9
r_ind.cx_k14 = 1.92e-8
;r_ind.cx_k14 = 4.3e-9  ;Barbosa94, Johnson82
;r_ind.cx_k15 = 9e-9    ;Schreier for L=6
r_ind.cx_k15 = 9e-10    ;McGrath for L=6
;r_ind.cx_k15 = 2.84e-9 ;Johnson82
r_ind.cx_k16 = 3.6e-10  ;McGrath

end
;-----------------------------------------------------------------------

;====================Mass source/loss functions=========================

;-----------------------------------------------------------------------
function F_s,n,h,r_ind,r_dep, ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.s 

S = r_ind.S_production
L = (ftint.is + ftint.ish + ftint.cx_k1 + ftint.cx_k2 + ftint.cx_k3 + ftint.cx_k4 + $
     ftint.cx_k9 + ftint.cx_k10 + ftint.cx_k11)/rootpi_h

F_s = S-L

return,F_s
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function F_sp,n,h,r_ind,r_dep,ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.sp

S = ftint.is + ftint.ish + ftint.rs2p + 2. * ftint.cx_k2 + ftint.cx_k4 + ftint.cx_k9 + $
     ftint.cx_k10 + ftint.cx_k12
L = ftint.isp + ftint.isph + ftint.rsp + ftint.cx_k8 + ftint.cx_k13 + ftint.cx_k16

F_sp = (S - L)/rootpi_h - n.sp*r_dep.transport

return,F_sp
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function F_s2p,n,h,r_ind,r_dep,ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.s2p

S = ftint.isp + ftint.isph + ftint.rs3p + ftint.cx_k4 + ftint.cx_k11 + ftint.cx_k13 + $
     ftint.cx_k14 + 2 * ftint.cx_k16
L = ftint.is2p + ftint.is2ph + ftint.rs2p + ftint.cx_k2 + ftint.cx_k12 + ftint.cx_k15

F_s2p = (S-L) / rootpi_h - n.s2p * r_dep.transport

return,F_s2p
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function F_s3p,n,h,r_ind,r_dep,ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.s3p

S = ftint.is2p + ftint.is2ph + ftint.cx_k15
L = ftint.is3p + ftint.is3ph + ftint.rs3p + ftint.cx_k4 + ftint.cx_k14 + ftint.cx_k16

F_s3p = (S - L) / rootpi_h - n.s3p*r_dep.transport

return,F_s3p
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
;function F_s4p,n,h,r_ind,r_dep,ftint
;;-----------------------------------------------------------------------
;S = r_dep.is3p*ftavg_el(n.s3p,n,h.s3p,h,h.s3p) + $
;    r_ind.is3ph*ftavg_elh(n.s3p,n,h.s3p,h,h.s3p)

;L = n.s4p*r_dep.transport + $
;    r_dep.rs4p*ftavg_el(n.s4p,n,h.s4p,h,h.s4p)

;F_s4p = S - L

;return,F_s4p
;end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function F_o,n,h,r_ind,r_dep,ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.o

S = r_ind.O_production 
L = ftint.io + ftint.ioh + ftint.cx_k5 + ftint.cx_k6 + ftint.cx_k7 + ftint.cx_k8 + $
     ftint.cx_k12 + ftint.cx_k14

F_o = S - L/rootpi_h

return,F_o
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function F_op,n,h,r_ind,r_dep,ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.op

S = ftint.io + ftint.ioh + ftint.ro2p + 2 * ftint.cx_k6 + ftint.cx_k8 + ftint.cx_k10 + $
     ftint.cx_k11+ ftint.cx_k12+ ftint.cx_k13+ ftint.cx_k14+ ftint.cx_k15
L = ftint.iop + ftint.ioph + ftint.rop + ftint.cx_k9

F_op = (S-L)/rootpi_h-n.op*r_dep.transport

return,F_op
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function F_o2p,n,h,r_ind,r_dep,ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.o2p

S = ftint.iop + ftint.ioph
L = ftint.io2p + ftint.io2ph + ftint.ro2p + ftint.cx_k6 + ftint.cx_k10 + ftint.cx_k11 + $
     ftint.cx_k13 + ftint.cx_k15

F_o2p = (S-L)/rootpi_h - n.o2p*r_dep.transport

return,F_o2p
end
;-----------------------------------------------------------------------

;=======================Energy source/loss functions=====================


;-----------------------------------------------------------------------
function Tpu,m,L
;input mass as atomic number
;-----------------------------------------------------------------------
mp = 1.67262158d-27    ;kg

; 71492.*2*!pi*L/9.925/3600.=12.57
vrel = 12.57*L - 42.0/sqrt(L)   ;km/s
vrel *= 1e3  ;m/s

Tpu = (2./3.)*0.5*m*mp*vrel^2   ;J
Tpu /= 1.60217646e-19       ;eV

return,Tpu
end
;-----------------------------------------------------------------------

;----------------------------------------------------------------------- 
function coulomb_logarithm,type,n1,t1,z1,mu1,n2,t2,z2,mu2
; Returns the Coulomb logarithm according to NRL plasma formulary Rev
; 2007 p. 34.
; Type must be set to one of three values:
;    'e_e' for electron electron collisions
;    'e_i' for electron ion collisions
;    'i_i' for ion ion collions
;
; n -- number density
; t -- temperature of species (in eV)
; z -- charge of species
; mu -- mass (in terms of proton mass)
;----------------------------------------------------------------------- 
lambda=dblarr(size(n1,/dim) > 1)
case type of
    'e_e':BEGIN
       lambda = 23.5 - alog(sqrt(n1) * T1^(-5./4)) - sqrt(1e-5 + (alog(T1)-2.)^2 / 16)
;        under=where(t1 le 10,undercnt,complement=over,ncomplement=overcnt)
;        if undercnt gt 0 then lambda[under] = 23. - alog(sqrt(n1)*T1^(-3./2))
;        if overcnt gt 0 then lambda[over] = 24. - alog(sqrt(n1)/T1)
    end
    'e_i':begin
        e_to_p=5.4461702d-4         ; electron mass / proton mass
        under=where(T2*e_to_p/mu2 lt T1,undercnt,complement=over,ncomplement=overcnt)
        underunder=where(t1[under] le 10.*z2^2,underundercnt,complement=overunder,$
          ncomplement=overundercnt)
        if underundercnt gt 0 then lambda[under[underunder]] = 23.0 - $
          alog(sqrt(n1)*z2*T1^(-3./2.))
        if overundercnt gt 0 then lambda[under[overunder]] = 24.0 - alog(sqrt(n1)/T1)
        if overcnt gt 0 then lambda[over] = 30.0 - alog(sqrt(n2)*T2^(-3./2.)*Z2^2/mu2)
    end
    'i_i':lambda=23.-alog(z1*z2*(mu1+mu2)/(mu1*T2+mu2*T1)*sqrt(n1*z1^2/T1+n2*z2^2/T2))
endcase

return, lambda
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function nu_ie,mu,z,ni,Ti,nel,Tel
; mu is atomic mass in amu
; z is charge number
; for collisions of an ion and the thermal electron population
;-----------------------------------------------------------------------
if check_math(mask=211,/noclear) gt 0 then stop
me = 9.10938188d-28    ;g
mp = 1.67262158d-24    ;g
mi = mu*mp

unity = replicate(1., (size(ni))[1]) 
ti = unity#ti
tel = unity#tel

lambda=coulomb_logarithm('e_i',nel,tel,1,0,ni,ti,z,mu)

nuarr= 1.8e-19*sqrt(me*mi)*Z^2*nel*ni*lambda/(mi*tel+me*ti)^(3./2.)

nu_ie = total(nuarr*ni*nel,1)/total(ni*nel,1)
if check_math(mask=211,/noclear) gt 0 then stop
return,nu_ie
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------  
function nu_ieh,mu,z,ni,Ti,neh,Teh
; mu is atomic mass in amu
; z is charge number
; for collisions of an ion species and the hot electron population
;-----------------------------------------------------------------------  
if check_math(mask=211,/noclear) gt 0 then stop
me = 9.10938188d-28    ;g
mp = 1.67262158d-24    ;g

mi = mu*mp

unity = replicate(1., (size(ni))[1]) 
neh = unity#neh
ti = unity#ti
teh = unity#teh

lambda=coulomb_logarithm('e_i',neh,teh,1,0,ni,ti,z,mu)

nuarr = 1.8e-19*sqrt(me*mi)*(Z^2)*neh*ni*lambda/(mi*Teh + me*Ti)^(3./2.)

nu_ieh = total(nuarr*ni*neh,1)/total(ni*neh,1)
if check_math(mask=211,/noclear) gt 0 then stop
return,nu_ieh
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
function nu_ee, nel, Tel, neh, Teh
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------
if check_math(mask=211,/noclear) gt 0 then stop
me = 9.10938188d-28    ;g

unity = replicate(1., (size(nel))[1]) 
neh = unity#neh
Teh = unity#Teh
Tel = unity#Tel

; We are assuming that the hot electrons have a scale height that is
; effectively infinite such that the density of hot electrons remains
; constant with latitude above the torus equator. This would violate
; the quasi-neutrality condition, so we must also postulate an equal
; number of protons are along for the ride.

lambda=coulomb_logarithm('e_e', nel, Tel)

nuarr = 1.8e-19 * me * nel * neh * lambda / (me * Tel + me * Teh)^(3./2.)
nu_ee = total(nuarr * nel * neh, 1)/total(nel * neh, 1)

if check_math(mask=211,/noclear) gt 0 then stop
return,nu_ee
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function nu_ii,mu1,z1,n1,T1,mu2,z2,n2,T2
; m is atomic mass
; z is atomic number
; for collision of given ion on ion population
;-----------------------------------------------------------------------
if check_math(mask=211,/noclear) gt 0 then stop
mp = 1.67262158d-24    ;g

m1 = mu1*mp
m2 = mu2*mp

unity = replicate(1., (size(n1))[1]) 
T1 = unity#T1
T2 = unity#T2

lambda=coulomb_logarithm('i_i',n1,t1,z1,mu1,n2,t2,z2,mu2)

nuarr = 1.8e-19*sqrt(m1*m2)*(Z1^2)*(Z2^2)*n1*n2*lambda/(m1*T2 + m2*T1)^(3./2.)

nu_ii = total(nuarr * n1 * n2, 1)/total(n1 * n2, 1)
if check_math(mask=211,/noclear) gt 0 then stop
return,nu_ii
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
pro update_temp,n,nT,T
;-----------------------------------------------------------------------
T.sp  = nT.sp  / n.sp
T.s2p = nT.s2p / n.s2p
T.s3p = nT.s3p / n.s3p
;T.s4p = nT.s4p/ n.s4p
T.op  = nT.op  / n.op
T.o2p = nT.o2p / n.o2p
T.el  = nT.el  / n.el
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------  
FUNCTION ft_rad,lat_dist,T,r_ind
;-----------------------------------------------------------------------  
log5000=Alog10(5000.)

unity = replicate(1., n_elements(lat_dist[0].z))

; Only the core electron population radiates energy away, since the
; hot electron temperature is held constant. Technically, the presence
; of a hot electron population will affect the population states of
; the ions, thus altering the radiation. However, for a hot electron
; fraction of 0.0023 and a Te_hot of 50-500 eV, the hot electrons
; cause <2% increase in total radiated power.
nsp  = lat_dist.sp
ns2p = lat_dist.s2p
ns3p = lat_dist.s3p
ns4p = lat_dist.s4p
nop  = lat_dist.op
no2p = lat_dist.o2p
nel  = lat_dist.el

; The electron temperature along a field line is constant for a
; Maxwellian distribution. A non-thermal distribution would require
; modifying this line to reflect the changing temperature with latitude.
tempindex=unity#(100.*(1.+alog10(T.el))/log5000)
densindex=100.*alog10(nel)/log5000

; radiation rates
emisSp  = interpolate(r_ind.emisSp,tempindex,densindex)
emisS2p = interpolate(r_ind.emisS2p,tempindex,densindex)
emisS3p = interpolate(r_ind.emisS3p,tempindex,densindex)
emisOp  = interpolate(r_ind.emisOp,tempindex,densindex)
emisO2p = interpolate(r_ind.emisO2p,tempindex,densindex)

rad_sp  = emisSp*nSp
rad_s2p = emisS2p*nS2p
rad_s3p = emisS3p*nS3p
rad_Op  = emisOp*nOp
rad_O2p = emisO2p*nO2p

rad_tot=rad_sp+rad_s2p+rad_s3p+rad_op+rad_o2p

Lavg = total(rad_tot*nel,1)/total(nel,1)

return, Lavg
END
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function EF_el,n, T, h, r_ind, r_dep, lat_dist, nu, ftint
;-----------------------------------------------------------------------
Teq = nu.sp_el* (T.sp- T.el) + $
      nu.s2p_el*(T.s2p-T.el) + $
      nu.s3p_el*(T.s3p-T.el) + $
;      nu.s4p_el*(T.s4p-T.el) + $
      nu.op_el* (T.op- T.el) + $
      nu.o2p_el*(T.o2p-T.el) + $
      nu.el_elh*(T.elh - T.el)

; Radiative losses in ev
rad = ft_rad(lat_dist, T, r_ind)

; In addition to Coulomb collisions and radiative losses, when
; ionization occurs, the electron population will lose at least as
; much energy as the ionization potential for a given ion species. 
ip_o = 13.61806       ; ionization potential in eV
ip_op = 35.11730
ip_o2p = 54.9355
ip_s = 10.36001
ip_sp = 23.3379
ip_s2p = 34.79
ip_s3p = 47.222

; total flux tube integrated amount of energy lost from the thermal
; electron population due to the ionization potential
en_ip = ftint.io * ip_o + ftint.iop * ip_op + ftint.io2p * ip_o2p + ftint.is * ip_s + $
        ftint.isp * ip_sp + ftint.is2p * ip_s2p + ftint.is3p * ip_s3p 

; Total number of electrons in the flux tube
Ne_tot = sqrt(!pi) * (1 * n.sp  * h.sp  + $
                      2 * n.s2p * h.s2p + $
                      3 * n.s3p * h.s3p + $
                      1 * n.op  * h.op  + $
                      2 * n.o2p * h.o2p) / (1.-n.protons)

; Ionization potential lost per electon 
ip_per_e = en_ip / Ne_tot

ip_loss = (2./3.) * ip_per_e * n.el
; The factor of 2/3 is neccecssary since these quantities are in
; energy, but we're really tracking temperature in this code. 
EF_el = Teq - (2./3.) * rad - r_dep.transport * n.el  * T.el-ip_loss

return,EF_el
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function EF_sp,n,T,r_ind,r_dep,h,nu, ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.sp
S = T.pu_s * (ftint.is + ftint.ish + ftint.cx_k1 + ftint.cx_k2 + ftint.cx_k4 + $
              ftint.cx_k9 + ftint.cx_k10) + $
    T.s2p * (ftint.rs2p + ftint.cx_k0 + ftint.cx_k2 + ftint.cx_k12)

L = T.sp * (ftint.isp + ftint.isph + ftint.rsp + ftint.cx_k0 + ftint.cx_k1 + $
            ftint.cx_k8 + ftint.cx_k13 + ftint.cx_k16 + $
            r_dep.transport * n.sp * rootpi_h)

Teq= nu.sp_s2p*(T.s2p-T.sp) + $
     nu.sp_s3p*(T.s3p-T.sp) + $
;     nu.sp_s4p*(T.s4p-T.sp) + $
     nu.sp_op* (T.op -T.sp) + $
     nu.sp_o2p*(T.o2p-T.sp) + $
     nu.sp_el* (T.el- T.sp) + $
     nu.sp_elh*(T.elh-T.sp)
  
EF_sp = (S - L)/rootpi_h + Teq

return,EF_sp
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function EF_s2p,n,T,r_ind,r_dep,h,nu, ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.s2p

S1 =  T.pu_s * (ftint.cx_k3 + ftint.cx_k11) + $
      T.sp * (ftint.isp + ftint.isph + ftint.cx_k0 + ftint.cx_k13 + ftint.cx_k16) + $
      T.s3p * (ftint.rs3p + ftint.cx_k4 + ftint.cx_k14 + ftint.cx_k16)

L1 = T.s2p * (ftint.is2p + ftint.is2ph + ftint.rs2p + ftint.cx_k0 + ftint.cx_k2 + $
              ftint.cx_k3 + ftint.cx_k12 + ftint.cx_k15 + $
              r_dep.transport * n.s2p * rootpi_h)

Teq= nu.sp_s2p* (T.sp- T.s2p) + $
     nu.s2p_s3p*(T.s3p-T.s2p) + $
;     nu.s2p_s4p*(T.s4p-T.s2p) + $
     nu.s2p_op* (T.op- T.s2p) + $
     nu.s2p_o2p*(T.o2p-T.s2p) + $
     nu.s2p_el* (T.el- T.s2p) + $
     nu.s2p_elh*(T.elh-T.s2p)

EF_s2p = (S1 - L1) / rootpi_h + Teq

return,EF_s2p
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function EF_s3p,n,T,r_ind,r_dep,h,nu, ftint
;-----------------------------------------------------------------------  
rootpi_h = sqrt(!dpi) * h.s3p

S1 = T.s2p * (ftint.is2p + ftint.is2ph + ftint.cx_k15)

L1 = T.s3p * (ftint.is3p + ftint.is3ph + ftint.rs3p + ftint.cx_k14 + ftint.cx_k16 + $
              r_dep.transport * n.s3p * rootpi_h)

Teq= nu.sp_s3p* (T.sp- T.s3p) + $
     nu.s2p_s3p*(T.s2p-T.s3p) + $
;     nu.s3p_s4p*(T.s4p-T.s3p) + $
     nu.s3p_op* (T.op- T.s3p) + $
     nu.s3p_o2p*(T.o2p-T.s3p) + $
     nu.s3p_el* (T.el- T.s3p) + $
     nu.s3p_elh*(T.elh-T.s3p)



EF_s3p = (S1 - L1)/rootpi_h + Teq

return,EF_s3p
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
;function EF_s4p,n,T,r_ind,r_dep,h,nu, ftint
;;-----------------------------------------------------------------------
;S = r_dep.is3p*ftavg_el(n.s3p,n,h.s3p,h,h.s3p)*T.s3p + $ ;cold e ionization 
;    r_ind.is3ph*ftavg_elh(n.s3p,n,h.s3p,h,h.s3p)*T.s3p   ;hot e ionization

;L = r_dep.transport*n.s4p*T.s4p

;Teq= nu.sp_s4p* (T.sp- T.s4p) + $
;     nu.s2p_s4p*(T.s2p-T.s4p) + $
;     nu.s3p_s4p*(T.s3p-T.s4p) + $
;     nu.s4p_op* (T.op- T.s4p) + $
;     nu.s4p_o2p*(T.o2p-T.s4p) + $
;     nu.s4p_el* (T.el- T.s4p) + $
;     nu.s4p_elh*(T.elh-T.s4p)

;EF_s4p = S - L + Teq

;return,EF_s4p
;end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function EF_op,n,T,r_ind,r_dep,h,nu, ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.op

S1 = T.pu_o * (ftint.io + ftint.ioh + ftint.cx_k5 + ftint.cx_k6 + ftint.cx_k8 + $
               ftint.cx_k12 + ftint.cx_k14) + $
     T.o2p * (ftint.ro2p + ftint.cx_k6 + ftint.cx_k10 + ftint.cx_k11 + ftint.cx_k13 + $
              ftint.cx_k15)

L1 = T.op * (ftint.iop + ftint.ioph + ftint.rop + ftint.cx_k5 + ftint.cx_k9 + $
             r_dep.transport * n.op * rootpi_h)

Teq= nu.sp_op* (T.sp- T.op) + $
     nu.s2p_op*(T.s2p-T.op) + $
     nu.s3p_op*(T.s3p-T.op) + $
;     nu.s4p_op*(T.s4p-T.op) + $
     nu.op_o2p*(T.o2p-T.op) + $
     nu.op_el* (T.el- T.op) + $
     nu.op_elh*(T.elh-T.op)

EF_op = (S1 - L1)/rootpi_h + Teq

return,EF_op
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function EF_o2p,n,T,r_ind,r_dep,h,nu, ftint
;-----------------------------------------------------------------------
rootpi_h = sqrt(!dpi) * h.o2p

S1 = T.pu_o * ftint.cx_k7 + $
     T.op * (ftint.iop + ftint.ioph)

L1 = T.o2p * (ftint.io2p + ftint.io2ph + ftint.ro2p + ftint.cx_k6 + ftint.cx_k7 + $
              ftint.cx_k10 + ftint.cx_k11 + ftint.cx_k13 + ftint.cx_k15 + $
              r_dep.transport * n.o2p * rootpi_h)

Teq= nu.sp_o2p* (T.sp- T.o2p) + $
     nu.s2p_o2p*(T.s2p-T.o2p) + $
     nu.s3p_o2p*(T.s3p-T.o2p) + $
;     nu.s4p_o2p*(T.s4p-T.o2p) + $
     nu.op_o2p* (T.op- T.o2p) + $
     nu.o2p_el* (T.el- T.o2p) + $
     nu.o2p_elh*(T.elh-T.o2p)

EF_o2p = (S1 - L1) / rootpi_h + Teq

return,EF_o2p
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
pro cm3_model,n,T,nT,h,nu,n1,T1,nT1,h1,nu1,np,Tp,nTp,r_ind,r_dep,r1_dep, lat_dist, $
              lat_dist1, ftint
;main program
;cubic centimeter torus chemistry model
;-----------------------------------------------------------------------
@cm3_model_common

if check_math(mask=211,/noclear) gt 0 then stop

; Update the temperature dependent reaction rates
get_dependent_rates,r_dep,r_ind,n,T,h

; Calculate the flux tube integrated number of reactions for all
; reactions
cm3_reactions, r_ind, r_dep, h, n, ftint

; Calulate the latitudinal distribution of densities
lat_distribution, n, h, lat_dist

; Calculate the nu values for all combinations of species at the 
; initial time.

nu.sp_s2p  = nu_ii(32, 1, lat_dist.sp, T.sp, 32, 2, lat_dist.s2p, T.s2p)
nu.sp_s3p  = nu_ii(32, 1, lat_dist.sp, T.sp, 32, 3, lat_dist.s3p, T.s3p)
;nu.sp_s4p  = nu_ii(32, 1, lat_dist.sp, T.sp, 32, 4, lat_dist.s4p, T.s4p)
nu.sp_op   = nu_ii(32, 1, lat_dist.sp, T.sp, 16, 1, lat_dist.op,  T.op)
nu.sp_o2p  = nu_ii(32, 1, lat_dist.sp, T.sp, 16, 2, lat_dist.o2p, T.o2p)

nu.s2p_s3p = nu_ii(32, 2, lat_dist.s2p, T.s2p, 32, 3, lat_dist.s3p, T.s3p)
;nu.s2p_s4p = nu_ii(32, 2, lat_dist.s2p, T.s2p, 32, 4, lat_dist.s4p, T.s4p)
nu.s2p_op  = nu_ii(32, 2, lat_dist.s2p, T.s2p, 16, 1, lat_dist.op,  T.op)
nu.s2p_o2p = nu_ii(32, 2, lat_dist.s2p, T.s2p, 16, 2, lat_dist.o2p, T.o2p)

;nu.s3p_s4p = nu_ii(32, 3, lat_dist.s3p, T.s3p, 32, 4, lat_dist.s4p, T.s4p)
nu.s3p_op  = nu_ii(32, 3, lat_dist.s3p, T.s3p, 16, 1, lat_dist.op,  T.op)
nu.s3p_o2p = nu_ii(32, 3, lat_dist.s3p, T.s3p, 16, 2, lat_dist.o2p, T.o2p)

;nu.s4p_op  = nu_ii(32, 4, lat_dist.s4p, T.s4p, 16, 1, lat_dist.op,  T.op)
;nu.s4p_o2p = nu_ii(32, 4, lat_dist.s4p, T.s4p, 16, 2, lat_dist.o2p, T.o2p)

nu.op_o2p  = nu_ii(16, 1, lat_dist.op, T.op, 16, 2, lat_dist.o2p, T.o2p)

nu.sp_el   = nu_ie(32, 1, lat_dist.sp,  T.sp,  lat_dist.el, T.el)
nu.s2p_el  = nu_ie(32, 2, lat_dist.s2p, T.s2p, lat_dist.el, T.el)
nu.s3p_el  = nu_ie(32, 3, lat_dist.s3p, T.s3p, lat_dist.el, T.el)
;nu.s4p_el  = nu_ie(32, 4, lat_dist.s4p, T.s4p, lat_dist.el, T.el)
nu.op_el   = nu_ie(16, 1, lat_dist.op,  T.op,  lat_dist.el, T.el)
nu.o2p_el  = nu_ie(16, 2, lat_dist.o2p, T.o2p, lat_dist.el, T.el)

nu.sp_elh  = nu_ieh(32, 1, lat_dist.sp,  T.sp,  n.elh, T.elh)
nu.s2p_elh = nu_ieh(32, 2, lat_dist.s2p, T.s2p, n.elh, T.elh)
nu.s3p_elh = nu_ieh(32, 3, lat_dist.s3p, T.s3p, n.elh, T.elh)
;nu.s4p_elh = nu_ieh(32, 4, lat_dist.s4p, T.s4p, n.elh, T.elh)
nu.op_elh  = nu_ieh(16, 1, lat_dist.op,  T.op,  n.elh, T.elh)
nu.o2p_elh = nu_ieh(16, 2, lat_dist.o2p, T.o2p, n.elh, T.elh)

nu.el_elh  = nu_ee(lat_dist.el, T.el, n.elh, T.elh)

; Calculate the ftavg values for all parent species

Fs   = F_s(n,h,r_ind,r_dep, ftint)
Fsp  = F_sp(n,h,r_ind,r_dep, ftint)
Fs2p = F_s2p(n,h,r_ind,r_dep, ftint)
Fs3p = F_s3p(n,h,r_ind,r_dep, ftint)
;Fs4p = F_s4p(n,h,r_ind,r_dep,ftint)
Fo   = F_o(n,h,r_ind,r_dep, ftint)
Fop  = F_op(n,h,r_ind,r_dep, ftint)
Fo2p = F_o2p(n,h,r_ind,r_dep, ftint)
 
;full time step advance
n1.s   = n.s   + dt * Fs
n1.sp  = n.sp  + dt * Fsp
n1.s2p = n.s2p + dt * Fs2p
n1.s3p = n.s3p + dt * Fs3p
;n1.s4p = n.s4p + dt * Fs4p
n1.o   = n.o   + dt * Fo
n1.op  = n.op  + dt * Fop
n1.o2p = n.o2p + dt * Fo2p

n1.el = (n1.sp + 2.*n1.s2p + 3.*n1.s3p + 4.*n1.s4p + n1.op + 2.*n1.o2p)/(1.-n1.protons)
n1.elh = n1.fh/(1.-n1.fh) * n1.el

;full time step advance for energy

EFsp  = EF_sp(n,T,r_ind,r_dep,h,nu, ftint)
EFs2p = EF_s2p(n,T,r_ind,r_dep,h,nu, ftint)
EFs3p = EF_s3p(n,T,r_ind,r_dep,h,nu, ftint)
;EFs4p = EF_s4p(n,T,r_ind,r_dep,h,nu, ftint)
EFop  = EF_op(n,T,r_ind,r_dep,h,nu, ftint)
EFo2p = EF_o2p(n,T,r_ind,r_dep,h,nu, ftint)
EFel  = EF_el(n,T,h, r_ind,r_dep,lat_dist,nu, ftint)

nT1.sp =  nT.sp  + dt * EFsp
nT1.s2p = nT.s2p + dt * EFs2p
nT1.s3p = nT.s3p + dt * EFs3p
;nT1.s4p = nT.s4p + dt * EFs4p
nT1.op =  nT.op  + dt * EFop
nT1.o2p = nT.o2p + dt * EFo2p
nT1.el =  nT.el  + dt * EFel

update_temp,n1,nT1,T1
get_scale_heights,h1,T1,n1
get_dependent_rates,r1_dep,r_ind,n1,T1,h1
lat_distribution, n1, h1, lat_dist1

if check_math(mask=211,/noclear) gt 0 then stop
; Calculate the nu values for all combinations of species at the
; full time step.
nu1.sp_s2p  = nu_ii(32, 1, n1.sp,  T1.sp, 32, 2, n1.s2p, T1.s2p)
nu1.sp_s3p  = nu_ii(32, 1, n1.sp,  T1.sp, 32, 3, n1.s3p, T1.s3p)
;nu1.sp_s4p  = nu_ii(32, 1, n1.sp,  T1.sp, 32, 4, n1.s4p, T1.s4p)
nu1.sp_op   = nu_ii(32, 1, n1.sp,  T1.sp, 16, 1, n1.op,  T1.op)
nu1.sp_o2p  = nu_ii(32, 1, n1.sp,  T1.sp, 16, 2, n1.o2p, T1.o2p)

nu1.s2p_s3p = nu_ii(32, 2, n1.s2p, T1.s2p, 32, 3, n1.s3p, T1.s3p)
;nu1.s2p_s4p = nu_ii(32, 2, n1.s2p, T1.s2p, 32, 4, n1.s4p, T1.s4p)
nu1.s2p_op  = nu_ii(32, 2, n1.s2p, T1.s2p, 16, 1, n1.op,  T1.op)
nu1.s2p_o2p = nu_ii(32, 2, n1.s2p, T1.s2p, 16, 2, n1.o2p, T1.o2p)

;nu1.s3p_s4p = nu_ii(32, 3, n1.s3p, T1.s3p, 32, 4, n1.s4p, T1.s4p)
nu1.s3p_op  = nu_ii(32, 3, n1.s3p, T1.s3p, 16, 1, n1.op,  T1.op)
nu1.s3p_o2p = nu_ii(32, 3, n1.s3p, T1.s3p, 16, 2, n1.o2p, T1.o2p)

;nu1.s4p_op  = nu_ii(32, 4, n1.s4p, T1.s4p, 16, 1, n1.op,  T1.op)
;nu1.s4p_o2p = nu_ii(32, 4, n1.s4p, T1.s4p, 16, 2, n1.o2p, T1.o2p)

nu1.op_o2p  = nu_ii(16, 1, n1.op,  T1.op, 16, 2, n1.o2p, T1.o2p)

nu1.sp_el  = nu_ie(32, 1, lat_dist1.sp,  T1.sp,  lat_dist1.el, T1.el)
nu1.s2p_el = nu_ie(32, 2, lat_dist1.s2p, T1.s2p, lat_dist1.el, T1.el)
nu1.s3p_el = nu_ie(32, 3, lat_dist1.s3p, T1.s3p, lat_dist1.el, T1.el)
;nu1.s4p_el = nu_ie(32, 4, lat_dist1.s4p, T1.s4p, lat_dist1.el, T1.el)
nu1.op_el  = nu_ie(16, 1, lat_dist1.op,  T1.op,  lat_dist1.el, T1.el)
nu1.o2p_el = nu_ie(16, 2, lat_dist1.o2p, T1.o2p, lat_dist1.el, T1.el)

nu1.sp_elh  = nu_ieh(32, 1, lat_dist1.sp,  T1.sp,  n1.elh, T1.elh)
nu1.s2p_elh = nu_ieh(32, 2, lat_dist1.s2p, T1.s2p, n1.elh, T1.elh)
nu1.s3p_elh = nu_ieh(32, 3, lat_dist1.s3p, T1.s3p, n1.elh, T1.elh)
;nu1.s4p_elh = nu_ieh(32, 4, lat_dist1.s4p, T1.s4p, n1.elh, T1.elh)
nu1.op_elh  = nu_ieh(16, 1, lat_dist1.op,  T1.op,  n1.elh, T1.elh)
nu1.o2p_elh = nu_ieh(16, 2, lat_dist1.o2p, T1.o2p, n1.elh, T1.elh)

nu1.el_elh = nu_ee(lat_dist1.el, T1.el, n1.elh, T1.elh)

;Improved Euler advance
np.s   = n.s   + dt * 0.5 * (Fs   + F_s(n1,h1,r_ind,r1_dep, ftint))
np.sp  = n.sp  + dt * 0.5 * (Fsp  + F_sp(n1,h1,r_ind,r1_dep, ftint))
np.s2p = n.s2p + dt * 0.5 * (Fs2p + F_s2p(n1,h1,r_ind,r1_dep, ftint))
np.s3p = n.s3p + dt * 0.5 * (Fs3p + F_s3p(n1,h1,r_ind,r1_dep, ftint))
;np.s4p = n.s4p + dt * 0.5 * (Fs4p + F_s4p(n1,h1,r_ind,r1_dep, ftint))
np.o   = n.o   + dt * 0.5 * (Fo   + F_o(n1,h1,r_ind,r1_dep, ftint))
np.op  = n.op  + dt * 0.5 * (Fop  + F_op(n1,h1,r_ind,r1_dep, ftint))
np.o2p = n.o2p + dt * 0.5 * (Fo2p + F_o2p(n1,h1,r_ind,r1_dep, ftint))

np.el = (np.sp + 2.*np.s2p + 3.*np.s3p + 4.*np.s4p + np.op + 2.*np.o2p)/(1.-np.protons)
np.elh = np.fh/(1.-np.fh) * np.el

nTp.sp =  nT.sp  + dt * 0.5 * (EFsp  + EF_sp(n1,T1,r_ind,r1_dep,h1,nu1, ftint))
nTp.s2p = nT.s2p + dt * 0.5 * (EFs2p + EF_s2p(n1,T1,r_ind,r1_dep,h1,nu1, ftint))
nTp.s3p = nT.s3p + dt * 0.5 * (EFs3p + EF_s3p(n1,T1,r_ind,r1_dep,h1,nu1, ftint))
;nTp.s4p = nT.s4p + dt * 0.5 * (EFs4p + EF_s4p(n1,T1,r_ind,r1_dep,h1,nu1, ftint))
nTp.op =  nT.op  + dt * 0.5 * (EFop  + EF_op(n1,T1,r_ind,r1_dep,h1,nu1, ftint))
nTp.o2p = nT.o2p + dt * 0.5 * (EFo2p + EF_o2p(n1,T1,r_ind,r1_dep,h1,nu1, ftint))
nTp.el =  nT.el  + dt * 0.5 * (EFel  + EF_el(n1,T1,h1, r_ind,r1_dep,lat_dist1,nu1, ftint))

update_temp,np,nTp,Tp

n = np
nT = nTp
T = Tp
if check_math(mask=211,/noclear) gt 0 then stop

return
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
PRO cm3_phaseplot,time,phase,_extra=extra
for i=0,n_elements(time)-2 do if phase[i+1]-phase[i] gt -270 and $
  phase[i+1]-phase[i] lt 270 then plots,time[i:i+1],phase[i:i+1],_extra=extra

end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
FUNCTION cm3_londist,lon,par
return,par[0]*cos((2*!pi/par[3])*(lon-par[1])*!dtor)+par[2]
end
;-----------------------------------------------------------------------


;------------------------------------------------------------------
pro lag_trans,u,dtheta,dt_lag, v = v, halfv = v_half

@cm3_model_common

; The greater the velocity, the slower the plasma (i.e. greater
; subcorotation).
IF NOT keyword_set(v) THEN v = lag_const+lag_amp*cos((l3-lag_phase)*!dtor)
IF NOT keyword_set(halfv) THEN v_half=lag_const+lag_amp*cos((l3+dtheta/2.-lag_phase)*!dtor)

rj = 7.1492e4 ;km

omega = v/(rdist*rj)  ;radians per second at 6.0 Rj
omega_half = v_half/(rdist*rj)

;theta_rad = l3*!dtor
dtheta_rad = dtheta*!dtor

flux = omega*u

; Upwind scheme (Too dissipative)
;flx = omega*(shift(u,-1)+u)/2.
;u = u - (dt_lag/dtheta_rad)*(flx - shift(flx,1))

; Two step Lax-Wendroff scheme
; Value of u at t=n+1/2 and x=j+1/2
u_half =.5 * (shift(u,-1) + u) - dt_lag/(2. * dtheta_rad) * (shift(flux,-1) - flux)

; The value of u at t=n+1/2 and x=j-1/2 is just the same as u_half,
; but shifted by +1, i.e. shift(u_half,1)

u = u - dt_lag/dtheta_rad * (u_half * omega_half - shift(u_half * omega_half, 1))

end
;------------------------------------------------------------------


;-------------------------------------------------------------------
PRO cm3_latavg_lag3,n_out,t_out,spfitp_out,s2pfitp_out,s3pfitp_out,opfitp_out,$
            lag_const1=lag_const1,lag_phase1=lag_phase1,lag_amp1=lag_amp1,runt1=runt1,$
            fehot_const1=fehot_const1,fehot_amp1=fehot_amp1,fehot_phase1=fehot_phase1,$
            s4fehot_amp1=s4fehot_amp1,s4fehot_phase1=s4fehot_phase1,tm1 = tm1, $
            Tehot=Tehot,no_initial_lonvar=no_initial_lonvar,plot=plot, $
            filename=filename,onebox=onebox,transport=transport,source=source,o_to_s=o_to_s,$
            contfile=contfile,protons=protons,lon=lon,neutral_amp=neutral_amp,$
            neutral_t0=neutral_t0,neutral_width=neutral_width,hote_amp=hote_amp,$
            hote_t0=hote_t0,hote_width=hote_width,trans_type1=trans_type1,$
            trans_exp1=trans_exp1, o2s_spike1 = o2s_spike1, rdist1 = rdist1, $
            domega = domega, n_height = n_height, printiter = printiter, $
            xwinsize = xwinsize, ywinsize = ywinsize

;-------------------------------------------------------------------
runtime=systime(1)
@cm3_model_common
getsysvariables

!p.multi=[0,1,1]
restore,file='/Users/steffl/iotorus/torus_model/latavg_lag/savefiles/avgmixratios.sav'

; Open new window(s) if needed to plot the results
IF keyword_set(plot) THEN BEGIN
   loadct,39
   oldwinid=!d.window


   IF NOT keyword_set(xwinsize) THEN xwinsize = 640
   IF NOT keyword_set(ywinsize) THEN ywinsize = 512

   window, /free, /pixmap, xs = xwinsize, ys = ywinsize
   comp_pixmap_id = !d.window
   window, /free, xs = xwinsize, ys = ywinsize
   comp_window_id = !d.window

   IF keyword_set(onebox) NE 1 THEN BEGIN
      window, /free, /pixmap, xs = xwinsize, ys = ywinsize
      lon_pixmap_id = !d.window
      window, /free, xs = xwinsize, ys = ywinsize
      lon_window_id = !d.window

      window, /free, /pixmap, xs = xwinsize, ys = ywinsize
      ehot_pixmap_id = !d.window
      window, /free, xs = xwinsize, ys = ywinsize
      ehot_window_id = !d.window
   ENDIF
ENDIF

; What is the transport exponent and is the transport rate dependednt
; on the neutral source or the total ion density
if keyword_set(trans_exp1) then trans_exp=trans_exp1 else trans_exp=1.
if keyword_set(trans_type1) then trans_type=1 else trans_type=0

; Model time step in seconds (slightly longer than the time for
; neutrals to cross a 15 degree azimuthal bin)
dt = 2e3     

if n_elements(runt1) EQ 0 THEN runt = 220.*8.64e4 ELSE runt=runt1*8.64e4; Total run time
nit = fix(runt/dt)+1     ; Number of iterations

IF keyword_set(onebox) THEN nbox = 1 ELSE nbox = 24 ; # of longitudinal boxes
   
dtheta = 360./nbox

; Radial distance of model, default is 6 R_J
IF keyword_set(rdist1) THEN rdist = rdist1 ELSE rdist = 6.

; we want the middle of the box, not the beginning
;l3 = (findgen(nbox)+.5)*dtheta
if keyword_set(lon) then l3 = lon ELSE l3 = findgen(nbox)*dtheta + dtheta/2.
lontemp = l3

; Neutral Cloud scale height in km
theta_offset = 6.4 * cos((l3 - 200.) * !dtor) * !dtor
zoff = abs(theta_offset * rdist * 71492.) ; in km, not cm! Always positive
IF NOT keyword_set(n_height) THEN n_height = 0.2 * 71492.

; System IV angular velocity (relative to System III) in degrees/day
IF NOT keyword_set(domega) THEN domega = 12.5

; Primary structures with values as a function of longitude. N.B. the
; tag "el" refers to the thermal electron population only. Similarly,
; "elh" is only for the hot electron population.
n = replicate({dens_lon,sp:0.0d,s2p:0.0d, s3p: 0.0d, s4p: 0.0d, op: 0.0d, $
  o2p: 0.0d, el: 0.0d, elh: 0.0d, s: 0.0d, o: 0.0d, fc:0d, fh:0d, protons:0d}, nbox)
nT = replicate({energy_lon, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p: 0.0d, op: 0.0d, $
                 o2p: 0.0d, el: 0.0d, elh: 0.0d}, nbox)
T = replicate({temp_lon,sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p: 0.0d,  $
                op: 0.0d, o2p: 0.0d, el: 0.0d, elh: 0.0d, pu_s: 0d, pu_o: 0d}, nbox)
h = replicate({scale_heights_lon, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
                o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}, nbox)

; Structures used to store intermediate values in the Euler advance subroutine
n1=n
nT1=nT
T1=T
h1=h

np=n
nTp=nT
Tp=T

; Structure containing volume coefficients for energy transfer via coulomb collisons
nu = replicate({nu, sp_s2p:0d, sp_s3p:0d, sp_s4p:0d, sp_op:0d, sp_o2p:0d, s2p_s3p:0d, $
      s2p_s4p:0d, s2p_op:0d, s2p_o2p:0d, s3p_s4p: 0d, s3p_op:0d, s3p_o2p:0d, $
      s4p_op: 0d, s4p_o2p:0d, op_o2p:0d, sp_el:0d, s2p_el:0d, s3p_el:0d, s4p_el:0d,$
      op_el:0d, o2p_el:0d, sp_elh:0d, s2p_elh:0d, s3p_elh:0d, s4p_elh:0d,$
      op_elh:0d, o2p_elh:0d, el_elh: 0d},nbox)

nu1 = nu

; The independent rates stored in this structure are constant with
; respect to electron temperature, however they may vary with time.
r_ind = {independent_rates,ish:0.0,isph:0.0,is2ph:0.0,is3ph:0.0,ioh:0.0,ioph:0.0,$
     io2ph:0.0, S_production: 0.0, O_production: 0.0, o2s_spike: 0.0, $
     o_to_s: 0.0, cx_k0: 0.0, cx_k1: 0.0, cx_k2: 0.0, cx_k3: 0.0, cx_k4: 0.0,$
     cx_k5: 0.0, cx_k6: 0.0, cx_k7: 0.0, cx_k8: 0.0, cx_k9: 0.0, $
     cx_k10: 0.0, cx_k11: 0.0, cx_k12: 0.0, cx_k13: 0.0, cx_k14: 0.0,$
     cx_k15: 0.0, cx_k16: 0.0, emisSp: fltarr(101,101), emisS2p:fltarr(101,101),$
     emisS3p: fltarr(101,101), emisOp: fltarr(101,101), emisO2p: fltarr(101,101), $
     emistemp: fltarr(101), emisden: fltarr(101)}

; Dependent rates are dependent on (hot) electron temperature and/or
; (hot) electron density and therefore will vary as a function of time
; and longitude
r_dep=replicate({dependent_rates, is: 0.0, isp: 0.0, is2p: 0.0, is3p: 0.0, $
  io: 0.0, iop:0.0, io2p: 0.0, rsp: 0.0, rs2p: 0.0, rs3p: 0.0, rop: 0.0, ro2p:0.0,$
  Transport:0.0}, nbox)
r1_dep=r_dep

; Latitudinal Density structure
lat_dist = replicate({lat_dist, z:dindgen(21)/5. * 71492., sp:dblarr(21), s2p:dblarr(21), $
             s3p:dblarr(21), s4p:dblarr(21), op:dblarr(21), o2p:dblarr(21), $
             el:dblarr(21)}, nbox)
lat_dist1 = lat_dist

; flux tube integrated number of reactions per second.
ftint = replicate({ftint, cx_k0:0d, cx_k1:0d, cx_k2:0d, cx_k3:0d, cx_k4:0d, cx_k5:0d, $
                  cx_k6:0d, cx_k7:0d, cx_k8:0d, cx_k9:0d, cx_k10:0d, cx_k11:0d, cx_k12:0d, $
                  cx_k13:0d, cx_k14:0d, cx_k15:0d, cx_k16:0d, is:0d, isp:0d, is2p:0d, $
                  is3p:0d, io:0d, iop:0d, io2p:0d, ish:0d, isph:0d, is2ph:0d, is3ph:0d, $
                  ioh:0d, ioph:0d, io2ph:0d, rsp:0d, rs2p:0d, rs3p:0d, rop:0d, ro2p:0d}, $
                  nbox)

; Output quantites, that will contain the values of density and
; temperature in every longitudinal box at every time step
n_out=n
t_out=t
h_out=h

spfitp_out=fltarr(4)
s2pfitp_out=fltarr(4)
s3pfitp_out=fltarr(4)
opfitp_out=fltarr(4)

spfitperr_out=fltarr(4)
s2pfitperr_out=fltarr(4)
s3pfitperr_out=fltarr(4)
opfitperr_out=fltarr(4)

; Initialize time, phase, and amplitude variables
IF n_elements(tm1) EQ 0 THEN tm = 170. ELSE tm = tm1
tm0=tm

; Initialize the parinfo structure, which is used in fitting the data
; with a sine wave.
parinfo=replicate({value:0d,limited:[0,0],limits:[0.d,0.d],parname:'',fixed:0},4)
parinfo.value=[.05,180.,.15,2*!pi]
parinfo[1].limited=[1,1]
parinfo[1].limits=[-1080.,1080.]
parinfo[0].parname='Amplitude'
parinfo[1].parname='Phase'
parinfo[2].parname='Constant Offset'
parinfo[3].parname='Period'
parinfo[3].fixed=1

; Is this a continuation of a previous run?
if keyword_set(contfile) then begin

    restore,file=contfile

    n_out_dim=size(n_out,/dim)

    n=n_out[*,n_out_dim[1]-1]
    T=T_out[*,n_out_dim[1]-1]
    h=h_out[*,n_out_dim[1]-1]

    IF n_elements(n) EQ 1 AND NOT keyword_set(onebox) THEN BEGIN
       n = cmreplicate(n, nbox)
       T = cmreplicate(T, nbox)
       h = cmreplicate(h, nbox)
       l3 = lontemp
    ENDIF

    n_out=n
    t_out=t
    h_out=h
    if size(trans_out,/type) eq 0 then trans_out=trans else $
      trans_out=trans_out[n_out_dim[1]-1]
    o2sout = otos

    source_out=source_out[n_out_dim[1]-1]
    net_source = source_out

    IF NOT keyword_set(onebox) THEN BEGIN 
       IF n_elements(spfitp_out) GT 4 THEN BEGIN
          spfitp_out=spfitp_out[*,n_out_dim[1]-1]
          s2pfitp_out=s2pfitp_out[*,n_out_dim[1]-1]
          s3pfitp_out=s3pfitp_out[*,n_out_dim[1]-1]
          opfitp_out=opfitp_out[*,n_out_dim[1]-1]

          spfitperr_out=spfitperr_out[*,n_out_dim[1]-1]
          s2pfitperr_out=s2pfitperr_out[*,n_out_dim[1]-1]
          s3pfitperr_out=s3pfitperr_out[*,n_out_dim[1]-1]
          opfitperr_out=opfitperr_out[*,n_out_dim[1]-1]
       ENDIF

       psp=spfitp_out[1]
       ps2p=s2pfitp_out[1]
       ps3p=s3pfitp_out[1]
       pop=opfitp_out[1]

       ampsp=spfitp_out[0]/spfitp_out[2]*100.
       amps2p=s2pfitp_out[0]/s2pfitp_out[2]*100.
       amps3p=s3pfitp_out[0]/s3pfitp_out[2]*100.
       ampop=opfitp_out[0]/opfitp_out[2]*100.
    ENDIF

; Start At DOY 170, which is June 18, in a leap year (as 2000 was)
    IF n_elements(tm1) EQ 0 THEN tm = 170. ELSE tm = tm1

endif else begin
    if keyword_set(no_initial_lonvar) then begin

; Peter's default density values
        n.sp=250.
        n.s2p=400.
        n.s3p=50.
;        n.s4p=1e-5
        n.op=700.
        n.o2p=50.
        n.s=10.d
        n.o=50.
    endif else begin
        n.sp=  cm3_londist(l3,[0.00997,190.,0.0939,2*!dpi])*1800.
        n.s2p= cm3_londist(l3,[0.00476,201.,0.2270,2*!dpi])*1800.
        n.s3p= cm3_londist(l3,[0.00351,355.,0.0223,2*!dpi])*1800.
        n.op=  cm3_londist(l3,[0.00793,40.2,0.2290,2*!dpi])*1800.
        n.o2p=0.123*n.op 
;       n.o2p=0.063*n.op ;Assume initial o++ is at a level of 6.3% o+ (Steffl et al.2004b)
        n.s=50.
        n.o=100.
    endelse

; January Latavg Conditions
    Te0 = 5.0d
    Ti0 = 70.0d
    Teh0 = 49.0d
    fehot_const = 0.0022d
    trans = 1.0d/(70.28*8.64e4)
    net_source = 19.95d27
    otos = 1.7d

; Total electron density
    n.el = (n.sp + 2 * n.s2p + 3 * n.s3p + 4 * n.s4p + n.op + 2 * n.o2p) * (1.-n.protons)
    n.elh = n.fh/(1.-n.fh) * n.el               ; Hot electron density
    n.fc = 1.-n.fh

; Assign Temp values
    T[*].sp = Ti0
    T[*].s2p = Ti0
    T[*].s3p = Ti0
;    T[*].s4p = Ti0
    T[*].op = Ti0
    T[*].o2p = Ti0
    T[*].el = Te0
    T[*].elh = Teh0

; Determine initial scale heights
    get_scale_heights,h,t,n

; Fit the initial conditions with a cosine function
    ftmix=ftint_mix(n,h)
    spfitp=mpfitfun('fitlonfunc',l3,ftmix.sp,replicate(1.,n_elements(n.sp)),$
                    parinfo=parinfo,yfit=s2yfit,perror=spfitperr,/quiet)
    s2pfitp=mpfitfun('fitlonfunc',l3,ftmix.s2p,replicate(1.,n_elements(n.s2p)),$
                     parinfo=parinfo,yfit=s3yfit,perror=s2pfitperr,/quiet)
    s3pfitp=mpfitfun('fitlonfunc',l3,ftmix.s3p,replicate(1.,n_elements(n.s3p)),$
                     parinfo=parinfo,yfit=s3pfitp,perror=s3pfitperr,/quiet)
    opfitp=mpfitfun('fitlonfunc',l3,ftmix.op,replicate(1.,n_elements(n.op)),$
                    parinfo=parinfo,yfit=opyfit,perror=opfitperr,/quiet)

    spfitp[1]  =  (spfitp[1]+1080.) mod 360.
    s2pfitp[1] = (s2pfitp[1]+1080.) mod 360.
    s3pfitp[1] = (s3pfitp[1]+1080.) mod 360.
    opfitp[1]  =  (opfitp[1]+1080.) mod 360.
    
    psp =spfitp[1]
    ps2p=s2pfitp[1]
    ps3p=s3pfitp[1]
    pop =opfitp[1]
    
    ampsp= spfitp[0]/spfitp[2]*100.
    amps2p=s2pfitp[0]/s2pfitp[2]*100.
    amps3p=s3pfitp[0]/s3pfitp[2]*100.
    ampop= opfitp[0]/opfitp[2]*100.  
    
    spfitp_out= spfitp
    s2pfitp_out=s2pfitp
    s3pfitp_out=s3pfitp
    opfitp_out= opfitp
    
    spfitperr_out= spfitperr
    s2pfitperr_out=s2pfitperr
    s3pfitperr_out=s3pfitperr
    opfitperr_out= opfitperr

    source_out=net_source
    n_out=n
    t_out=t
    h_out=h
    trans_out=avg(r_dep.transport)
 ENDELSE

; Allow for protons to be included in the model. Right now, any added
; protons do not take part in the chemistry i.e. no p + O -> H + O+
; type of reactions. Therefore, their only real effect is to increase
; the electron denisty while keeping the ion densities
; constant. Default is not to include protons
IF keyword_set(protons) THEN BEGIN
    n.protons  = protons
    n1.protons = protons
    np.protons = protons
 ENDIF

; Initialize System III varying corotational lag and hot electron source 
if n_elements(lag_const1) gt 0 then lag_const=lag_const1 else lag_const=1.2 ; km/s
if n_elements(lag_phase1) gt 0 then lag_phase=lag_phase1 else lag_phase=20. 
if n_elements(lag_amp1) gt 0 then lag_amp=lag_amp1 else lag_amp=.2    ; km/s

if n_elements(fehot_const1) gt 0 then fehot_const=fehot_const1
if n_elements(fehot_amp1) gt 0 then fehot_amp=fehot_amp1 else fehot_amp=  0.2
if n_elements(fehot_phase1) gt 0 then fehot_phase=fehot_phase1 else fehot_phase=20.
if n_elements(s4fehot_amp1) gt 0 then s4fehot_amp=s4fehot_amp1 else s4fehot_amp=  0.2
if n_elements(s4fehot_phase1) gt 0 then s4fehot_phase=s4fehot_phase1 else s4fehot_phase=20.

if keyword_set(Tehot) then Teh0=Tehot
if keyword_set(transport) then trans = 1./(transport*8.64e4)
if keyword_set(source) then net_source = source
if keyword_set(o_to_s) then otos= o_to_s
IF keyword_set(o2s_spike1) THEN o2s_spike = o2s_spike1 ELSE o2s_spike = otos
r_ind.o_to_s = otos
r_ind.o2s_spike = o2s_spike

; Initialize time variable neutral source and hot electron fraction parameters
tau0=1./trans/8.64e4
net_source0 = net_source
fh0 = fehot_const

if not keyword_set(neutral_amp) then neutral_amp=0.
if not keyword_set(neutral_t0) then neutral_t0=249. ; 5-Sept-2000 as per Delamere et al 2004
if not keyword_set(neutral_width) then neutral_width=25.0

if not keyword_set(hote_amp) then hote_amp=0.0000
if not keyword_set(hote_t0) then hote_t0= 269.
if not keyword_set(hote_width) then hote_width=45.5
IF keyword_set(n_height) THEN BEGIN 
   h.s = n_height
   h.o = n_height
ENDIF

; Read in total emitted power arrays.
read_emission_tables,'home/dcoffin/1D_Model/Eddiecode/iotorus/torus_model/tepSII_v7.1.5',temp,den,emisSp
read_emission_tables,'home/dcoffin/1D_Model/Eddiecode/iotorus/torus_model/tepSIII_v7.1.5',temp,den,emisS2p
read_emission_tables,'home/dcoffin/1D_Model/Eddiecode/iotorus/torus_model/tepSIV_v7.1.5',temp,den,emisS3p
read_emission_tables,'home/dcoffin/1D_Model/Eddiecode/iotorus/torus_model/tepOII_v7.1.5',temp,den,emisOp
read_emission_tables,'home/dcoffin/1D_Model/Eddiecode/iotorus/torus_model/tepOIII_v7.1.5',temp,den,emisO2p

; Store info into r_ind structure. These fields should not get modified by
; any other subroutine
r_ind.emistemp = temp
r_ind.emisden  = den
r_ind.emisSp   = emisSp
r_ind.emisS2p  = emisS2p
r_ind.emisS3p  = emisS3p
r_ind.emisOp   = emisOp
r_ind.emisO2p  = emisO2p

;set initial pickup temperatures
Lshl = rdist
T.pu_s = Tpu(32.,Lshl)
T.pu_o = Tpu(16.,Lshl)
T.elh = Teh0

T1.pu_s = T.pu_s
T1.pu_o = T.pu_o
T1.elh = Teh0

Tp.pu_s = T.pu_s
Tp.pu_o = T.pu_o
Tp.elh = Teh0

; Call get_independent_rates to get the reaction rates for charge
; exchange, hot electron ionization, and neutral source production. As
; the name implies, these rates do not change over the run.
get_independent_rates,r_ind,T,h

net_source0 = net_source
iondens0 = avg(n.sp + n.s2p + n.s3p + n.s4p + n.op + n.o2p)

; Angular velocity (rad/sec) for circular orbit about Jupiter at rdist
omega_neutrals = 1./sqrt((rdist * !rjkm * 1e3)^3/(6.673e-11 * 1.8987e27))
v_neutrals = omega_neutrals * rdist * !rjkm ; km/s
v_corotation = 2. * !pi * rdist * !rjkm / (9.925 * 3600.) ; km/s
rel_v_neutrals = v_corotation - v_neutrals

; Final consistency checks
; Check that the sum of the thermal & hot electron fractions equals one
n.fc = 1. - n.fh

;  Ensure charge neutrality in the thermal electron population. We
;  assume that the hot electrons are linked to a separate hot proton population.
n.el = (n.sp + 2.*n.s2p + 3.*n.s3p + 4.*n.s4p + n.op + 2.*n.o2p)/(1.-n.protons)
n.elh = n.fh/(1.-n.fh) * n.el

nT.el  = n.el  * T.el
nT.elh  = n.elh * T.elh
nT.sp  = n.sp  * T.sp
nT.s2p = n.s2p * T.s2p
nT.s3p = n.s3p * T.s3p
;nT.s4p = n.s4p * T.s4p
nT.op  = n.op  * T.op
nT.o2p = n.o2p * T.o2p

n1=n
nT1=nT
T1=T
h1=h

np=n
nTp=nT
Tp=T

get_scale_heights,h,t,n

IF keyword_set(printiter) THEN format="('% CM3_LATAVG_LAG3: Iteration ', i5, ' of ', i5)"
badtime = 0l
;=================================================================
; Loop over # of interations
;=================================================================
for j=0,nit-1 do begin 

   IF keyword_set(printiter) AND (j MOD 100) EQ 0 THEN print, $
     strtrim(j, 2), strtrim(nit, 2), format = format

   if check_math(mask=211,/noclear) gt 1 then begin
      print,'Math error detected! Exiting loop.'
      break
   endif

   time = tm0 + j * dt / 86400.

; Time variable neutral source
   net_source = net_source0*(1. + neutral_amp*exp(-(time-neutral_t0)^2/neutral_width^2))

; Let the neutral source be the sum of two independent sources: a
; baseline source with an O/S of o_to_s and a Gaussian increase with a
; potentially different O/S of o2s_spike
   r_ind.o_to_s = (o_to_s + $
                   o2s_spike * neutral_amp*exp(-(time-neutral_t0)^2/neutral_width^2))/$
                  (1+neutral_amp*exp(-(time-neutral_t0)^2/neutral_width^2))

; Time variable hot electron fraction
   n.fh = fehot_const * (1. + hote_amp * exp(-(time - hote_t0)^2 / hote_width^2)) * $
          (1. + fehot_amp * cos((l3 - fehot_phase) * !dtor) + $
           s4fehot_amp * cos(((l3 - s4fehot_phase - time * domega) mod 360) * !dtor))

; This was an attempt at revision
;   n.fh = fehot_const * (1. + hote_amp * exp(-(time - hote_t0)^2 / hote_width^2) + $
;            fehot_amp * cos((l3 - fehot_phase) * !dtor) + $
;            s4fehot_amp * cos(((l3 - s4fehot_phase - time * domega) mod 360) * !dtor))
   
   negative_fh = where(n.fh LT 0, nbadfh)
   IF nbadfh GT 0 THEN BEGIN
      badtime = [badtime, j]
      n.fh >= 0.
   ENDIF

   n1.fh = n.fh
   np.fh = n.fh

   n.fc  = 1d0-n.fh
   n1.fc = 1d0-n.fh
   np.fc = 1d0-n.fh

; Update the hot electron densities
   n.elh = n.fh/(1.-n.fh) * n.el
   nT.elh  = n.elh * T.elh

; Call cm3_model to advance the model one time step
   cm3_model, n,  T,  nT,  h,  nu,  $
              n1, T1, nT1, h1, nu1, $
              np, Tp, nTp, $
              r_ind, r_dep, r1_dep, lat_dist, lat_dist1, ftint

;-------------------------------------------------------------
; Azimuthal transport loop. While a time step of 1e4 seconds works
; well for the chemistry timescales, the angular velocity of the
; neutrals at 6 Rj is such that they will cross several azimuthal bins
; in a single timestep. This is problematic, since the neutral density
; in an azimuthal bin at the start of the time step will be
; significantly different from the neutral density at the end of the
; timestep, purely due to the azimuthal transport. The best solution,
; I believe is to decrease the model time step to something
; approaching the time for a neutral to cross 1 azimuthal bin. For 15
; degree bins, this is 1,928 seconds, so round up to 2e3. This change
; is reflected in the new value of dt above.
   IF NOT keyword_set(onebox) THEN BEGIN 
      ns = n.s
      lag_trans, ns, dtheta, dt, v = rel_v_neutrals, halfv = rel_v_neutrals
      n.s = ns
    
      no = n.o
      lag_trans, no, dtheta, dt, v = rel_v_neutrals, halfv = rel_v_neutrals
      n.o = no

      IF lag_const GT 0 THEN BEGIN
         nsp = n.sp
         lag_trans,nsp,dtheta, dt
         n.sp = nsp
         nTsp = nT.sp
         lag_trans,nTsp,dtheta, dt
         nT.sp = nTsp
         
         ns2p = n.s2p
         lag_trans,ns2p,dtheta, dt
         n.s2p = ns2p
         nTs2p = nT.s2p
         lag_trans,nTs2p,dtheta, dt
         nT.s2p = nTs2p
         
         ns3p = n.s3p
         lag_trans,ns3p,dtheta, dt
         n.s3p = ns3p
         nTs3p = nT.s3p
         lag_trans,nTs3p,dtheta, dt
         nT.s3p = nTs3p
         
;         ns4p = n.s4p
;         lag_trans,ns4p,dtheta, dt
;         n.s4p = ns4p
;         nTs4p = nT.s4p
;         lag_trans,nTs4p,dtheta, dt
;         nT.s4p = nts4p
       
         nop = n.op
         lag_trans,nop,dtheta, dt
         n.op = nop
         nTop = nT.op
         lag_trans,nTop,dtheta, dt
         nT.op = nTop
         
         no2p = n.o2p
         lag_trans,no2p,dtheta, dt
         n.o2p = no2p
         nTo2p = nT.o2p
         lag_trans,nTo2p,dtheta, dt
         nT.o2p = nTo2p
      
         nel = n.el
         lag_trans,nel,dtheta, dt
         n.el = nel
         nTel = nT.el
         lag_trans,nTel,dtheta, dt
         nT.el = nTel
      ENDIF
   ENDIF
; End Azimuthal Transport Loop    
;-------------------------------------------------------------

; Update the temperature and scale height structures
   update_temp,n,nT,T
   get_scale_heights,h,t,n

; Fitting things every 2000 seconds is overkill, so only do this on
; every 5th timestep (i.e. every 1e4 s)
   IF j MOD 5 EQ 0 THEN BEGIN 
      tm = [tm,time]
      
      ftmix=ftint_mix(n,h)
      n_out=[[n_out],[n]]
      t_out=[[t_out],[t]]
      h_out=[[h_out],[h]]
      trans_out=[trans_out,avg(r_dep.transport)]
      o2sout = [o2sout, r_ind.o_to_s]
      source_out=[source_out,net_source]

      IF keyword_set(plot) THEN BEGIN 
; Plot longitudunally averaged mixing ratios vs. time

; Open a new window if needed then plot the results
         device,window_state=openwindows
         IF openwindows[comp_pixmap_id] EQ 0 THEN BEGIN
            window, /free, /pixmap, xs = xwinsize, ys = ywinsize
            comp_pixmap_id = !d.window
         ENDIF ELSE wset, comp_pixmap_id

         !p.multi = 0
         IF j EQ 0 THEN ftmix_out = [ftmix, ftmix] ELSE ftmix_out=[ftmix_out, ftmix]
         plot,tm,(NOT keyword_set(onebox)) ? avg(ftmix_out.sp,0) : ftmix_out.sp,/xsty,$
              /nodata,yr=[.01,.4],/ysty,xtit='Doy 2000',ytit='Avg. FT int. mix ratio',$
              /ylog,xr=[min(tm),390]
         oplot,tm,(NOT keyword_set(onebox)) ? avg(ftmix_out.sp,0)  : ftmix_out.sp, $
               col=57,thick=3
         oplot,tm,(NOT keyword_set(onebox)) ? avg(ftmix_out.s2p,0) : ftmix_out.s2p, $
               col=147,thick=3
         oplot,tm,(NOT keyword_set(onebox)) ? avg(ftmix_out.s3p,0) : ftmix_out.s3p, $
               col=210,thick=3
         oplot,tm,(NOT keyword_set(onebox)) ? avg(ftmix_out.op,0)  : ftmix_out.op, $
               col=254,thick=3
         oplot,tm,(NOT keyword_set(onebox)) ? avg(ftmix_out.o2p,0) : ftmix_out.o2p, $
               col=100,thick=3
         
         oplot,avgtime,avgs2mix,linestyle=1
         oplot,avgtime,avgs3mix,linestyle=1
         oplot,avgtime,avgs4mix,linestyle=1
         oplot,avgtime,avgo2mix,linestyle=1
         
         plots,avgtime[0],avgs2mix[0],psym=1,symsize=1.2
         plots,avgtime[0],avgs3mix[0],psym=2,symsize=1.2
         plots,avgtime[0],avgs4mix[0],psym=4,symsize=1.2
         plots,avgtime[0],avgo2mix[0],psym=5,symsize=1.2
;        plots,avgtime[0],avgo2mix[0]*.123,psym=6,symsize=1.2

         plots,avgtime[41],avgs2mix[41],psym=1,symsize=1.2
         plots,avgtime[41],avgs3mix[41],psym=2,symsize=1.2
         plots,avgtime[41],avgs4mix[41],psym=4,symsize=1.2
         plots,avgtime[41],avgo2mix[41],psym=5,symsize=1.2
;        plots,avgtime[41],avgo2mix[41]*.123,psym=6,symsize=1.2

         jan14=366.+14
         plots,jan14,0.060,psym=1,symsize=1.2
         plots,jan14,0.212,psym=2,symsize=1.2
         plots,jan14,0.034,psym=4,symsize=1.2
         plots,jan14,0.242,psym=5,symsize=1.2
         plots,jan14,0.030,psym=6,symsize=1.2

         IF openwindows[comp_window_id] EQ 0 THEN BEGIN
            window, /free, xs = xwinsize, ys = ywinsize
            comp_window_id = !d.window
         ENDIF ELSE wset, comp_window_id

         device,copy=[0,0,xwinsize,ywinsize,0,0,comp_pixmap_id]

      ENDIF

      IF n_elements(onebox) EQ 0 THEN BEGIN
         
; Fit the mixing ratios with a sine wave. Extract the amplitude,
; phase, and constant offset
         spfitp=mpfitfun('fitlonfunc',l3,ftmix.sp,replicate(1.,n_elements(n.sp)),$
                         parinfo=parinfo,yfit=s2yfit,perror=spfitperr,/quiet)
         s2pfitp=mpfitfun('fitlonfunc',l3,ftmix.s2p,replicate(1.,n_elements(n.s2p)),$
                          parinfo=parinfo,yfit=s3yfit,perror=s2pfitperr,/quiet)
         s3pfitp=mpfitfun('fitlonfunc',l3,ftmix.s3p,replicate(1.,n_elements(n.s3p)),$
                          parinfo=parinfo,yfit=s3pfitp,perror=s3pfitperr,/quiet)
         opfitp=mpfitfun('fitlonfunc',l3,ftmix.op,replicate(1.,n_elements(n.op)),$
                         parinfo=parinfo,yfit=opyfit,perror=opfitperr,/quiet)
         
         spfitp[1]  =  (spfitp[1]+1080.) mod 360.
         s2pfitp[1] = (s2pfitp[1]+1080.) mod 360.
         s3pfitp[1] = (s3pfitp[1]+1080.) mod 360.
         opfitp[1]  =  (opfitp[1]+1080.) mod 360.
         
         psp = [psp,spfitp[1]]
         ps2p = [ps2p,s2pfitp[1]]
         ps3p = [ps3p,s3pfitp[1]]
         pop = [pop,opfitp[1]]
         
         ampsp= [ampsp,spfitp[0]/spfitp[2]*100.]
         amps2p= [amps2p,s2pfitp[0]/s2pfitp[2]*100.]
         amps3p= [amps3p,s3pfitp[0]/s3pfitp[2]*100.]
         ampop= [ampop,opfitp[0]/opfitp[2]*100.]  
         
         spfitp_out=[[spfitp_out],[spfitp]]
         s2pfitp_out=[[s2pfitp_out],[s2pfitp]]
         s3pfitp_out=[[s3pfitp_out],[s3pfitp]]
         opfitp_out=[[opfitp_out],[opfitp]]
         
         spfitperr_out=[[spfitperr_out],[spfitperr]]
         s2pfitperr_out=[[s2pfitperr_out],[s2pfitperr]]
         s3pfitperr_out=[[s3pfitperr_out],[s3pfitperr]]
         opfitperr_out=[[opfitperr_out],[opfitperr]]
         
; Open a new window if needed then plot the results
         IF keyword_set(plot) THEN BEGIN 
            device,window_state=openwindows
            IF openwindows[lon_pixmap_id] EQ 0 THEN BEGIN
               window, /free, /pixmap, xs = xwinsize, ys = ywinsize
               lon_pixmap_id = !d.window
            ENDIF ELSE wset, lon_pixmap_id
            
            !p.multi=[0,1,2]
            plot,tm,[0,360],yrange=[0,360],/ysty,/nodata,xtit='Day of Year 2000',$
                 ytit='System III longitude of peak',ytickint=60,yminor=6,/xsty
            cm3_phaseplot,tm,psp,col=57,thick=3
            cm3_phaseplot,tm,ps2p,col=147,thick=3
            cm3_phaseplot,tm,ps3p,col=210,thick=3
            cm3_phaseplot,tm,pop,col=254,thick=3
            
; Plot the amplitude (As a percentage of the contant offset)
            plot,tm,[0,100],yrange=[0,40],/ysty,/nodata,xtit='Day of Year 2000',$
                 ytit='Amplitude of variation (%)',ytickint=10,/xsty
            plots,tm,ampsp,col=57,thick=3
            plots,tm,amps2p,col=147,thick=3
            plots,tm,amps3p,col=210,thick=3
            plots,tm,ampop,col=254,thick=3
            
            IF openwindows[lon_window_id] EQ 0 THEN BEGIN
               window, /free, xs = xwinsize, ys = ywinsize
               lon_window_id = !d.window
            ENDIF
            device,copy=[0,0,xwinsize,ywinsize,0,0,lon_window_id]

; Plot the ion mixing ratios and ion temperatures vs. longitude    
            device,window_state=openwindows
            IF openwindows[ehot_pixmap_id] EQ 0 THEN BEGIN
               window, /free, /pixmap, xs = xwinsize, ys = ywinsize
               ehot_pixmap_id = !d.window
            ENDIF ELSE wset, ehot_pixmap_id

            !p.multi=[0,1,1]
            plot,l3,fehotvariation,/ynozero,xr=[0,360],/xsty,xtit='System III longitude',$
                 ytit='Fraction of hot electrons',yr=[.5,1.5],psym=-1

            IF openwindows[lon_window_id] EQ 0 THEN BEGIN
               window, /free, xs = xwinsize, ys = ywinsize
               lon_window_id = !d.window
            ENDIF
            device,copy=[0,0,xwinsize,ywinsize,0,0,lon_window_id]
         ENDIF
      ENDIF
   ENDIF
ENDFOR

IF n_elements(badtime) GT 1 THEN BEGIN
   badtime = badtime[1:*]
   print, 'CM3_LATAVG_LAG3: n.fh went negative!'
   stop
ENDIF

save,file='cm3_lag'+filename+'.sav',n_out,t_out,h_out,spfitp_out,s2pfitp_out,s3pfitp_out,$
     opfitp_out,spfitperr_out,s2pfitperr_out,s3pfitperr_out,opfitperr_out,tm,Teh0,$
     source_out,otos,lag_const,lag_amp,lag_phase,fehot_const,fehot_amp,fehot_phase,$
     trans_out,l3

IF n_elements(onebox) EQ 0 THEN BEGIN   
   set_plot,'ps'
   device,/land,/color,bits=8,/helvetica,file='cm3_lag'+filename+'.ps'
   loadct,39
   !p.multi=[0,1,3]
   plot,tm,[0,360],yrange=[0,360],/ysty,/nodata,xtit='Day of year 2000',$
        ytit='System III longitude of peak',ytickint=60,yminor=6,/xsty
   cm3_phaseplot,tm,psp,col=57,thick=3
   cm3_phaseplot,tm,ps2p,col=147,thick=3
   cm3_phaseplot,tm,ps3p,col=210,thick=3
   cm3_phaseplot,tm,pop,col=254,thick=3
   plots,[avgtime[0],avgtime[0]],[!y.crange[0],!y.crange[1]]
   plots,[avgtime[41],avgtime[41]],[!y.crange[0],!y.crange[1]]
   
; Plot the amplitude (As a percentage of the contant offset)
   plot,tm,[0,100],yrange=[0,40],/ysty,/nodata,xtit='Day of year 2000',$
        ytit='Amplitude of variation (%)',ytickint=10,/xsty
   plots,tm,ampsp,col=57,thick=3
   plots,tm,amps2p,col=147,thick=3
   plots,tm,amps3p,col=210,thick=3
   plots,tm,ampop,col=254,thick=3
   plots,[avgtime[0],avgtime[0]],[!y.crange[0],!y.crange[1]]
   plots,[avgtime[41],avgtime[41]],[!y.crange[0],!y.crange[1]]

; Plot composition vs longitude
   if not n_elements(onebox) eq 1 then begin
      plot,l3,ftmix.op,psym=4,yrange=[0.01,0.4],/ysty,xrange=[0,360],/xsty,/nodata,xminor=4,$
           xtickint=60,/ylog,xtit='System III longitude',ytit='FT-integrated Ion mixing ratio'
      oplot,l3,ftmix.sp,psym=-1,col=57
      oplot,l3,ftmix.s2p,psym=-2,col=147
      oplot,l3,ftmix.s3p,psym=-5,col=210
      oplot,l3,ftmix.op,psym=-4,col=254
      oplot,l3,ftmix.o2p,psym=-3,col=100
      oplot,l3,n.el/1e5   
   endif


; Plot composition vs time
   ftmix_out=ftint_mix(n_out,h_out)
   plot,tm,avg(ftmix_out.sp,0),/xsty,/nodata,yr=[.01,.4],/ysty,xtit='Day of year 2000',$
        ytit='Avg. FT int. mix ratio',/ylog
   oplot,tm,avg(ftmix_out.sp,0),col=57,thick=3
   oplot,tm,avg(ftmix_out.s2p,0),col=147,thick=3
   oplot,tm,avg(ftmix_out.s3p,0),col=210,thick=3
   oplot,tm,avg(ftmix_out.op,0),col=254,thick=3
   oplot,tm,avg(ftmix_out.o2p,0),col=100,thick=3
   oplot,avgtime,avgs2mix,linestyle=1
   oplot,avgtime,avgs3mix,linestyle=1
   oplot,avgtime,avgs4mix,linestyle=1
   oplot,avgtime,avgo2mix,linestyle=1
   
   plots,avgtime[0],avgs2mix[0],psym=1,symsize=1.2
   plots,avgtime[0],avgs3mix[0],psym=2,symsize=1.2
   plots,avgtime[0],avgs4mix[0],psym=4,symsize=1.2
   plots,avgtime[0],avgo2mix[0],psym=5,symsize=1.2
   plots,avgtime[0],avgo2mix[0]*.123,psym=6,symsize=1.2
   
   plots,avgtime[41],avgs2mix[41],psym=1,symsize=1.2
   plots,avgtime[41],avgs3mix[41],psym=2,symsize=1.2
   plots,avgtime[41],avgs4mix[41],psym=4,symsize=1.2
   plots,avgtime[41],avgo2mix[41],psym=5,symsize=1.2
   plots,avgtime[41],avgo2mix[41]*.123,psym=6,symsize=1.2
   
   plots,[avgtime[0],avgtime[0]],[!y.crange[0],!y.crange[1]]
   plots,[avgtime[41],avgtime[41]],[!y.crange[0],!y.crange[1]]
   
   jan14=366.+14
   plots,jan14,0.060,psym=1,symsize=1.2
   plots,jan14,0.212,psym=2,symsize=1.2
   plots,jan14,0.034,psym=4,symsize=1.2
   plots,jan14,0.242,psym=5,symsize=1.2
   plots,jan14,0.030,psym=6,symsize=1.2
   
   sigma=textoidl('\sigma')
   xyouts,.0,.99,/normal,'A!DN!N = '+strsigfig(neutral_amp,4),charsize=.8
   xyouts,.0,.97,/normal,'t!D0!LN!N = '+strsigfig(neutral_t0,4),charsize=.8
   xyouts,.0,.95,/normal,sigma+'!DN!N = '+strsigfig(neutral_width,4),charsize=.8
   
   xyouts,.15,.99,/normal,'A!Dhot e!N = '+strsigfig(hote_amp,4),charsize=.8
   xyouts,.15,.97,/normal,'t!D0!Lhot e!N = '+strsigfig(hote_t0,4),charsize=.8
   xyouts,.15,.95,/normal,sigma+'!Dhot e!N = '+strsigfig(hote_width,4),charsize=.8
   
   xyouts,.3,.99,/normal,'T!De hot!N = '+strsigfig(teh0,4),charsize=.8
   xyouts,.3,.97,/normal,'O / S = '+strsigfig(otos,4),charsize=.8
   xyouts,.3,.95,/normal,'source = '+strsigfig(net_source,4),charsize=.8
   
   xyouts,.45,.99,/normal,'lag amp = '+strsigfig(lag_amp,4),charsize=.8
   xyouts,.45,.97,/normal,'lag const = '+strsigfig(lag_const,4),charsize=.8
   xyouts,.45,.95,/normal,'lag phase = '+strsigfig(lag_phase,4),charsize=.8
   
   xyouts,.6,.99,/normal,'f!De hot!N amp = '+strsigfig(fehot_amp,4),charsize=.8
   xyouts,.6,.97,/normal,'f!De hot!N const = '+strsigfig(fehot_const,4),charsize=.8
   xyouts,.6,.95,/normal,'f!De hot!N phase = '+strsigfig(fehot_phase,4),charsize=.8
   
   xyouts,.75,.99,/normal,'Sys4 f!De hot!N amp = '+strsigfig(s4fehot_amp,4),charsize=.8
   xyouts,.75,.97,/normal,'Sys4 f!De hot!N phase = '+strsigfig(s4fehot_phase,4),charsize=.8
   
   !p.multi = [0, 2, 2, 0, 1]

   restore,'~/iotorus/torus_model/latavg_lag/savefiles/fitparvslon_o3.sav'
   
; Phase plots
   plot,hourfit/24.,dawns2lonfitp[1,*],psym=1,xr=[275,320],/xsty,$
        yr=[0,360],/ysty,/nodata,chars=1.,xtickint=5., yminor = 6,$
        ytickint=60,ymargin=[0,4],xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
   oplot,hourfit/24.,dawns2lonfitp[1,*],psym=1,col=57
   
   for i=0,(size(spfitp_out,/dim))[1]-2 do $
         if abs(spfitp_out[1,i+1]-spfitp_out[1,i]) lt 100 then $
         oplot,tm[i:i+1],spfitp_out[1,i:i+1],thick=3 ;,col=100
   xyouts,277,300,charsize=1.5,'S II'
   
   al_legend,/bottom,/right,['UVIS Data','Model'],psym=[1,0],number=2

   plot,hourfit/24.,dawns4lonfitp[1,*],psym=1,xr=[275,320],/xsty,xtit='Day of Year 2000',$
        yr=[0,360],/ysty,/nodata,chars=1.,xtickint=5.,ytickint=60,ymargin=[4,0], yminor = 6
   oplot,hourfit/24.,dawns4lonfitp[1,*],psym=1,col=210
   for i=0,(size(s3pfitp_out,/dim))[1]-2 do $
         if abs(s3pfitp_out[1,i+1]-s3pfitp_out[1,i]) lt 100 then $
         oplot,tm[i:i+1],s3pfitp_out[1,i:i+1],thick=3 ;,col=100
   xyouts,.06,.18,/normal,orient=90,'Longitude of Peak Mixing ratio',charsize=1.5
   xyouts,277,300,charsize=1.5,'S IV'
   
   
; Amplitude plots
   plot,hourfit/24.,dawns2lonfitp[0,*]/dawns2lonfitp[2,*],psym=1,xr=[275,320],/xsty,$
        yr=[0,.3],/ysty,/nodata,chars=1.,xtickint=5.,yminor = 5, $
        ymargin=[0,4],xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
   oplot,hourfit/24.,dawns2lonfitp[0,*]/dawns2lonfitp[2,*],psym=1,col=57
   oplot,tm,spfitp_out[0,*]/spfitp_out[2,*],thick=3 ;,col=100
   xyouts,315,.25,charsize=1.5,'S II'
   al_legend,pos=[276,.29],['UVIS Data','Model'],psym=[1,0],number=2
   
   plot,hourfit/24.,dawns4lonfitp[0,*]/dawns4lonfitp[2,*],psym=1,xr=[275,320],/xsty,$
        xtit='Day of Year 2000',$
        yr=[0,.3],/ysty,/nodata,chars=1.,xtickint=5.,ymargin=[4,0], yminor = 5
   oplot,hourfit/24.,dawns4lonfitp[0,*]/dawns4lonfitp[2,*],psym=1,col=210
   oplot,tm,s3pfitp_out[0,*]/s3pfitp_out[2,*],thick=3
   xyouts,315,.25,charsize=1.5,'S IV'
   
   xyouts,.06,.13,/normal,orient=90,'Amplitude of Azimuthal Variation',charsize=1.5
   
   device,/close
ENDIF

IF keyword_set(plot) THEN BEGIN
   wdelete, comp_pixmap_id
   IF n_elements(onebox) EQ 0 THEN BEGIN
      wdelete, lon_pixmap_id
      wdelete, ehot_pixmap_id
   ENDIF
ENDIF

getsysvariables,/restore

print,'Runtime of: '+strtrim(systime(1)-runtime,2)+' seconds'

end
;-------------------------------------------------------------------

