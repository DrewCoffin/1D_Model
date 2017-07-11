;+
;  IO TORUS CUBIC CENTIMETER ONEBOX MODEL, SUBROUTINES
;     
;  Main routine: cm3_latavg_onebox.pro
;-

;-----------------------------------------------------------------------
FUNCTION ion_ion, h1, h2

h_prime = h1 * h2 / sqrt(h1*h1 + h2*h2)

return, h_prime
END
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------

FUNCTION ion_neutral, hi, hn, z0

hi2 = hi*hi
hn2 = hn*hn

a = (hi2 + hn2)/(hi2 * hn2)
b = -2 * z0/hn2
c = z0*z0/hn2

return, sqrt(1./a) * exp((b*b - 4.*a*c)/(4.*a))
END
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
PRO cm3_reactions, r_ind, r_dep, h, n, ftint
common info,runt,dt,zoff,rdist
;rootpi = !rootpi

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
ftint.cx_k0  = !rootpi * r_ind.cx_k0  * n.sp  * n.s2p * sp_s2p
ftint.cx_k1  = !rootpi * r_ind.cx_k1  * n.s   * n.sp  * s_sp
ftint.cx_k2  = !rootpi * r_ind.cx_k2  * n.s   * n.s2p * s_s2p
ftint.cx_k3  = !rootpi * r_ind.cx_k3  * n.s   * n.s2p * s_s2p
ftint.cx_k4  = !rootpi * r_ind.cx_k4  * n.s   * n.s3p * s_s3p
ftint.cx_k5  = !rootpi * r_ind.cx_k5  * n.o   * n.op  * o_op
ftint.cx_k6  = !rootpi * r_ind.cx_k6  * n.o   * n.o2p * o_o2p
ftint.cx_k7  = !rootpi * r_ind.cx_k7  * n.o   * n.o2p * o_o2p
ftint.cx_k8  = !rootpi * r_ind.cx_k8  * n.o   * n.sp  * o_sp
ftint.cx_k9  = !rootpi * r_ind.cx_k9  * n.s   * n.op  * s_op
ftint.cx_k10 = !rootpi * r_ind.cx_k10 * n.s   * n.o2p * s_o2p
ftint.cx_k11 = !rootpi * r_ind.cx_k11 * n.s   * n.o2p * s_o2p
ftint.cx_k12 = !rootpi * r_ind.cx_k12 * n.o   * n.s2p * o_s2p
ftint.cx_k13 = !rootpi * r_ind.cx_k13 * n.o2p * n.sp  * sp_o2p
ftint.cx_k14 = !rootpi * r_ind.cx_k14 * n.o   * n.s3p * o_s3p
ftint.cx_k15 = !rootpi * r_ind.cx_k15 * n.o2p * n.s2p * s2p_o2p
ftint.cx_k16 = !rootpi * r_ind.cx_k16 * n.s3p * n.sp  * sp_s3p

; Electron Impact Ionization (Thermal)
ftint.is   = !rootpi * r_dep.is   * n.s   * e_s
ftint.isp  = !rootpi * r_dep.isp  * n.sp  * e_sp
ftint.is2p = !rootpi * r_dep.is2p * n.s2p * e_s2p
ftint.is3p = !rootpi * r_dep.is3p * n.s3p * e_s3p

ftint.io   = !rootpi * r_dep.io   * n.o   * e_o
ftint.iop  = !rootpi * r_dep.iop  * n.op  * e_op
ftint.io2p = !rootpi * r_dep.io2p * n.o2p * e_o2p

; Electron Impact Ionization (Hot Population)
ftint.ish   = !rootpi * r_ind.ish   * n.elh * n.s   * h.s
ftint.isph  = !rootpi * r_ind.isph  * n.elh * n.sp  * h.sp
ftint.is2ph = !rootpi * r_ind.is2ph * n.elh * n.s2p * h.s2p
ftint.is3ph = !rootpi * r_ind.is3ph * n.elh * n.s3p * h.s3p

ftint.ioh   = !rootpi * r_ind.ioh   * n.elh * n.o   * h.o
ftint.ioph  = !rootpi * r_ind.ioph  * n.elh * n.op  * h.op
ftint.io2ph = !rootpi * r_ind.io2ph * n.elh * n.o2p * h.o2p

; Recombination (thermal only)
ftint.rsp = !rootpi  * r_dep.rsp  * n.sp  * e_sp
ftint.rs2p = !rootpi * r_dep.rs2p * n.s2p * e_s2p
ftint.rs3p = !rootpi * r_dep.rs3p * n.s3p * e_s3p

ftint.rop = !rootpi  * r_dep.rop  * n.op  * e_op
ftint.ro2p = !rootpi * r_dep.ro2p * n.o2p * e_o2p

END
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
pro get_scale_heights,h,T,n
;all h in km
;
;
print,h,T,n
;-----------------------------------------------------------------------
if check_math(mask=211,/noclear) gt 0 then stop
@cm3_model_common

; The plasma isn't rigidly corotating anymore. It seems that instead
; of the corotation angular velocity, we should use the actual angular
; velocity as a function of System III longitude.
v_corot = 2. * !pi * rdist * !rjkm / (9.925 * 3600.)
v = v_corot - lag_const + lag_amp * cos((l3-lag_phase)*!dtor)
omega = v/(rdist*!rjkm) 

mp = 1.67262e-27
ms = 32.*mp 
mo = 16.*mp

; Scale height from equation on p. 11,049 of Bagenal 1994
h.sp  = sqrt(2.0*T.sp *1.6e-19*(1+1.*T.el/T.sp)/(3.0*ms*omega^2))/1e3
h.s2p = sqrt(2.0*T.s2p*1.6e-19*(1+2.*T.el/T.s2p)/(3.0*ms*omega^2))/1e3
h.s3p = sqrt(2.0*T.s3p*1.6e-19*(1+3.*T.el/T.s3p)/(3.0*ms*omega^2))/1e3
;h.s4p = sqrt(2.0*T.s4p*1.6e-19*(1+4.*T.el/T.s4p)/(3.0*ms*omega^2))/1e3
h.op  = sqrt(2.0*T.op *1.6e-19*(1+1.*T.el/T.op)/(3.0*mo*omega^2))/1e3
print, "T.op, T.el = ", T.op, T.el
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
lat_dist.el = lat_dist.sp + 2. * lat_dist.s2p + 3. * lat_dist.s3p + 4. * lat_dist.s4p + $
              lat_dist.op + 2. * lat_dist.o2p

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
r_dep.is   = cfit(16,16, T.el)
r_dep.isp  = cfit(16,15, T.el)
r_dep.is2p = cfit(16,14, T.el)
r_dep.is3p = cfit(16,13, T.el)
r_dep.io   = cfit( 8, 8, T.el)
r_dep.iop  = cfit( 8, 7, T.el)
r_dep.io2p = cfit( 8, 6, T.el)

; Total electron recombination rates from the work of Sultana Nahar
;     S II, S III, O III: (ApJS, 101, 423 [1995]; ApJS 106, 213 [1996])
;     O I, O II:          (ApJS, 120, 131 [1999])
; Recombination rates for S I and S IV come from Mazzotta et al 1998, which
; provides formulae for the dielectronic recombination rate. The
; radiative recombination rate comes from Dima Verner's rrfit.f code,
; which uses: Shull & Van Steenberg, (1982, ApJS, 48, 95) for Sulfur 
;             Pequignot et al. (1991, A&A, 251, 680) for Oxygen
recombination_rates, T.el, rsp, rs2p, rs3p, rop, ro2p
r_dep.rsp  = rsp 
r_dep.rs2p = rs2p
r_dep.rs3p = rs3p
r_dep.rop  = rop
r_dep.ro2p = ro2p

fs = 1.0/(1.0 + r_ind.o_to_s)
fo = r_ind.o_to_s/(1.0 + r_ind.o_to_s)

; If transport_type=0, then the transport rate is inversely
; proprtional to the neutral source to the trans_exp power. If the
; trans_type keyword is set, then the transport rate will be inversely
; proportional to the ion density to the trans_exp power.
if keyword_set(trans_type) then BEGIN
   iondens = n.sp + n.s2p + n.s3p + n.s4p + n.op + n.o2p
   tau = tau0 * (iondens0/iondens)^trans_exp 
ENDIF ELSE tau = tau0 * (net_source0/net_source)^trans_exp
r_dep.transport = 1/(tau*8.64e4)

; The neutral source rate at the Equator is given by:
;  n(0) = N_tot / (sqrt(pi) x H)
; i.e. if the neutrals are more spread out along the field line, less
; will be produced at the equator, and more at higher latitudes. 
; Since scale heights are in km, and densities in cm need to multiply
; by 1e5.
r_ind.S_production = fs*net_source/(!rootpi*(h[0].s*1e5))
r_ind.O_production = fo*net_source/(!rootpi*(h[0].o*1e5))
end

;-----------------------------------------------------------------------
pro get_independent_rates,r_ind,T,h
;-----------------------------------------------------------------------
@cm3_model_common
; The hot electron temperature is fixed, so this should be done
; outside of any loops.
r_ind.ish   = cfit(16,16,T[0].elh,c)
r_ind.isph  = cfit(16,15,T[0].elh,c)
r_ind.is2ph = cfit(16,14,T[0].elh,c)
r_ind.is3ph = cfit(16,13,T[0].elh,c)
r_ind.ioh   = cfit(8,8,T[0].elh,c)
r_ind.ioph  = cfit(8,7,T[0].elh,c)
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
rootpi_h = !rootpi * h.s 

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
rootpi_h = !rootpi * h.sp

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
rootpi_h = !rootpi * h.s2p

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
rootpi_h = !rootpi * h.s3p

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
rootpi_h = !rootpi * h.o

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
rootpi_h = !rootpi * h.op

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
rootpi_h = !rootpi * h.o2p

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

printf,21, '------------Inputs for interpolation------------------'
; The electron temperature along a field line is constant for a
; Maxwellian distribution. A non-thermal distribution would require
; modifying this line to reflect the changing temperature with latitude.
tempindex=unity#(100.*(1.+alog10(T.el))/log5000)
densindex=100.*alog10(nel)/log5000

; radiation rates
emisSp  = interpolate(r_ind.emisSp,tempindex,densindex)
printf,21,FORMAT = ' ("r_ind.emisSp: ",F0)', r_ind.emisSp(0,0)*1e15
printf,21,FORMAT = ' ("T.el: ",F0)', T.el
printf,21,FORMAT = ' ("lat%elec: ",F0)', lat_dist.el(0)
printf,21,FORMAT = ' ("lat%sp: ",F0)', lat_dist.sp(0)
printf,21, " --------Outputs of interpolation---------------"
emisS2p = interpolate(r_ind.emisS2p,tempindex,densindex)
emisS3p = interpolate(r_ind.emisS3p,tempindex,densindex)
emisOp  = interpolate(r_ind.emisOp,tempindex,densindex)
emisO2p = interpolate(r_ind.emisO2p,tempindex,densindex)

; total power per cc
rad_sp  = emisSp*nSp
rad_s2p = emisS2p*nS2p
rad_s3p = emisS3p*nS3p
rad_Op  = emisOp*nOp
rad_O2p = emisO2p*nO2p

; print, "rad_sp = ", rad_sp
; total power emitted by cc at equator
rad_tot = rad_sp + rad_s2p + rad_s3p + rad_op + rad_o2p

Lavg = total(rad_tot*nel,1)/total(nel,1)

printf,21,FORMAT = '("rad_sp: ",F0)', rad_sp(0)
printf,21,FORMAT = '("rad_s2p: ",F0)', rad_s2p(0)
printf,21,FORMAT = '("rad_s3p: ",F0)', rad_s3p(0)
printf,21,FORMAT = '("rad_op: ",F0)', rad_op(0)
printf,21,FORMAT = '("rad_o2p: ",F0)', rad_o2p(0)
printf,21,FORMAT = '("ftrad: ",F0)', total(rad_tot*nel,1)
printf,21,FORMAT = '("elec_tot: ",F0)', total(nel,1)

return, Lavg
END
;-----------------------------------------------------------------------



;-----------------------------------------------------------------------  
FUNCTION ft_rad_tot,lat_dist,T,r_ind
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

; total power per cc
rad_sp  = emisSp*nSp
rad_s2p = emisS2p*nS2p
rad_s3p = emisS3p*nS3p
rad_Op  = emisOp*nOp
rad_O2p = emisO2p*nO2p

; total power emitted by cc at equator
rad_tot=rad_sp+rad_s2p+rad_s3p+rad_op+rad_o2p

rad_tot = total(rad_tot,1)

Lavg = total(rad_tot*nel,1)/total(nel,1)

;return, Lavg
return, rad_tot
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
Ne_tot = !rootpi * (1 * n.sp  * h.sp  + $
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

printf,21,FORMAT = '("----------------EF_elec components: --------------------------- ",F0)', EF_el
printf,21,FORMAT = '("Teq: ",F0)', Teq
printf,21,FORMAT = '("rad: ",F0)', rad
printf,21,FORMAT = '("r_dep.transport: ",F11)', r_dep.transport*1e7
printf,21,FORMAT = '("n.el: ",F0)', n.el
printf,21,FORMAT = '("T.el: ",F0)', T.el
printf,21,FORMAT = '("ip_loss: ",F0)', ip_loss

return,EF_el
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function EF_sp,n,T,r_ind,r_dep,h,nu, ftint
;-----------------------------------------------------------------------
rootpi_h = !rootpi * h.sp
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
rootpi_h = !rootpi * h.s2p

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
rootpi_h = !rootpi * h.s3p

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
rootpi_h = !rootpi * h.op

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
rootpi_h = !rootpi * h.o2p

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


;+
;  PROCEDURE
;    printEF
;
;  PURPOSE
;    Take Andrew Steffl's energy code and breakdown the sources
;    and losses of individual processes
;
;  created by Adam Shinn
;-
pro printEF, n, T, h, r_ind, r_dep, lat_dist, nu, ftint
;-----------------------------------------------------------------------
;function EF_el,n, T, h, r_ind, r_dep, lat_dist, nu, ftint
;-----------------------------------------------------------------------
  Teq = nu.sp_el *(T.sp  - T.el) + $
        nu.s2p_el*(T.s2p - T.el) + $
        nu.s3p_el*(T.s3p - T.el) + $
        ;nu.s4p_el*(T.s4p-T.el) + $
        nu.op_el * (T.op - T.el) + $
        nu.o2p_el*(T.o2p - T.el) + $
        nu.el_elh*(T.elh - T.el)

  ; Radiative losses in ev
  rad = ft_rad(lat_dist, T, r_ind)

  ; In addition to Coulomb collisions and radiative losses, when
  ; ionization occurs, the electron population will lose at least as
  ; much energy as the ionization potential for a given ion species. 
  ip_o = 13.61806               ; ionization potential in eV
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
  Ne_tot = !rootpi * (1 * n.sp  * h.sp  + $
                      2 * n.s2p * h.s2p + $
                      3 * n.s3p * h.s3p + $
                      1 * n.op  * h.op  + $
                      2 * n.o2p * h.o2p) / (1.-n.protons)

  ; Ionization potential lost per electon 
  ip_per_e = en_ip / Ne_tot

  ip_loss = (2./3.) * ip_per_e * n.el
  ; The factor of 2/3 is necessary since these quantities are in
  ; energy, but we're really tracking temperature in this code. 
  EF_el = Teq - (2./3.) * rad - r_dep.transport * n.el  * T.el-ip_loss

  ;return, EF_el
  ;end
  ;-----------------------------------------------------------------------

  ;-----------------------------------------------------------------------
  ;function EF_sp,n,T,r_ind,r_dep,h,nu, ftint
  ;-----------------------------------------------------------------------
  rootpi_h = !rootpi * h.sp
  S_sp = T.pu_s * (ftint.is + ftint.ish + ftint.cx_k1 + ftint.cx_k2 + ftint.cx_k4 + $
                   ftint.cx_k9 + ftint.cx_k10) + $
         T.s2p * (ftint.rs2p + ftint.cx_k0 + ftint.cx_k2 + ftint.cx_k12)
  
  L_sp = T.sp * (ftint.isp + ftint.isph + ftint.rsp + ftint.cx_k0 + ftint.cx_k1 + $
                 ftint.cx_k8 + ftint.cx_k13 + ftint.cx_k16 + $
                 r_dep.transport * n.sp * rootpi_h)

  Teq_sp = nu.sp_s2p*(T.s2p-T.sp) + $
           nu.sp_s3p*(T.s3p-T.sp) + $
           ;nu.sp_s4p*(T.s4p-T.sp) + $
           nu.sp_op* (T.op -T.sp) + $
           nu.sp_o2p*(T.o2p-T.sp) + $
           nu.sp_el* (T.el- T.sp) + $
           nu.sp_elh*(T.elh-T.sp)
  
  EF_sp = (S_sp - L_sp)/rootpi_h + Teq_sp

  src_sp = ( S_sp)/rootpi_h + Teq_sp
  lss_sp = (-L_sp)/rootpi_h + Teq_sp

  frtpi = (3./2.)/rootpi_h
  print, '%% SII'
  print, '%% ionized SI by e-.............', frtpi*T.pu_s*ftint.is
  print, '%% ionized SI by hot e-.........', frtpi*T.pu_s*ftint.ish
  print, '%% recombination of SIII........', frtpi*T.s2p*ftint.rs2p
  print, '%% k0...........................', frtpi*T.s2p*ftint.cx_k0
  print, '%% k1...........................', frtpi*T.pu_s*ftint.cx_k1
  print, '%% k2...........................', frtpi*(T.pu_s*ftint.cx_k2 + T.s2p*ftint.cx_k2)
  print, '%% k4...........................', frtpi*T.pu_s*ftint.cx_k4
  print, '%% k9...........................', frtpi*T.pu_s*ftint.cx_k9
  print, '%% k10..........................', frtpi*T.pu_s*ftint.cx_k10
  print, '%% k12..........................', frtpi*T.s2p*ftint.cx_k12
  print, '%% total source.................', frtpi*S_sp + (3./2.)*Teq_sp
  print, '%% .............................'
  print, '%% ionized SII by e-............', frtpi*T.sp*ftint.isp
  print, '%% ionized SII by hot e-........', frtpi*T.sp*ftint.isph
  print, '%% recombination of SII.........', frtpi*T.sp*ftint.rsp
  print, '%% k0...........................', frtpi*T.sp*ftint.cx_k0
  print, '%% k1...........................', frtpi*T.sp*ftint.cx_k1
  print, '%% k8...........................', frtpi*T.sp*ftint.cx_k8
  print, '%% k13..........................', frtpi*T.sp*ftint.cx_k13
  print, '%% k16..........................', frtpi*T.sp*ftint.cx_k16
  print, '%% thermal equil SIII...........', (3./2.)*nu.sp_s2p*(T.s2p-T.sp)
  print, '%% thermal equil SIV............', (3./2.)*nu.sp_s3p*(T.s3p-T.sp)
  print, '%% thermal equil OII............', (3./2.)*nu.sp_op* (T.op -T.sp)
  print, '%% thermal equil OIII...........', (3./2.)*nu.sp_o2p*(T.o2p-T.sp)
  print, '%% thermal equil e-.............', (3./2.)*nu.sp_el* (T.el- T.sp)
  print, '%% thermal equil hot e-.........', (3./2.)*nu.sp_elh*(T.elh-T.sp)
  print, '%% radial transport.............', frtpi*T.sp*r_dep.transport*n.sp*rootpi_h
  print, '%% total loss...................', frtpi*L_sp + (-3./2.)*Teq_sp

  ;return,EF_sp
  ;end
  ;-----------------------------------------------------------------------

  ;-----------------------------------------------------------------------
  ;function EF_s2p,n,T,r_ind,r_dep,h,nu, ftint
  ;-----------------------------------------------------------------------
  rootpi_h = !rootpi * h.s2p

  S_s2p =  T.pu_s * (ftint.cx_k3 + ftint.cx_k11) + $
           T.sp * (ftint.isp + ftint.isph + ftint.cx_k0 + ftint.cx_k13 + ftint.cx_k16) + $
           T.s3p * (ftint.rs3p + ftint.cx_k4 + ftint.cx_k14 + ftint.cx_k16)

  L_s2p = T.s2p * (ftint.is2p + ftint.is2ph + ftint.rs2p + ftint.cx_k0 + ftint.cx_k2 + $
                   ftint.cx_k3 + ftint.cx_k12 + ftint.cx_k15 + $
                   r_dep.transport * n.s2p * rootpi_h)

  Teq_s2p = nu.sp_s2p* (T.sp- T.s2p) + $
            nu.s2p_s3p*(T.s3p-T.s2p) + $
            ;nu.s2p_s4p*(T.s4p-T.s2p) + $
            nu.s2p_op* (T.op- T.s2p) + $
            nu.s2p_o2p*(T.o2p-T.s2p) + $
            nu.s2p_el* (T.el- T.s2p) + $
            nu.s2p_elh*(T.elh-T.s2p)

  EF_s2p = (S_s2p - L_s2p) / rootpi_h + Teq_s2p

  src_s2p = ( S_s2p) / rootpi_h + Teq_s2p
  lss_s2p = (-L_s2p) / rootpi_h + Teq_s2p

  ;return,EF_s2p
  ;end
  ;-----------------------------------------------------------------------


  ;-----------------------------------------------------------------------
  ;function EF_s3p,n,T,r_ind,r_dep,h,nu, ftint
  ;-----------------------------------------------------------------------  
  rootpi_h = !rootpi * h.s3p

  S_s3p = T.s2p * (ftint.is2p + ftint.is2ph + ftint.cx_k15)

  L_s3p = T.s3p * (ftint.is3p + ftint.is3ph + ftint.rs3p + ftint.cx_k14 + ftint.cx_k16 + $
                   r_dep.transport * n.s3p * rootpi_h)

  Teq_s3p = nu.sp_s3p* (T.sp- T.s3p) + $
            nu.s2p_s3p*(T.s2p-T.s3p) + $
            ;nu.s3p_s4p*(T.s4p-T.s3p) + $
            nu.s3p_op* (T.op- T.s3p) + $
            nu.s3p_o2p*(T.o2p-T.s3p) + $
            nu.s3p_el* (T.el- T.s3p) + $
            nu.s3p_elh*(T.elh-T.s3p)

  EF_s3p = (S_s3p - L_s3p)/rootpi_h + Teq_s3p

  src_s3p = ( S_s3p)/rootpi_h + Teq_s3p
  lss_s3p = (-L_s3p)/rootpi_h + Teq_s3p

  ;return,EF_s3p
  ;end
  ;-----------------------------------------------------------------------


  ;-----------------------------------------------------------------------
  ;function EF_s4p,n,T,r_ind,r_dep,h,nu, ftint
  ;-----------------------------------------------------------------------
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
  ;function EF_op,n,T,r_ind,r_dep,h,nu, ftint
  ;-----------------------------------------------------------------------
  rootpi_h = !rootpi * h.op

  S_op = T.pu_o * (ftint.io + ftint.ioh + ftint.cx_k5 + ftint.cx_k6 + ftint.cx_k8 + $
                   ftint.cx_k12 + ftint.cx_k14) + $
         T.o2p * (ftint.ro2p + ftint.cx_k6 + ftint.cx_k10 + ftint.cx_k11 + ftint.cx_k13 + $
                  ftint.cx_k15)

  L_op = T.op * (ftint.iop + ftint.ioph + ftint.rop + ftint.cx_k5 + ftint.cx_k9 + $
                 r_dep.transport * n.op * rootpi_h)

  Teq_op = nu.sp_op* (T.sp- T.op) + $
           nu.s2p_op*(T.s2p-T.op) + $
           nu.s3p_op*(T.s3p-T.op) + $
           ;nu.s4p_op*(T.s4p-T.op) + $
           nu.op_o2p*(T.o2p-T.op) + $
           nu.op_el* (T.el- T.op) + $
           nu.op_elh*(T.elh-T.op)

  EF_op = (S_op - L_op)/rootpi_h + Teq_op

  ;return,EF_op
  ;end
  ;-----------------------------------------------------------------------

  ;-----------------------------------------------------------------------
  ;function EF_o2p,n,T,r_ind,r_dep,h,nu, ftint
  ;-----------------------------------------------------------------------
  rootpi_h = !rootpi * h.o2p

  S_o2p = T.pu_o * ftint.cx_k7 + $
          T.op * (ftint.iop + ftint.ioph)
  
  L_o2p = T.o2p * (ftint.io2p + ftint.io2ph + ftint.ro2p + ftint.cx_k6 + ftint.cx_k7 + $
                   ftint.cx_k10 + ftint.cx_k11 + ftint.cx_k13 + ftint.cx_k15 + $
                   r_dep.transport * n.o2p * rootpi_h)

  Teq_o2p = nu.sp_o2p* (T.sp- T.o2p) + $
            nu.s2p_o2p*(T.s2p-T.o2p) + $
            nu.s3p_o2p*(T.s3p-T.o2p) + $
            ;nu.s4p_o2p*(T.s4p-T.o2p) + $
            nu.op_o2p* (T.op- T.o2p) + $
            nu.o2p_el* (T.el- T.o2p) + $
            nu.o2p_elh*(T.elh-T.o2p)

  EF_o2p = (S_o2p - L_o2p) / rootpi_h + Teq_o2p

  ;return,EF_o2p
  ;end
  ;-----------------------------------------------------------------------

  f = 3./2.

  print
  print, '%% ENERGY INFORMATION'
  print, '%% sp  source...........', f*src_sp
  print, '%% sp  loss.............', f*lss_sp
  print, '%% s2p source...........', f*src_s2p
  print, '%% s2p loss.............', f*lss_s2p
  print, '%% s3p source...........', f*src_s3p
  print, '%% s3p loss.............', f*lss_s3p
  print, '%% op..................', f*EF_op
  print, '%% o2p.................', f*EF_o2p
  print

end




;-----------------------------------------------------------------------
pro cm3_latavg_model,n,T,nT,h,nu,n1,T1,nT1,h1,nu1,np,Tp,nTp,r_ind,r_dep,r1_dep, lat_dist, $
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
 
;full time step advance, ensure no values go negative
n1.s   = (n.s   + dt * Fs)   > 0.
n1.sp  = (n.sp  + dt * Fsp)  > 0.
n1.s2p = (n.s2p + dt * Fs2p) > 0.
n1.s3p = (n.s3p + dt * Fs3p) > 0.
;n1.s4p = (n.s4p + dt * Fs4p) > 0.
n1.o   = (n.o   + dt * Fo)   > 0.
n1.op  = (n.op  + dt * Fop)  > 0.
n1.o2p = (n.o2p + dt * Fo2p) > 0.

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

nT1.sp  = (nT.sp  + dt * EFsp)  > 0.
nT1.s2p = (nT.s2p + dt * EFs2p) > 0.
nT1.s3p = (nT.s3p + dt * EFs3p) > 0.
;nT1.s4p = (nT.s4p + dt * EFs4p) > 0.
nT1.op  = (nT.op  + dt * EFop)  > 0.
nT1.o2p = (nT.o2p + dt * EFo2p) > 0.
nT1.el  = (nT.el  + dt * EFel)  > 0.
;print, "EFel = ", EFel

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
np.s   = (n.s   + dt * 0.5 * (Fs   + F_s(n1,h1,r_ind,r1_dep, ftint)))   > 0.
np.sp  = (n.sp  + dt * 0.5 * (Fsp  + F_sp(n1,h1,r_ind,r1_dep, ftint)))  > 0.
np.s2p = (n.s2p + dt * 0.5 * (Fs2p + F_s2p(n1,h1,r_ind,r1_dep, ftint))) > 0.
np.s3p = (n.s3p + dt * 0.5 * (Fs3p + F_s3p(n1,h1,r_ind,r1_dep, ftint))) > 0.
;np.s4p = (n.s4p + dt * 0.5 * (Fs4p + F_s4p(n1,h1,r_ind,r1_dep, ftint))) > 0.
np.o   = (n.o   + dt * 0.5 * (Fo   + F_o(n1,h1,r_ind,r1_dep, ftint)))   > 0.
np.op  = (n.op  + dt * 0.5 * (Fop  + F_op(n1,h1,r_ind,r1_dep, ftint)))  > 0.
np.o2p = (n.o2p + dt * 0.5 * (Fo2p + F_o2p(n1,h1,r_ind,r1_dep, ftint))) > 0.

np.el = (np.sp + 2.*np.s2p + 3.*np.s3p + 4.*np.s4p + np.op + 2.*np.o2p)/(1.-np.protons)
np.elh = np.fh/(1.-np.fh) * np.el

nTp.sp  = (nT.sp  + dt * 0.5 * (EFsp  + EF_sp(n1,T1,r_ind,r1_dep,h1,nu1, ftint)))  > 0.
nTp.s2p = (nT.s2p + dt * 0.5 * (EFs2p + EF_s2p(n1,T1,r_ind,r1_dep,h1,nu1, ftint))) > 0.
nTp.s3p = (nT.s3p + dt * 0.5 * (EFs3p + EF_s3p(n1,T1,r_ind,r1_dep,h1,nu1, ftint))) > 0.
;nTp.s4p = (nT.s4p + dt * 0.5 * (EFs4p + EF_s4p(n1,T1,r_ind,r1_dep,h1,nu1, ftint))) > 0.
nTp.op  = (nT.op  + dt * 0.5 * (EFop  + EF_op(n1,T1,r_ind,r1_dep,h1,nu1, ftint)))  > 0.
nTp.o2p = (nT.o2p + dt * 0.5 * (EFo2p + EF_o2p(n1,T1,r_ind,r1_dep,h1,nu1, ftint))) > 0.
nTp.el  = (nT.el  + dt * 0.5 * (EFel  + EF_el(n1,T1,h1, r_ind,r1_dep,lat_dist1,nu1, ftint))) > 0.
print, "EFel, EFel(stuff) = ", EFel, EF_el(n1,T1,h1, r_ind,r1_dep,lat_dist1,nu1, ftint)

update_temp,np,nTp,Tp

n = np
nT = nTp
T = Tp
if check_math(mask=211,/noclear) gt 0 then stop

return
end

;+
;  SUBROUTINE
;    energy_budget
;
;  ABSTRACT
;    This code was written to combine the overall energy budget (coded by
;    Delamere for his cubic centimeter model) as well as include the
;    individual reactions (coded by Steffl) which are necessary for this
;    particular cubic centimeter, latitudinally averaged, model.
;
;  PUROPOSE
;    Produce a breakdown of the total energy budget for a given run of
;    the model. The breakdown is supplied as a structure.
;
;  Note: output 'energy' structure needs to be previously defined
;
;  Written by Adam Shinn
;-
pro energy_budget, n, h, T, r_dep, r_ind, ftint, lat_dist, nu, energy;, print = print

  ; compute energy budget
  ;rootpi = sqrt(!dpi)

  ; thermal equilibrium of sp
  Teq_sp = nu.sp_s2p*(T.s2p-T.sp) + $
           nu.sp_s3p*(T.s3p-T.sp) + $
           nu.sp_op* (T.op -T.sp) + $
           nu.sp_o2p*(T.o2p-T.sp) + $
           nu.sp_el* (T.el- T.sp)

  ; ionization, charge exchange, fast neutrals and transport of sp
  ion_sp   = T.pu_s*ftint.is /(!rootpi*h.sp)
  ionh_sp  = T.pu_s*ftint.ish/(!rootpi*h.sp)
  cx_sp    = T.pu_s*(ftint.cx_k1 + ftint.cx_k2 + ftint.cx_k4 + ftint.cx_k9 + ftint.cx_k10)/(!rootpi*h.sp)
  fast_sp  = T.sp*(ftint.cx_k1 + ftint.cx_k8)/(!rootpi*h.sp)
  trans_sp = T.sp*(r_dep.transport*n.sp)

  ; thermal equilibrium of s2p
  Teq_s2p = nu.sp_s2p* (T.sp- T.s2p) + $
            nu.s2p_s3p*(T.s3p-T.s2p) + $
            nu.s2p_op* (T.op- T.s2p) + $
            nu.s2p_o2p*(T.o2p-T.s2p) + $
            nu.s2p_el* (T.el- T.s2p)

  ; ionization, charge exchange, fast neutrals and transport of s2p
  ion_s2p   = T.sp*ftint.isp /(!rootpi*h.s2p)
  ionh_s2p  = T.sp*ftint.isph/(!rootpi*h.s2p)
  cx_s2p    = T.pu_s*(ftint.cx_k3 + ftint.cx_k11)/(!rootpi*h.s2p)
  fast_s2p  = T.s2p*(ftint.cx_k3)/(!rootpi*h.s2p)
  trans_s2p = T.s2p*(r_dep.transport*n.s2p)

  ; thermal equilibrium of s3p
  Teq_s3p = nu.sp_s3p* (T.sp- T.s3p) + $
            nu.s2p_s3p*(T.s2p-T.s3p) + $
            nu.s3p_op* (T.op- T.s3p) + $
            nu.s3p_o2p*(T.o2p-T.s3p) + $
            nu.s3p_el* (T.el- T.s3p)

  ; ionization, fast neutrals and transport of s3p
  ion_s3p   = T.s2p*ftint.is2p /(!rootpi*h.s3p)
  ionh_s3p  = T.s2p*ftint.is2ph/(!rootpi*h.s3p)
  ;cx_s3p    = T.s2p*ftint.cx_k15/(!rootpi*h.s3p)
  trans_s3p = T.s3p*(r_dep.transport*n.s3p)

  ; thermal equilibrium of op
  Teq_op = nu.sp_op* (T.sp- T.op) + $
           nu.s2p_op*(T.s2p-T.op) + $
           nu.s3p_op*(T.s3p-T.op) + $
           nu.op_o2p*(T.o2p-T.op) + $
           nu.op_el* (T.el- T.op)

  ; ionization, charge exchange, fast neutrals and transport of op
  ion_op   = T.pu_o*ftint.io /(!rootpi*h.op)
  print, "T.pu_o, ftint.io, h.op = ", T.pu_o, ftint.io, h.op
  ionh_op  = T.pu_o*ftint.ioh/(!rootpi*h.op)
  cx_op    = T.pu_o*(ftint.cx_k5 + ftint.cx_k6 + ftint.cx_k8 + ftint.cx_k12 + ftint.cx_k14)/(!rootpi*h.op)
  fast_op  = T.op*(ftint.cx_k5 + ftint.cx_k9)/(!rootpi*h.op)
  print, "T%op, ft%cx(7), ft%cx(10) = ", T.op, ftint.cx_k5, ftint.cx_k9
  trans_op = T.op*(r_dep.transport*n.op)
  
  ; thermal equilibrium of o2p
  Teq_o2p = nu.sp_o2p* (T.sp- T.o2p) + $
            nu.s2p_o2p*(T.s2p-T.o2p) + $
            nu.s3p_o2p*(T.s3p-T.o2p) + $
            nu.op_o2p* (T.op- T.o2p) + $
            nu.o2p_el* (T.el- T.o2p)

  ; ionization, charge exchange, fast neutrals and transport of o2p
  ion_o2p   = T.op*ftint.iop /(!rootpi*h.o2p)
  ionh_o2p  = T.op*ftint.ioph/(!rootpi*h.o2p)
  cx_o2p    = T.pu_o*(ftint.cx_k7)/(!rootpi*h.o2p)
  fast_o2p  = T.o2p*(ftint.cx_k7)/(!rootpi*h.o2p)
  trans_o2p = T.o2p*(r_dep.transport*n.o2p)

  ; thermal equilibrium of hot electrons
  Teq_elh = nu.sp_elh *(T.elh - T.sp ) + $
            nu.s2p_elh*(T.elh - T.s2p) + $
            nu.s3p_elh*(T.elh - T.s3p) + $
            nu.op_elh *(T.elh - T.op ) + $
            nu.o2p_elh*(T.elh - T.o2p) + $
            nu.el_elh *(T.elh - T.el )
 
  ; conversion from temperature to energy
  f = 3./2.

  ; energy of thermal equilibration of hot electrons and 
  ; the total thermal equilibration
  energy.eh_eq  = f*Teq_elh
  energy.tot_eq = f*(Teq_sp + Teq_s2p + Teq_s3p + Teq_op + Teq_o2p)

  ; total energy from ionization and charge exchange for both sulfur and oxygen
  energy.s_ion = f*(ion_sp + ionh_sp + ion_s2p + ionh_s2p + ion_s3p + ionh_s3p)
  ;print, "energy.sp, h_sp, s2p, h_s2p, s3p, h_s3p = ", ion_sp,  ionh_sp,  ion_s2p,  ionh_s2p,  ion_s3p,  ionh_s3p
  energy.o_ion = f*(ion_op + ionh_op + ion_o2p + ionh_o2p)
  energy.s_cx  = f*(cx_sp + cx_s2p); + cx_s3p)
  energy.o_cx  = f*(cx_op + cx_o2p)

  ; total power in
  energy.P_in = energy.s_ion + energy.o_ion + energy.s_cx + energy.o_cx + energy.eh_eq; + energy.tot_eq

  ; total power from fast neutrals and uv radiation
  energy.Pfast = f*(fast_sp + fast_s2p + fast_op + fast_o2p)
  energy.Puv = ft_rad(lat_dist, T, r_ind)

  ; total power from transport
  energy.Ptrans = f*r_dep.transport*(n.sp*T.sp + n.s2p*T.s2p + n.s3p*T.s3p + n.op*T.op + n.o2p*T.o2p + n.el*T.el)
  energy.Ptrans_eh = f*r_dep.transport*n.el*n.fh*T.elh

  ; total power out
  energy.P_out = energy.puv + energy.pfast + energy.ptrans + energy.ptrans_eh; - energy.tot_eq

  ;energy = {energy, s_ion: 0.0, s_cx: 0.0, o_ion: 0.0, o_cx: 0.0, eh_eq: 0.0, tot_eq: 0.0, P_in: 0.0, $
  ;          Puv: 0.0, Pfast: 0.0, Ptrans: 0.0, Ptrans_eh: 0.0, P_out: 0.0}

  ;; if keyword_set(print) then begin
  ;;    print, '$$--------------------------------'
  ;;    print, '$$ IN-CODE ENERGY BUDGET'
  ;;    print, '$$--------------------------------'
  ;;    print, '$$ ionized S............', energy.s_ion
  ;;    print, '$$ ionized O............', energy.o_ion
  ;;    print, '$$ charge exchange S....', energy.s_cx
  ;;    print, '$$ charge exchange O....', energy.o_cx
  ;;    print, '$$ equil with ehot......', energy.eh_eq + energy.tot_eq
  ;;    print, '$$ total in.............', energy.P_in + energy.tot_eq
  ;;    print, '$$ puv..................', energy.Puv
  ;;    print, '$$ fast/ena.............', energy.pfast - energy.tot_eq
  ;;    print, '$$ transport............', energy.ptrans + energy.ptrans_eh
  ;;    print, '$$ total out............', energy.P_out - energy.tot_eq
  ;;    print, '$$ in/out...............', (energy.P_in + energy.tot_eq)/(energy.P_out - energy.tot_eq)
  ;;    print
  ;; endif

end

