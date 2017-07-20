;+
;  IO TORUS CUBIC CENTIMETER ONEBOX MODEL, MAIN PROCEDURE
;     
;  Supporting subroutines can be found in cm3_latavg_model.pro
;
;  Written by Andrew Steffl, tailored by Adam Shinn to in order to
;  create a widget that runs this model for a wide range of parameter
;  space.
;-
;-------------------------------------------------------------------
PRO cm3_latavg_onebox, n_out, t_out, h_out, f_out, r_ind, r_dep, energy, $
                       runt1=runt1, $
                       lag_const1=lag_const1, lag_phase1=lag_phase1, lag_amp1=lag_amp1, $
                       fehot_const1=fehot_const1, $
                       tm1 = tm1, $
                       Tehot=Tehot, plot=plot, $
                       filename=filename, nosave = nosave, $
                       transport=transport, source=source, o_to_s=o_to_s, $
                       contfile=contfile, protons=protons, lon=lon, $
                       neutral_amp=neutral_amp, neutral_t0=neutral_t0, neutral_width=neutral_width, $
                       hote_amp=hote_amp, hote_t0=hote_t0, hote_width=hote_width, $
                       trans_type1=trans_type1, trans_exp1=trans_exp1, o2s_spike1 = o2s_spike1, rdist1 = rdist1, $
                       n_height = n_height, printiter = printiter, $
                       xwinsize = xwinsize, ywinsize = ywinsize
;-------------------------------------------------------------------
runtime=systime(1)
@cm3_model_common
getsysvariables

model_dir = 'model/cm3_latavg_onebox/'
fname = '/home/dcoffin/1D_Model/Eddiecode/iotorus/model/cm3_latavg_onebox/theta.dat'
openw,21,fname
IF NOT keyword_set(filename) THEN filename = '_output'

tepdir = model_dir
!p.multi=[0,1,1]

print, " source = ", source

; Open new window(s) if needed to plot the results
;;IF keyword_set(plot) THEN BEGIN
  ;; loadct,39
  ;; device, decomposed = 0
  ;; oldwinid=!d.window
 
   ;restore, file = model_dir + 'avgmixratios.sav'

  ;; IF NOT keyword_set(xwinsize) THEN xwinsize = 640
 ;;  IF NOT keyword_set(ywinsize) THEN ywinsize = 512

;;   window, /free, /pixmap, xs = xwinsize, ys = ywinsize
;;   comp_pixmap_id = !d.window
;;   window, /free, xs = xwinsize, ys = ywinsize
;;   comp_window_id = !d.window

;;ENDIF

; The plasma is not corotating with magnetic field
if n_elements(lag_const1) gt 0 then lag_const=lag_const1 else lag_const=0.0 ; km/s
if n_elements(lag_phase1) gt 0 then lag_phase=lag_phase1 else lag_phase=0. 
if n_elements(lag_amp1) gt 0 then lag_amp=lag_amp1 else lag_amp=0.    ; km/s

; What is the transport exponent and is the transport rate dependednt
; on the neutral source or the total ion density
if keyword_set(trans_exp1) then trans_exp=trans_exp1 else trans_exp=1.
if keyword_set(trans_type1) then trans_type=1 else trans_type=0

; Model time step in seconds (slightly longer than the time for
; neutrals to cross a 15 degree azimuthal bin)
dt = 750.     
printf, 21, FORMAT = '("dt: ",I0)', dt ;;;;;

if n_elements(runt1) EQ 0 THEN runt = 220.*8.64e4 ELSE runt=runt1*8.64e4; Total run time
printf,21,FORMAT = '("runt: ",I0)', runt ;;;;;

nit = fix(runt/dt)+1     ; Number of iterations
printf,21,FORMAT = '("nit: ",I0)', nit ;;;;;

; Radial distance of model, default is 6 R_J
IF keyword_set(rdist1) THEN rdist = rdist1 ELSE rdist = 6.

printf,21,FORMAT = '("rdist: ",I0)', rdist ;;;;; 

; System III longitude of box
if keyword_set(lon) then l3 = lon ELSE l3 = 110.
lontemp = l3

printf,21,FORMAT = '("lontemp: ",I0)', lontemp ;;;;;

; Neutral Cloud scale height in km
theta_offset = 6.4 * cos((l3 - 200.) * !dtor) * !dtor

printf,21,FORMAT = '("theta_offset: ",F0)', theta_offset ;;;;;

zoff = abs(theta_offset * rdist * 71492.) ; in km, not cm! Always positive

printf,21,FORMAT = '("zoff: ",F0)', zoff ;;;;;

IF NOT keyword_set(n_height) THEN n_height = 0.5 * 71492.

printf,21,FORMAT = '("n_height: ",I0)', n_height ;;;;;
; N.B. the tag "el" refers to the thermal electron population only. Similarly,
; "elh" is only for the hot electron population.
n = {dens_lon,sp:0.0d,s2p:0.0d, s3p: 0.0d, s4p: 0.0d, op: 0.0d, $
  o2p: 0.0d, el: 0.0d, elh: 0.0d, s: 0.0d, o: 0.0d, fc:0d, fh:0d, protons:0d}

printf,21,FORMAT = '("n: ",I0)', n ;;;;;

nT = {energy_lon, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p: 0.0d, op: 0.0d, $
                 o2p: 0.0d, el: 0.0d, elh: 0.0d}

printf,21,FORMAT = '("nT: ",I0)', nT ;;;;;

T = {temp_lon,sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p: 0.0d,  $
                op: 0.0d, o2p: 0.0d, el: 0.0d, elh: 0.0d, pu_s: 0d, pu_o: 0d}

printf,21,FORMAT = '("T: ",I0)', T ;;;;;

h = {scale_heights_lon, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
                o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}

printf,21,FORMAT = '("h: ",I0)', h ;;;;;

; Structures used to store intermediate values in the Euler advance subroutine
n1=n
nT1=nT
T1=T
h1=h

np=n
nTp=nT
Tp=T

; Structure containing volume coefficients for energy transfer via coulomb collisons
nu = {nu, sp_s2p:0d, sp_s3p:0d, sp_s4p:0d, sp_op:0d, sp_o2p:0d, s2p_s3p:0d, $
      s2p_s4p:0d, s2p_op:0d, s2p_o2p:0d, s3p_s4p: 0d, s3p_op:0d, s3p_o2p:0d, $
      s4p_op: 0d, s4p_o2p:0d, op_o2p:0d, sp_el:0d, s2p_el:0d, s3p_el:0d, s4p_el:0d,$
      op_el:0d, o2p_el:0d, sp_elh:0d, s2p_elh:0d, s3p_elh:0d, s4p_elh:0d,$
      op_elh:0d, o2p_elh:0d, el_elh: 0d}

printf,21,FORMAT = '("nu: ",I0)', nu ;;;;;

; structure containing energy values
;; energy = {energy, s_ion: 0.0, s_cx: 0.0, o_ion: 0.0, o_cx: 0.0, eh_eq: 0.0, tot_eq: 0.0, P_in: 0.0, $
;;           Puv: 0.0, Pfast: 0.0, Ptrans: 0.0, Ptrans_eh: 0.0, P_out: 0.0}

; **to be deleted**
; neutral loss structure, values generated by get_neutral_loss_rate
;nl = {neutral_loss, sei: 0.0, seih: 0.0, oei: 0.0, oeih: 0.0, $
;      scx: 0.0, ocx: 0.0, k0: 0.0, k1: 0.0, k2: 0.0, k3: 0.0, k4: 0.0, k5: 0.0, $
;      k6: 0.0, k7: 0.0, k8: 0.0, k9: 0.0, k10:0.0, k11: 0.0, k12: 0.0, $
;      k14: 0.0, s_tot: 0.0, o_tot: 0.0, scx_tot: 0.0, ocx_tot: 0.0, $
;      sion_tot: 0.0, oion_tot: 0.0, fast_O_k5: 0.0, fast_O_k7: 0.0, $
;      fast_S_k1: 0.0, fast_S_k3: 0.0, fast_O_k9: 0.0, fast_S_k8: 0.0}

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
r_dep = {dependent_rates, is: 0.0, isp: 0.0, is2p: 0.0, is3p: 0.0, $
  io: 0.0, iop:0.0, io2p: 0.0, rsp: 0.0, rs2p: 0.0, rs3p: 0.0, rop: 0.0, ro2p:0.0,$
  Transport:0.0}

printf,21,FORMAT = '("r_dep: ",I0)', r_dep ;;;;;
r1_dep=r_dep

; Latitudinal Density structure
lat_dist = {lat_dist, z:dindgen(21)/5. * 71492., sp:dblarr(21), s2p:dblarr(21), $
             s3p:dblarr(21), s4p:dblarr(21), op:dblarr(21), o2p:dblarr(21), $
             el:dblarr(21)}

printf,21,FORMAT = '("lat_dist: ",I0)', lat_dist ;;;;;
lat_dist1 = lat_dist

; flux tube integrated number of reactions per second.
ftint = {ftint, cx_k0:0d, cx_k1:0d, cx_k2:0d, cx_k3:0d, cx_k4:0d, cx_k5:0d, $
                  cx_k6:0d, cx_k7:0d, cx_k8:0d, cx_k9:0d, cx_k10:0d, cx_k11:0d, cx_k12:0d, $
                  cx_k13:0d, cx_k14:0d, cx_k15:0d, cx_k16:0d, is:0d, isp:0d, is2p:0d, $
                  is3p:0d, io:0d, iop:0d, io2p:0d, ish:0d, isph:0d, is2ph:0d, is3ph:0d, $
                  ioh:0d, ioph:0d, io2ph:0d, rsp:0d, rs2p:0d, rs3p:0d, rop:0d, ro2p:0d}

printf,21,FORMAT = '("ftint: ",I0)', ftint ;;;;;

; Output quantites, that will contain the values of density and
; temperature in every longitudinal box at every time step
n_out=n

printf,21,FORMAT = '("n_out: ",I0)', n_out ;;;;;

t_out=t

printf,21,FORMAT = '("t_out: ",I0)', t_out ;;;;;

h_out=h

printf,21,FORMAT = '("h_out: ",I0)', h_out ;;;;;

; Initialize time, phase, and amplitude variables
IF n_elements(tm1) EQ 0 THEN tm = 170. ELSE tm = tm1
tm0=tm

printf,21,FORMAT = '("tm: ",I0)', tm ;;;;;

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

    printf,21,FORMAT = '("trans_out: ",I0)', trans_out ;;;;;

    o2s_out = otos

    printf,21,FORMAT = '("o2s_out: ",I0)', o2s_out ;;;;;

    source_out=source_out[n_out_dim[1]-1]

    printf,21,FORMAT = '("source_out: ",I0)', source_out ;;;;;

    net_source = source_out

    printf,21,FORMAT = '("net_source: ",I0)', net_source ;;;;;

; Start At DOY 170, which is June 18, in a leap year (as 2000 was)
    IF n_elements(tm1) EQ 0 THEN tm = 170. ELSE tm = tm1

endif else begin

   n.sp=  0.060 * 1800.
   printf,21,FORMAT = '("n.sp: ",I0)', n.sp ;;;;;
   n.s2p= 0.212 * 1800.
   printf,21,FORMAT = '("n.s2p: ",I0)', n.s2p ;;;;;
   n.s3p= 0.034 * 1800.
   printf,21,FORMAT = '("n.s3p: ",I0)', n.s3p ;;;;;
   n.op=  0.242 * 1800.
   printf,21,FORMAT = '("n.op: ",I0)', n.op ;;;;;
   n.o2p= 0.123 * n.op 
   printf,21,FORMAT = '("n.o2p: ",I0)', n.o2p ;;;;;

   n.s=25.
   printf,21,FORMAT = '("n.s: ",I0)', n.s ;;;;;
   n.o=50.
   printf,21,FORMAT = '("n.o: ",I0)', n.o ;;;;;
   n.protons = 0.1

; January Latavg Conditions
    Te0 = 5.0d
    printf,21,FORMAT = '("Te0: ",I0)', Te0 ;;;;;
    Ti0 = 70.0d
    printf,21,FORMAT = '("Ti0: ",I0)', Ti0 ;;;;;
    Teh0 = 40.0d
    printf,21,FORMAT = '("Teh0: ",I0)', Teh0 ;;;;;
    n.fh = fehot_const1 ; = 0.003d 
    printf,21,FORMAT = '("n.fh: ",F0)', n.fh ;;;;;
    trans = 1.0d/(70.28*8.64e4)
    printf,21,FORMAT = '("trans: ",I0)', trans ;;;;;
; with a 0.5 Rj scale height for neutrals, this comes out to 9.9e-4 cm-3 s-1
; at the equator
    net_source = source ;6.3e6
    printf,21,FORMAT = '("net_source: ",I0)', net_source ;;;;;
    otos = 1.7d
    printf,21,FORMAT = '("otos: ",I0)', otos ;;;;;

; Total electron density
    n.el = (n.sp + 2 * n.s2p + 3 * n.s3p + n.op + 2 * n.o2p) / (1.-n.protons)
    ;printf,21,FORMAT = '("n.protons: ",F0)', n.protons ;;;;
    printf,21,FORMAT = '("n.el: ",I0)', n.el ;;;;;
    n.elh = n.fh/(1.-n.fh) * n.el               ; Hot electron density
    printf,21,FORMAT = '("n.elh: ",F0)', n.elh ;;;;;
    n.fc = 1.-n.fh
    printf,21,FORMAT = '("n.fh: ",F0)', n.fh ;;;;;

; Assign Temp values
    T[*].sp = Ti0
    T[*].s2p = Ti0
    T[*].s3p = Ti0
    ;T[*].s4p = Ti0
    T[*].op = Ti0
    T[*].o2p = Ti0
    T[*].el = Te0
    T[*].elh = Teh0

; Determine initial scale heights
    get_scale_heights,h,t,n

    source_out=net_source
    n_out=n
    t_out=t
    h_out=h
    trans_out=avg(r_dep.transport)
    o2s_out = otos
 ENDELSE

; Allow for protons to be included in the model. Right now, any added
; protons do not take part in the chemistry i.e. no p + O -> H + O+
; type of reactions. Therefore, their only real effect is to increase
; the electron denisty while keeping the ion densities
; constant. Default is not to include protons
IF keyword_set(protons) THEN BEGIN
    n.protons  = protons
    printf,21,FORMAT = '("n.protons: ",I0)', n.protons ;;;;;
    n1.protons = protons
    np.protons = protons
 ENDIF


if n_elements(fehot_const1) gt 0 then fehot_const=fehot_const1
if keyword_set(Tehot) then Teh0=Tehot
if keyword_set(transport) then trans = 1./(transport*8.64e4)
if keyword_set(source) then net_source = source
if keyword_set(o_to_s) then otos= o_to_s
IF keyword_set(o2s_spike1) THEN o2s_spike = o2s_spike1 ELSE o2s_spike = otos
r_ind.o_to_s = otos
r_ind.o2s_spike = o2s_spike

; Initialize time variable neutral source and hot electron fraction parameters
tau0=1./trans/8.64e4
printf,21,FORMAT = '("tau0: ",I0)', tau0 ;;;;;
net_source0 = net_source
printf,21,FORMAT = '("net_source0: ",I0)', net_source0 ;;;;;
fh0 = fehot_const
printf,21,FORMAT = '("fh0: ",I0)', fh0 ;;;;;

if not keyword_set(neutral_amp) then neutral_amp=0.
printf,21,FORMAT = '("neutral_amp: ",I0)', neutral_amp ;;;;;
if not keyword_set(neutral_t0) then neutral_t0=249. ; 5-Sept-2000 as per Delamere et al 2004
printf,21,FORMAT = '("neutral_t0: ",I0)', neutral_t0 ;;;;;
if not keyword_set(neutral_width) then neutral_width=25.0
printf,21,FORMAT = '("neutral_width: ",I0)', neutral_width ;;;;;

if not keyword_set(hote_amp) then hote_amp=0.0000
printf,21,FORMAT = '("hote_amp: ",I0)', hote_amp ;;;;;
if not keyword_set(hote_t0) then hote_t0= 269.
printf,21,FORMAT = '("hote_t0: ",I0)', hote_t0 ;;;;;
if not keyword_set(hote_width) then hote_width=45.5
printf,21,FORMAT = '("hote_width: ",I0)', hote_width ;;;;;
IF keyword_set(n_height) THEN BEGIN 
   h.s = n_height
   printf,21,FORMAT = '("h.s: ",I0)', h.s ;;;;;
   h.o = n_height
   printf,21,FORMAT = '("h.0: ",I0)', h.o ;;;;;
ENDIF

; Read in total emitted power arrays.
read_emission_tables, tepdir + 'tepSII_v7.0.2',temp,den,emisSp
read_emission_tables, tepdir + 'tepSIII_v7.0.2',temp,den,emisS2p
read_emission_tables, tepdir + 'tepSIV_v7.0.2',temp,den,emisS3p
read_emission_tables, tepdir + 'tepOII_v7.0.2',temp,den,emisOp
read_emission_tables, tepdir + 'tepOIII_v7.0.2',temp,den,emisO2p

; Store info into r_ind structure. These fields should not get modified by
; any other subroutine
r_ind.emistemp = temp
;rintf,21,FORMAT = '("r_ind.emistemp: ",I0)', r_ind.emistemp ;;;;;
r_ind.emisden  = den
;rintf,21,FORMAT = '("r_ind.emisden: ",I0)', r_ind.emisden ;;;;;
r_ind.emisSp   = emisSp
;rintf,21,FORMAT = '("r_ind.emisSp: ",F0)', r_ind.emisSp ;;;;;
r_ind.emisS2p  = emisS2p
;printf,21,FORMAT = '("r_ind.emisS2p: ",I0)', r_ind.emisS2p ;;;;;
r_ind.emisS3p  = emisS3p
;rintf,21,FORMAT = '("r_ind.emisS3p: ",I0)', r_ind.emisS3p ;;;;;
r_ind.emisOp   = emisOp
;rintf,21,FORMAT = '("r_ind.emisOp: ",I0)', r_ind.emisOp ;;;;;
r_ind.emisO2p  = emisO2p
;rintf,21,FORMAT = '("r_ind.emisO2p: ",I0)', r_ind.emisO2p ;;;;;

;set initial pickup temperatures
Lshl = rdist
T.pu_s = Tpu(32.,Lshl)
printf,21,FORMAT = '("T.pu_s: ",I0)', T.pu_s ;;;;;
T.pu_o = Tpu(16.,Lshl)
printf,21,FORMAT = '("T.pu_o: ",I0)', T.pu_o ;;;;;
T.elh = Teh0
printf,21,FORMAT = '("T.elh: ",I0)', T.elh ;;;;;

T1.pu_s = T.pu_s
printf,21,FORMAT = '("T1.pu_s: ",I0)', T1.pu_s ;;;;;
T1.pu_o = T.pu_o
printf,21,FORMAT = '("T1.pu_o: ",I0)', T1.pu_o ;;;;;
T1.elh = Teh0
printf,21,FORMAT = '("T1.elh: ",I0)', T1.elh ;;;;;

Tp.pu_s = T.pu_s
printf,21,FORMAT = '("Tp.pu_s: ",I0)', Tp.pu_s ;;;;;
Tp.pu_o = T.pu_o
printf,21,FORMAT = '("Tp.pu_o: ",I0)', Tp.pu_o ;;;;;
Tp.elh = Teh0
printf,21,FORMAT = '("Tp.elh: ",I0)', Tp.elh ;;;;;

; Call get_independent_rates to get the reaction rates for charge
; exchange, hot electron ionization, and neutral source production. As
; the name implies, these rates do not change over the run.
get_independent_rates,r_ind,T,h

net_source0 = net_source

iondens0 = avg(n.sp + n.s2p + n.s3p + n.s4p + n.op + n.o2p)

printf,21,FORMAT = '("iondens0: ",I0)', iondens0 ;;;;;
; Angular velocity (rad/sec) for circular orbit about Jupiter at rdist
;omega_neutrals = 1./sqrt((rdist * !rjkm * 1e3)^3/(6.673e-11 * 1.8987e27))
rdistkm = rdist*!rjkm
printf,21,FORMAT = '("rdistkm: ",I0)', rdistkm ;;;;;
omega_neutrals = 1./sqrt( (1e9)*(rdistkm*rdistkm*rdistkm) /(6.673e-11 * 1.8987e27))
printf,21,FORMAT = '("omega_neutrals: ",I0)', omega_neutrals ;;;;;
v_neutrals = omega_neutrals * rdistkm ; km/s
printf,21,FORMAT = '("v_neutrals: ",I0)', v_neutrals ;;;;;
v_corotation = 2. * !pi * rdistkm / (9.925 * 3600.) ; km/s
printf,21,FORMAT = '("v_corotation: ",I0)', v_corotation ;;;;;
rel_v_neutrals = v_corotation - v_neutrals
printf,21,FORMAT = '("rel_v_neutrals: ",I0)', rel_v_neutrals ;;;;;

; Final consistency checks
; Check that the sum of the thermal & hot electron fractions equals one
n.fc = 1. - n.fh
printf,21,FORMAT = '("n.fc: ",I0)', n.fc ;;;;;

;  Ensure charge neutrality in the thermal electron population. We
;  assume that the hot electrons are linked to a separate hot proton population.
;n.el = (n.sp + 2.*n.s2p + 3.*n.s3p + 4.*n.s4p + n.op + 2.*n.o2p)/(1.-n.protons)
n.el = (n.sp + 2.*n.s2p + 3.*n.s3p + n.op + 2.*n.o2p)/(1.-n.protons)
printf,21,FORMAT = '("n.el: ",I0)', n.el ;;;;;
n.elh = n.fh/(1.-n.fh) * n.el
printf,21,FORMAT = '("n.elh: ",I0)', n.elh ;;;;;

nT.el  = n.el  * T.el
printf,21,FORMAT = '("nT.el: ",I0)', nT.el ;;;;;
nT.elh  = n.elh * T.elh
printf,21,FORMAT = '("nT.elh: ",I0)', nT.elh ;;;;;
nT.sp  = n.sp  * T.sp
printf,21,FORMAT = '("nT.sp: ",I0)', nT.sp ;;;;;
nT.s2p = n.s2p * T.s2p
printf,21,FORMAT = '("nT.s2p: ",I0)', nT.s2p ;;;;;
nT.s3p = n.s3p * T.s3p
printf,21,FORMAT = '("nT.s3p: ",I0)', nT.s3p ;;;;;
;nT.s4p = n.s4p * T.s4p
nT.op  = n.op  * T.op
printf,21,FORMAT = '("nT.op: ",I0)', nT.op ;;;;;
nT.o2p = n.o2p * T.o2p
printf,21,FORMAT = '("nT.o2p: ",I0)', nT.o2p ;;;;;

n1=n
nT1=nT
T1=T
h1=h

np=n
nTp=nT
Tp=T

get_scale_heights,h,t,n

;IF keyword_set(printiter) THEN format="('% CM3_LATAVG_ONEBOX: Iteration ', i5, ' of ', i5)"
badtime = 0l
;=================================================================
; Loop over # of interations
;=================================================================
for j=0,nit-1 do begin 

   ;IF keyword_set(printiter) AND (j MOD 100) EQ 0 THEN print, $
   ;  strtrim(j, 2), strtrim(nit, 2), format = format

   if check_math(mask=211,/noclear) gt 1 then begin
      print,'Math error detected! Exiting loop.'
      break
   endif

   time = tm0 + j * dt / 86400.
   printf,21,FORMAT = '("time: ",I0)', time ;;;;;
; Time variable neutral source
   net_source = net_source0*(1. + neutral_amp*exp(-(time-neutral_t0)^2/neutral_width^2))
   printf,21,FORMAT = '("net_source: ",I0)', net_source ;;;;;
; Let the neutral source be the sum of two independent sources: a
; baseline source with an O/S of o_to_s and a Gaussian increase with a
; potentially different O/S of o2s_spike
   r_ind.o_to_s = (otos + $
                   o2s_spike * neutral_amp*exp(-(time-neutral_t0)^2/neutral_width^2))/$
                  (1+neutral_amp*exp(-(time-neutral_t0)^2/neutral_width^2))
  ;printf,21,FORMAT = '("r_ind.o_to_s: ",I0)', r_ind.o_to_s ;;;;;
; Time variable hot electron fraction
   n.fh = fehot_const * (1. + hote_amp * exp(-(time - hote_t0)^2 / hote_width^2))
   printf,21,FORMAT = '("n.fh: ",I0)', n.fh ;;;;;
   
   negative_fh = where(n.fh LT 0, nbadfh)
   IF nbadfh GT 0 THEN BEGIN
      badtime = [badtime, j]
      n.fh >= 0.
   ENDIF

   n1.fh = n.fh
   np.fh = n.fh

   n.fc  = 1d - n.fh
   n1.fc = 1d - n.fh
   np.fc = 1d - n.fh

; Update the hot electron densities
   n.elh = n.fh/(1.-n.fh) * n.el ; ** could be re-written n.fh/n.fc **
   printf,21,FORMAT = '("n.elh: ",I0)', n.elh ;;;;;
   nT.elh  = n.elh * T.elh
   printf,21,FORMAT = '("nT.elh: ",I0)', nT.elh ;;;;;

; Call cm3_latavg_model to advance the model one time step
   cm3_latavg_model, n,  T,  nT,  h,  nu,  $
                     n1, T1, nT1, h1, nu1, $
                     np, Tp, nTp, $
                     r_ind, r_dep, r1_dep, lat_dist, lat_dist1, ftint

; Update the temperature and scale height structures
   update_temp,n,nT,T
   get_scale_heights,h,t,n

  if not keyword_set(plot) then begin

  if j mod 5 eq 0 then begin
    ftmix = ftint_mix(n,h)
    n_out = n
    t_out = t
    h_out = h
    f_out = ftmix
    trans_out = avg(r_dep.transport)
    o2s_out = r_ind.o_to_s
    source_out = net_source
    spr = ftmix.s2p/ftmix.s3p
    if j eq 0 then dspr = 1 else dspr = abs((spr_old - spr)/.01)
    spr_old = spr

    energy_budget, n, h, T, r_dep, r_ind, ftint, lat_dist, nu, energy
    P_IO = (energy.P_in + energy.tot_eq)/(energy.P_out - energy.tot_eq)

  endif
  
  endif else begin
  
  ; Fitting things every 2000 seconds is overkill, so only do this on
  ; every 5th timestep (i.e. every 1e4 s)

   IF j MOD 5 EQ 0 THEN BEGIN 
      tm = [tm,time]
      
      ftmix=ftint_mix(n,h)
      printf,21,FORMAT = '("ftmix: ",F0)', ftmix ;;;;;
      n_out = n
      t_out = t
      h_out = h
      f_out = ftmix
      trans_out=[trans_out,avg(r_dep.transport)]
      printf,21,FORMAT = '("trans_out: ",I0)', trans_out ;;;;;
      o2s_out = [o2s_out, r_ind.o_to_s]
      printf,21,FORMAT = '("o2s_out: ",I0)', o2s_out ;;;;;
      source_out = [source_out, net_source]
      printf,21,FORMAT = '("source_out: ",I0)', source_out ;;;;;

      energy_budget, n, h, T, r_dep, r_ind, ftint, lat_dist, nu, energy
      P_IO = (energy.P_in + energy.tot_eq)/(energy.P_out - energy.tot_eq)
      printf,21,FORMAT = '("P_IO: ",I0)', P_IO ;;;;;
      spr = ftmix.s2p/ftmix.s3p
      printf,21,FORMAT = '("spr: ",I0)', spr ;;;;;
     
      ;if j EQ 0 then ftmix_out = [ftmix, ftmix] else ftmix_out=[ftmix_out, ftmix]

      if j eq 0 then begin
        ; initialize plot
        ;loadct, 39
        ;device, decomposed = 0
   ;;     plot, [tm,runt1+tm], [0.001,1.0], /nodata, $
;;              ytitle = 'ft int ratio', xtitle = 'time [doy]', /xstyle, $ ;, ystyle = 8, $
;;              xmargin = [8, 8], ymargin = [4, 4], charsize = 1.25, /ylog
 ;;       xyouts, .15, .95, /normal,  'sp', col =  57, charsize = 1.25
 ;;       xyouts, .20, .95, /normal, 's2p', col = 147, charsize = 1.25
 ;;       xyouts, .25, .95, /normal, 's3p', col = 210, charsize = 1.25
 ;;       xyouts, .30, .95, /normal,  'op', col = 254, charsize = 1.25
  ;;      xyouts, .35, .95, /normal, 'o2p', col = 100, charsize = 1.25
 ;;       xyouts, .40, .95, /normal, 'd/dt s2p/s3p', col = 180, charsize = 1.25
        dspr = 1 ; initialize dspr for break check
      endif else begin
        dspr = abs((spr_old - spr)/.01)
        printf,21,FORMAT = '("dspr: ",I0)', dspr ;;;;;
        ;print, dspr

 ;;       oplot, [time], [ftmix.sp ], col =  57, psym = 3
 ;;       oplot, [time], [ftmix.s2p], col = 147, psym = 3
 ;;       oplot, [time], [ftmix.s3p], col = 210, psym = 3
 ;;       oplot, [time], [ftmix.op ], col = 254, psym = 3
 ;;       oplot, [time], [ftmix.o2p], col = 100, psym = 3
 ;;       oplot, [time], [     dspr], col = 180, psym = 3
      endelse

      spr_old = spr
    endif ; mod 5

    endelse ; plot
    
    ; **IMPORTANT**
    ; this set of conditional statements decides when the model stops
    if (abs(1. - P_IO) le 0.2 && dspr le 0.1) && time gt 200. then break

ENDFOR



; energy information
if keyword_set(printiter) then begin
   ;printEF,  n, T, h, r_ind, r_dep, lat_dist, nu, ftint
     print, '$$--------------------------------'
     print, '$$ IN-CODE ENERGY BUDGET'
     print, '$$--------------------------------'
     print, '$$ ionized S............', energy.s_ion
     print, '$$ ionized O............', energy.o_ion
     print, '$$ charge exchange S....', energy.s_cx
     print, '$$ charge exchange O....', energy.o_cx
     print, '$$ equil with ehot......', energy.eh_eq + energy.tot_eq > 0.
     print, '$$ total in.............', energy.P_in + energy.tot_eq
     print, '$$ puv..................', energy.Puv
     print, '$$ fast/ena.............', energy.pfast - energy.tot_eq
     print, '$$ transport............', energy.ptrans + energy.ptrans_eh > 0.
     print, '$$ total out............', energy.P_out - energy.tot_eq 
     print, '$$ in/out...............', (energy.P_in + energy.tot_eq > 0.)/(energy.P_out - energy.tot_eq > 0.)
     print
endif

IF n_elements(badtime) GT 1 THEN BEGIN
   badtime = badtime[1:*]
   print, 'CM3_LATAVG_ONEBOX: n.fh went negative!'
   stop
ENDIF

if not keyword_set(nosave) then begin
save,file='cm3_latavg_onebox'+filename+'.sav', n_out, t_out, h_out, f_out, $
     ;tm, $
     time, $
     Teh0, source_out, otos, fehot_const, trans_out, l3
end

IF keyword_set(plot) THEN BEGIN
   wdelete, comp_pixmap_id
ENDIF

getsysvariables,/restore
close, 21
;if keyword_set(printiter) then print,'Runtime of: '+strtrim(systime(1)-runtime,2)+' seconds'

end
