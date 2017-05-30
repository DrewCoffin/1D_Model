Pro runmodel

  ; convert from total # neutrals in flux tube to cm^-3 s^-1
  ; = sqrt(pi)*.5*RJcm
  defsysv, '!rjkm', 71492. ; Jovian Radius in km
  tot2cms = sqrt(!dpi)*0.5*!rjkm*1e5

  ; create structure to hold all input parameters
  input = create_struct( 'source',7.2e-4*tot2cms, $ ; [# cm^-3 s^-1]
                         'o_to_s',     1.7    , $ ; ratio
                    'fehot_const',     0.003  , $ ; fraction of hot e
                      'transport',    50.0     , $ ; [days]
                          'tehot',    40.0     , $ ; [ev]
  ; note: the following perturbation parameters aren't used by widget
                      'lag_const',     0.0     , $ ; [km/s] subcorotational plasma speed
                    'neutral_amp',     0.0     , $ ; neutral perturb, amplitude
                     'neutral_t0',     0.0     , $ ; neutral perturb, initial time [DOY 2000]
                  'neutral_width',     0.0     , $ ; width of neutral perturb [days]
                       'hote_amp',     0.0     , $ ; hot e fraction perturb
                        'hote_t0',     0.0     , $ ; hot e fraction perturb, inital time [DOY 2000]
                     'hote_width',     0.0       ) ; width of hot e fraction perturb [days]
  
  ; structure containing energy values
  energy = {energy, s_ion: 0.0, s_cx: 0.0, o_ion: 0.0, o_cx: 0.0, eh_eq: 0.0, tot_eq: 0.0, P_in: 0.0, $
            Puv: 0.0, Pfast: 0.0, Ptrans: 0.0, Ptrans_eh: 0.0, P_out: 0.0}

  ; record start time
  start = systime(/julian)

  cm3_latavg_onebox, n_out, t_out, h_out, f_out, r_ind, r_dep, energy, $
                     /plot, /printiter, $
                     tm   = 0.01, $
                     ;runt = 500., $ latest dataset
                     runt = 500., $
                     /Nosave, $
                     ;filename = '_dim10_5-5', $
                     source = input.source, o_to_s = input.o_to_s, fehot_const = input.fehot_const, $
                     tehot = input.tehot, transport = input.transport, lag_const = input.lag_const, $
                     neutral_amp = input.neutral_amp, neutral_t0 = input.neutral_t0, neutral_width = input.neutral_width, $
                     hote_amp = input.hote_amp, hote_t0 = input.hote_t0, hote_width = input.hote_width, $
                     lag_amp = 0, lag_phase = 0., $
                     protons = 0.1, n_height = 0.5*!rjkm
  stop = systime(/julian)
  mrun = (stop - start)*84600d

  ;; reminder of contents of n, h, and t structures:
  ;; n = {dens_lon,sp:0.0d,s2p:0.0d, s3p: 0.0d, s4p: 0.0d, op: 0.0d, $
  ;;      o2p: 0.0d, el: 0.0d, elh: 0.0d, s: 0.0d, o: 0.0d, fc:0d, fh:0d, protons:0d}
  ;; T = {temp_lon,sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p: 0.0d,  $
  ;;      op: 0.0d, o2p: 0.0d, el: 0.0d, elh: 0.0d, pu_s: 0d, pu_o: 0d}
  ;; h = {scale_heights_lon, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
  ;;      o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}

  n = n_out
  h = h_out
  t = t_out


  ; the following requires chianti to be set up...

  ;; Here are the lines to actually call "emiss_calc" and calculate the emission coefficients. 
  ;; Sulfur ions
  ;s2em = emiss_calc(16, 2, temp = log(t.el/!kev), dens = log(n.el), /no_de, radt = 1., /quiet);, _extra = extra)
  ;s3em = emiss_calc(16, 3, temp = log(t.el/!kev), dens = log(n.el), /no_de, radt = 1., /quiet);, _extra = extra)
  ;s4em = emiss_calc(16, 4, temp = log(t.el/!kev), dens = log(n.el), /no_de, radt = 1., /quiet);, _extra = extra)

  ;; Oxygen ions
  ;o2em = emiss_calc(8, 2, temp = log(t.el/!kev), dens = log(n.el), /no_de, radt = 1., /quiet);, _extra = extra)
  ;o3em = emiss_calc(8, 3, temp = log(t.el/!kev), dens = log(n.el), /no_de, radt = 1., /quiet);, _extra = extra)
  
  ;h = 6.626068d-34 ; [m2 kg / s]
  ;c = 2.998d8      ; [m / s] 
  ;A = 1d-10        ; [m]
  ;cc = 1d31        ; [cm3] mass of torus

  ;; P = N*hc /lambda
  ;Ps2em = (reform(s2em.em)*n.sp )*(h*c) / (s2em.lambda*A)
  ;Ps3em = (reform(s3em.em)*n.s2p)*(h*c) / (s3em.lambda*A)
  ;Ps4em = (reform(s4em.em)*n.s3p)*(h*c) / (s4em.lambda*A)
  ;Po2em = (reform(o2em.em)*n.op )*(h*c) / (o2em.lambda*A)
  ;Po3em = (reform(o3em.em)*n.o2p)*(h*c) / (o3em.lambda*A)

  ;; only use EUV range of wavelength
  ;mineuv =  500.
  ;maxeuv = 1200.
  ;avgeuv = (mineuv + maxeuv)/2.

  ;Ps2em = Ps2em[where(abs(s2em.lambda - avgeuv) le avgeuv - mineuv)]
  ;Ps3em = Ps3em[where(abs(s3em.lambda - avgeuv) le avgeuv - mineuv)]
  ;Ps4em = Ps4em[where(abs(s4em.lambda - avgeuv) le avgeuv - mineuv)]
  ;Po2em = Po2em[where(abs(o2em.lambda - avgeuv) le avgeuv - mineuv)]
  ;Po3em = Po3em[where(abs(o3em.lambda - avgeuv) le avgeuv - mineuv)]

  ;Ptot = total([Ps2em, Ps3em, Ps4em, Po2em, Po3em])*cc
  ;print
  ;print, '%  Total Power in UV......', Ptot

  print
  print,'% Runtime of: ' + string(mrun)  + ' seconds'
  print
  
  End
