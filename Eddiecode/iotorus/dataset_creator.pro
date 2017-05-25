;+
;  PROCEDURE
;   dataset_creator
;
;  Using Andrew Steffl's io plasma torus chemistry model,
;  vary four input parameters in order to create a block of data to be
;  displayed by a slider widget.
;
;  VARIED QUANTITIES
;  - source
;  - o_to_s
;  - transport (radial transport timescale)
;  - fehot_const
;
;-
pro dataset_creator

  ; start time
  start = string(systime(/julian), format = '(c())')
  print, 'start time:' + start

  ; name file after the start time of the dataset, i.e., month-day_hour
  startfilestr = string(systime(/julian), format = '(C(CMOI02,"-",CDI02,"_",CHI2.2))')

  ; variation per parameter
  ;dim = 10

  ; convert from total # neutrals in flux tube to cm^-3 s^-1
  ; = sqrt(pi)*.5*RJcm
  tot2cms = (1d-4)*sqrt(!dpi)*0.5*!rjkm*1d5

  ; determine parameter span for neutral source rate [10^-4 cm^-3 s^-1]
  incsrc =  2.0 ; step size
  minsrc = 15.2 ; min value
  maxsrc = 21.2 ; max value
  dimsrc = (maxsrc - minsrc)/incsrc ; dimension of array
  source = (maxsrc - minsrc)*findgen(dimsrc)/(dimsrc - 1.) + minsrc
  source = source*tot2cms ; convert to units used by code [Sn-tot]

  ; determine parameter span for o to s ratio
  incots = 0.1 ; step size
  minots = 1.5  ; min value
  maxots = 2.0  ; max value
  dimots = (maxots - minots)/incots ; dimension of array
  o_to_s = (maxots - minots)*findgen(dimots)/(dimots - 1.) + minots

  ; determine parameter span for fraction of hot electrons
  incfeh = 0.05 ; step size
  minfeh = 0.2   ; min value
  maxfeh = 0.3   ; max value
  dimfeh = (maxfeh - minfeh)/incfeh ; dimension of array
  fehotc = (maxfeh - minfeh)*findgen(dimfeh)/(dimfeh - 1.) + minfeh
  fehotc = fehotc/100. ; convert to units used by code

  ; determine parameter span for radial transport time [days]
  inctsp = 2.5  ; step size
  mintsp = 24.5 ; min value
  maxtsp = 32.0 ; max value
  dimtsp = (maxtsp - mintsp)/inctsp ; dimension of array
  trnspt = (maxtsp - mintsp)*findgen(dimtsp)/(dimtsp - 1.) + mintsp

  ; input data
  datainput = create_struct(  'source', source, $
                              'o_to_s', o_to_s, $
                         'fehot_const', fehotc, $
                           'transport', trnspt, $
                               'tehot',     80. )


  ; single value n, T, and h structures
  n = {dens_lon, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p: 0.0d, op: 0.0d, $
       o2p: 0.0d, el: 0.0d, elh: 0.0d, s: 0.0d, o: 0.0d, fc: 0d, fh: 0d, protons: 0d}
  T = {temp_lon, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p: 0.0d,  $
       op: 0.0d, o2p: 0.0d, el: 0.0d, elh: 0.0d, pu_s: 0d, pu_o: 0d}
  h = {scale_heights_lon, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
       o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}

  ntags = n_tags(n)
  htags = n_tags(h)
  ttags = n_tags(T)

  ; structure containing energy values
  energy_budget = {energy, s_ion: 0.0, s_cx: 0.0, o_ion: 0.0, o_cx: 0.0, eh_eq: 0.0, tot_eq: 0.0, P_in: 0.0, $
                   Puv: 0.0, Pfast: 0.0, Ptrans: 0.0, Ptrans_eh: 0.0, P_out: 0.0}

  ; energy value structures
  ;; Ec = {energy_conserve, s_ion: 0.0, s_ion_h: 0.0, s_cx: 0.0, o_ion: 0.0, o_ion_h: 0.0, $
  ;;       o_cx: 0.0, s_tot_in: 0.0, o_tot_in: 0.0, P_pu: 0.0, eh_eq: 0.0, P_in: 0.0, $
  ;;       Puv: 0.0, Pfast: 0.0, Pfast_O: 0.0, Pfast_S: 0.0, Ptrans: 0.0, $
  ;;       Ptrans_O: 0.0, Ptrans_S: 0.0, Ptrans_e: 0.0, Ptrans_eh: 0.0, $
  ;;       P_out: 0.0, in_out: 0.0, ion_e_eq: 0.0}
  
  ;; nl = {neutral_loss, sei: 0.0, seih: 0.0, oei: 0.0, oeih: 0.0, $
  ;;       scx: 0.0, ocx: 0.0, k1: 0.0, k2: 0.0, k3: 0.0, k4: 0.0, k5: 0.0, $
  ;;       k6: 0.0, k7: 0.0, k8: 0.0, k9: 0.0, k10:0.0, k11: 0.0, k12: 0.0, $
  ;;       k14: 0.0, s_tot: 0.0, o_tot: 0.0, scx_tot: 0.0, ocx_tot: 0.0, $
  ;;       sion_tot: 0.0, oion_tot: 0.0, fast_O_k5: 0.0, fast_O_k7: 0.0, $
  ;;       fast_S_k1: 0.0, fast_S_k3: 0.0, fast_O_k9: 0.0, fast_S_k8: 0.0}
  
  energy = replicate(energy_budget, dimsrc, dimots, dimfeh, dimtsp)

  ; output for one run
  output = create_struct( 'sp', 0., $
                         's2p', 0., $
                         's3p', 0., $
                         's4p', 0., $
                          'op', 0., $
                         'o2p', 0., $
                         'Tsp', 0., $
                        'Ts2p', 0., $
                        'Ts3p', 0., $
                        'Ts4p', 0., $
                         'Top', 0., $
                        'To2p', 0., $
                         'Tel', 0., $
                           'n',  n  )

  ; entire data output
  dataoutput = replicate(output, dimsrc, dimots, dimfeh, dimtsp)
  
  ;sqrtpi = sqrt(!pi)

  ; loop through every possibility - dim^4 iterations
  for aa = 0, dimsrc-1 do begin ; source
    for bb = 0, dimots-1 do begin ; o_to_s
      for cc = 0, dimfeh-1 do begin ; fehot_const
        for dd = 0, dimtsp-1 do begin ; transport
          print,aa
          ; note: code is designed to run until equilibrium, or until 'runt'
          cm3_latavg_onebox, n_out, t_out, h_out, f_out, r_ind, r_dep, energy_budget, $
                             runt = 500., $
                             ;runt = 2.*datainput.transport[dd]
                             /nosave, $
                             source      = datainput.source[aa], $
                             o_to_s      = datainput.o_to_s[bb], $
                             fehot_const = datainput.fehot_const[cc], $
                             transport   = datainput.transport[dd], $
                             tehot       = datainput.tehot, $
                             lag_const = 0., lag_amp = 0, lag_phase = 0., $
                             neutral_amp = 0., hote_amp = 0., $
                             tm = .01, protons = 0.1, n_height = .5*!rjkm
          
          ; store energy
          energy[aa,bb,cc,dd] = energy_budget

          ; store ft mix ratios
          dataoutput[aa,bb,cc,dd].sp  = f_out.sp
          dataoutput[aa,bb,cc,dd].s2p = f_out.s2p
          dataoutput[aa,bb,cc,dd].s3p = f_out.s3p
          dataoutput[aa,bb,cc,dd].s4p = f_out.s4p
          dataoutput[aa,bb,cc,dd].op  = f_out.op
          dataoutput[aa,bb,cc,dd].o2p = f_out.o2p

          ; store electron data
          dataoutput[aa,bb,cc,dd].Tel  = t_out.el
          dataoutput[aa,bb,cc,dd].n    = n_out

          ; store temperatures
          dataoutput[aa,bb,cc,dd].Tsp  = t_out.sp
          dataoutput[aa,bb,cc,dd].Ts2p = t_out.s2p
          dataoutput[aa,bb,cc,dd].Ts3p = t_out.s3p
          dataoutput[aa,bb,cc,dd].Ts4p = t_out.s4p
          dataoutput[aa,bb,cc,dd].Top  = t_out.op
          dataoutput[aa,bb,cc,dd].To2p = t_out.o2p

        endfor ; dd
      endfor ; cc
    endfor ; bb
  endfor ; aa

  ; save structures to IDL save file in the widget directory
  save, datainput, dataoutput, energy, filename = 'widget/torusmodel-data_' + startfilestr + '.sav'

  ; end time
  stop = string(systime(/julian), format = '(c())')
  print, '  end time:' + stop

end

;+
;  print info of supplied data set, if no data set supplied, find
;  latest hypercube file
;-
pro dataset_info, filename = filename

  ; if file not specified, retrieve the newest hypercube file
  if not keyword_set(filename) then filename = (file_search('torusmodel-data*.sav'))[-1]

  ; convert from total # neutrals in flux tube to cm^-3 s^-1
  ; = sqrt(pi)*.5*RJcm
  tot2cms = (1d-4)*sqrt(!dpi)*0.5*!rjkm*1d5

  print
  print, '***************************************'
  print, 'PROCEDURE DATASET_INFO'
  print, 'FILE: ' + filename
  print, '***************************************'
  restore, filename, /verbose
  print
  
  tag_names = tag_names(datainput)
  n_tags = n_tags(datainput)

  ; **should be in above procedure**
  input = {units:strarr(n_tags), conv:dblarr(n_tags)}
  input.units = ['10^-4 cm^-3 s^-1', 'NONE', '%', 'DAYS', 'eV']
  input.conv = [tot2cms, 1, 0.01, 1, 1]

  for ii = 0, n_tags-1 do begin

    min = min(datainput.(ii), max = max)
    n = n_elements(datainput.(ii))

    print, tag_names[ii]
    print, 'UNITS = ' + input.units[ii]
    print, '  MIN = ' + string(          min/input.conv[ii], format = '(F7.3)')
    print, '  MAX = ' + string(          max/input.conv[ii], format = '(F7.3)')
    print, '    N = ' + string(                           n, format =   '(I5)')
    print, ' STEP = ' + string(((max-min)/input.conv[ii])/n, format = '(F7.3)')
    print
    
  endfor

end
