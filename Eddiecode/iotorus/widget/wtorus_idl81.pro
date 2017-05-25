;+
;  PROCEDURE
;    wtorus_idl81
;    (widget to cubic centimeter io torus model)
;
;  PURPOSE
;    Create a graphical user interface to explore the chemical makeup
;    and energy flow of the Io plasma torus, using a cubic centimeter
;    chemistry model of the torus. The basic plan is to vary four of
;    the model variables: neutral source rate, oxygen to sulfur
;    ratio, fraction of hot electrons, and the radial transport time
;    while keeping the temperature of hot electrons constant. With the
;    chosen input parameters, the user can then see the output mixing
;    ratios of the ion species, energy flow, and a simulated spectra.
;    For more information about the model, see:
;     - Delamere, Steffl, and Bagenal (2004), Modeling temporal
;       variability of plasma conditions in the Io torus during the Cassini era
;
;  SETUP AND CALL SEQUENCE
;    Note: this widget is meant to run on IDL 8.1 and later versions.
;    A handling procedure, wtorus.pro, compiles and runs this
;    widget. Wtorus also sets up the Chianti database, which is
;    required to simulate spectra.
;    Note: Before running, be sure to edit wtorus.pro, setting the 
;    chiantipath variable to the main folder of the chianti database.
;    If you do not have Chianti installed on your system, visit
;    http://www.chiantidatabase.org/chianti_download.html and download
;    the two stand-alone files. Create the file structure as follows:
;    ~/chianti/idl/
;    ~/chianti/dbase/
;    Unpack the zipped files in their respective folders and Chianti
;    is ready to go. The handling procedure, given the appropriate 
;    path to the folder 'chianti', will set up Chianti and add its
;    library to IDL's !path variable.
;      Note: This procedure does not need to be compiled, wtorus.pro
;            will compile this procedure using IDL's resolve_routine.
;    With Chianti installed, run the widget by calling the handling
;    procedure as follows:
;    IDL> wtorus
;
;  MOD HISTORY
;    Created by Adam Shinn May 2012
;-

function energytextbox, aa, bb, cc, dd, energystruct = energystruct, outputstruct = outputstruct
; given subscripits and energy structure, return string array for "power breakdown"

  ; obtain power data from energy structure
  PSion   = energystruct[aa,bb,cc,dd].s_ion
  POion   = energystruct[aa,bb,cc,dd].o_ion
  PSexch  = energystruct[aa,bb,cc,dd].s_cx
  POexch  = energystruct[aa,bb,cc,dd].o_cx
  Phot    = energystruct[aa,bb,cc,dd].eh_eq
  Puv     = energystruct[aa,bb,cc,dd].Puv
  Pneut   = energystruct[aa,bb,cc,dd].Pfast
  Ptrans  = energystruct[aa,bb,cc,dd].Ptrans + energystruct[aa,bb,cc,dd].Ptrans_eh
  Peq     = energystruct[aa,bb,cc,dd].tot_eq

  ; total power in and power out
  PinTOT = energystruct[aa,bb,cc,dd].P_in
  PouTOT = energystruct[aa,bb,cc,dd].P_out
  
  ; equilibration
  Phot   += Peq > 0.
  PinTOT += Peq > 0.
  Pneut  -= Peq > 0.
  PouTOT -= Peq > 0.

  ; mixing ratios
  sp  = outputstruct[aa,bb,cc,dd].sp
  s2p = outputstruct[aa,bb,cc,dd].s2p
  s3p = outputstruct[aa,bb,cc,dd].s3p
  op  = outputstruct[aa,bb,cc,dd].op
  o2p = outputstruct[aa,bb,cc,dd].o2p

  format = '(F7.4)'

  text = ['Power/CC [eV s^-1]/cm^-3       ' , $
          '                               ' , $
          '   Input                       ' , $
          '    S ionization...............' + string(PSion , format = format), $
          '    O ionization...............' + string(POion , format = format), $
          '    S charge exchange..........' + string(PSexch, format = format), $
          '    O charge exchange..........' + string(POexch, format = format), $
          '    Hot electrons..............' + string(Phot  , format = format), $
          '    Total Input (100%).........' + string(PinTOT, format = format), $
          '                               ' , $
          '   Output                      ' , $
          '    UV.........................' + string(Puv   , format = format), $
          '    Fast neutrals..............' + string(Pneut , format = format), $
          '    Transport..................' + string(Ptrans, format = format), $
          '    Total Output (100%)........' + string(PouTOT, format = format), $
          '                               ' , $
          '   Power In / Power Out........' + string(PinTOT/PouTOT, format = format), $
          '                               ' , $
          'Mixing Ratios [unitless]       ' , $
          '                               ' , $
          '    S+ ........................' + string(sp , format = format), $
          '    S++ .......................' + string(s2p, format = format), $
          '    S+++ ......................' + string(s3p, format = format), $
          '    O+ ........................' + string(op , format = format), $
          '    O++ .......................' + string(o2p, format = format), $
          '                               ' ]

  return, text
end

;+
; given the density structure, compute the total electron density
;-
function nel_tot, n

  nel = (n.sp + 2.*n.s2p + 3.*n.s3p + $
         n.op + 2.*n.o2p)/(1. - n.protons)

return, nel
end

;+
; given the density structure, compute the total electron density
;-
function avg_ionch, n

  tot_ion = n.sp + n.s2p + n.s3p + $
            n.op + n.o2p

  avgion = (n.sp + 2.*n.s2p + 3.*n.s3p + $
            n.op + 2.*n.o2p)/tot_ion

return, avgion
end

;+
; function used for obtaining spectra using chianti's emiss_calc
;-
function spectra_get, aa, bb, cc, dd, outputstruct = outputstruct, $
                 min = min, max = max, fwhm = fwhm, margin = margin, $
                 xwav = xwav

  ; sample spectra
  if keyword_set(max) then maxwav = max else maxwav = 1800d
  if keyword_set(min) then minwav = min else minwav = 1100d
  if keyword_set(fwhm) then fwhm = fwhm else fwhm = 2.
  if keyword_set(margin) then margin = margin else margin = .15

  ; emission line calculation: sulfur ions
  s2em = emiss_calc(16, 2, temp = alog10(outputstruct[aa,bb,cc,dd].tel/!kev), dens = alog10(outputstruct[aa,bb,cc,dd].n.el), /no_de, radt = 1., /quiet)
  s3em = emiss_calc(16, 3, temp = alog10(outputstruct[aa,bb,cc,dd].tel/!kev), dens = alog10(outputstruct[aa,bb,cc,dd].n.el), /no_de, radt = 1., /quiet)
  s4em = emiss_calc(16, 4, temp = alog10(outputstruct[aa,bb,cc,dd].tel/!kev), dens = alog10(outputstruct[aa,bb,cc,dd].n.el), /no_de, radt = 1., /quiet)

  ; emission line calculation: oxygen ions
  o2em = emiss_calc(8, 2, temp = alog10(outputstruct[aa,bb,cc,dd].tel/!kev), dens = alog10(outputstruct[aa,bb,cc,dd].n.el), /no_de, radt = 1., /quiet)
  o3em = emiss_calc(8, 3, temp = alog10(outputstruct[aa,bb,cc,dd].tel/!kev), dens = alog10(outputstruct[aa,bb,cc,dd].n.el), /no_de, radt = 1., /quiet)

  ; units of em are currently [(photon/sec)/ion]/cc
  ; multiply by ion density to get => [photon/sec]/cc
  s2emiss = reform(s2em.em)*outputstruct[aa,bb,cc,dd].n.sp 
  s3emiss = reform(s3em.em)*outputstruct[aa,bb,cc,dd].n.s2p
  s4emiss = reform(s4em.em)*outputstruct[aa,bb,cc,dd].n.s3p
  o2emiss = reform(o2em.em)*outputstruct[aa,bb,cc,dd].n.op 
  o3emiss = reform(o3em.em)*outputstruct[aa,bb,cc,dd].n.o2p

  ; discrete
  xwavi = [s2em.lambda, s3em.lambda, s4em.lambda, o2em.lambda, o3em.lambda]
  wsort = sort(xwavi) ; sort by wavelength
  xwavi = xwavi[wsort]
  
  yptsi = [s2emiss, s3emiss, s4emiss, o2emiss, o3emiss]
  yptsi = yptsi[wsort]

  yname = [s2em.ion_name, s3em.ion_name, s4em.ion_name, o2em.ion_name, o3em.ion_name]
  yname = yname[wsort]

  ; cut down arrays to necessary range
  avgwav = (minwav + maxwav)/2.
  wrange = where(abs(xwavi - avgwav) le avgwav - minwav)
  xwavi = xwavi[wrange]
  yptsi = yptsi[wrange]
  yname = yname[wrange]

  ; non discrete
  xwav = (maxwav - minwav)*dindgen(1000)/999d + minwav
  ypts = dblarr(1000)

  ; relation between sigma and fwhm:
  ; => fwhm = 2*sigma*sqrt(2*ln(2))
  sigma = fwhm/(2d*sqrt(2d*alog(2d)))
  sigma2 = sigma*sigma
  a = 1d/(sigma*sqrt(2d*!dpi))
  for ii = 0, 999 do begin
    ; treat every emission as a normalized gaussian, add up all gaussians to create simulated spectra
    ; f(x) = (1/sigma*sqrt(2pi)) * exp( -.5(x - x0)^2 / sigma^2 )
    y = a*yptsi[0]*exp( (-.5*(xwav[ii]-xwavi[0])*(xwav[ii]-xwavi[0]))/sigma2 )
    for jj = 1, n_elements(xwavi)-1 do y += a*yptsi[jj]*exp( (-.5*(xwav[ii]-xwavi[jj])*(xwav[ii]-xwavi[jj]))/sigma2 )

    ypts[ii] = y
  endfor
  
  return, ypts
  
  ; note: xwav is also an output, placed as a keyword
  
end

;+
; event handler for the io torus widget
;-
pro wtorus_event, event
  ; use common block to carry data and window variables
  common drawblock, mixtempwinvalue, datainput, dataoutput, energy

  ; Use WIDGET_CONTROL to get the user value of any widget touched and put
  ; that value into 'eventval':
  widget_control, event.id, get_uvalue = eventval

  ; get the widget info structure
  widget_control, event.top, get_uvalue = widgetinfo;, /no_copy

  ; The movement of sliders is easily handled with a CASE statement.
  ; When the slider is moved, the value of 'eventval' becomes 'SLIDE':

  ; obtain values from input sliders
  widget_control, widgetinfo.slidersource   , get_value = aa
  widget_control, widgetinfo.slidero_to_s   , get_value = bb
  widget_control, widgetinfo.sliderfehot    , get_value = cc
  widget_control, widgetinfo.slidertransport, get_value = dd

  ; retreive widget window by uname
  mixtempwindow = widget_info(event.top, find_by_uname = 'drawmixtemp')

  ; case statement broken down into widget actions:
  case tag_names(event, /structure_name) of
    'WIDGET_SLIDER':begin
      widget_control, event.id, get_uvalue = event_UV
      case event_UV of
        'slide-input':begin

          ; update slider labels
          source = datainput.source[aa]/((1e-4)*sqrt(!dpi)*0.5*!rjkm*1e5)
          widget_control, widgetinfo.labelsource, set_value = string(source, format = '(F5.2)')
          widget_control, widgetinfo.labelo_to_s, set_value = string(datainput.o_to_s[bb], format = '(F5.2)')
          fehot = datainput.fehot_const[cc]*100.
          widget_control, widgetinfo.labelfehot, set_value = string(fehot, format = '(F5.2)')
          widget_control, widgetinfo.labeltransport, set_value = string(datainput.transport[dd], format = '(F5.2)')

          ; plot
          wtorus_plotdraw, aa, bb, cc, dd, outputstruct = dataoutput, energystruct = energy, window = mixtempwindow

          ; refresh electron total value
          nel_tot = nel_tot(dataoutput[aa,bb,cc,dd].n)
          widget_control, widgetinfo.edensval, set_value = string(nel_tot, format = '(F8.2)') + ' cm^-3'

          ; refresh average ion charge label
          avg_ion = avg_ionch(dataoutput[aa,bb,cc,dd].n)
          widget_control, widgetinfo.avgionlbl, set_value = '(average ion charge: ' + string(avg_ion, format = '(F4.2)') + ')'

          ; refresh cold electron value
          widget_control, widgetinfo.ectmpval, set_value = string(dataoutput[aa,bb,cc,dd].Tel, format = '(F7.2)') + ' eV'

          ; refresh power breakdown
          if widgetinfo.powerinfo ne 0 then widget_control, widgetinfo.powerinfo, set_value = energytextbox(aa,bb,cc,dd, energystruct = energy, outputstruct = dataoutput), bad_id = bad

        end
        'slide-gaussian':begin

          ; obtain values from spectra knobs
          widget_control, widgetinfo.minfield      , get_value = min
          widget_control, widgetinfo.maxfield      , get_value = max
          ;widget_control, widgetinfo.slidergaussian, get_value = fwhm

          widget_control, event.id, get_value = fwhm
          spectra_draw, aa, bb, cc, dd, outputstruct = dataoutput, window = mixtempwindow, min = min, max = max, fwhm = fwhm
        end
        else: ; do nothing
      endcase
    end ; widget_slider
    'WIDGET_BUTTON':begin
       widget_control, event.id, get_uvalue = event_UV
       case event_UV of
          'powerinfobutton':begin
             ; if it already exists, destroy the old box and create a new one
             if widgetinfo.powerinfo ne 0 then widget_control, widgetinfo.hidebox, /destroy, bad_id = bad

             ; create a new base with text widget inside
             popbase = widget_base(/column)
             powerinfo = widget_text(popbase, $
                                     ysize = 30, $
                                     xsize = 50, $
                                     value = energytextbox(0, 0, 0, 0, energystruct = energy, outputstruct = dataoutput), $
                                     uvalue = 'textpowerinfo')

             ; store text box identifiers in identifier structure
             widgetinfo.powerinfo = powerinfo
             widgetinfo.hidebox = popbase ; this identifier is used to destroy the base if needed

             ; save identifier structure
             widget_control, event.top, set_uvalue = widgetinfo, /no_copy

             ; realize the new widget base and hand off control to xmanager
             widget_control, popbase, /realize
             xmanager, 'popup', popbase, group_leader = event.top, /no_block
          end
          'plotspectra':begin

             ; obtain values from spectra knobs
             widget_control, widgetinfo.minfield      , get_value = min
             widget_control, widgetinfo.maxfield      , get_value = max
             widget_control, widgetinfo.slidergaussian, /sensitive, get_value = fwhm

             spectra_draw, aa, bb, cc, dd, outputstruct = dataoutput, window = mixtempwindow, min = min, max = max, fwhm = fwhm
          end
          'savebutton':begin
             
             ; obtain input values
             source = datainput.source[aa]/((1e-4)*sqrt(!dpi)*0.5*!rjkm*1e5)
             o_to_s = datainput.o_to_s[bb]
             fehot = datainput.fehot_const[cc]*100.
             transport = datainput.transport[dd]

             ; create structure to save input parameters this instance
             input = create_struct('source',    source, $
                                   'o_to_s',    o_to_s, $
                                    'fehot',     fehot, $
                                'transport',  transport )

             ; obtain energy information
             energytext = energytextbox(aa,bb,cc,dd, energystruct = energy, outputstruct = dataoutput)

             ; obtain power data from energy structure
             PSion   = energy[aa,bb,cc,dd].s_ion
             POion   = energy[aa,bb,cc,dd].o_ion
             PSexch  = energy[aa,bb,cc,dd].s_cx
             POexch  = energy[aa,bb,cc,dd].o_cx
             Phot    = energy[aa,bb,cc,dd].eh_eq
             Puv     = energy[aa,bb,cc,dd].Puv
             Pneut   = energy[aa,bb,cc,dd].Pfast
             Ptrans  = energy[aa,bb,cc,dd].Ptrans + energy[aa,bb,cc,dd].Ptrans_eh
             Peq     = energy[aa,bb,cc,dd].tot_eq

             ; total power in and power out
             PinTOT = energy[aa,bb,cc,dd].P_in
             PouTOT = energy[aa,bb,cc,dd].P_out
  
             ; equilibration
             Phot   += Peq > 0.
             PinTOT += Peq > 0.
             Pneut  -= Peq > 0.
             PouTOT -= Peq > 0.

             ; create a structure to save the energy budget at this instance
             energybudget = create_struct( 'PSion',  PSion, $
                                           'POion',  POion, $   
                                          'PSexch', PSexch, $
                                          'POexch', POexch, $  
                                            'Phot',   Phot, $    
                                             'Puv',    Puv, $     
                                           'Pneut',  Pneut, $   
                                          'Ptrans', Ptrans, $  
                                             'Peq',    Peq, $
                                          'PinTOT', PinTOT, $
                                          'PouTOT', PouTOT, $
                                           'units', '[eV s^-1]/cm^-3')

             ; mixing ratios
             sp  = dataoutput[aa,bb,cc,dd].sp
             s2p = dataoutput[aa,bb,cc,dd].s2p
             s3p = dataoutput[aa,bb,cc,dd].s3p
             op  = dataoutput[aa,bb,cc,dd].op
             o2p = dataoutput[aa,bb,cc,dd].o2p

             ; temperature [eV] for each neutral species
             Tsp  = dataoutput[aa,bb,cc,dd].Tsp
             Ts2p = dataoutput[aa,bb,cc,dd].Ts2p
             Ts3p = dataoutput[aa,bb,cc,dd].Ts3p
             Top  = dataoutput[aa,bb,cc,dd].Top
             To2p = dataoutput[aa,bb,cc,dd].To2p

             ; total equitorial electron density [cm^-3]
             nel_tot = nel_tot(dataoutput[aa,bb,cc,dd].n)

             ; average ion charge
             avg_ion = avg_ionch(dataoutput[aa,bb,cc,dd].n)

             ; cold electron value in eV
             Tel = dataoutput[aa,bb,cc,dd].Tel

             ; create info structure to save model output units
             ; note: flux tube integrated mixing ratios are unitless
             ounits = create_struct( 'sp',      '', $
                                    's2p',      '', $
                                    's3p',      '', $
                                     'op',      '', $
                                    'o2p',      '', $
                                    'Tsp',    'eV', $
                                   'Ts2p',    'eV', $
                                   'Ts3p',    'eV', $
                                    'Top',    'eV', $
                                   'To2p',    'eV', $
                                'nel_tot', 'cm^-3', $
                                'avg_ion',      '', $
                                    'Tel',    'eV'  )
          
             ; create structure to save model output at this instance
             output = create_struct( 'sp',      sp, $
                                    's2p',     s2p, $
                                    's3p',     s3p, $
                                     'op',      op, $
                                    'o2p',     o2p, $
                                    'Tsp',     Tsp, $
                                   'Ts2p',    Ts2p, $
                                   'Ts3p',    Ts3p, $
                                    'Top',     Top, $
                                   'To2p',    To2p, $
                                'nel_tot', nel_tot, $
                                'avg_ion', avg_ion, $
                                    'Tel',     Tel, $
                                  'units',  ounits  )


             ; obtain values from spectra knobs
             widget_control, widgetinfo.minfield      , get_value = min
             widget_control, widgetinfo.maxfield      , get_value = max
             widget_control, widgetinfo.slidergaussian, get_value = fwhm

             ; obtain spectra in units of [photons sec^-1 cm^-3]
             ypts = spectra_get(aa, bb, cc, dd, outputstruct = dataoutput, min = min, max = max, fwhm = fwhm, xwav = xwav)

             sunits = create_struct('yphotons', '[photons sec^-1 cm^-3]', $
                                      'minwav', 'Angstrom', $
                                      'maxwav', 'Angstrom', $
                                        'fwhm', 'Angstrom', $
                                        'xwav', 'Angstrom'  )

             ; create structure to save chianti spectra at this instance
             spectra = create_struct('yphotons', ypts, $
                                       'minwav',  min, $
                                       'maxwav',  max, $
                                         'fwhm', fwhm, $
                                         'xwav', xwav, $
                                        'units', sunits)

             ; name of save file
             savfile = 'wtorus_' + string(aa, bb, cc, dd, fwhm, format = "(I02,'-',I02,'-',I02,'-',I02,'_fw',I02)") + '.sav'

             ; save structures into IDL save file
             save, input, energybudget, output, spectra, filename = savfile

             print, '# ----INSTANCE SAVED----'
             print, '# filename: ' + savfile
             print, '#'
             print, '# STRUCTURES:'
             print, '# input, energybudget, output, spectra'
             print, '#'
             print, '# SPECTRA:'
             print, '# fwhm:  ' + string(fwhm, format = '(I2)')
             print, '# minwav:' + string( min, format = '(I4)')
             print, '# maxwav:' + string( max, format = '(I4)')
             print, '# ----------------------'

          end
       endcase
    end ; widget_button
    'CW_FIELD':begin

       ; obtain values from spectra knobs
       widget_control, widgetinfo.minfield      , get_value = min
       widget_control, widgetinfo.maxfield      , get_value = max
       widget_control, widgetinfo.slidergaussian, get_value = fwhm

       spectra_draw, aa, bb, cc, dd, outputstruct = dataoutput, window = mixtempwindow, min = min, max = max, fwhm = fwhm
    end
    else: ; do nothing
  endcase

end ; event handler

pro wtorus_plotdraw, aa, bb, cc, dd, outputstruct = outputstruct, energystruct = energystruct, $
                     window = window, init = init
  ; retrieve object window
  widget_control, window, get_value = windowobj

  ; select window object
  windowobj.select

  ; retrieve plot by name
  mix   = windowobj['mixplot']
  tmp   = windowobj['templot']
  yaxis = windowobj['tempaxis']
  pwrin = windowobj['pwrinplot']
  pwrou = windowobj['pwrouplot']

  ; erase spectra plot
  spectra = windowobj['spectraplot']
  if obj_valid(spectra) then spectra.setdata, 0

  ; check to make sure plot objects are present
  if keyword_set(init) then init = 1 else if obj_valid(mix) || obj_valid(tmp) then init = 0 else init = 1
  
  ; mixing ratios for each neutral species
  sp  = outputstruct[aa,bb,cc,dd].sp
  s2p = outputstruct[aa,bb,cc,dd].s2p
  s3p = outputstruct[aa,bb,cc,dd].s3p
  s4p = outputstruct[aa,bb,cc,dd].s4p
  op  = outputstruct[aa,bb,cc,dd].op
  o2p = outputstruct[aa,bb,cc,dd].o2p

  ; temperature [eV] for each neutral species
  Tsp  = outputstruct[aa,bb,cc,dd].Tsp
  Ts2p = outputstruct[aa,bb,cc,dd].Ts2p
  Ts3p = outputstruct[aa,bb,cc,dd].Ts3p
  Ts4p = outputstruct[aa,bb,cc,dd].Ts4p
  Top  = outputstruct[aa,bb,cc,dd].Top
  To2p = outputstruct[aa,bb,cc,dd].To2p

  ; vector of bar graph quantities
  ymix  = [ sp,  s2p,  s3p,  op,  o2p]
  ytemp = [Tsp, Ts2p, Ts3p, Top, To2p]

  ; obtain power data from energy structure
  PSion   = energystruct[aa,bb,cc,dd].s_ion
  POion   = energystruct[aa,bb,cc,dd].o_ion
  PSexch  = energystruct[aa,bb,cc,dd].s_cx
  POexch  = energystruct[aa,bb,cc,dd].o_cx
  Phot    = energystruct[aa,bb,cc,dd].eh_eq
  Puv     = energystruct[aa,bb,cc,dd].Puv
  Pneut   = energystruct[aa,bb,cc,dd].Pfast
  Ptrans  = energystruct[aa,bb,cc,dd].Ptrans + energystruct[aa,bb,cc,dd].Ptrans_eh
  Peq     = energystruct[aa,bb,cc,dd].tot_eq

  ; total power in and power out
  PinTOT = energystruct[aa,bb,cc,dd].P_in
  PouTOT = energystruct[aa,bb,cc,dd].P_out
  
  ; equilibration
  Phot   += Peq > 0.
  PinTOT += Peq > 0.
  Pneut  -= Peq > 0.
  PouTOT -= Peq > 0.

  ; turn power quantities into percentages
  ;PinTOT = PSion + POion + PSexch + POexch + Phot
  ;PouTOT = PUV + Pneut + Ptrans
  Pin = 100.*[PSion, POion, PSexch, POexch, Phot]/PinTOT
  Pou = 100.*[PUV, Pneut, Ptrans]/PouTOT

  ; plot temperature on same scale as mixing ratio
  ; note: a proper temperature scale will replace the right-side y-axis
  maxtemp = 500. ; eV
  maxmix = 1.
  ytemp = ytemp*(maxmix/maxtemp)

  ; initialize plots
  if init then begin
    margin = .15

    ; color code mixing ratios and temperature
    mixcolor  = 'cyan'
    tempcolor = 'orange'

    ; vector of bar graph quantities
    names = ['', 's+', 's++', 's+++', 'o+', 'o++', '']

    ; mixing ratio plot
    mix = barplot(ymix, nbars = 2, index = 0, /current, $
                  layout = [1, 3, 1], $
                  fill_color =  mixcolor, $
                  axis_style = 1, $
                  font_size = 12, $
                  /font_style, $
                  yrange = [0, maxmix], $
                  yticklen = 0.08, $
                  yminor = 3, $
                  xtickname = names, xminor = 0, $
                  margin = margin, $
                  title = 'ion species mix ratios and temperatures', $
                  ytitle = 'flux tube integrated mix ratios', $
                  name = 'mixplot')

    ; add temperature bars
    tmp = barplot(ytemp, nbars = 2, index = 1, fill_color = tempcolor, /overplot, $
                  axis_style = 1, $
                  name = 'templot')

    ; insert another y axis for temperature
    tmpaxis = axis('y', location = [n_elements(ytemp), 0], $
                   tickdir = 1, $
                   textpos = 1, $
                   tickfont_size = 12, $
                   /tickfont_style, $
                   title = 'Temp [eV]', $
                   tickvalues = (indgen(6)/5.)*maxmix, $
                   tickname = string(maxtemp*indgen(6)/5., format = '(I3)'), $
                   color = tempcolor, $
                   name = 'tempaxis')

    ; colors of I/O power bar plots
    inputcolor  = 'lavender'
    outputcolor = 'light_pink'

    ; power input plot
    pwrin = barplot(Pin, /current, $
                    fill_color =  inputcolor, $
                    axis_style = 1, $
                    font_size = 12, $
                    /font_style, $
                    yrange = [0, 100], $
                    xtickname = ['', 'P-S!Lion!N', 'P-O!Lion!N', 'P-S!LCHEX!N', 'P-O!LCHEX!N', 'P!Lhot E!N', ''], $
                    xmajor = 7, $
                    xminor = 0, $
                    margin = [.090, .375, .500, .375], $
                    ytitle = 'power input [%]', $
                    name = 'pwrinplot')

    ; power output plot
    pwrou = barplot(Pou, /current, $
                    fill_color =  outputcolor, $
                    axis_style = 1, $
                    font_size = 12, $
                    /font_style, $
                    yrange = [0, 100], $
                    xtickname = ['', 'P!LUV!N', 'P!LENA!N', 'P!Ltransport!N', ''], $
                    xmajor = 5, $
                    xminor = 0, $
                    margin = [.500, .375, .090, .375], $
                    width = .56, $
                    ytitle = 'power output [%]', $
                    name = 'pwrouplot')

    ; move power output plot to the right
    pwrou.translate, 50, 0

  endif else begin ; endif initialize, else use setdata method
    ; note: use 'refresh' object method for a smooth transition between slider choices

    ; disable refresh
    mix.refresh, /disable
    tmp.refresh, /disable
    pwrin.refresh, /disable
    pwrou.refresh, /disable

    ; set new data
    mix.setdata, ymix
    tmp.setdata, ytemp
    pwrin.setdata, Pin
    pwrou.setdata, Pou

    ; refresh object
    mix.refresh
    tmp.refresh
    pwrin.refresh
    pwrou.refresh
  endelse
end

pro spectra_draw, aa, bb, cc, dd, outputstruct = outputstruct, $
                  min = min, max = max, fwhm = fwhm, init = init, window = window, margin = margin

  ; sample spectra
  if keyword_set(max) then maxwav = max else maxwav = 1800d
  if keyword_set(min) then minwav = min else minwav = 1100d

  ; retrieve object window
  widget_control, window, get_value = windowobj

  ; select window object
  windowobj.select

  ; retrieve plot by name
  spectra = windowobj['spectraplot']

  if obj_valid(spectra) then init = 0 else init = 1

  ; obtain spectra
  ypts = spectra_get(aa, bb, cc, dd, outputstruct = outputstruct, min = min, max = max, fwhm = fwhm, xwav = xwav)

  max = max(ypts)
  if init then begin
    spectra = plot(xwav, ypts, /current, $
                   layout = [1, 3, 3], $
                   axis_style = 1, $
                   font_size = 12, $
                   /font_style, $
                   yrange = [-.05*max, 1.1*max], $
                   xrange = [minwav, maxwav], $
                   margin = [.15, .15, .05, .15], $
                   ytitle = '[photons sec!U-1!N cm!U-3!N]', $
                   xtitle = 'Angstrom', $
                   name = 'spectraplot')
    ; move spectra plot up a bit
    spectra.translate, 0, 10
  endif else begin
    spectra.refresh, /disable
    spectra.setdata, xwav, ypts
    spectra.xrange = [minwav, maxwav]
    spectra.yrange = [-.05*max, 1.1*max]
    spectra.refresh
  endelse

end

pro wtorus_idl81, group = group
  ; main procedure

  ; attempt to size widget appropriately depending on screen size
  device, get_screen_size = screen_size
  ;screen_size = [1366, 768]
  if screen_size[0] lt 1400 then begin
     ; if less than default size, set total widget size
     ; to the total screen size minus 30 pixels for padding
     x_scroll = 1040
     y_scroll = 700
     xsize = 1050
     ysize =  710
     scroll = 1
  endif else begin
     ; default size of widget
     x_scroll = 1155
     y_scroll = 875
     xsize = 1152
     ysize = 850
     scroll = 0
  endelse

  ; top-level base widget to hold torus widget
  base = widget_base(title = 'Io Plasma Torus CC Model: Widget Interface', $
                     x_scroll_size = x_scroll, $
                     y_scroll_size = y_scroll, $
                     xsize = xsize, $
                     ysize = ysize, $
                     scroll = scroll, $
                     /column)

  ; in/out base (buttons on left, plots on right)
  iobase = widget_base(base, /row)

  ; create a base for the sliders
  ; create a base for the plot(s)
  inputbase = widget_base(iobase, /column)
  plotbase  = widget_base(iobase, /column)

  ; if file not specified, retrieve the newest model data file
  filename = (file_search('torusmodel-data*.sav'))[-1]

  ; torus model data save file contains 'datainput', 'dataoutput' and 'energy' structures
  restore, filename, /verbose

  ; use common block to carry data and window variables
  common drawblock, mixtempwinvalue, datainput, dataoutput, energy

  ; define slider maximum and minumum
  ; note: slider values are subscripts to data arrays
  minvalue = 0
  maxvalue = size(dataoutput)-1

  ; slider size
  xsize = 275

  ; label size and format
  labelxsize = 37
  valueformat = '(F5.2)'

  ; input title
  inputtitle = widget_label(inputbase, value = '- INPUT PARAMETERS -', ysize = 30)

  ; sliders
  sliderbase   = widget_base  ( inputbase, /frame, /row)
  slidersource = widget_slider(sliderbase, /suppress_value, $
                               minimum = 0, maximum = maxvalue[1], $
                               title = 'Neutral source rate [10^-4 cm^-3 s^-1]', $
                               uvalue = 'slide-input', $
                               value = 0, $
                               xsize = xsize)
  labelsource  = widget_label (sliderbase, /sunken_frame, xsize = labelxsize, ysize = 30, $
                               value = string(datainput.source[0]/((1e-4)*sqrt(!dpi)*0.5*!rjkm*1e5), format = valueformat))

  sliderbase   = widget_base  (inputbase, /frame, /row)
  slidero_to_s = widget_slider(sliderbase, /suppress_value, $
                               minimum = 0, maximum = maxvalue[2], $
                               title = 'O to S ratio of neutral source', $
                               uvalue = 'slide-input', $
                               value = 0, $
                               xsize = xsize)
  labelo_to_s  = widget_label (sliderbase, /sunken_frame, xsize = labelxsize, ysize = 30, $
                               value = string(datainput.o_to_s[0], format = valueformat))  

  sliderbase  = widget_base  (inputbase, /frame, /row)
  sliderfehot = widget_slider(sliderbase, /suppress_value, $
                              minimum = 0, maximum = maxvalue[3], $
                              title = 'Fraction of hot electrons [%]', $
                              uvalue = 'slide-input', $
                              value = 0, $
                              xsize = xsize)
  labelfehot = widget_label  (sliderbase, /sunken_frame, xsize = labelxsize, ysize = 30, $
                              value = string(datainput.fehot_const[0]*100., format = valueformat))

  sliderbase      = widget_base  (inputbase, /frame, /row)
  slidertransport = widget_slider(sliderbase, /suppress_value, $
                                  minimum = 0, maximum = maxvalue[4], $
                                  title = 'Radial transport timescale [days]', $
                                  uvalue = 'slide-input', $
                                  value = 0, $
                                  xsize = xsize)
  labeltransport = widget_label  (sliderbase, /sunken_frame, xsize = labelxsize, ysize = 30, $
                               value = string(datainput.transport[0], format = valueformat))


  ; extra parameters are going to be displayed using two label widgets
  ; therefore, you need a row base for each extra parameter
  valuexsize = 70
  
  ; two different sizes for plot window depending on screen size
  if screen_size[0] lt 1400 then begin
     xsize = 700
     ysize = 700
  endif else begin
     xsize = 800
     ysize = 850
  endelse

  ; plot window
  plotrow = widget_base(plotbase, /column)
  mixtempwindow = widget_window(plotrow, sensitive = 1, $
                                uvalue = 'drawmixtemp', $
                                uname  = 'drawmixtemp', $
                                xsize = xsize, $
                                ysize = ysize)

  ; electron title
  electtitle = widget_label(inputbase, value = '- ELECTRONS -', ysize = 30)
  
  ; total density
  nel_tot = nel_tot(dataoutput[0,0,0,0].n)
  parambase  = widget_base (inputbase, /row, /align_center)
  edenslbl   = widget_label(parambase, $
                            value = '   Total equatorial electron density:')
  edensval   = widget_label(parambase, /frame, xsize = 95, $
                            value = string(nel_tot, format = '(F8.2)') + ' cm^-3')

  ; average ion charge
  avg_ion = avg_ionch(dataoutput[0,0,0,0].n)
  parambase   = widget_base (inputbase, /row, /align_right)
  avgionlbl   = widget_label(parambase, value = '(average ion charge: ' + string(avg_ion, format = '(F4.2)') + ')')

  ; temperature of cold electrons
  parambase  = widget_base (inputbase, /row, /align_center)
  ectmplbl   = widget_label(parambase, $
                            value = '       Temperature of cold electrons:')
  ectmpval   = widget_label(parambase, /frame, xsize = 95, $
                            value = string(dataoutput[0,0,0,0].Tel, format = '(F7.2)') + ' eV')

  ; temperature of hot electrons
  extraparam = widget_base (inputbase, /row, /align_center)
  tehotlbl   = widget_label(extraparam, $
                            value = 'Temperature of hot electrons (const):')
  tehotval   = widget_label(extraparam, frame = 2, xsize = 95, $
                            value = string(datainput.tehot, format = '(F4.1)') + ' eV')

  ; power breakdown title
  powertitle = widget_label(inputbase, value = '- POWER BREAKDOWN -', ysize = 30)

  ; power breakdown textbox
  powerinfobutton = widget_button(inputbase, $
                                  uvalue = 'powerinfobutton', $
                                  value = 'View Power Breakdown')

  ; spectra input title
  spectratitle = widget_label(inputbase, value = '- SIMULATE SPECTRA -', ysize = 30)

  ; spectra inputs
  spectrabase = widget_base(inputbase, /column, /align_bottom, xsize = 325, /frame, title = 'spectra knobs')
  plotspectra = widget_button(spectrabase, value = 'plot spectra', uvalue = 'plotspectra')
  minfield = cw_field(spectrabase, /floating, /all_events, $
                      title = 'Min wavelength', $
                      xsize  =     11, $
                      value  =     1100., $
                      uvalue = 'minlambda')
  maxfield = cw_field(spectrabase, /floating, /all_events, $
                      title = 'Max wavelength', $
                      xsize  =     11, $
                      value  =     1800., $
                      uvalue = 'maxlambda')
  slidergaussian = widget_slider(spectrabase, /frame, $
                                 minimum = 1, maximum = 30, $
                                 title = 'FWHM [Angstroms]', $
                                 uvalue = 'slide-gaussian', $
                                 value = 2., $ ; set initial (default) fwhm here
                                 sensitive = 0.)
                                 
  savebutton = widget_button(inputbase, $
                             value = 'Save Current Selection to File', $
                             uvalue = 'savebutton')                                 

  ; store slider + label objects in structure in order to refer to them all in one event
  ; (necessary for looking up all slider values when one slider is moved)
  widgetinfo = create_struct(    'slidersource', slidersource   , $
                                 'slidero_to_s', slidero_to_s   , $
                                  'sliderfehot', sliderfehot    , $
                              'slidertransport', slidertransport, $
                                  'labelsource', labelsource    , $
                                  'labelo_to_s', labelo_to_s    , $
                                   'labelfehot', labelfehot     , $
                               'labeltransport', labeltransport , $
                               'slidergaussian', slidergaussian , $
                                     'minfield', minfield       , $
                                     'maxfield', maxfield       , $
                              'powerinfobutton', powerinfobutton, $
                                    'powerinfo', 0l             , $
                                      'hidebox', 0l             , $  
                                     'edensval', edensval       , $
                                     'ectmpval', ectmpval       , $
                                    'avgionlbl', avgionlbl      , $
                                   'savebutton', savebutton       )

  ; realize the widgets (i.e. display it on screen)
  widget_control, base, /realize

  ; give widget info (IDs) to widget_control
  widget_control, base, set_uvalue = widgetinfo, /no_copy

  ; hand off control of the widget to the XMANAGER:
  xmanager, 'wtorus_idl81', base, group_leader = group, event_handler = 'wtorus_event', /no_block

  ; initialize mix ratio and temp plot with bottom slider values
  wtorus_plotdraw, 0, 0, 0, 0, outputstruct = dataoutput, energystruct = energy, window = mixtempwindow, /init

end
