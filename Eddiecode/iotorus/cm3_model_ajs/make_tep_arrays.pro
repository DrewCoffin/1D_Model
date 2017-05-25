function make_tep_values, Tel, nel,ns,nsp,ns2p,ns3p,ns4p,no,nop,no2p, $
  min = min, max = max
  ;Tel in eV, nel in #/cm^3, output is total emitted power in ergs/s 
  
  if keyword_set(max) then maxwav = max else maxwav = 1150d
  if keyword_set(min) then minwav = min else minwav = 550d



  conversion=11604.5221 ; conversion between eV and Kelvin from NIST
  ; emission line calculation: sulfur sions
  s1em = emiss_calc_ajs(16, 1, temp = alog10(Tel*conversion), dens = alog10(nel), radt = 1., /quiet)
  s2em = emiss_calc_ajs(16, 2, temp = alog10(Tel*conversion), dens = alog10(nel), radt = 1., /quiet)
  s3em = emiss_calc_ajs(16, 3, temp = alog10(Tel*conversion), dens = alog10(nel), radt = 1., /quiet)
  s4em = emiss_calc_ajs(16, 4, temp = alog10(Tel*conversion), dens = alog10(nel), radt = 1., /quiet)
  s5em = emiss_calc_ajs(16,5,temp = alog10(Tel*conversion), dens = alog10(nel), radt = 1., /quiet)
  ; emission line calculation: oxygen ions
  o1em = emiss_calc_ajs(8, 1, temp = alog10(Tel*conversion), dens = alog10(nel), radt = 1., /quiet)
  o2em = emiss_calc_ajs(8, 2, temp = alog10(Tel*conversion), dens = alog10(nel), radt = 1., /quiet)
  o3em = emiss_calc_ajs(8, 3, temp = alog10(Tel*conversion), dens = alog10(nel), radt = 1., /quiet)




s1emiss=(reform(s1em.em))*ns
  s2emiss = (reform(s2em.em))*nsp
  s3emiss = (reform(s3em.em))*ns2p
  s4emiss = (reform(s4em.em))*ns3p
  s5emiss=  (reform(s5em.em))*ns4p
  o1emiss = (reform(o1em.em))*no
  o2emiss = (reform(o2em.em))*nop
  o3emiss = (reform(o3em.em))*no2p

  ; discrete
  xwavi = [s1em.lambda,s2em.lambda, s3em.lambda, s4em.lambda, s5em.lambda, o1em.lambda, o2em.lambda, o3em.lambda]
  wsort = sort(xwavi) ; sort by wavelength
  xwavi = xwavi[wsort]

  yptsi = [s1emiss, s2emiss, s3emiss, s4emiss, s5emiss, o1emiss, o2emiss, o3emiss]
  yptsi = yptsi[wsort]

  yname = [s1em.ion_name, s2em.ion_name, s3em.ion_name, s4em.ion_name, s5em.ion_name, o1em.ion_name, o2em.ion_name, o3em.ion_name]
  yname = yname[wsort]

  ; cut down arrays to necessary range
  avgwav = (minwav + maxwav)/2.
  wrange = where(abs(xwavi - avgwav) le avgwav - minwav)
  xwavi = xwavi[wrange]
  yptsi = yptsi[wrange]
  yname = yname[wrange]


  output=total(yptsi)



  return, output

  ; note: xwav is also an output, placed as a keyword

end




pro make_tep_arrays

  ; start time
  tinitial = systime(/julian)
  start = string(systime(/julian), format = '(c())')
  print, 'start time:' + start
  
  openr,1,'temperatures_for_tep.txt'
  tel_array=fltarr(101)
  readf,1,tel_array
  close,1
  
  
  openr,1,'densities_for_tep.txt'
  nel_array=fltarr(101)
  readf,1,nel_array
  close,1


;tel=5. ;eV nominal
;nel=2000. ; #/cc nominal
output_array1=fltarr(101,101)
output_array2=fltarr(101,101)
output_array3=fltarr(101,101)
output_array4=fltarr(101,101)
output_array5=fltarr(101,101)
output_array6=fltarr(101,101)
output_array7=fltarr(101,101)
output_array8=fltarr(101,101)
minvalue=1.
maxvalue=600000.
ns=1
nsp=0
ns2p=0
ns3p=0
ns4p=0
no=0
nop=0
no2p=0
;a=0
;b=0
for a=0,100 do begin
for b=0,100 do begin
  tel=tel_array[a]
  nel=nel_array[b]
output=make_tep_values(tel,nel,ns,nsp,ns2p,ns3p,ns4p,no,nop,no2p,min = minvalue, max = maxvalue)
output_array1(a,b)=output
endfor
endfor

;;;



ns=0
nsp=1.
ns2p=0
ns3p=0
ns4p=0
no=0
nop=0
no2p=0
;a=0
;b=0
for a=0,100 do begin
  for b=0,100 do begin
    tel=tel_array[a]
    nel=nel_array[b]
    output=make_tep_values(tel,nel,ns,nsp,ns2p,ns3p,ns4p,no,nop,no2p,min = minvalue, max = maxvalue)
    output_array2(a,b)=output
  endfor
endfor

;;;


ns=0
nsp=0
ns2p=1.
ns3p=0
ns4p=0
no=0
nop=0
no2p=0
;a=0
;b=0
for a=0,100 do begin
  for b=0,100 do begin
    tel=tel_array[a]
    nel=nel_array[b]
    output=make_tep_values(tel,nel,ns,nsp,ns2p,ns3p,ns4p,no,nop,no2p,min = minvalue, max = maxvalue)
    output_array3(a,b)=output
  endfor
endfor

;;


ns=0
nsp=0
ns2p=0
ns3p=1
ns4p=0
no=0
nop=0
no2p=0
;a=0
;b=0
for a=0,100 do begin
  for b=0,100 do begin
    tel=tel_array[a]
    nel=nel_array[b]
    output=make_tep_values(tel,nel,ns,nsp,ns2p,ns3p,ns4p,no,nop,no2p,min = minvalue, max = maxvalue)
    output_array4(a,b)=output
  endfor
endfor


ns=0
nsp=0
ns2p=0
ns3p=0
ns4p=1
no=0
nop=0
no2p=0
;a=0
;b=0
for a=0,100 do begin
  for b=0,100 do begin
    tel=tel_array[a]
    nel=nel_array[b]
    output=make_tep_values(tel,nel,ns,nsp,ns2p,ns3p,ns4p,no,nop,no2p,min = minvalue, max = maxvalue)
    output_array5(a,b)=output
  endfor
endfor

ns=0
nsp=0
ns2p=0
ns3p=0
ns4p=0
no=1.
nop=0
no2p=0
;a=0
;b=0
for a=0,100 do begin
  for b=0,100 do begin
    tel=tel_array[a]
    nel=nel_array[b]
    output=make_tep_values(tel,nel,ns,nsp,ns2p,ns3p,ns4p,no,nop,no2p,min = minvalue, max = maxvalue)
    output_array6(a,b)=output
  endfor
endfor

ns=0
nsp=0
ns2p=0
ns3p=0
ns4p=0
no=0
nop=1.
no2p=0
;a=0
;b=0
for a=0,100 do begin
  for b=0,100 do begin
    tel=tel_array[a]
    nel=nel_array[b]
    output=make_tep_values(tel,nel,ns,nsp,ns2p,ns3p,ns4p,no,nop,no2p,min = minvalue, max = maxvalue)
    output_array7(a,b)=output
  endfor
endfor

ns=0
nsp=0
ns2p=0
ns3p=0
ns4p=0
no=0
nop=0
no2p=1.
;a=0
;b=0
for a=0,100 do begin
  for b=0,100 do begin
    tel=tel_array[a]
    nel=nel_array[b]
    output=make_tep_values(tel,nel,ns,nsp,ns2p,ns3p,ns4p,no,nop,no2p,min = minvalue, max = maxvalue)
    output_array8(a,b)=output
  endfor
endfor

; end time
tfinal=systime(/julian)
stop = string(systime(/julian), format = '(c())')
print, '  end time: ' + stop
delta_t = 24*60*(tfinal-tinitial)
print,'  time elapsed : ' , delta_t


write_csv,'ns_output_tep_values.txt',output_array1 ; 101 by 101 output array
write_csv,'nsp_output_tep_values.txt',output_array2 ; 101 by 101 output array
write_csv,'ns2p_output_tep_values.txt',output_array3 ; 101 by 101 output array
write_csv,'ns3p_output_tep_values.txt',output_array4 ; 101 by 101 output array
write_csv,'ns4p_output_tep_values.txt',output_array5 ; 101 by 101 output array
write_csv,'no_output_tep_values.txt',output_array6 ; 101 by 101 output array
write_csv,'nop_output_tep_values.txt',output_array7 ; 101 by 101 output array
write_csv,'no2p_output_tep_values.txt',output_array8 ; 101 by 101 output array
stop

end