; te = electron temperature in eV

pro recombination_rates,T_e,sp,s2p,s3p,op,o2p,s4p

common recombination_rates_common,rrtemp,rrop,rro2p,rrs2p,rrs3p

if n_elements(rrtemp) ne 81 then begin
    readcol,'/Users/shinn/pro/iotorus/recombination_O.dat',rrtemp,rrop,$
            rro2p,rro3p,rrov,rrovi,rrovii,rrov2p,/silent,format='F,F,F,F,F,F,F,F'
    
    readcol,'/Users/shinn/pro/iotorus/recombination_S.dat',rrs2p,rrs3p,$
            /silent,format='F,F'
endif

k_eV=8.61739e-5
logtemp=alog10(T_e/k_eV)

; s2p = rate for S(III) + e -> S(II) + [gamma]
sp=rrfit(16,16,T_e/k_eV)+dielectronic_recombination(16,16,T_e)
s2p=interpol(rrs2p,rrtemp,logtemp)
s3p=interpol(rrs3p,rrtemp,logtemp)
if n_params() eq 7 then s4p=rrfit(16,13,T_e/k_eV)+$
  dielectronic_recombination(16,13,T_e)

op=interpol(rrop,rrtemp,logtemp)
o2p=interpol(rro2p,rrtemp,logtemp)

end


