; Modified 9/1/2004 by AJS to accept CHIANTI v. 4.2 TEP files. These
; new files use the data in emissions_chianti_4.2.var and are 101x101 arrays.

pro read_emission_tables,file,temp,den,emis

ntemp = 101
nden = 101

temp = fltarr(ntemp)
den = fltarr(nden)
emis = dblarr(ntemp,nden)

openr,lun,file,/get

readf,lun,temp
readf,lun,den
readf,lun,emis
free_lun,lun

emis /= 1.60217646e-12   ;ergs/sec to eV/s
;print,temp
;print,den
;surface,emis

return
end
