;+
; Project     : STEREO
;
; Name        : IS_RICE_COMP
;
; Purpose     : Check if a RICE-compressed file
;
; Category    : system utility 
;
; Syntax      : IDL> chk=is_rice_comp(file)
;
; Inputs      : FILE = input file to check
;
; Outputs     : CHK = 1 or 0
;
; Keywords    : ERR= error string
;
; History     : 9-Apr-2012, Zarro (ADNET) - written
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-
                                                                                         
function is_rice_comp,file,_ref_extra=extra,err=err

;-- search extension header for RICE compression keyword 

i=0
repeat begin
 err=''
 mrd_head,file,header,ext=i,err=err,_extra=extra,/no_check
 if is_blank(terr) then begin
  chk=where(stregex(header,'cmp.+rice.+compression',/bool,/fold),count)
  if count ne 0 then return,1b
  i=i+1
 endif
endrep until is_string(err)

err='Not a RICE-compressed file.'

return,0b
end
