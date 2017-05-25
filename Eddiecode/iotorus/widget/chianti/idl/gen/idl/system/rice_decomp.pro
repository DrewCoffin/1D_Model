;+
; Project     : STEREO
;
; Name        : RICE_DECOMP
;
; Purpose     : Decompress RICE compressed file
;
; Category    : system utility 
;
; Syntax      : IDL> rfile=rice_decomp(file)
;
; Inputs      : FILE = RICE compressed file name
;
; Outputs     : RFILE = decompressed file name
;
; Keywords    : ERR= error string
;
; History     : 21-Nov-2011, Zarro (ADNET) - written
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-
                                                                                         
function rice_decomp,file,err=err,_ref_extra=extra

if ~is_rice_comp(file,err=err) then return,''

if ~have_proc('mreadfits_tilecomp') then begin
 epath=local_name('$SSW/vobs/ontology/idl/jsoc')
 if is_dir(epath) then add_path,epath,/quiet,/append
 if ~have_proc('mreadfits_tilecomp') then begin
  err='Missing RICE decompressor function - mreadfits_tilecomp.'
  message,err,/info
  return,''
 endif
endif

;-- always return to current directory in case we switched or had errors

cd,current=cdir
error=0
catch,error
if error ne 0 then begin
 catch,/cancel
 cd,cdir
 return,''
endif

rdir=file_dirname(file)
rfile=file_basename(file)

;-- kluge for Windows

if os_family() eq 'Windows' then begin
 cd,rdir & hide=1
endif else rfile=file

mreadfits_tilecomp,rfile,index,/nodata,fnames_uncomp=fnames_uncomp,$
 /silent,/noshell,/only_uncompress,_extra=extra,hide=hide

cd,cdir
chk=loc_file(fnames_uncomp,count=count)
if count eq 0 then begin
 err='Decompression failed.'
 message,err,/info
 return,''
endif

return,fnames_uncomp
end
