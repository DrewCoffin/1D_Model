;+
; Project     : SDO
;
; Name        : AIA__DEFINE
;
; Purpose     : Class definition for SDO/AIA
;
; Category    : Objects
;
; History     : Written 15 June 2010, D. Zarro (ADNET)
;
; Contact     : dzarro@solar.stanford.edu
;-

function aia::init,_ref_extra=extra

if ~self->fits::init(_extra=extra) then return,0

;-- setup environment

self->setenv,_extra=extra

return,1 & end

;---------------------------------------------------

function aia::search,tstart,tend,_ref_extra=extra

return,vso_files(tstart,tend,inst='aia',_extra=extra,window=30)
end

;----------------------------------------------------
pro aia::read,file,_ref_extra=extra

err=''
self->getfile,file,local_file=afile,err=err,_extra=extra
if is_blank(afile) or is_string(err) then return
self->empty

;-- check if RICE-compressed

rice=0b
if is_rice_comp(afile) then begin
 rfile=rice_decomp(afile,err=err)
 if is_string(err) then begin
  message,err,/info
  return
 endif 
 rice=1b
endif else rfile=afile
  
chk=get_fits_det(rfile)
if ~stregex(chk,'AIA',/bool,/fold) then begin
 err='Input file not an AIA image.'
 message,err,/info
 return
endif

mrd_head,rfile,header
index=fitshead2struct(header)
if index.lvl_num lt 1.5 then begin
 if ~self->have_path(err=err,/verbose) then return
 read_sdo,rfile,index,data
 aia_prep,index,data,oindex,odata,_extra=extra
 data=temporary(odata)
 index=oindex
endif else mreadfits,rfile,index,data,_extra=extra
index2map,index,data,map,/no_copy
self->set,index=index,map=map,/no_copy
self->set,/log_scale,grid=30,/limb

if rice then file_delete,rfile,/quiet

return & end

;-----------------------------------------------------

pro aia::write,file,k,err=err,out_dir=out_dir,$
             verbose=verbose,_ref_extra=extra

if is_blank(file) then begin
 message,'Output filename not entered.',/info
 return
endif

odir=curdir()
ofile=file_break(file,path=path)
if is_string(out_dir) then odir=out_dir else if is_string(path) then odir=path
if ~file_test(odir,/directory,/write) then begin
 message,'Output directory name is invalid or does not have write access.',/info
 return
endif

pmap=self->get(/map,/pointer)
index=self->get(/index)
if ~ptr_exist(pmap) or ~is_struct(index) then begin
 message,'No data saved in map object.',/info
 return
endif

mwritefits,index,(*pmap).data,outfile=ofile,outdir=odir,loud=verbose,_extra=extra

return & end

;-----------------------------------------------------------------------------
;-- check for AIA branch in !path

function aia::have_path,err=err,verbose=verbose

err=''
if ~have_proc('read_sdo') then begin
 ssw_path,/ontology,/quiet
 if ~have_proc('read_sdo') then begin
  err='VOBS/Ontology branch of SSW not installed.'
  if keyword_set(verbose) then message,err,/info
  return,0b
 endif
endif

if ~have_proc('aia_prep') then begin
 ssw_path,/aia,/quiet
 if ~have_proc('aia_prep') then begin
  err='SDO/AIA branch of SSW not installed.'
  if keyword_set(verbose) then message,err,/info
  return,0b
 endif
endif

return,1b
end

;------------------------------------------------------------------------
;-- setup AIA environment variables

pro aia::setenv,_extra=extra

if is_string(chklog('AIA_CALIBRATION')) then return
mklog,'$SSW_AIA','$SSW/sdo/aia',/local
file_env=local_name('$SSW/sdo/aia/setup/setup.aia_env')
file_setenv,file_env,_extra=extra
return & end

;------------------------------------------------------
pro aia__define,void                 

void={aia, inherits fits}

return & end
