;+
; Project     : VSO
;
; Name        : THREAD
;
; Purpose     : Wrapper around IDL-IDLBridge object to run any procedure
;               in a background thread
;
; Category    : utility objects
;
; Example     : IDL> thread,'proc_name',arg1,arg2,...arg10,key1=key1,key2=key2...
;
; Inputs      : PROC = procedure name
;               ARGi = arguments accepted by PROC (up to 10)
;
; Keywords    : KEYi = keywords accepted by PROC
;
; Restrictions: Minimal error checking, can't handle functions
;               or return modified argument or keyword values (yet)
;
; History     : 22-Feb-2012, Zarro - Written
;-

;--- call back routine to notify when thread is complete

pro thread_callback, status, error, oBridge, userdata

;-- check for modified return values (can't return these to
;   caller level yet)

if status eq 2 then begin
 ndata=n_elements(userdata)
 if ndata gt 0 then begin
  for i=0,ndata-1 do begin
   var_val=obridge->getvar(userdata[i])
  endfor
 endif
 message,'Completed.',/info
endif

if is_string(error) then message,error,/info

return & end

;---------------------------------------------------------------------------------

pro thread,proc,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,_ref_extra=extra,reset=reset,verbose=verbose,err=err

if ~have_proc(proc) then begin
 pr_syntax,"thread,'procedure name',arg0,arg1,...arg9,keyword=keyword"
 return
endif

common thread,obridge

if keyword_set(reset) then if obj_valid(obridge) then obj_destroy,obridge
verbose=keyword_set(verbose)

err=''
if ~since_version('6.3') then begin
 err='Requires at least IDL version 6.3.'
 message,err,/info
 return
endif

;-- create IDL-IDL bridge object
;-- make sure thread object has same SSW IDL environment/path as parent

if verbose then output=''

if ~obj_valid(obridge) then begin
 obridge = Obj_New('IDL_IDLBridge',callback='thread_callback',output=output)
 chk=obridge->getvar("getenv('SSW_GEN')")
 if is_blank(chk) then oBridge->execute, '@' + pref_get('IDL_STARTUP')
endif

if obridge->status() then begin
 message,'Thread busy. Come back again later.',/info
 return
endif

;-- pass argument values into bridge object

cmd=proc
for i=1,n_params()-1 do begin
 var_val=''
 var_name='p'+strtrim(string(i),2)
 if exist(scope_varfetch(var_name)) then var_val=scope_varfetch(var_name)
 obridge->setvar,var_name,var_val
 cmd=cmd+','+var_name
 lname=arg_present(var_name)? scope_varname(var_name,level=1):''
 userdata=append_arr(userdata,var_name)
endfor

;-- pass keyword values into bridge object

extra_cmd=''
if is_string(extra) then begin
 ntags=n_elements(extra)
 for i=0,ntags-1 do begin
  var_val=''
  var_name=extra[i]
  if exist(scope_varfetch(var_name,/ref))then var_val=scope_varfetch(var_name,/ref)
  obridge->setvar,var_name,var_val
  extra_cmd=extra_cmd+var_name+'='+var_name
  if i lt (ntags-1) then extra_cmd=extra_cmd+','
  lname=arg_present(var_name)? scope_varname(var_name,level=1):''
  userdata=append_arr(userdata,var_name)
 endfor
 cmd=cmd+','+extra_cmd
endif

;-- send command to thread

if verbose then message,cmd,/info
if exist(userdata) then obridge->setproperty,userdata=userdata
obridge->execute,cmd,/nowait

;-- check status

case obridge->status(err=err) of
 1: message,'Submitted',/info
 2: message,'Completed.',/info
 3: message,'Failed - '+err,/info
 4: message,'Aborted - '+err,/info
 else: message,'Idle',/info
endcase

if is_string(err) then message,err,/info

return & end
