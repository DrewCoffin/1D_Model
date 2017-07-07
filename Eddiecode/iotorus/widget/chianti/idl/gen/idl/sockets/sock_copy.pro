;+                                                                              
; Project     : HESSI                                                           
;                                                                               
; Name        : SOCK_COPY                                                       
;                                                                               
; Purpose     : copy file via sockets                                      
;                                                                               
; Category    : utility system sockets                                          
;                                                                               
; Syntax      : IDL> sock_copy,url,out_name,out_dir=out_dir                           
;                                                                               
; Inputs      : URL = remote file name to copy with URL path               
;               OUT_NAME = optional output name for copied file
;                                                                               
; Outputs     : None                                                           
;                                                                               
; Keywords    : OUT_DIR = output directory to copy file                         
;               ERR   = string error message                                    
;               USE_NETWORK= force using sock_get (usesIDL network object)
;               LOCAL_FILE = local name of copied file
;               BACKGROUND = download in background
;                                                                  
; History     : 27-Dec-2001,  D.M. Zarro (EITI/GSFC) - Written                  
;               23-Dec-2005, Zarro (L-3Com/GSFC) - removed COMMON               
;               26-Dec-2005, Zarro (L-3Com/GSFC) 
;                - added /HEAD_FIRST to HTTP->COPY to ensure checking for              
;                  file before copying                               
;               18-May-2006, Zarro (L-3Com/GSFC) 
;                - added IDL-IDL bridge capability for background copying                  
;               10-Nov-2006, Zarro (ADNET/GSFC) 
;                - replaced read_ftp call by ftp object
;                1-Feb-2007, Zarro (ADNET/GSFC)
;                - allowed for vector OUT_DIR
;               4-June-2009, Zarro (ADNET)
;                - improved FTP support
;               27-Dec-2009, Zarro (ADNET)
;                - piped FTP copy thru sock_get
;               18-March-2010, Zarro (ADNET)
;                - moved err and out_dir keywords into _ref_extra
;                8-Oct-2010, Zarro (ADNET)
;                - dropped support for COPY_FILE. Use LOCAL_FILE to
;                  capture name of downloaded file.
;               30-March-2011, Zarro (ADNET)
;                - restored capability to download asynchronously 
;                  using /background
;               16-December-2011, Zarro (ADNET)
;                - force using sock_get if proxy server is being used
;                                                                               
; Contact     : DZARRO@SOLAR.STANFORD.EDU                                       
;-                                                                              
                                                                                
pro sock_copy_main,url,out_name,_ref_extra=extra,local_file=local_file,$
                     use_network=use_network
                                         
if is_blank(url) then begin                                                     
 pr_syntax,'sock_copy,url'                                                      
 return                                                                         
endif                                                                           
                         
network_object=keyword_set(use_network) or (use_proxy() and since_version(6.4))
n_url=n_elements(url)                                                           
local_file=strarr(n_url)    
for i=0,n_url-1 do begin                                                        
 is_ftp=stregex(url[i],'ftp://',/bool)
 do_get=is_ftp or network_object
 if do_get then begin
  sock_get,url[i],out_name,_extra=extra,local_file=temp
 endif else begin
  if ~obj_valid(sock) then sock=obj_new('http',_extra=extra,read_timeout=10,/swap_endian)             
  sock->copy,url[i],out_name,_extra=extra,local_file=temp
 endelse
 if file_test(temp) then local_file[i]=temp
endfor                                                 
                         
if n_url eq 1 then local_file=local_file[0]
                                                                                
if obj_valid(sock) then obj_destroy,sock
                                                                                
return                                                                          
end                                                                             

;-----------------------------------------------------------------------

pro sock_copy_thread,url,out_name,reset=reset,local_file=local_file,$
                    _extra=extra,verbose=verbose,out_dir=out_dir,$
                     new_thread=new_thread

if is_blank(url) then begin                                                    
 pr_syntax,'sock_copy,url'                                                     
 return                                                                        
endif 

;-- stash variable name for returning local_file names into userdata
;   so that callback function can access it. Also store new_thread
;   keyword so that callback knows to clean up bridge object when
;   done.

lname=arg_present(local_file)? scope_varname(local_file,level=1):''
new_thread=keyword_set(new_thread)
userdata=[lname[0],trim(new_thread)]

if is_blank(out_dir) then out_dir=curdir()

;-- create IDL-IDL bridge object
;-- recycle last saved bridge object unless /new_thread
;-- make sure thread object has same SSW IDL environment/path as parent

common sock_copy,obridge_save
if keyword_set(reset) then if obj_valid(obridge_sav) then obj_destroy,obridge_sav
verbose=keyword_set(verbose)
if verbose then output=''
if ~new_thread and obj_valid(obridge_save) then obridge=obridge_save

if ~obj_valid(obridge) then begin
 obridge = Obj_New('IDL_IDLBridge',callback='sock_copy_callback',output=output)
 chk=obridge->getvar("getenv('SSW_GEN')")
 if is_blank(chk) then oBridge->execute, '@' + pref_get('IDL_STARTUP')
endif

if ~new_thread then obridge_save=obridge

if obridge->status() then begin
 message,'Thread busy. Come back again later.',/info
 return
endif

;-- command arguments

local_file=''
obridge->setvar,"url",url
obridge->setvar,"err",''
obridge->setvar,"out_name",''
obridge->setvar,"local_file",local_file
obridge->setvar,"out_dir",out_dir
obridge->setvar,"userdata",userdata
if is_string(out_name) then obridge->setvar,"out_name",out_name
obridge->setvar,"verbose",verbose

cmd='sock_copy,url,out_name,verbose=verbose,out_dir=out_dir,local_file=local_file'
if is_struct(extra) then begin
 tags=tag_names(extra)
 extra_cmd=''
 ntags=n_elements(tags)
 for i=0,ntags-1 do begin
  obridge->setvar,tags[i],extra.(i)
  extra_cmd=extra_cmd+tags[i]+'='+tags[i]
  if i lt (ntags-1) then extra_cmd=extra_cmd+','
 endfor
 cmd=cmd+','+extra_cmd
endif

;-- send copy command to thread

obridge->execute,cmd,/nowait

;-- check status

status=obridge->status(error=error)
case status of
 1: message,'Request sent to SOCK_COPY...',/info
 2: message,'Completed.',/info
 3: message,'Failed - '+error,/info
 4: message,'Aborted - '+error,/info
 else: message,'Idle',/info
endcase

if is_string(error) then message,error,/info

return & end

;-----------------------------------------------------------------------
;--- callback routine to notify when thread is complete

pro sock_copy_callback, status, error,obridge,userdata

err=obridge->getvar('err')
if is_string(err) then message,err,/info

userdata=obridge->getvar('userdata')
lname=userdata[0]
new_thread=fix(userdata[1])
verbose=obridge->getvar('verbose')
local_file=obridge->getvar('local_file')

;-- pass local_file names into lname variable passed on command line

if status eq 2 then begin
 message,'Completed.',/info
 if is_string(lname) then (scope_varfetch(lname,level=1))=local_file
endif

if is_string(error) then message,error,/info
if new_thread then obj_destroy,obridge

return & end

;--------------------------------------------------------------------------
pro sock_copy,url,out_name,_ref_extra=extra,background=background

if keyword_set(background) and since_version('6.3') then $
 sock_copy_thread,url,out_name,_extra=extra else $
  sock_copy_main,url,out_name,_extra=extra

return & end
