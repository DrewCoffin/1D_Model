;+
; Project     : VSO
;
; Name        : SOCK_RESPONSE
;
; Purpose     : Wrapper around IDLnetURL object to check response
;               headers
;
; Category    : utility system sockets
;
; Syntax      : IDL> header=sock_response(url)
;
; Inputs      : URL = remote URL file name to check
;
; Outputs     : HEADER = response header 
;
; Keywords    : CODE = response code
;             : HOST_ONLY = only check host (not full path)
;             : SIZE = number of bytes in return content
;
; History     : 24-Aug-2011, Zarro (ADNET) - Written
;              
;-

;-----------------------------------------------------------------  
function sock_response_callback, status, progress, data  

;-- since we only need the response header, we just read
;   the first set of bytes and return

if (progress[0] eq 1) then begin
 if ptr_valid(data) then begin
  ourl=(*data).ourl
  if obj_valid(ourl) then begin
   ourl->GetProperty, RESPONSE_CODE=rspCode, $
         RESPONSE_HEADER=rspHdr, RESPONSE_FILENAME=rspFn
   if rspcode eq 0 then return,1
   (*data).rsphdr=rsphdr
   (*data).rspcode=rspcode
   (*data).rspfn=rspfn
   return,0
  endif
 endif
endif
return,0

end

;-----------------------------------------------------------------------------
  
function sock_response,url,err=err,verbose=verbose,code=code,passive=passive,$
                       _extra=extra,host_only=host_only,size=bsize

err='' 

if ~since_version('6.4') then begin
 err='Requires IDL version 6.4 or greater.'
 message,err,/info
 return,''
endif

if is_blank(url) then begin
 pr_syntax,'header=sock_response(url)'
 return,''
endif

durl=url
if ~has_url_scheme(durl) then durl='http://'+durl

stc=parse_url(durl)
url_path=stc.path
if is_blank(stc.path) then url_path='/'

user='anonymous'
pass='nobody'
if is_string(stc.username) then user=stc.username
if is_string(stc.password) then pass=stc.password
if is_string(username) then user=username
if is_string(password) then pass=password
url_scheme=stc.scheme

;-- check for PROXY environment variables

new_extra=sock_proxy(url,_extra=extra,verbose=verbose)

;-- create pointer to pass to callback function 

callback_function='sock_response_callback'
if is_number(passive) then passive= (0 > passive < 1) else passive=1

ourl = obj_new('IDLnetUrl')
callback_data=ptr_new({ourl:ourl,rspcode:0l,rsphdr:'', rspfn:''})

if is_string(stc.query) then url_path=url_path+'?'+stc.query
if keyword_set(host_only) then url_path='/'

ourl->setproperty,callback_function=callback_function,url_scheme=url_scheme,$
                  callback_data=callback_data,$
                  url_host=stc.host,url_username=user,url_password=pass,$
                  url_path=url_path,_extra=new_extra,ftp_connection_mode=1-passive

;-- have to use a catch since canceling the callback triggers it

error=0
catch, error
IF (error ne 0) then begin
 catch,/cancel
 message,/reset
 goto, bail
endif

result = oUrl->Get(/buffer)  

bail:

;-- retrieve response values from pointer since network object resets
;   them

rsphdr=(*callback_data).rsphdr
code=(*callback_data).rspcode
heap_free,callback_data

;-- determine content size

if arg_present(bsize) then begin
 bsize=0l
 h=stregex(rsphdr,'Content-Length: *([0-9]+) *.*',/sub,/extract,/fold)
 if is_number(h[1]) then bsize=long(h[1])
endif


return,rsphdr & end  
