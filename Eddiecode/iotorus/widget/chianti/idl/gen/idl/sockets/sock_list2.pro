;+
; Project     : VSO
;
; Name        : SOCK_LIST2
;
; Purpose     : Wrapper around IDLnetURL object to list URL.
;
; Category    : utility system sockets
;
; Syntax      : IDL> sock_list2,url,output
;
; Inputs      : URL = URL to list
;
; Outputs     : OUTPUT = string array (optional) 
;
; Keywords    : see IDLnetURL doc
;
; History     : 20-July-2011, Zarro (ADNET) - Written
;              
;-

pro sock_list2,url,output,_extra=extra,err=err,verbose=verbose
err=''

if ~since_version('6.4') then begin
 err='Requires IDL version 6.4 or greater.'
 message,err,/info
 return
endif

if is_blank(url) then begin
 pr_syntax,'sock_list2,url,[output]'
 return
endif


verbose=keyword_set(verbose)
CATCH, errorStatus
IF (errorStatus NE 0) THEN BEGIN  
      CATCH, /CANCEL  
  
      ; Display the error msg in a dialog and at the IDL  
      ; command line.  
;      r = DIALOG_MESSAGE(!ERROR_STATE.msg, TITLE='URL Error', $  
;         /ERROR)  
      err=!error_state.msg
      if verbose then message,err,/info  
      err='Remote listing failed.' 
      ; Get  more details about the error and  
      ; display at the IDL command line.
      if obj_valid(ourl) then begin  
       ourl->GetProperty, RESPONSE_CODE=rspCode, $  
         RESPONSE_HEADER=rspHdr, RESPONSE_FILENAME=rspFn  
       if verbose then begin
        PRINT, 'rspCode = ', rspCode  
        PRINT, 'rspHdr= ', rspHdr  
        PRINT, 'rspFn= ', rspFn  
       endif
      endif
      return
ENDIF  
  
;-- set properties

curl=strtrim(url,2)
if ~has_url_scheme(curl) then curl='http://'+curl
stc=parse_url(curl)
file=file_break(stc.path)
path=file_break(stc.path,/path)+'/'

user='anonymous'
pass='nobody'
if is_string(stc.username) then user=stc.username
if is_string(stc.password) then pass=stc.password
if is_string(username) then user=username
if is_string(password) then pass=password
url_scheme=stc.scheme
active=0
if is_number(passive) then active=1-passive

;-- check for HTTP_PROXY and NO_PROXY environment variables

new_extra=sock_proxy(curl,_extra=extra,verbose=verbose)

ourl = obj_new('IDLnetUrl')
url_path=path+file
if is_string(stc.query) then url_path=url_path+'?'+stc.query
ourl->setproperty,url_scheme=url_scheme,$
                  verbose=verbose,$
                  url_host=stc.host,url_username=user,url_password=pass,$
                  url_path=url_path,_extra=new_extra,ftp_connection_mode=active

;-- read URL 

output = oUrl->Get(_extra=extra,/string_array)
  
obj_destroy,oUrl
 
return & end  
