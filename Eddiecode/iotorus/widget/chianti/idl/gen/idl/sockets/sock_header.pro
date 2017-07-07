;+
; Project     : HESSI
;
; Name        : SOCK_HEADER
;
; Purpose     : return the HTTP header of a remote URL
;
; Category    : utility sockets 
;
; Syntax      : IDL> header=sock_header(url,code=code)
;                   
; Inputs      : URL = remote URL 
;
; Outputs     : HEADER = string header
;
; Keywords    : ERR   = string error message
;               HOST_ONLY = only check host name (without full path)
;
; History     : 28-Feb-2012, Zarro (ADNET) - written
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-

function sock_header,url,_ref_extra=extra,host_only=host_only

err=''

if is_blank(url) then begin
 pr_syntax,'header=sock_header(url)'
 return,''
endif

durl=url
if ~has_url_scheme(durl) then durl='http://'+durl

stc=parse_url(durl)
if stc.scheme eq 'ftp' then begin
 message,'Not applicable to FTP.',/info
 return,''
endif

if keyword_set(host_only) then durl=stc.host

http=obj_new('http',_extra=extra)
http->head,durl,response,_extra=extra
http->extract_content,response,_extra=extra

obj_destroy,http

return,response

end


