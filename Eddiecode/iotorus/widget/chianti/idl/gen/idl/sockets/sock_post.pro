;+
; Project     : VSO
;
; Name        : SOCK_POST
;
; Purpose     : Wrapper around IDLnetURL object to issue POST request
;
; Category    : utility system sockets
;
; Syntax      : IDL> output=sock_post(url,content)
;
; Inputs      : URL = remote URL file to send content
;               CONTENT = string content to post
;
; Outputs     : OUTPUT = server response
;
; Keywords    : HEADERS = optional string array with headers 
;                         For example: ['Accept: text/xml']
;
; History     : 23-November-2011, Zarro (ADNET) - Written
;              
;-

function sock_post,url,content,err=err,verbose=verbose,$
                       _extra=extra

err='' & output=''

if ~since_version('6.4') then begin
 err='Requires IDL version 6.4 or greater.'
 message,err,/info
 return,output
endif

if is_blank(url) or is_blank(content) then begin
 pr_syntax,'output=sock_post(url,content,headers=headers'
 return,output
endif

;-- parse out URL

durl=url
if ~has_url_scheme(durl) then durl='http://'+durl
stc=parse_url(durl)
if is_blank(stc.host) then begin
 err='Host name missing from URL.'
 message,err,/info
 return,output
endif

;-- check for PROXY environment variables

new_extra=sock_proxy(durl,_extra=extra,verbose=verbose)

ourl = obj_new('IDLnetUrl')

ourl->setproperty,url_host=stc.host,url_scheme='http',url_path=stc.path,$
                  _extra=new_extra,url_password=stc.password,url_username=stc.username,$
                   url_query=stc.query

cdir=curdir()
error=0
catch, error
IF (error ne 0) then begin
 catch,/cancel
 goto,bail
endif

;-- have to send output to writeable temp directory

tdir=get_temp_dir()
sdir=concat_dir(tdir,'temp'+get_rid())
file_mkdir,sdir
cd,sdir
result = ourl->put(content,/buffer,/post)

;-- clean up

bail: cd,cdir
obj_destroy,ourl
sresult=file_search(sdir,'*.*',count=count)
if count eq 1 then result=sresult[0]
if file_test(result) then output=rd_ascii(result) else $
 message,'POST request failed.',/verb

file_delete,sdir,/quiet,/recursive

return,output & end  
