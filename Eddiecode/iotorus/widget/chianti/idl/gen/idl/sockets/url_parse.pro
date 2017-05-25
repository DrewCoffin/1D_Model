;+
; Project     : Hinode/EIS
;
; Name        : url_parse
;
; Purpose     : parse URL into its components
;
; Category    : utility system sockets
;
; Example     : IDL> r=url_parse('http://host.domain:port/path/filename?query#fragment')
;
; Inputs      : URL = url to parse
;
; Outputs     : R = structure with tags:
;               SCHEME = URI (e.g. http:// or ftp://)
;               HOST = host or server name (inc. domain)
;               PORT = port number (e.g. 80)
;               PATH = path/filename
;               QUERY = query (e.g. a=1&b=c)
;               FRAGMENT = #fragment
;
; History     : Written 10-Feb-2007, D.M. Zarro (ADNET/GSFC)
;               Modified 21-Jul-2007, Zarro (ADNET)
;                - renamed to avoid conflict with IDL 6.4 version
;               Modified 20-Aug-2008, Zarro (ADNET) 
;                - added support for single name hosts (e.g. localhost:0)
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-

function url_parse,url

res={scheme:'',host:'',port:'',path:'',query:'',fragment:''}
if is_string(url) then begin
 url=strtrim(url,2)
 url=str_replace(url,'\','/')
 streg='([a-z]+://)?([^/\:?]+\.?[^/\:?]+)?(:[0-9]+)?(/?[^#?]+)?(\?[^#]+)?(#.+)?'
 np=n_elements(url)
 res=replicate(res,np)
 r=stregex(strcompress(url,/rem),streg,/fold,/sub,/extra)
 maxlen=max(strlen(r[3,*]))
 r[3,*]=strmid(r[3,*],1,maxlen)
 ntags=n_tags(res)
 for i=0,ntags-1 do res.(i)=reform(r[i+1,*])
endif

if ~stregex(res.host,'[a-z0-9]+',/bool,/fold) then res.host=''
;if is_blank(res.path) and is_string(res.host) then res.host=''

return,res


end
