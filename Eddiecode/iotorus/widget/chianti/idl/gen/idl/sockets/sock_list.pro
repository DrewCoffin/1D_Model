;+
; Project     : HESSI
;
; Name        : SOCK_LIST
;
; Purpose     : list remote WWW page via sockets
;
; Category    : utility system sockets
;
; Syntax      : IDL> sock_list,url,page
;                  
; Inputs      : URL = URL path to list [e.g. www.cnn.com]
;
; Opt. Outputs: PAGE= captured HTML 
;
; Keywords    : ERR   = string error message
;               USE_NETWORK = set to use IDL network object (via sock_list2)
;
; History     : 27-Dec-2001,  D.M. Zarro (EITI/GSFC)  Written
;               26-Dec-2003, Zarro (L-3Com/GSFC) - added FTP capability
;               23-Dec-2005, Zarro (L-3Com/GSFC) - removed COMMON
;               27-Dec-2009, Zarro (ADNET)
;                - piped FTP list thru sock_list2
;               16-Dec-2011, Zarro (ADNET)
;                - use sock_list2 if using a PROXY server
;               13-May-2012, Zarro (ADNET)
;                - added USE_NETWORK
;
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-

pro sock_list,url,page,_ref_extra=extra,use_network=use_network

;-- check if using FTP or PROXY

page=''

use_net=keyword_set(use_network)
is_ftp=stregex(url,'ftp://',/bool) 
if use_net or is_ftp or (use_proxy() and since_version(6.4)) then begin
 sock_list2,url,page,_extra=extra

endif else begin

;-- else use HTTP object

 http=obj_new('http',_extra=extra)
 http->list,url,page,_extra=extra 
 obj_destroy,http

endelse

if n_params(0) eq 1 then print,page

return

end


