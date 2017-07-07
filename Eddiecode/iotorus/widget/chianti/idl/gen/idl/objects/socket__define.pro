;+
; Project     : RHESSI
;
; Name        : SOCKET__DEFINE
;
; Purpose     : Object wrapper around SOCK_COPY 
;
; Category    : Objects
;
; History     : Written 22 March 2011, D. Zarro (ADNET)
;
; Contact     : dzarro@solar.stanford.edu
;-

;------------------------------------------------------
pro socket::copy,url,out_name,_ref_extra=extra

sock_copy,url,out_name,_extra=extra

return

end

;------------------------------------------------------
pro socket__define

temp =  {socket,inherits dotprop}

return & end
