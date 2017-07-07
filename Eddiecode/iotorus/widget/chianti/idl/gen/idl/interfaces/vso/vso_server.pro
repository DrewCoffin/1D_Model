;+
; Project     : VSO
;
; Name        : VSO_SERVER
;
; Purpose     : check and return available VSO proxy server and URI
;
; Category    : vso sockets
;
; Inputs      : None
;
; Outputs     : SERVER = VSO proxy server name
;
; Optional
;  Outputs    : URI = URI name
;
; Keywords    : NETWORK = 1 if network is up
;               NO_CHECK = return server and URI names without checking network status,
;
; History     : 1-Dec-2005,  D.M. Zarro (L-3Com/GSFC), Written
;               22-Dec-2005, J.Hourcle  (L-3Com/GSFC).  changed server; URI need not resolve
;               19-Nov-2010, J.Hourcle  (Wyle/GSFC).  Failover to SAO
;               or NAO, 'no_check' ignored
;               22-Mar-2010, Zarro (ADNET). Modified to 
;               return blank string instead of stopping on error.
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-

function vso_server,uri,_ref_extra=extra,network=network,$
          no_check=no_check

    ; check=1b-keyword_set(no_check)
    network=0b

    uri='http://virtualsolar.org/VSO/VSOi'

;-- check endpoint (soap proxy)

    proxies = [ $
	 'http://sdac.virtualsolar.org/cgi-bin/vsoi_tabdelim' $
        ,'http://vso.tuc.noao.edu/cgi-bin/drms_test/vsoi_tabdelim' $
	$ , 'http://kurasuta.cfa.harvard.edu/cgi-bin/vso_tabdelim' $
     ]

    proxy = '';
    for i = 0, sizeof(proxies)-1 do begin
        proxy = proxies[i]
        network=have_network(proxy, _extra=extra)
        if ( network ) then break
    endfor

    if ( ~network ) then begin
 ;       on_error, 2
        message, 'No VSO servers available.  Please check your network connection',/cont
        return,''
    endif
    

;- return URI


    return,proxy
end
