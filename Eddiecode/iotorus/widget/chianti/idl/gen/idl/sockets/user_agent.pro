;+
; Project     : VSO
;
; Name        : USER_AGENT
;
; Purpose     : Fake HTTP user-agent string to trick server into
;               thinking that a valid browser client is being used.
;
; Inputs      : None
;
; Outputs     : User-Agent string
;
; Keywords    : None
;
; History     : 15-January-2012, Zarro (ADNET) - written
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-


function user_agent

agentStr= 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.9) Gecko/2009041408 Red Hat/3.0.9-1.el5 Firefox/3.0.9'

return,agentStr

end
