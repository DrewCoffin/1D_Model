;
;  p -- vector containing the parameters of the fit. 
;       p[0]  =  amplitude of sine wave
;       p[1]  =  phase
;       p[2]  =  offset
;       p[3]  =  Period
;       p[4]  = 
;       p[5]  = 
;       p[6]  = 
;
;
;


function fitlonfunc,lon,p

if p[0] lt 0 then begin
    p[0]=-p[0]
    p[1]=p[1]+180.
endif
if p[1] lt 0 then p[1]=p[1]+360.
if p[1] gt 360 then p[1]=p[1]-360.

return,p[0]*cos((2*!pi/p[3])*(lon-p[1])*!dtor)+p[2]

end



