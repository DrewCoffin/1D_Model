;+
; NAME: dielectronic_recombination
;
;
;
; PURPOSE: ; This code uses the formula for calculating the total
; dielectronic recombination rate given by Mazzotta et al (1998 A&AS
; 133, 403-409) 
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE: dielectronic_recombination(iz,in,t)
;
;
;
; INPUTS:
;   iz -- atomic number
;   in -- number of electrons in final state
;   t  -- electron temperature in eV
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;   returns dielectronic recombination coefficient in cm^3 s^-1
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;   right now I have only entered the coefficients for the O and S ion
;   species present in the Io torus
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;   AJS -- 11/30/2004
;-

function dielectronic_recombination,iz,in,t

case iz of
    8:case in of
        8: begin
            c=0.969000e-9
            e=15.60
        end
        7: begin
            c=0.382000e-8
            e=18.23
        end
        6: begin
            c=[0.55713e-8,0.35028e-10,0.54359e-8]
            e=[38.17,1.88,17.81]
        end 
    endcase
    16:case in of
        16: begin
            c=0.137e-8
            e=14.95
        end
        15: begin
            c=[0.80729e-8,0.11012e-9]
            e=[17.56,7.07]
        end
        14: begin
            c=[0.18172e-7,0.59195e-10]
            e=[16.62,2.40]
        end
        13: begin
            c=0.17710e-7
            e=13.46
        end
    endcase
    else: begin
        print,'DIELCTRONIC_RECOMBINATION: the ion species requested has not been'+$
              'added to this code. Returning'
        return,0
    end
endcase

alpha_d=0.

for i=0,n_elements(c)-1 do alpha_d += t^(-1.5)*c[i]*exp(-e[i]/t)

return,alpha_d
end
