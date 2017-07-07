;+
; Project     : RHESSI
;
; Name        : MK_COLOR_MAP
;
; Purpose     : Add RGB color information to a map
;
; Category    : imaging
;
; Syntax      : map=mk_color_map(map,red,green,blue)
;
; Inputs      : MAP = image 
;               RED, GREEN, BLUE = color table arrays
;               (if not entered, read from color table)
;
; Outputs     : OMAP= map with RGB color arrays appended as fields
; 
; Keywords    : NO_COPY = throw away original map
;
; Written     : Zarro (ADNET) 23-Dec-2010
;                 
; Contact     : dzarro@solar.stanford.edu
;-

function mk_color_map,map,red,green,blue,no_copy=no_copy

if ~valid_map(map,have_colors=have_colors) then return,map
if have_colors then begin
 message,'Input map already colorized.',/cont
 return,map
endif
if ~exist(red)*exist(green)*exist(blue) then tvlct,red,green,blue,/get
if keyword_set(no_copy) then $ 
 return,create_struct(temporary(map),'red',red,'green',green,'blue',blue) else $
  return,create_struct(map,'red',red,'green',green,'blue',blue) 

end
