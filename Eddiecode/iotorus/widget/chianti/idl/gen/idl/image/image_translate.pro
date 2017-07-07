function image_translate, incube, offsets, interp=interp, cubic=cubic
;
;   Name: image_translate
;
;   Purpose: tranlate/align an image or cube via poly2d in user supplied offsets
;
;   Input Parameters:
;      incube  - the data cube 
;      offsets - relative x/y pixel offsets for each image (2xNIMAGES)
;                (often generated by a cross correlation such as
;                 korrel.pro, get_off.pro etc)
;
;   Keyword Parameters:
;      interp - bilinear interpolation (see poly2D doc)
;      cubic  - cubic convolution interp (see poly2D doc)
;
;   Output:
;      function returns translated/warped/aligned cube
;
;   Calling Sequence:
;      outcube=image_translate(incube, xyoffsets [,/interp] [,/cubic])
;
;   History:
;      15-October - S.L.Freeland - based upon 'translate.pro' (slater/morrison)
;                   standardized interface, documented, modernize
;
;   Category:
;      2D , 3D, CUBE, Alignment, Correlation, Image, POLY2D
;-

cubic=keyword_set(cubic)
if n_elements(interp) eq 0 then interp = 0

p = [0.,0,1,0] & q = [0.,1,0,0]            ; init P & Q  (see poly2D)

nimages=data_chk(incube,/nimage)          ; number of images in
case nimages  of          
  0: begin 
        box_message,['Need 2D or 3D image input', $
            'IDL> outcube=image_translate(incube, offsets [,/interp] [/cubic] )']
     endcase
  1: begin                                                  ; 2D case
        p(0) = -offsets(0) & q(0) = -offsets(1) 
        cube_out = poly_2d(incube,p,q,interp,cubic=cubic)
     endcase
  else: begin                                               ; 3D case
     cube_out = incube
     for i=0,nimages-1 do begin
        p(0) = -offsets(0,i) & q(0) = -offsets(1,i)
        cube_out(0,0,i) = poly_2d(incube(*,*,i),p,q,interp)
     endfor
  endcase
endcase

return,cube_out
end
