;+
; Project     : HESSI
;
; Name        : FTOTAL
;
; Purpose     : Wrapper for IDL's TOTAL function that is less restrictive
;
; Category    : Utility
;
; Explanation : IDL's TOTAL function crashes if the dimension you want to sum over
;               has only one element, or if the dimension isn't there.  In the first case,
;               FTOTAL returns the array with the single dimension removed.  In the second
;               case, FTOTAL returns the array unchanged. Otherwise, FTOTAL acts just
;               like TOTAL.
;
; Syntax      : IDL> result = ftotal (array, dimension)
;
; Inputs      : array - array to be summed
;
; Opt. Inputs : dimension - dimension over which to sum (starting at 1)
;
; Outputs     : Returns sum of elements in array, or sum over dimension if dimension is specified
;
; Keywords    : Any keywords that total takes
;				PRODUCT - Use the function Product instead of Total as it has the same restrictions
;
; Common      : None
;
; Restrictions: None
;
; Side effects: None
;
; History     : 25-Aug-2001, Kim Tolbert
;				30-aug-2010, richard.schwartz@nasa.gov, added PRODUCT support
;				also array is unchanged
;				also protect against removing a 0 from a list starting with 1
;
; Contact     : kim.tolbert@gsfc.nasa.gov
;-

function ftotal, array, dimension, product=product, _extra=_extra

if not exist(dimension) then return, $
	~keyword_set(product) ? total(array, _extra=_extra) : product(array, _extra=_extra)

dim = size(array, /dim)

; case where dimension want to total over doesn't exist
if dimension gt n_elements(dim) then return, array

; normal case - just call total
if dim[dimension-1] gt 1 then $
	return, ~keyword_set(product) ? total(array, dimension, _extra=_extra) : product(array, dimension, _extra=_extra)

; case where the dimension is there, but is 1
; make sure we don't eliminate other dimensions that are 1 though
; protect against removing a 0 from a list starting with 1, ras, 30-aug-2010
if dim gt 1 then remove, dimension-1, dim
temp = array
return, reform (temp, dim)


end
