function rfits2,files,index=fnum,date_obs=date_obs,time_obs=time_obs, $
         head=head,scale=scale, xsize=xsize, ysize=ysize, qdebug=qdebug, $
	 nodata=nodata
;+
; NAME:
;        RFITS2
; PURPOSE:
;        Reads multiple standard fits or compressed fits files into array
; CALLING SEQUENCE:
;        result = rfits2(filename)
; INPUTS:
;        files= string or string array containing the file name(s)
; OPTIONAL (KEYWORD) INPUT PARAMETERS:
;        index = nonnegative integer. If this parameter is present, a period
;              and the index number are appended to the filename (e.g., '.34').
;              This option makes handling of data in the MCCD file naming
;              convention easier.
;        scale = if set, floating point array is returned, reflecting
;                that true value = read value * bscale + bzero.
;	 xsize = the images can be resized to this value as each image is read.
;	 ysize = the images can be resized to this value as each image is read.
;		If xsize is passed, but ysize is not, then ysize=xsize.
;	 nodata= If set, then just read the headers
; OUTPUTS:
;        result = byte, integer, or long array, containing the FITS data array.
;               The dimensionality of result reflects the structure of the FITS
;               data.  If keyword scale is set, then floating point array.
; OPTIONAL (KEYWORD) OUTPUT PARAMETERS:
;        date_obs = date of observation (string).
;        time_obs = time of observation (string).
;        head     = string vector, containing the full FITS header (each element
;                 of the vector contains one FITS keyword parameter).
; COMMON BLOCKS:
;        None.
; SIDE EFFECTS:
;        None.
; RESTRICTIONS:
;        Only simple FITS files are read. FITS extensions (e.g., groups and
;        tables) are not supported.
; MODIFICATION HISTORY:
;	SLF - 7-Jul-1993 - handle compressed fits data
;	21-Jul-93 (MDM) - The FITS header was not being assembled properly for
;			  multiple file extraction
;	 7-Oct-93 (SLF) - Use home directory for temporary (uncompressed) files
;	 8-Oct-93 (SLF) - try DIR_GBO_SC first, HOME if it does not exist
;	28-Jul-94 (MDM) - Added XSIZE and YSIZE
;       15-sep-94 (SLF) - past date_obs and time_obs out correctly
;	30-Apr-96 (MDM) - Modified to allocate the memory on the first
;			  at the beginning and to insert the data into it
;			  (rather than to constantly append)
;			- Added /qdebug
;	22-Aug-96 (MDM) - Added /nodata
;-
;
qresize = 0 
if (keyword_set(xsize)) then begin
    if (n_elements(ysize) eq 0) then ysize = xsize
    qresize = 1
end
;
ifiles=files
comp=where(strpos(ifiles,'.Z') ne -1, ccount)
scomps=ccount gt 0					; some compressed

; first, decompress all compressed files
if scomps then begin
   message,/info,'Decompressing files...'
   if file_exist(getenv('DIR_GBO_SC'),/dir) then $
	outdir=getenv('DIR_GBO_SC') else outdir = getenv('HOME')
   file_uncompress,ifiles(comp),cfiles,/nor, outdir=outdir
   ifiles(comp) = cfiles				; update fits name
endif

; read the first file
if (keyword_set(qdebug)) then print, 'RFITS2 Now Reading: ' + ifiles(0)
img = rfits(ifiles(0), head=head, date_obs=date_obs, time_obs=time_obs, $
	 	scale=scale, nodata=nodata)
if (qresize) then begin
    img = congrid(img, xsize, ysize)
    sxaddpar, head, 'NAXIS1', xsize
    sxaddpar, head, 'NAXIS2', ysize
end
;
if (n_elements(ifiles) gt 1) then begin
    img0 = temporary(img)
    img = make_array(n_elements(img0(*,0)), n_elements(img0(0,*)), n_elements(ifiles), type=data_type(img0))
    img(0,0,0) = img0
end
;
; if applicable, read the rest (and append output)
for i=1,n_elements(ifiles)  - 1 do begin
    if (keyword_set(qdebug)) then print, 'RFITS2 Now Reading: ' + ifiles(i)
    img0=rfits(ifiles(i), head=head0, date_obs=date_obs0, time_obs=time_obs0, $
	scale=scale, nodata=nodata)
    if (qresize) then begin
	img0 = congrid(img0, xsize, ysize)
	sxaddpar, head0, 'NAXIS1', xsize
	sxaddpar, head0, 'NAXIS2', ysize
    end
   ;img= [[[img]], [[img0]]]
   img(0,0,i) = img0
   head = [[head], [head0]]
   date_obs = [date_obs, date_obs0]
   time_obs = [time_obs, time_obs0]
endfor 

; delete uncompressed files generated by this routine

if scomps then file_delete,cfiles

return, img

end

