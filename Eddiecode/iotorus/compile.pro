;+
; IDL script file for compiling nescessary procedures
; type @compile into the IDL command prompt to execute
;-

; if compile script has already run, set boolean 'compile' equal to false
if ~keyword_set(compile) then compile = 1 else compile = 0

; get full path of current directory
cd, current = current

if compile then defsysv, '!rjkm', 71492.
if compile then defsysv, '!kev', 8.617385d-5 ; Boltzmann's constant in eV
if compile then defsysv, '!rootpi', sqrt(!dpi)

; jhpl astron library
astrolib = ':' + expand_path('+' + current + '/model/astrolib/')
if compile then !path = !path + astrolib

; create a string of all directories under cc model folder separated by ':'s
ccmodelpros = ':' + expand_path(current + '/model/cm3_latavg_onebox/')
if compile then !path = !path + ccmodelpros

; compile the model's supporting procedures
.r cm3_latavg_model

; CHIANTI
; create a string of all directories under chianti folder separated by ':'s
;chiantipros = ':' + expand_path('+' + current + '/chianti/idl/')
;if compile then !path = !path + chiantipros
;use_chianti, current + '/chianti/dbase'
