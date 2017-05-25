;+
; IDL procedure for setting up the Chianti database as well as
; compiling and running the io torus widget.
; To run the widget
; IDL> wtorus
;
; ASSUMPTIONS
;  - user is using a unix-based operating system (can be modified to work on a windows platform)
;  - user has IDL 8.1 or later
;  - Chianti is set up and its file structure is all lower case
;-
pro wtorus

  ; set variable to path of main chianti directory, assuming chianti is in the current working directory
  cd, current = current
  chiantipath = current + '/chianti/'

  ; or, if chianti is located elsewhere, edit the following variable
  ;chiantipath = '/Users/shinn/pro/chianti/' ; <---- user specific path

  ; use 'exists' keyword to determine whether to execute the following
  defsysv, '!chianti', exists = compile
  if ~compile then defsysv, '!chianti', 1
 
  if ~compile then defsysv, '!rjkm', 71492. ; RJ in km
  if ~compile then defsysv, '!kev', 8.617385d-5 ; Boltzmann's constant in eV
  if ~compile then defsysv, '!rootpi', sqrt(!dpi)

  ; create a string of all directories under chianti folder separated by ':'s
  if ~compile then chiantipros = ':' + expand_path('+' + chiantipath + 'idl/')

  ; if compile script has not been run before, add chianti path
  if ~compile then !path = !path + chiantipros
  if ~compile then use_chianti, chiantipath + 'dbase'

  ; compile the widget
  ; note: with this command, it is unnecessary to compile the widget seperately
  resolve_routine, 'wtorus_idl81'

  ; run the widget
  wtorus_idl81

end
