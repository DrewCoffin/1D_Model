
pro getsysvariables,restore=restore

common setsysvariables_common,devname,bangp

if keyword_set(restore) then begin
    if keyword_set(devname) then set_plot,devname
    if keyword_set(bangp) then !p=bangp
    return
endif

devname=!d.name
bangp=!p

end
