


phase = (180. + 180.)*findgen(11)/10. - 180.
time  = findgen(11)

plot, time, phase

t0 = systime(1)
for i = 0, n_elements(time) - 2 do if phase[i+1] - phase[i] gt -270  $
   and phase[i+1] - phase[i] lt 270 then plots, time[i:i+1], phase[i:i+1], linestyle = 3, color = 200
t1 = systime(1)
print, t1 - t0


t0 = systime(1)
pshift = shift(phase, 1) - phase
ind = where(pshift lt abs(270.))
if ind[0] ne -1 then ind = ind[1:n_elements(ind)-1] & $
   plots, time[ind], phase[ind]+15, linestyle = 4, color = 150
t1 = systime(1)
print, t1 - t0
