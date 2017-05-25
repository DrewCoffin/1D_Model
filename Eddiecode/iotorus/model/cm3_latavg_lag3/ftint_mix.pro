; Simple function to convert density structures output from the
; latitude-averaged cm3_model and convert them into mixing ratio
; structures.

function ftint_mix,n,h

mix=replicate(create_struct(name='ftint_mix_ratios',n[0]),size(n,/dim))

mix.el=n.el
mix.elh=n.elh
mix.s = n.s
mix.o = n.o
mix.fc = n.fc
mix.fh = n.fh

N_sp=n.sp*sqrt(!pi)*h.sp
N_s2p=n.s2p*sqrt(!pi)*h.s2p
N_s3p=n.s3p*sqrt(!pi)*h.s3p
N_s4p=n.s4p*sqrt(!pi)*h.s4p
N_op=n.op*sqrt(!pi)*h.op
N_o2p=n.o2p*sqrt(!pi)*h.o2p

protons = n.protons
nel_tot=(n_sp+2*n_s2p+3*n_s3p+4*n_s4p+n_op+2*n_o2p)/(1.-protons)
;nel_tot=(n_sp+2*n_s2p+3*n_s3p+4*n_s4p+n_op+2*n_o2p)

mix.sp=N_sp/Nel_tot
mix.s2p=N_s2p/Nel_tot
mix.s3p=N_s3p/Nel_tot
mix.s4p=N_s4p/Nel_tot
mix.op=N_op/Nel_tot
mix.o2p=N_o2p/Nel_tot

return,mix
end
