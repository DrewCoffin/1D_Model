;-----------------------------------------------------------------------
PRO lat_distribution, n, h, lat_dist
;-----------------------------------------------------------------------
common info,runt,dt,zoff,rdist

; Calculate the densities as a function of z in increments of .2 Rj.
unity = replicate(1., n_elements(lat_dist[0].z))

ns0   = unity # n.s
nsp0  = unity # n.sp
ns2p0 = unity # n.s2p
ns3p0 = unity # n.s3p
;ns4p0 = unity # n.s4p
no0   = unity # n.o
nop0  = unity # n.op
no2p0 = unity # n.o2p

hs   = unity # h.s
hsp  = unity # h.sp
hs2p = unity # h.s2p
hs3p = unity # h.s3p
;hs4p = unity # h.s4p
ho   = unity # h.o
hop  = unity # h.op
ho2p = unity # h.o2p

; Neutrals
lat_dist.s   = ns0   * exp(-lat_dist.z^2/hs^2)
lat_dist.o   = no0   * exp(-lat_dist.z^2/ho^2)

; Ions
lat_dist.sp  = nsp0  * exp(-lat_dist.z^2/hsp^2)
lat_dist.s2p = ns2p0 * exp(-lat_dist.z^2/hs2p^2)
lat_dist.s3p = ns3p0 * exp(-lat_dist.z^2/hs3p^2)
;lat_dist.s4p = ns4p0 * exp(-lat_dist.z^2/hs4p^2)
lat_dist.op  = nop0  * exp(-lat_dist.z^2/hop^2)
lat_dist.o2p = no2p0 * exp(-lat_dist.z^2/ho2p^2)

; Electrons
lat_dist.el = lat_dist.sp + 2 * lat_dist.s2p + 3 * lat_dist.s3p + 4 * lat_dist.s4p + $
              lat_dist.op + 2 * lat_dist.o2p

END
