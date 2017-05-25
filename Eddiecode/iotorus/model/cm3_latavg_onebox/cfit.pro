; cfit.pro 
;
; Adapted from cfit.f
;
;* Version 2, March 24, 1997
;*****************************************************************************
;*** This subroutine calculates rates of direct collisional ionization 
;*** for all ionization stages of all elements from H to Ni (Z=28)
;*** by use of the fits from G. S. Voronov, 1997, ADNDT, 65, 1
;*** Input parameters:  iz - atomic number 
;***                    in - number of electrons from 1 to iz 
;***                    t  - temperature, eV
;*** Output parameter:  c  - rate coefficient, cm^3 s^(-1)
;*****************************************************************************
;
function cfit,iz,in,te,c
@cm3_model_common

if n_elements(cf) ne 5046l then begin 
    cf=fltarr(6,29,29)

    cf[*, 1, 1]=[0.,13.6,0.,2.91E-08,0.2320,0.39]
    cf[*, 2, 2]=[0., 24.6,0.,1.75E-08,0.1800,0.35]
    cf[*, 2, 1]=[0., 54.4,1.,2.05E-09,0.2650,0.25]
    cf[*, 3, 3]=[0.,  5.4,0.,1.39E-07,0.4380,0.41]
    cf[*, 3, 2]=[0., 75.6,1.,2.01E-09,0.2090,0.23]
    cf[*, 3, 1]=[0.,122.4,1.,9.60E-10,0.5820,0.17]
    cf[*, 4, 4]=[0.,  9.3,0.,1.02E-07,0.3750,0.27]
    cf[*, 4, 3]=[0., 18.2,1.,2.08E-08,0.4390,0.21]
    cf[*, 4, 2]=[0.,153.9,0.,2.67E-09,0.6120,0.27]
    cf[*, 4, 1]=[0.,217.7,1.,4.27E-10,0.6580,0.15]
    cf[*, 5, 5]=[0.,  8.3,0.,6.49E-08,0.2000,0.26]
    cf[*, 5, 4]=[0., 25.2,1.,1.24E-08,0.2670,0.22]
    cf[*, 5, 3]=[0., 37.9,1.,3.27E-09,0.2950,0.23]
    cf[*, 5, 2]=[0.,259.4,1.,4.95E-10,0.4890,0.09]
    cf[*, 5, 1]=[0.,340.2,1.,2.19E-10,0.6570,0.15]
    cf[*, 6, 6]=[0., 11.3,0.,6.85E-08,0.1930,0.25]
    cf[*, 6, 5]=[0., 24.4,1.,1.86E-08,0.2860,0.24]
    cf[*, 6, 4]=[0., 47.9,1.,6.35E-09,0.4270,0.21]
    cf[*, 6, 3]=[0., 64.5,1.,1.50E-09,0.4160,0.13]
    cf[*, 6, 2]=[0.,392.1,1.,2.99E-10,0.6660,0.02]
    cf[*, 6, 1]=[0.,490.0,1.,1.23E-10,0.6200,0.16]
    cf[*, 7, 7]=[0., 14.5,0.,4.82E-08,0.0652,0.42]
    cf[*, 7, 6]=[0., 29.6,0.,2.98E-08,0.3100,0.30]
    cf[*, 7, 5]=[0., 47.5,1.,8.10E-09,0.3500,0.24]
    cf[*, 7, 4]=[0., 77.5,1.,3.71E-09,0.5490,0.18]
    cf[*, 7, 3]=[0., 97.9,0.,1.51E-09,0.0167,0.74]
    cf[*, 7, 2]=[0.,552.1,0.,3.71E-10,0.5460,0.29]
    cf[*, 7, 1]=[0.,667.0,1.,7.77E-11,0.6240,0.16]
    cf[*, 8, 8]=[0., 13.6,0.,3.59E-08,0.0730,0.34]
    cf[*, 8, 7]=[0., 35.1,1.,1.39E-08,0.2120,0.22]
    cf[*, 8, 6]=[0., 54.9,1.,9.31E-09,0.2700,0.27]
    cf[*, 8, 5]=[0., 77.4,0.,1.02E-08,0.6140,0.27]
    cf[*, 8, 4]=[0.,113.9,1.,2.19E-09,0.6300,0.17]
    cf[*, 8, 3]=[0.,138.1,0.,1.95E-09,0.3600,0.54]
    cf[*, 8, 2]=[0.,739.3,0.,2.12E-10,0.3960,0.35]
    cf[*, 8, 1]=[0.,871.4,1.,5.21E-11,0.6290,0.16]
    cf[*, 9, 9]=[0., 17.4,1.,7.00E-08,0.1780,0.29]
    cf[*, 9, 8]=[0., 35.0,0.,5.41E-08,0.5710,0.27]
    cf[*, 9, 7]=[0., 62.7,1.,9.37E-09,0.3190,0.20]
    cf[*, 9, 6]=[0., 87.1,1.,4.92E-09,0.3230,0.24]
    cf[*, 9, 5]=[0.,114.2,0.,7.06E-09,0.6840,0.27]
    cf[*, 9, 4]=[0.,157.2,1.,1.28E-09,0.6480,0.16]
    cf[*, 9, 3]=[0.,185.2,1.,5.61E-10,0.7380,0.16]
    cf[*, 9, 2]=[0.,953.9,0.,1.66E-10,0.5420,0.29]
    cf[*, 9, 1]=[0.,1103.1,1.,3.74E-11,0.6590,0.15]
    cf[*,10,10]=[0., 21.6,1.,1.50E-08,0.0329,0.43]
    cf[*,10, 9]=[0., 41.0,0.,1.98E-08,0.2950,0.20]
    cf[*,10, 8]=[0., 63.5,1.,7.03E-09,0.0677,0.39]
    cf[*,10, 7]=[0., 97.1,1.,4.24E-09,0.0482,0.58]
    cf[*,10, 6]=[0.,126.2,1.,2.79E-09,0.3050,0.25]
    cf[*,10, 5]=[0.,157.9,0.,3.45E-09,0.5810,0.28]
    cf[*,10, 4]=[0.,207.3,1.,9.56E-10,0.7490,0.14]
    cf[*,10, 3]=[0.,239.1,1.,4.73E-10,0.9920,0.04]
    cf[*,10, 2]=[0.,1196.0,1.,3.92E-11,0.2620,0.20]
    cf[*,10, 1]=[0.,1360.6,1.,2.77E-11,0.6610,0.13]
    cf[*,11,11]=[0.,  5.1,1.,1.01E-07,0.2750,0.23]
    cf[*,11,10]=[0., 47.3,1.,7.35E-09,0.0560,0.35]
    cf[*,11, 9]=[0., 71.6,1.,8.10E-09,0.1480,0.32]
    cf[*,11, 8]=[0., 98.9,0.,1.14E-08,0.5530,0.28]
    cf[*,11, 7]=[0.,138.4,1.,2.63E-09,0.2300,0.29]
    cf[*,11, 6]=[0.,172.2,1.,1.85E-09,0.3630,0.22]
    cf[*,11, 5]=[0.,208.5,0.,2.82E-09,0.6740,0.27]
    cf[*,11, 4]=[0.,264.2,1.,6.72E-10,0.7520,0.14]
    cf[*,11, 3]=[0.,299.9,1.,2.80E-10,0.7810,0.15]
    cf[*,11, 2]=[0.,1465.1,1.,4.63E-11,0.5580,0.16]
    cf[*,11, 1]=[0.,1648.7,1.,2.16E-11,0.7430,0.13]
    cf[*,12,12]=[0.,  7.6,0.,6.21E-07,0.5920,0.39]
    cf[*,12,11]=[0., 15.2,0.,1.92E-08,0.0027,0.85]
    cf[*,12,10]=[0., 80.1,1.,5.56E-09,0.1070,0.30]
    cf[*,12, 9]=[0.,109.3,1.,4.35E-09,0.1590,0.31]
    cf[*,12, 8]=[0.,141.3,0.,7.10E-09,0.6580,0.25]
    cf[*,12, 7]=[0.,186.5,1.,1.70E-09,0.2420,0.28]
    cf[*,12, 6]=[0.,224.9,1.,1.22E-09,0.3430,0.23]
    cf[*,12, 5]=[0.,266.0,0.,2.20E-09,0.8970,0.22]
    cf[*,12, 4]=[0.,328.2,1.,4.86E-10,0.7510,0.14]
    cf[*,12, 3]=[0.,367.5,1.,2.35E-10,1.0300,0.10]
    cf[*,12, 2]=[0.,1761.8,1.,2.06E-11,0.1960,0.25]
    cf[*,12, 1]=[0.,1962.7,1.,1.75E-11,0.8350,0.11]
    cf[*,13,13]=[0.,  6.0,1.,2.28E-07,0.3870,0.25]
    cf[*,13,12]=[0., 18.8,0.,1.18E-07,2.2100,0.25]
    cf[*,13,11]=[0., 28.5,1.,4.40E-09,0.1060,0.24]
    cf[*,13,10]=[0.,120.0,0.,1.75E-08,0.8720,0.22]
    cf[*,13, 9]=[0.,153.8,1.,2.61E-09,0.1590,0.31]
    cf[*,13, 8]=[0.,198.5,1.,1.85E-09,0.1520,0.36]
    cf[*,13, 7]=[0.,241.4,1.,1.14E-09,0.2280,0.29]
    cf[*,13, 6]=[0.,284.6,1.,8.00E-10,0.4170,0.16]
    cf[*,13, 5]=[0.,390.2,1.,5.83E-10,0.4970,0.23]
    cf[*,13, 4]=[0.,399.4,0.,4.93E-10,0.7060,0.16]
    cf[*,13, 3]=[0.,442.0,1.,9.77E-11,0.2780,0.17]
    cf[*,13, 2]=[0.,2086.6,0.,3.94E-11,0.2860,0.36]
    cf[*,13, 1]=[0.,2304.1,1.,1.38E-11,0.8350,0.11]
    cf[*,14,14]=[0.,  8.2,1.,1.88E-07,0.3760,0.25]
    cf[*,14,13]=[0., 16.4,1.,6.43E-08,0.6320,0.20]
    cf[*,14,12]=[0., 33.5,1.,2.01E-08,0.4730,0.22]
    cf[*,14,11]=[0., 54.0,1.,4.94E-09,0.1720,0.23]
    cf[*,14,10]=[0.,166.8,1.,1.76E-09,0.1020,0.31]
    cf[*,14, 9]=[0.,205.3,1.,1.74E-09,0.1800,0.29]
    cf[*,14, 8]=[0.,246.5,1.,1.23E-09,0.5180,0.07]
    cf[*,14, 7]=[0.,303.5,1.,8.27E-10,0.2390,0.28]
    cf[*,14, 6]=[0.,351.1,1.,6.01E-10,0.3050,0.25]
    cf[*,14, 5]=[0.,401.4,1.,4.65E-10,0.6660,0.04]
    cf[*,14, 4]=[0.,476.4,1.,2.63E-10,0.6660,0.16]
    cf[*,14, 3]=[0.,523.5,1.,1.18E-10,0.7340,0.16]
    cf[*,14, 2]=[0.,2437.7,0.,3.36E-11,0.3360,0.37]
    cf[*,14, 1]=[0.,2673.2,1.,1.19E-11,0.9890,0.08]
    cf[*,15,15]=[0., 10.5,1.,1.99E-07,0.5350,0.24]
    cf[*,15,14]=[0., 19.8,1.,5.88E-08,0.5370,0.21]
    cf[*,15,13]=[0., 30.2,1.,2.96E-08,0.8650,0.16]
    cf[*,15,12]=[0., 51.4,1.,1.01E-08,0.5460,0.20]
    cf[*,15,11]=[0., 65.0,1.,2.36E-09,0.1920,0.17]
    cf[*,15,10]=[0.,220.4,0.,6.66E-09,1.0000,0.18]
    cf[*,15, 9]=[0.,263.2,1.,1.24E-09,0.2150,0.26]
    cf[*,15, 8]=[0.,309.4,0.,2.27E-09,0.7340,0.23]
    cf[*,15, 7]=[0.,371.7,1.,6.14E-10,0.2560,0.27]
    cf[*,15, 6]=[0.,424.5,1.,4.69E-10,0.3420,0.23]
    cf[*,15, 5]=[0.,479.6,0.,6.14E-10,0.3340,0.39]
    cf[*,15, 4]=[0.,560.4,0.,3.22E-10,0.8500,0.12]
    cf[*,15, 3]=[0.,611.9,1.,9.32E-11,0.7340,0.16]
    cf[*,15, 2]=[0.,2816.9,0.,3.79E-11,0.8050,0.22]
    cf[*,15, 1]=[0.,3069.9,1.,9.73E-12,0.9910,0.08]
    cf[*,16,16]=[0., 10.4,1.,5.49E-08,0.1000,0.25]
    cf[*,16,15]=[0., 23.3,1.,6.81E-08,0.6930,0.21]
    cf[*,16,14]=[0., 34.8,1.,2.14E-08,0.3530,0.24]
    cf[*,16,13]=[0., 47.3,1.,1.66E-08,1.0300,0.14]
    cf[*,16,12]=[0., 72.6,1.,6.12E-09,0.5800,0.19]
    cf[*,16,11]=[0., 88.1,1.,1.33E-09,0.0688,0.35]
    cf[*,16,10]=[0.,280.9,0.,4.93E-09,1.1300,0.16]
    cf[*,16, 9]=[0.,328.2,1.,8.73E-10,0.1930,0.28]
    cf[*,16, 8]=[0.,379.1,0.,1.35E-09,0.4310,0.32]
    cf[*,16, 7]=[0.,447.1,1.,4.59E-10,0.2420,0.28]
    cf[*,16, 6]=[0.,504.8,1.,3.49E-10,0.3050,0.25]
    cf[*,16, 5]=[0.,564.7,0.,5.23E-10,0.4280,0.35]
    cf[*,16, 4]=[0.,651.6,0.,2.59E-10,0.8540,0.12]
    cf[*,16, 3]=[0.,707.2,1.,7.50E-11,0.7340,0.16]
    cf[*,16, 2]=[0.,3223.9,0.,2.67E-11,0.5720,0.28]
    cf[*,16, 1]=[0.,3494.2,1.,6.32E-12,0.5850,0.17]
    cf[*,17,17]=[0., 13.0,1.,1.69E-07,0.4300,0.24]
    cf[*,17,16]=[0., 23.8,1.,6.96E-08,0.6700,0.20]
    cf[*,17,15]=[0., 39.6,1.,3.40E-08,0.8650,0.18]
    cf[*,17,14]=[0., 53.5,1.,1.10E-08,0.3280,0.25]
    cf[*,17,13]=[0., 67.8,1.,1.11E-08,1.3700,0.10]
    cf[*,17,12]=[0., 97.0,1.,3.17E-09,0.3300,0.24]
    cf[*,17,11]=[0.,114.2,1.,1.01E-09,0.1960,0.16]
    cf[*,17,10]=[0.,348.3,0.,2.11E-09,0.3130,0.37]
    cf[*,17, 9]=[0.,400.1,1.,6.32E-10,0.1730,0.30]
    cf[*,17, 8]=[0.,455.6,0.,9.48E-10,0.3440,0.36]
    cf[*,17, 7]=[0.,529.3,1.,3.69E-10,0.2730,0.26]
    cf[*,17, 6]=[0.,592.0,1.,2.85E-10,0.3430,0.23]
    cf[*,17, 5]=[0.,656.7,0.,4.81E-10,0.6580,0.27]
    cf[*,17, 4]=[0.,749.8,1.,1.31E-10,0.6230,0.16]
    cf[*,17, 3]=[0.,809.4,1.,6.13E-11,0.7360,0.16]
    cf[*,17, 2]=[0.,3658.4,0.,1.90E-11,0.3790,0.36]
    cf[*,17, 1]=[0.,3946.3,1.,5.14E-12,0.5530,0.18]
    cf[*,18,18]=[0., 15.8,1.,5.99E-08,0.1360,0.26]
    cf[*,18,17]=[0., 27.6,1.,6.07E-08,0.5440,0.21]
    cf[*,18,16]=[0., 40.9,1.,3.43E-08,0.8340,0.17]
    cf[*,18,15]=[0., 52.3,0.,3.00E-08,1.0300,0.25]
    cf[*,18,14]=[0., 75.0,1.,8.73E-09,0.3660,0.31]
    cf[*,18,13]=[0., 91.0,1.,5.78E-09,0.3140,0.34]
    cf[*,18,12]=[0.,124.3,1.,2.98E-09,0.7030,0.16]
    cf[*,18,11]=[0.,143.5,1.,7.25E-10,0.2070,0.15]
    cf[*,18,10]=[0.,422.4,1.,1.40E-09,0.6960,0.13]
    cf[*,18, 9]=[0.,478.7,1.,4.78E-10,0.1640,0.31]
    cf[*,18, 8]=[0.,539.0,0.,8.02E-10,0.4390,0.32]
    cf[*,18, 7]=[0.,618.3,1.,2.88E-10,0.2590,0.27]
    cf[*,18, 6]=[0.,686.1,1.,2.32E-10,0.3620,0.22]
    cf[*,18, 5]=[0.,755.7,0.,3.33E-10,0.4120,0.36]
    cf[*,18, 4]=[0.,854.8,1.,1.27E-10,0.9100,0.13]
    cf[*,18, 3]=[0.,918.0,1.,5.21E-11,0.7810,0.15]
    cf[*,18, 2]=[0.,4120.7,0.,1.66E-11,0.4350,0.33]
    cf[*,18, 1]=[0.,4426.2,1.,4.32E-12,0.5540,0.18]
    cf[*,19,19]=[0.,  4.3,1.,2.02E-07,0.2720,0.31]
    cf[*,19,18]=[0., 31.6,1.,4.01E-08,0.3710,0.22]
    cf[*,19,17]=[0., 45.8,1.,1.50E-08,0.4330,0.21]
    cf[*,19,16]=[0., 60.9,1.,1.94E-08,0.8890,0.16]
    cf[*,19,15]=[0., 82.7,1.,6.95E-09,0.4940,0.18]
    cf[*,19,14]=[0., 99.4,1.,4.11E-09,0.5400,0.17]
    cf[*,19,13]=[0.,117.6,1.,2.23E-09,0.5190,0.16]
    cf[*,19,12]=[0.,154.7,1.,2.15E-09,0.8280,0.14]
    cf[*,19,11]=[0.,175.8,0.,1.61E-09,0.6420,0.13]
    cf[*,19,10]=[0.,504.0,1.,1.07E-09,0.6950,0.13]
    cf[*,19, 9]=[0.,564.7,1.,3.78E-10,0.1730,0.30]
    cf[*,19, 8]=[0.,629.4,0.,6.24E-10,0.4180,0.33]
    cf[*,19, 7]=[0.,714.6,1.,2.29E-10,0.2450,0.28]
    cf[*,19, 6]=[0.,786.6,1.,1.86E-10,0.3440,0.23]
    cf[*,19, 5]=[0.,861.1,0.,2.69E-10,0.3960,0.37]
    cf[*,19, 4]=[0.,968.0,1.,1.06E-10,0.9120,0.13]
    cf[*,19, 3]=[0.,1053.4,1.,4.24E-11,0.7370,0.16]
    cf[*,19, 2]=[0.,4610.9,0.,1.38E-11,0.4160,0.34]
    cf[*,19, 1]=[0.,4934.1,1.,3.67E-12,0.5550,0.18]
    cf[*,20,20]=[0.,  6.1,0.,4.40E-07,0.8480,0.33]
    cf[*,20,19]=[0., 11.9,0.,5.22E-08,0.1510,0.34]
    cf[*,20,18]=[0., 50.9,1.,2.06E-08,0.4180,0.20]
    cf[*,20,17]=[0., 67.3,1.,1.72E-08,0.6380,0.19]
    cf[*,20,16]=[0., 84.5,1.,1.26E-08,1.0100,0.14]
    cf[*,20,15]=[0.,108.8,1.,4.72E-09,0.5260,0.17]
    cf[*,20,14]=[0.,127.2,1.,2.89E-09,0.5480,0.17]
    cf[*,20,13]=[0.,147.2,1.,1.64E-09,0.5520,0.15]
    cf[*,20,12]=[0.,188.3,1.,1.57E-09,0.7990,0.14]
    cf[*,20,11]=[0.,211.3,1.,4.32E-10,0.2320,0.14]
    cf[*,20,10]=[0.,591.9,0.,9.47E-10,0.3110,0.38]
    cf[*,20, 9]=[0.,657.2,1.,2.98E-10,0.1630,0.31]
    cf[*,20, 8]=[0.,726.6,0.,4.78E-10,0.3590,0.36]
    cf[*,20, 7]=[0.,817.6,1.,1.86E-10,0.2440,0.28]
    cf[*,20, 6]=[0.,894.5,1.,1.56E-10,0.3640,0.22]
    cf[*,20, 5]=[0.,974.0,0.,2.16E-10,0.3570,0.39]
    cf[*,20, 4]=[0.,1087.0,1.,7.70E-11,0.6550,0.15]
    cf[*,20, 3]=[0.,1157.0,1.,3.58E-11,0.7360,0.16]
    cf[*,20, 2]=[0.,5128.9,0.,1.28E-11,0.5200,0.30]
    cf[*,20, 1]=[0.,5469.9,1.,3.08E-12,0.5280,0.19]
    cf[*,21,21]=[0.,  6.6,1.,3.16E-07,0.2040,0.28]
    cf[*,21,20]=[0., 12.8,1.,8.61E-08,0.1810,0.25]
    cf[*,21,19]=[0., 24.8,1.,5.08E-08,0.3570,0.24]
    cf[*,21,18]=[0., 73.5,1.,1.00E-08,0.4530,0.15]
    cf[*,21,17]=[0., 91.9,1.,6.76E-09,0.4600,0.15]
    cf[*,21,16]=[0.,110.7,1.,5.27E-09,0.5610,0.17]
    cf[*,21,15]=[0.,138.0,1.,3.40E-09,0.5600,0.16]
    cf[*,21,14]=[0.,158.1,1.,2.18E-09,0.6120,0.15]
    cf[*,21,13]=[0.,180.0,1.,1.26E-09,0.6100,0.14]
    cf[*,21,12]=[0.,225.1,1.,1.24E-09,0.8520,0.13]
    cf[*,21,11]=[0.,249.8,1.,3.62E-10,0.3490,0.05]
    cf[*,21,10]=[0.,687.4,1.,5.52E-10,0.3750,0.28]
    cf[*,21, 9]=[0.,756.7,1.,5.64E-10,0.8730,0.15]
    cf[*,21, 8]=[0.,830.8,1.,4.50E-10,1.0500,0.13]
    cf[*,21, 7]=[0.,927.5,1.,2.73E-10,0.8660,0.15]
    cf[*,21, 6]=[0.,1009.0,1.,1.56E-10,0.7150,0.17]
    cf[*,21, 5]=[0.,1094.0,0.,1.81E-10,1.1400,0.36]
    cf[*,21, 4]=[0.,1213.0,1.,4.29E-11,0.7840,0.15]
    cf[*,21, 3]=[0.,1288.0,0.,2.21E-11,0.0270,0.82]
    cf[*,21, 2]=[0.,5674.9,1.,4.51E-12,0.9180,0.04]
    cf[*,21, 1]=[0.,6033.8,0.,2.03E-12,0.0170,0.70]
    cf[*,22,22]=[0.,  6.8,1.,1.60E-07,0.3600,0.28]
    cf[*,22,21]=[0., 13.6,0.,2.14E-07,0.8800,0.28]
    cf[*,22,20]=[0., 27.5,1.,2.85E-08,0.2270,0.21]
    cf[*,22,19]=[0., 43.3,1.,3.48E-08,0.3900,0.23]
    cf[*,22,18]=[0., 99.3,1.,1.00E-08,0.5790,0.18]
    cf[*,22,17]=[0.,119.5,1.,7.01E-09,0.6380,0.17]
    cf[*,22,16]=[0.,140.8,1.,4.95E-09,0.7170,0.16]
    cf[*,22,15]=[0.,170.4,1.,2.99E-09,0.6930,0.17]
    cf[*,22,14]=[0.,192.1,1.,2.10E-09,0.7220,0.16]
    cf[*,22,13]=[0.,215.9,1.,1.62E-09,0.7650,0.14]
    cf[*,22,12]=[0.,265.0,1.,1.11E-09,0.8850,0.12]
    cf[*,22,11]=[0.,291.5,0.,9.09E-10,0.9720,0.06]
    cf[*,22,10]=[0.,787.8,1.,4.41E-10,0.3590,0.29]
    cf[*,22, 9]=[0.,863.1,1.,4.39E-10,0.7810,0.17]
    cf[*,22, 8]=[0.,941.9,1.,3.73E-10,1.0500,0.13]
    cf[*,22, 7]=[0.,1044.0,1.,2.28E-10,0.8580,0.15]
    cf[*,22, 6]=[0.,1131.0,1.,1.34E-10,0.7570,0.16]
    cf[*,22, 5]=[0.,1221.0,0.,1.55E-10,1.1500,0.36]
    cf[*,22, 4]=[0.,1346.0,1.,3.80E-11,0.8350,0.14]
    cf[*,22, 3]=[0.,1426.0,0.,1.89E-11,0.0280,0.82]
    cf[*,22, 2]=[0.,6249.1,1.,4.01E-12,0.9680,0.03]
    cf[*,22, 1]=[0.,6625.0,1.,1.62E-12,0.6570,0.14]
    cf[*,23,23]=[0.,  6.7,0.,8.82E-07,0.3590,0.32]
    cf[*,23,22]=[0., 14.7,0.,3.11E-07,0.4320,0.29]
    cf[*,23,21]=[0., 29.3,1.,3.50E-08,0.2470,0.25]
    cf[*,23,20]=[0., 46.7,0.,5.32E-08,1.1100,0.16]
    cf[*,23,19]=[0., 65.3,1.,8.98E-09,0.1400,0.37]
    cf[*,23,18]=[0.,128.1,1.,5.87E-09,0.5170,0.17]
    cf[*,23,17]=[0.,150.6,1.,5.11E-09,0.6790,0.16]
    cf[*,23,16]=[0.,173.4,1.,3.71E-09,0.7610,0.15]
    cf[*,23,15]=[0.,205.8,1.,2.24E-09,0.7110,0.17]
    cf[*,23,14]=[0.,230.5,1.,1.65E-09,0.7640,0.15]
    cf[*,23,13]=[0.,256.0,1.,1.26E-09,0.7620,0.14]
    cf[*,23,12]=[0.,308.0,1.,8.86E-10,0.8860,0.12]
    cf[*,23,11]=[0.,336.3,0.,3.89E-10,0.1420,0.39]
    cf[*,23,10]=[0.,896.0,1.,3.80E-10,0.4090,0.27]
    cf[*,23, 9]=[0.,976.0,0.,4.84E-10,0.1730,0.57]
    cf[*,23, 8]=[0.,1060.0,1.,2.49E-10,0.6500,0.14]
    cf[*,23, 7]=[0.,1168.0,0.,5.91E-10,1.6100,0.18]
    cf[*,23, 6]=[0.,1260.0,0.,5.02E-10,2.1200,0.15]
    cf[*,23, 5]=[0.,1355.0,1.,5.38E-11,0.1370,0.40]
    cf[*,23, 4]=[0.,1486.0,1.,5.56E-11,0.7080,0.10]
    cf[*,23, 3]=[0.,1571.0,0.,2.84E-11,0.0240,0.79]
    cf[*,23, 2]=[0.,6851.3,0.,2.54E-11,2.9200,0.09]
    cf[*,23, 1]=[0.,7246.1,0.,1.32E-11,3.5100,0.07]
    cf[*,24,24]=[0.,  6.8,1.,1.03E-07,0.2170,0.27]
    cf[*,24,23]=[0., 16.5,0.,2.45E-07,0.3810,0.32]
    cf[*,24,22]=[0., 31.0,0.,1.09E-07,0.5180,0.27]
    cf[*,24,21]=[0., 49.1,1.,1.52E-08,0.1820,0.30]
    cf[*,24,20]=[0., 69.5,0.,3.25E-08,1.3600,0.13]
    cf[*,24,19]=[0., 90.6,1.,5.50E-09,0.1430,0.37]
    cf[*,24,18]=[0.,160.2,1.,5.13E-09,0.6570,0.16]
    cf[*,24,17]=[0.,184.7,1.,3.85E-09,0.7220,0.15]
    cf[*,24,16]=[0.,209.3,1.,2.81E-09,0.7590,0.15]
    cf[*,24,15]=[0.,244.4,1.,1.76E-09,0.7320,0.16]
    cf[*,24,14]=[0.,271.0,1.,1.30E-09,0.7640,0.15]
    cf[*,24,13]=[0.,298.0,1.,1.02E-09,0.8100,0.13]
    cf[*,24,12]=[0.,354.8,1.,7.19E-10,0.8870,0.12]
    cf[*,24,11]=[0.,384.2,1.,1.61E-10,0.1500,0.22]
    cf[*,24,10]=[0.,1011.0,1.,4.64E-10,0.9710,0.12]
    cf[*,24, 9]=[0.,1097.0,1.,3.31E-10,0.9240,0.14]
    cf[*,24, 8]=[0.,1185.0,1.,2.49E-10,0.9310,0.15]
    cf[*,24, 7]=[0.,1299.0,1.,1.68E-10,0.9100,0.14]
    cf[*,24, 6]=[0.,1396.0,1.,1.01E-10,0.8050,0.15]
    cf[*,24, 5]=[0.,1496.0,0.,1.17E-10,1.2100,0.35]
    cf[*,24, 4]=[0.,1634.0,1.,2.91E-11,0.8840,0.13]
    cf[*,24, 3]=[0.,1721.0,0.,1.45E-11,0.0350,0.80]
    cf[*,24, 2]=[0.,7482.0,1.,3.07E-12,0.9670,0.03]
    cf[*,24, 1]=[0.,7894.8,1.,1.46E-12,0.1830,0.39]
    cf[*,25,25]=[0.,  7.4,1.,8.56E-08,0.1320,0.26]
    cf[*,25,24]=[0., 15.6,0.,1.18E-07,0.3590,0.19]
    cf[*,25,23]=[0., 33.7,0.,8.54E-08,0.3970,0.32]
    cf[*,25,22]=[0., 51.2,1.,1.80E-08,0.2720,0.18]
    cf[*,25,21]=[0., 72.4,1.,8.22E-09,0.1610,0.32]
    cf[*,25,20]=[0., 95.0,0.,2.15E-08,1.5400,0.11]
    cf[*,25,19]=[0.,119.3,1.,3.65E-09,0.1470,0.37]
    cf[*,25,18]=[0.,194.5,1.,3.91E-09,0.6990,0.15]
    cf[*,25,17]=[0.,221.8,1.,2.92E-09,0.7190,0.15]
    cf[*,25,16]=[0.,248.3,1.,2.23E-09,0.8060,0.14]
    cf[*,25,15]=[0.,286.0,1.,1.39E-09,0.7350,0.16]
    cf[*,25,14]=[0.,314.4,1.,1.04E-09,0.7610,0.15]
    cf[*,25,13]=[0.,343.6,1.,8.28E-10,0.8090,0.13]
    cf[*,25,12]=[0.,403.0,1.,5.60E-10,0.7870,0.14]
    cf[*,25,11]=[0.,435.2,1.,1.52E-10,0.2990,0.08]
    cf[*,25,10]=[0.,1133.0,1.,4.03E-10,1.0400,0.11]
    cf[*,25, 9]=[0.,1244.0,1.,2.74E-10,0.9230,0.14]
    cf[*,25, 8]=[0.,1317.0,1.,2.18E-10,0.9900,0.14]
    cf[*,25, 7]=[0.,1437.0,1.,1.49E-10,0.9680,0.13]
    cf[*,25, 6]=[0.,1539.0,1.,8.70E-11,0.8020,0.15]
    cf[*,25, 5]=[0.,1644.0,0.,1.02E-10,1.2200,0.35]
    cf[*,25, 4]=[0.,1788.0,1.,2.54E-11,0.8830,0.13]
    cf[*,25, 3]=[0.,1880.0,0.,1.28E-11,0.0330,0.81]
    cf[*,25, 2]=[0.,8141.0,1.,2.77E-12,1.0100,0.02]
    cf[*,25, 1]=[0.,8571.9,1.,1.32E-12,0.2190,0.37]
    cf[*,26,26]=[0.,  7.9,0.,2.52E-07,0.7010,0.25]
    cf[*,26,25]=[0., 16.2,1.,2.21E-08,0.0330,0.45]
    cf[*,26,24]=[0., 30.6,0.,4.10E-08,0.3660,0.17]
    cf[*,26,23]=[0., 54.8,0.,3.53E-08,0.2430,0.39]
    cf[*,26,22]=[0., 75.0,1.,1.04E-08,0.2850,0.17]
    cf[*,26,21]=[0., 99.0,1.,1.23E-08,0.4110,0.21]
    cf[*,26,20]=[0.,125.0,1.,9.47E-09,0.4580,0.21]
    cf[*,26,19]=[0.,151.1,1.,4.71E-09,0.2800,0.28]
    cf[*,26,18]=[0.,233.6,1.,3.02E-09,0.6970,0.15]
    cf[*,26,17]=[0.,262.1,1.,2.34E-09,0.7640,0.14]
    cf[*,26,16]=[0.,290.0,1.,1.76E-09,0.8050,0.14]
    cf[*,26,15]=[0.,331.0,1.,1.14E-09,0.7730,0.15]
    cf[*,26,14]=[0.,361.0,1.,8.66E-10,0.8050,0.14]
    cf[*,26,13]=[0.,392.0,1.,6.61E-10,0.7620,0.14]
    cf[*,26,12]=[0.,457.0,1.,4.41E-10,0.6980,0.16]
    cf[*,26,11]=[0.,489.3,1.,1.18E-10,0.2110,0.15]
    cf[*,26,10]=[0.,1262.0,1.,3.61E-10,1.1600,0.09]
    cf[*,26, 9]=[0.,1360.0,1.,2.45E-10,0.9780,0.13]
    cf[*,26, 8]=[0.,1470.0,1.,1.87E-10,0.9880,0.14]
    cf[*,26, 7]=[0.,1582.0,1.,1.33E-10,1.0300,0.12]
    cf[*,26, 6]=[0.,1690.0,1.,7.84E-11,0.8480,0.14]
    cf[*,26, 5]=[0.,1800.0,0.,8.90E-11,1.2000,0.35]
    cf[*,26, 4]=[0.,1960.0,1.,2.29E-11,0.9360,0.12]
    cf[*,26, 3]=[0.,2046.0,0.,1.12E-11,0.0340,0.81]
    cf[*,26, 2]=[0.,8828.0,1.,2.46E-12,1.0200,0.02]
    cf[*,26, 1]=[0.,9277.7,1.,9.79E-13,0.6640,0.14]
    cf[*,27,27]=[0.,  7.9,1.,8.89E-08,0.1270,0.24]
    cf[*,27,26]=[0., 17.1,1.,5.65E-08,0.1940,0.23]
    cf[*,27,25]=[0., 33.5,1.,3.06E-08,0.2010,0.22]
    cf[*,27,24]=[0., 51.3,0.,2.27E-08,0.5740,0.10]
    cf[*,27,23]=[0., 79.5,0.,1.93E-08,0.1950,0.42]
    cf[*,27,22]=[0.,102.0,0.,1.27E-08,0.1260,0.47]
    cf[*,27,21]=[0.,129.0,1.,3.58E-09,0.1940,0.29]
    cf[*,27,20]=[0.,158.0,0.,1.17E-08,1.9800,0.07]
    cf[*,27,19]=[0.,186.1,1.,1.78E-09,0.1120,0.42]
    cf[*,27,18]=[0.,275.0,1.,2.41E-09,0.7390,0.14]
    cf[*,27,17]=[0.,305.0,1.,1.86E-09,0.7620,0.14]
    cf[*,27,16]=[0.,336.0,1.,1.41E-09,0.8040,0.14]
    cf[*,27,15]=[0.,379.0,1.,9.54E-10,0.8130,0.14]
    cf[*,27,14]=[0.,411.0,1.,7.12E-10,0.8030,0.14]
    cf[*,27,13]=[0.,444.0,1.,5.34E-10,0.7180,0.15]
    cf[*,27,12]=[0.,512.0,1.,3.62E-10,0.6580,0.17]
    cf[*,27,11]=[0.,546.6,1.,1.05E-10,0.2520,0.12]
    cf[*,27,10]=[0.,1397.0,1.,3.10E-10,1.1700,0.09]
    cf[*,27, 9]=[0.,1486.0,1.,1.56E-10,0.5720,0.15]
    cf[*,27, 8]=[0.,1603.0,1.,1.32E-10,0.6820,0.13]
    cf[*,27, 7]=[0.,1735.0,1.,9.08E-11,0.5110,0.17]
    cf[*,27, 6]=[0.,1846.0,0.,3.45E-10,2.8400,0.11]
    cf[*,27, 5]=[0.,1962.0,1.,5.01E-11,0.7140,0.11]
    cf[*,27, 4]=[0.,2119.0,1.,1.92E-11,0.1170,0.42]
    cf[*,27, 3]=[0.,2219.0,1.,1.95E-11,1.2000,0.00]
    cf[*,27, 2]=[0.,9544.0,0.,1.68E-11,3.5200,0.05]
    cf[*,27, 1]=[0.,10012.1,1.,1.45E-12,0.6350,0.15]
    cf[*,28,28]=[0.,  7.6,0.,1.65E-07,0.4520,0.28]
    cf[*,28,27]=[0., 18.2,0.,8.42E-08,0.6190,0.16]
    cf[*,28,26]=[0., 35.3,1.,1.89E-08,0.2200,0.21]
    cf[*,28,25]=[0., 54.9,1.,1.48E-08,0.2160,0.21]
    cf[*,28,24]=[0., 76.0,0.,1.13E-08,0.5180,0.09]
    cf[*,28,23]=[0.,108.0,0.,1.16E-08,0.1850,0.44]
    cf[*,28,22]=[0.,133.0,0.,8.68E-09,0.1380,0.46]
    cf[*,28,21]=[0.,162.0,1.,2.45E-09,0.1630,0.32]
    cf[*,28,20]=[0.,193.0,0.,9.24E-09,2.2500,0.05]
    cf[*,28,19]=[0.,225.0,0.,2.41E-09,0.0270,0.79]
    cf[*,28,18]=[0.,321.0,1.,1.92E-09,0.7380,0.14]
    cf[*,28,17]=[0.,352.0,1.,1.50E-09,0.7610,0.14]
    cf[*,28,16]=[0.,384.0,1.,1.16E-09,0.8030,0.14]
    cf[*,28,15]=[0.,430.0,1.,8.08E-10,0.8560,0.13]
    cf[*,28,14]=[0.,464.0,1.,6.09E-10,0.8500,0.13]
    cf[*,28,13]=[0.,499.0,1.,4.48E-10,0.7180,0.15]
    cf[*,28,12]=[0.,571.3,1.,3.00E-10,0.6220,0.18]
    cf[*,28,11]=[0.,607.0,1.,7.90E-11,0.1600,0.19]
    cf[*,28,10]=[0.,1541.0,1.,2.78E-10,1.2500,0.08]
    cf[*,28, 9]=[0.,1648.0,1.,1.92E-10,1.0400,0.12]
    cf[*,28, 8]=[0.,1756.0,1.,1.51E-10,1.1100,0.12]
    cf[*,28, 7]=[0.,1894.0,1.,1.05E-10,1.0900,0.11]
    cf[*,28, 6]=[0.,2011.0,1.,6.04E-11,0.8490,0.14]
    cf[*,28, 5]=[0.,2131.0,0.,6.91E-11,1.2100,0.35]
    cf[*,28, 4]=[0.,2295.0,1.,1.84E-11,0.9910,0.11]
    cf[*,28, 3]=[0.,2399.0,0.,9.03E-12,0.0420,0.79]
    cf[*,28, 2]=[0.,10290.0,1.,2.61E-12,0.5680,0.16]
    cf[*,28, 1]=[0.,10775.3,1.,1.39E-12,0.7290,0.13]
endif

if iz lt 1 or iz gt 28 then return,-1
if in lt 1 or in gt iz then return,-1

u=cf[1,iz,in]/te

if max(u) gt 80 then return,-2
c=cf[3,iz,in]*(1.0+cf[2,iz,in]*sqrt(u))/(cf[4,iz,in]+u)*u^cf[5,iz,in]*exp(-u)

return,c
end




