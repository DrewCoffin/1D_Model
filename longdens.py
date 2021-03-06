#Read data files
import os
import sys
import numpy as np
import math as m
import matplotlib.pyplot as plt

def numzeros(i):
     if i < 1000:
          base =  m.log(i)/m.log(10)
          return 3 - int(base)
     else: 
          return 0

def getdens(output, days, spec, runloc, dim): #finds the max density per day
     denvals = []
     for i in days:
          dayden = []
          path = './' + runloc + '/' + spec + '/DENS/DENS' + spec + numzeros(i)*'0' + str(i) + '.dat'
          with open(path) as datafile:
               data = [next(datafile) for x in xrange(dim)]  
               datalines = [data[i].split() for i in range(len(data))]
               for j in range(len(datalines)):
                    dayden.append(float(datalines[j][1])) 
          denvals.append(max(dayden))
     denvals = np.array(denvals)
     return denvals  

lng = 24
maxday = 3000
#print lng, rad
days = range(1, maxday)
output = []
specs = ['sp', 's2p', 's3p', 'op', 'o2p', 'elec']
output = [getdens(output, days, specs[i], 'plots/data', lng) for i in range(len(specs))]

#Plot max density array

plt.subplot(231)
plt.title('S+ density')
#plt.xlabel('Days')
plt.ylabel('Peak density')
plt.semilogx(days, output[0])

plt.subplot(232)
plt.title('S++ density')
#plt.xlabel('Days')
#plt.ylabel('Peak density')
plt.semilogx(days, output[1])

plt.subplot(233)
plt.title('S+++ density')
#plt.xlabel('Days')
#plt.ylabel('Peak density')
plt.semilogx(days, output[2])

plt.subplot(234)
plt.title('O+ density')
plt.xlabel('Days')
plt.ylabel('Peak density')
plt.semilogx(days, output[3])

plt.subplot(235)
plt.title('O++ density')
plt.xlabel('Days')
#plt.ylabel('Peak density')
plt.semilogx(days, output[4])

plt.subplot(236)
plt.title('e- density')
plt.xlabel('Days')
#plt.ylabel('Peak density')
plt.semilogx(days, output[5])

plt.show()
