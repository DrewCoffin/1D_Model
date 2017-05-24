import os
import sys

def readval(sp, param):
   file = open('plots/data/' + sp + '/' + param + '/' + param + sp + '0200.dat')
   line = file.readline().strip()
   file.close()
   llist = line.split(' ')
   lval = llist.pop()
   return float(lval)

colval = [0.1126, 240.44, 0.2110, 188.3, 0.0123, 160.08, 0.2962, 160.63, 0.0161, 142.53, 721.23]
sparr = ['sp', 's2p', 's3p', 'op', 'o2p']
pararr = ['MIXR', 'TEMP']
val = []

for i in range(len(sparr)):
	for j in range(len(pararr)):
		val.append(str(sparr[i] + pararr[j]))
		val.append(readval(sparr[i], pararr[j]))

val.append('elecDENS')
val.append(readval('elec','DENS'))

labels = [val[2*k] for k in range(len(colval))]
diff = [(val[2*k+1] - colval[k])/colval[k] for k in range(len(colval))]

difflist = sum(zip(labels, diff), ())

print(val)
print(difflist)
