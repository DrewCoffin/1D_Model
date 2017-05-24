import os
import sys

quantity="DENS"

os.chdir(str(sys.argv[1]))
os.chdir(quantity)

output = str(sys.argv[1])+"Peak.dat"
out = open(output,'w')

for i in range(int(sys.argv[2])):
  if(i>99):  file = quantity+str(sys.argv[1])+"0"+str(i)+".dat" 
  if(i>9 and i<99):  file = quantity+str(sys.argv[1])+"00"+str(i)+".dat" 
  if(i<9):  file = quantity+str(sys.argv[1])+"000"+str(i)+".dat" 
  f = open(file)
  l = []
  l = [ line.split() for line in f]
  l=zip(*l)
  l[1]=[float(x) for x in l[1]]
  l[0]=[float(x) for x in l[0]]
  a=max(l[1])
  a=l[1].index(a)
  a=l[0][a]
  out.write(str(i)+" "+str(a)+'\n')
    

