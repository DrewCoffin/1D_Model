import os
import sys
import shutil

os.chdir(str(sys.argv[1]))
os.chdir("DENS")

output = str(sys.argv[1])+"PeakRatio.dat"
out = open(output,'w')

for i in range(int(sys.argv[2])):
  if(i>99):  file = "DENS"+str(sys.argv[1])+"0"+str(i)+".dat" 
  if(i>9 and i<99):  file = "DENS"+str(sys.argv[1])+"00"+str(i)+".dat" 
  if(i<9):  file = "DENS"+str(sys.argv[1])+"000"+str(i)+".dat" 
  f = open(file)
  l = []
  l = [ line.split() for line in f]
  l=zip(*l)
  l[1]=[float(x) for x in l[1]]
  l[0]=[float(x) for x in l[0]]
  a=max(l[1])
  bigMax=a
  a=l[1].index(a)
  dex=a
  a=l[0][a]
  concavity=1
  found=1
  old=bigMax
  smallMax=0
  smalldex=0
  oldChange=0
  for j in range(22):
    change=l[1][(dex+j+1)%23]-l[1][(dex+j)%23]
    if(concavity==1 and found==1 and change > oldChange): concavity=0
    if(concavity==0 and found==1 and l[1][(dex+j+1)%23] < old):
      smallMax=old
      smalldex=(dex+j)%23
      found=0   
    old = l[1][(dex+j+1)%23]
  place=l[0][smalldex]
  minimum=min(l[1])
  out.write(str(i)+" "+str((smallMax-minimum)/(bigMax-minimum))+" "+str(place)+'\n')
  
src= str(sys.argv[1])+"PeakRatio.dat"
dest= "../../."
os.remove("../../"+src)
shutil.move(src,dest)
