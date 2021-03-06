#!/bin/bash

npes=24
days=5000 

./changeDimension.sh $npes

rm FFT.dat

make all

if [ $? -eq 0 ] 
  then 

  time mpirun -n $npes ./torus > runlog

  if [ $? -ge 0 ] 
    then 

    mv DENS*.dat plots/.  
    mv MIXR*.dat plots/.  
    mv TEMP*.dat plots/.  
#    mv INTS*.dat plots/.  
    mv LOAD*.dat plots/.  
 #   mv intensity*.dat plots/.  

    cd plots

#      ./plotifyWload $days DENS
#      mv animated.avi ../dens.avi
#      ./plotify $days MIXR
#      mv animated.avi ../mixr.avi
#      ./plotify $days TEMP	
#      mv animated.avi ../temp.avi
#      ./plotify $days INTS
#      mv animated.avi ../intensity.avi
#      ./overlay $days 
#      mv overlay.jpeg ../.
#      mv animated.gif ../dens.gif
#      ./plotify 150 TEMP
#      mv animated.gif ../temp.gif

      ./plotSpecies.sh $days DENS s3p
      mv animated.avi ../SPonly.avi

#      ./plotify $days MIXR
#      mv animated.avi ../mixr.avi

      mv DENS*.dat data/.  
      mv MIXR*.dat data/.  
      mv TEMP*.dat data/.  
#      mv INTS*.dat data/.  
      mv LOAD*.dat data/.  
#      mv intensity*.dat data/intensity/.

      cd data
        ./organize.sh    
#        ./peakWload.sh $days
        ./peakPlot.sh $days
        ./plotRatio.sh $days
        cp peaks.jpeg ../../.
        cp peakRatio.jpeg ../../.
      cd ..

    cd ..

  fi

fi

make clean

#python colo.py

#gifview -a dens.gif &
#mplayer -loop -0 dens.avi &
#vlc dens.avi &
#vlc intensity.avi &
#mplayer dens.avi &
#display peaks.jpeg 
#display peakRatio.jpeg  
#display overlay.jpeg



