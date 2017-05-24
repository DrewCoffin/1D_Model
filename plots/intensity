#!/bin/bash

FILE='test'

echo 'set terminal jpeg'                                      >$FILE
echo 'set key off'                                           >>$FILE
#echo 'set logscale y'                                       >>$FILE
echo 'set size 1, 1'                                         >>$FILE
echo "set xlabel 'System III Longitude'"                     >>$FILE
echo "set ylabel 'Intenisty'"                                >>$FILE
echo "set yrange [0.9:1.1]"                                  >>$FILE
echo "set xrange [0:360]"                                    >>$FILE
echo "set xtics 45"                                          >>$FILE
echo "set output 'intensityOverlay.jpeg'" >> $FILE  

for i in $(seq 1 $1)
do
  
  if [ "$i" -lt 10 ]
  then   
#    echo "set title 'Flux Tube Emission intensity (000$i)'" >>$FILE
#    echo "set output 'intense000$i.jpeg'" >> $FILE  
  
    echo "plot 'intensity000$i.dat'  with lines"  >> $FILE
  else
    if [ "$i" -lt 100 ] 
    then
#      echo "set title 'Flux Tube Emission intensity (00$i)'" >>$FILE
#      echo "set output 'intense00$i.jpeg'" >> $FILE  
  
      echo "plot 'intensity00$i.dat'  with lines"  >> $FILE
    else
      if [ "$i" -lt 1000 ] 
      then
#        echo "set title 'Flux Tube Emission intensity (0$i)'" >>$FILE
#        echo "set output 'intense0$i.jpeg'" >> $FILE  
  
        echo "plot 'intensity0$i.dat'  with lines"  >> $FILE
      else
#        echo "set title 'Flux Tube Emission intensity ($i)'" >>$FILE
#        echo "set output 'intense$i.jpeg'" >> $FILE  
  
        echo "plot 'intensity$i.dat'  with lines"  >> $FILE
      fi
    fi
  fi


done

gnuplot $FILE

rm intense*.jpeg
rm intensity*.dat

#convert -delay 15 intense* animated.gif

