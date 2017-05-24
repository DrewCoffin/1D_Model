#!/bin/bash

FILE='test'

echo 'set terminal jpeg'                                      >$FILE
echo 'set key off'                                           >>$FILE
#echo 'set logscale y'                                       >>$FILE
echo 'set size 1, 1'                                         >>$FILE
echo "set xlabel 'System III Longitude'"                     >>$FILE
echo "set ylabel 'Electron Density'"                         >>$FILE
echo "set yrange [1000:2000]"                                  >>$FILE
echo "set xrange [0:360]"                                    >>$FILE
echo "set xtics 45"                                          >>$FILE

for i in $(seq 1 $1)
do
  
  if [ "$i" -lt 10 ]
  then   
    echo "set title 'Electron Density (000$i)'" >>$FILE
    echo "set output 'elecdens000$i.jpeg'" >> $FILE  
  
    echo "plot 'elec000$i.dat'  with lines"  >> $FILE
  else
    if [ "$i" -lt 100 ] 
    then
      echo "set title 'Electron Density (00$i)'" >>$FILE
      echo "set output 'elecdens00$i.jpeg'" >> $FILE  
  
      echo "plot 'elec00$i.dat'  with lines"  >> $FILE
    else
      if [ "$i" -lt 1000 ] 
      then
        echo "set title 'Elecron Density (0$i)'" >>$FILE
        echo "set output 'elecdens0$i.jpeg'" >> $FILE  
  
        echo "plot 'elec0$i.dat'  with lines"  >> $FILE
      else
        echo "set title 'Electron Density ($i)'" >>$FILE
        echo "set output 'elecdens$i.jpeg'" >> $FILE  
  
        echo "plot 'elec$i.dat'  with lines"  >> $FILE
      fi
    fi
  fi


done

gnuplot $FILE

convert -delay 15 elecdens* animated.gif

