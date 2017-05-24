#!/bin/bash

n=$( printf "%.0f" $1 )

FILE='test'

echo 'set terminal jpeg'                                      >$FILE
echo 'set key off'                                           >>$FILE
#echo 'set logscale y'                                       >>$FILE
echo 'set size 1, 1'                                         >>$FILE
echo "set xlabel 'System III Longitude'"                     >>$FILE
echo "set ylabel 'Intenisty'"                                >>$FILE
#echo "set yrange [0.5:1.5]"                                >>$FILE
echo "set xrange [0:360]"                                    >>$FILE
echo "set xtics 45"                                          >>$FILE
echo "set title 'Flux Tube Emission intensity'"              >>$FILE
echo "set output 'overlay.jpeg'" >> $FILE  

echo "plot 'intensity0000.dat'  with lines," '\'   >> $FILE

notlast()
{
  if [ "$i" -lt 10 ]
  then   
    echo "   'intensity000$i.dat'  with lines,"  '\'  >> $FILE
  else
    if [ "$i" -lt 100 ] 
    then
      echo "   'intensity00$i.dat'  with lines," '\'  >> $FILE
    else
      if [ "$i" -lt 1000 ] 
      then
        echo "   'intensity0$i.dat'  with lines," '\' >> $FILE
      else
        echo "   'intensity$i.dat'  with lines,"  '\' >> $FILE
      fi
    fi
  fi
}
  
last()
{
  if [ "$n" -lt 10 ]
  then   
    echo "   'intensity000$n.dat'  with lines"    >> $FILE
  else
    if [ "$n" -lt 100 ] 
    then
      echo "   'intensity00$n.dat'  with lines"   >> $FILE
    else
      if [ "$n" -lt 1000 ] 
      then
        echo "   'intensity0$n.dat'  with lines"  >> $FILE
      else
        echo "   'intensity$n.dat'  with lines"   >> $FILE
      fi
    fi
  fi
}

for i in $(seq 1 $n)
do
  if [ "$i" -eq "$n" ]
  then
    last
  else
    notlast
  fi

done

gnuplot $FILE

