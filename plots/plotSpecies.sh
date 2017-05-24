#!/bin/bash

FILE='test'
FILELIST='list.dat'
PREFIX=$2
species=$3

echo '' > $FILELIST

echo 'set terminal jpeg'                                      >$FILE
#echo 'set logscale y'                                        >>$FILE
echo 'set key below'                                         >>$FILE
echo 'set size 1, 1'                                         >>$FILE
echo "set xlabel 'System III Longitude'"                     >>$FILE
echo "set ylabel '$PREFIX'"                                  >>$FILE
#echo "set yrange [0:500]"                                    >>$FILE
echo "set xrange [0:360]"                                    >>$FILE
echo "set xtics 45"                                          >>$FILE
echo "set grid ytics"                                        >>$FILE
echo "set yrange [0:300]"                                    >>$FILE

#if [ $2 -eq "DENS" ]
#then
#  echo "This part works"
#  echo "set yrange [0:500]"                                    >>$FILE
#  echo "set ytics 50"                                         >>$FILE
#fi

for i in $(seq 1 $1)
do
  
  if [ "$i" -lt 10 ]
  then   
    echo "set title 'Temporal Variability of Flux Tube Content (00$i)'" >>$FILE
    echo "set output 'gplot00$i.jpeg'" >> $FILE  
    echo "gplot00$i.jpeg" >> $FILELIST  
  
    echo "plot '"$PREFIX"sp000$i.dat'  with lines title 'Sulfur + (II)',  " '\' >> $FILE
    echo "     '"$PREFIX"s3p000$i.dat'  with lines title 'Sulfur +++ (IV)',   " '\'  >> $FILE
    echo "     'LOAD000$i.dat'  with lines title 'Normalized Loading'   "  >> $FILE
  else
    if [ "$i" -lt 100 ] 
    then
      echo "set title 'Temporal Variability of Flux Tube Content (0$i)'" >>$FILE
      echo "set output 'gplot0$i.jpeg'" >> $FILE  
      echo "gplot0$i.jpeg" >> $FILELIST  
  
      echo "plot '"$PREFIX"sp00$i.dat'  with lines title 'Sulfur + (II)',  " '\'  >> $FILE
      echo "     '"$PREFIX"s3p00$i.dat'  with lines title 'Sulfur +++ (IV)',   " '\'  >> $FILE
      echo "     'LOAD00$i.dat'  with lines title 'Normalized Loading' ,  "   >> $FILE
    else
      echo "set title 'Temporal Variability of Flux Tube Content ($i)'" >>$FILE
      echo "set output 'gplot$i.jpeg'" >> $FILE  
      echo "gplot$i.jpeg" >> $FILELIST  
    
      echo "plot '"$PREFIX"sp0$i.dat'  with lines title 'Sulfur + (II)',  " '\'  >> $FILE
      echo "     '"$PREFIX"s3p0$i.dat'  with lines title 'Sulfur +++ (IV)',  " '\'  >> $FILE
      echo "     'LOAD0$i.dat'  with lines title 'Normalized Loading'   "  >> $FILE
    fi
  fi


done

gnuplot $FILE

#convert -delay 10 gplot* animated.gif
mencoder -nosound -really-quiet -ovc lavc -lavcopts \
  vcodec=mpeg4:mbd=2:trell:autoaspect:vqscale=3 \
  -vf scale=768:432 -mf type=jpeg:fps=10 \
  mf://@list.dat -o animated.avi

rm gplot*.jpeg
