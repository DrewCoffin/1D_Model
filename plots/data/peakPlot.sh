#!/bin/bash

python findPeaks.py sp $1
mv sp/DENS/spPeak.dat ./.
python findPeaks.py s3p $1
mv s3p/DENS/s3pPeak.dat ./.

old=200

sed -i "s/xrange\[0:$old\]/xrange\[0:$1\]/g" gplot
sed -i 's/xrange \[0:'$old'\]/xrange \[0:'$1'\]/g' gplotRatio
sed -i 's/old='$old'/old='$1'/g' peakPlot.sh

gnuplot gplot
