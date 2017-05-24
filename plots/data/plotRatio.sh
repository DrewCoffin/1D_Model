#!/bin/bash

python findTwoPeaks.py sp $1
python findTwoPeaks.py s3p $1

gnuplot gplotRatio
