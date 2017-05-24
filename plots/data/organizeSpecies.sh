#!/bin/bash

species=$1

mkdir "$species"
good=$?
if [ $good -ne 0 ]
  then
  echo "'$species' folder already exists. Please remove or change name to procede." 
fi

if [ $good -eq 0 ]
  then
  mkdir "$species"/DENS
  mv DENS"$species"*.dat "$species"/DENS/.
  mkdir "$species"/MIXR
  mv MIXR"$species"*.dat "$species"/MIXR/.
  mkdir "$species"/TEMP
  mv TEMP"$species"*.dat "$species"/TEMP/.
  mkdir "$species"/INTS
  mv INTS"$species"*.dat "$species"/INTS/.
fi

