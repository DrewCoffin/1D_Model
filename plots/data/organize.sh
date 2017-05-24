#!/bin/bash

rm -r sp/
rm -r s2p/
rm -r s3p/
rm -r op/
rm -r o2p/
rm -r elec/

./organizeSpecies.sh sp
./organizeSpecies.sh s2p
./organizeSpecies.sh s3p
./organizeSpecies.sh op
./organizeSpecies.sh o2p
./organizeSpecies.sh elec
