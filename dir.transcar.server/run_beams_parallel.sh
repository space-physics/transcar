#!/bin/bash
# upgraded to Bash 4 by Michael Hirsch 2014
# original by Matt Zettergren 2013
# this script for loops transcar, making a new precinput.dat each time for the 
# respective beam energies.
# this program is meant to be run from dir.transcar.server directory- cd there first

BeamEnergyTableFN=BT_E1E2prev.csv

teea ()
{
tee --append "$1"
}


export -f teea
# jobs is equal to number of CPU cores by default
parallel -S labHST0 -S labHST1 \
    --eta --progress --joblog parallellog --colsep ',' ./beamRunner.sh :::: $BeamEnergyTableFN
