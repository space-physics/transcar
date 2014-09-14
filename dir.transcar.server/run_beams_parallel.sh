#!/bin/bash
# Michael Hirsch 2014
# using GNU parallel 20130922
# this script for loops transcar, making a new precinput.dat each time for the 
# respective beam energies.
# this program is meant to be run from dir.transcar.server directory- cd there first
#
# USAGE HINTS:
# 0) be sure the remotes have the appropriate directory structure for dir.transcar.server (use git)
# 1) consider using ssh-add with -t more than long enough to complete your calculation, 
# since a new SSH login is emitted for each beam! The simulation process will not
# continue for new beams if your ssh-add -t has expired!
# 2) when simulation is completed, the "return" argument of parallel uses rsync to 
# bring back simulation output to your PC.

BeamEnergyTableFN=BT_E1E2prev.csv
RODIR=../matt2013local
exedir=transcar/dir.transcar.server

# purge output directory
[[ -d $RODIR ]] && \rm -r $RODIR
ssh labHST0 -t "[[ -d $exedir/$RODIR ]] && rm -r $exedir/$RODIR"
ssh labHST1 -t "[[ -d $exedir/$RODIR ]] && rm -r $exedir/$RODIR"

#-S labHST0,labHST1
# jobs is equal to number of CPU cores by default
# note --cleanup doesn't work with parallel 20130922 b/c we're returning a whole directory tree (no rm -r in --cleanup)
parallel -S labHST0,labHST1 --return $RODIR \
    --nice 18 --halt 2 --eta --progress --joblog parallellog --colsep ',' \
    --workdir $exedir \
    "./beamRunner.sh" $RODIR :::: $BeamEnergyTableFN
