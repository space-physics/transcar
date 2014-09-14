#!/bin/bash
# upgraded to Bash 4 by Michael Hirsch 2014
# original by Matt Zettergren 2013
# this script for loops transcar, making a new precinput.dat each time for the 
# respective beam energies.
# this program is meant to be run from dir.transcar.server directory- cd there first
#
# USAGE HINTS:
# 0) because I don't yet use the rsync --transfer, --cleanup option of GNU parallel
# you must ensure the remote nodes have the latest version of your code. I handle this
# with git
# 1) consider using ssh-add with -t more than long enough to complete your calculation, 
# since a new SSH login is emitted for each beam! The simulation process will not
# continue for new beams if your ssh-add -t has expired!
# 2) when simulation is completed, I use SFTP to gather the results to my main 
# analysis PC
# so you see that I could skip steps 0 and 2 by using parallel rsync, but I haven't
# tested this yet and wanted to avoid a mishap.

BeamEnergyTableFN=BT_E1E2prev.csv

# jobs is equal to number of CPU cores by default
parallel -S labHST0 -S labHST1 \
    --eta --progress --joblog parallellog --colsep ',' transcar/dir.transcar.server/beamRunner.sh :::: $BeamEnergyTableFN
