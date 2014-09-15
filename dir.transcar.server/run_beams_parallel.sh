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

remotes=(labHST0 labHST1)

#-------------
remjp=$(IFS=,; echo "${remotes[*]}") #puts array into comma separated string for GNU parallel

# purge output directory
[[ -d $RODIR ]] && \rm -r $RODIR
for remote in "${remotes[@]}"; do
    ssh $remote -t "[[ -d $exedir/$RODIR ]] && rm -r $exedir/$RODIR"
done

#-S labHST0,labHST1
# jobs is equal to number of CPU cores by default
# note --cleanup doesn't work with parallel 20130922 b/c we're returning a whole directory tree (no rm -r in --cleanup)
nice parallel \
    -S $remjp --return $RODIR \
    --nice 18 --halt 2 --eta --progress --joblog parallellog --colsep ',' \
    --workdir $exedir \
    "./beamRunner.sh" $RODIR :::: $BeamEnergyTableFN

#-- check results for proper simulation finish
# not possible using find/tail/grep without for loop -- would need gawk -- easier to do this
for f in $(find $RODIR -mindepth 2 -maxdepth 2 -type f -name "TranscarErrors.txt"); do
outcome=$(tail -n1 $f)
[[ $outcome != *fin\ normale ]] && echo "abnormal completion in $f"
done



