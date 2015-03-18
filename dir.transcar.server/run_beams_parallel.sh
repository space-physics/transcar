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
RODIR=../conttanh0
exedir=code/transcar/dir.transcar.server
localonly=1
remotes=(labHST0 labHST1)
hact=0   # 2 term all jobs on error now, 1 let existing jobs finish, 0 keep running/starting 

#------- start code --------------
if [[ $localonly -eq 0 ]]; then
  remjp=$(IFS=,; echo "${remotes[*]}") #puts array into comma separated string for GNU parallel

  # purge output directory
  [[ -d $RODIR ]] && \rm -r $RODIR #local
  for remote in "${remotes[@]}"; do #remote
      ssh $remote -t "[[ -d $exedir/$RODIR ]] && rm -r $exedir/$RODIR"
  done

# jobs is equal to number of CPU cores by default
# note --cleanup doesn't work with parallel 20130922 b/c we're returning a whole directory tree (no rm -r in --cleanup)
 
  parallel \
    -S $remjp --return $RODIR \
    --nice 18 --halt $hact --eta --progress --joblog parallellog --colsep ',' \
    --workdir $exedir \
    "./beamRunner.sh" $RODIR :::: $BeamEnergyTableFN

else #local only
  parallel \
    --nice 18 --halt $hact --eta --progress --joblog parallellog --colsep ',' \
    --workdir $exedir \
    "./beamRunner.sh" $RODIR :::: $BeamEnergyTableFN

fi
#-- check results for proper simulation finish
./checkoutcome.sh $RODIR



