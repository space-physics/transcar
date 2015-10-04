#!/bin/bash
# Michael Hirsch 2014
# using GNU parallel 20130922
# this script for loops transcar, making a new precinput.dat each time for the
# respective beam energies.
# this program is meant to be run from transcar/transcar directory- cd there first
#
# USAGE HINTS:
# 0) be sure the remotes have the appropriate directory structure for transcar/transcar (use git)
# 1) consider using ssh-add with -t more than long enough to complete your calculation,
# since a new SSH login is emitted for each beam. The simulation process will not
# continue for new beams if your ssh-add -t has expired.
# 2) when simulation is completed, the "return" argument of parallel uses rsync to
# bring back simulation output to your PC.

localonly=1
remotes=(irs4 irs3 swoboj)

BeamEnergyTableFN=BT_E1E2prev.csv
RODIR=$1
[[ -z $RODIR ]] && { echo "you must specify an output directory"; exit 1; }
exedir=code/transcar/transcar

flux0=70114000000.0

hact=0   # 2 term all jobs on error now, 1 let existing jobs finish, 0 keep running/starting

#------- start code --------------
[[ -d $RODIR ]] && \rm -r $RODIR  #cleanup local output

if [[ $localonly -eq 0 ]]; then
  remjp=$(IFS=,; echo "${remotes[*]}") #puts array into comma separated string for GNU parallel

  # purge remote output directory and setup ssh agent for duration of sims
  for remote in "${remotes[@]}"; do
    ssh-add -t 7200 "$HOME/.ssh/$remote" #use ssh agent so as to not have to retype password
    ssh $remote -t "(cd $exedir && git pull && cd dir.source && make -s && cd; cd code/transcar-utils && git pull && cd ../hist-utils && git pull)"
    ssh $remote -t "[[ -d $exedir/$RODIR ]] && rm -r $exedir/$RODIR"
  done

# jobs is equal to number of CPU cores by default
# note --cleanup doesn't work with parallel 20130922 b/c we're returning a whole directory tree (no rm -r in --cleanup)

  parallel \
    -S $remjp --return $RODIR \
    --nice 18 --halt $hact --eta --progress --joblog parallellog --colsep ',' \
    --workdir $exedir \
    "python3 transcar_run.py" $RODIR $flux0 :::: $BeamEnergyTableFN
    #./beamRunner.sh $RODIR :::: $BeamEnergyTableFN

  ssh-add -D #remove ssh keys from memory

else #local only
  parallel \
    --nice 18 --halt $hact --eta --progress --joblog parallellog --colsep ',' \
    --workdir $exedir \
    "python3 transcar_run.py" $RODIR $flux0 :::: $BeamEnergyTableFN
    #./beamRunner.sh $RODIR :::: $BeamEnergyTableFN

fi
#-- check results for proper simulation finish
./checkoutcome.sh $RODIR

