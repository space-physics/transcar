#!/bin/bash
# Michael Hirsch 2014
# prototype for parallel execution, by looping 

BeamLooper()
{
RTdir=$1
[[ ! -d $RTdir ]] && { echo "Error, root directory $RTdir does not exist. Aborting"; exit 1; }

Efile=$2;
[[ ! -a $Efile ]] && { echo "Error, energy CSV file $Efile does not exist. Aborting"; exit 2; }


cIFS=$IFS #we will restore this after the while loop
IFS=,
while read E1 E2 pr1 pr2 
do
	./run_beams.sh $RTdir $E1 $E2 $pr1 $pr2 # this is what parallel will do

     #echo $RTdir $E1 $E2 $pr1 $pr2 #debug
     #error handling
     LastErr=$?
     case $LastErr in 
     99) echo "loop_beams: Debug exit $E1, code $LastErr"; exit 99 ;;
     *) echo "loop_beams: $E1 returned $LastErr" ;;
     esac

done < $Efile
IFS=$cIFS
}

BeamLooper $1 $2
