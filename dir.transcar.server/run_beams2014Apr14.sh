#!/bin/bash
# upgraded to Bash 4 by Michael Hirsch 2014
# original by Matt Zettergren 2013
# this script for loops transcar, making a new precinput.dat each time for the 
# respective beam energies.

runBeams()
{

if [[ -z $1 ]]; then
 echo "error: must specify root directory e.g. ./run_beams ../BT"; exit 1
fi

#set -e #stop on any error
set +e #don't stop on error unless specifically directed

RODIR=$1
BMlog=$RODIR/Beams.log

flux0=70114000000.0

energies=(52.726	63.557	76.487	91.921	110.34	132.34	158.59	189.93)
energies+=( 227.34	272.	325.3	388.94	464.9	555.57)	
energies=(${energies[@]} 663.81 793.02	947.25	1131.4	1351.1	1613.5)
energies=(${energies[@]} 1926.7	2300.5  2746.8	3279.5	3915.4	4674.4)
energies=(${energies[@]} 5580.6	6662.2	7953.4  9494.7	11335.	13531.)
energies=(${energies[@]} 16152.	19282.)

prev=(49.521 59.731 71.92 86.47	103.84	124.57	149.32	178.86	214.13	256.23)
prev=(${prev[@]} 306.48	366.46	438.07	523.55	625.58	747.38	892.78	1066.3)
prev=(${prev[@]} 1273.5	1520.8	1816.1	2168.5	2589.2	3091.3	3690.8	4406.3)
prev=(${prev[@]} 5260.5	6280.2	7497.3	8950.3	10685.	12755.	15227.	18177.)

tstart=3000. #seconds, where precipitation starts?
tfin=3600. #seconds, where precipitation stops?

nEnergies=${#energies[@]}
PrecFN="dir.input/precinput.dat"
echo "Preparing to run for $nEnergies energies:  " ${energies[@]} | teea $BMlog
#echo "overwriting $PrecFN for TRANSCAR Input file"
for ((i=0;  i<$(( $nEnergies-1 )); i++))
do
CurrDir=$RODIR/beam${energies[$i]}
i1=$(( $i+1 ))

echo "computing energy # $i1" | teea $BMlog
  #calculate the beam characteristics
  E1=${energies[$i]}
  E2=${energies[$i1]}
  dE=$(echo "$E2 - $E1" | bc)
  Esum=$(echo "$E2 + $E1" | bc)

  # bc will not give any decimal without scale argument in this case
  flux=$(echo "scale=4; $flux0 / 0.5 / $Esum / $dE" | bc) 
  #get 'padding' so that beams are truncated correctly in energy
  pr=${prev[$i]};
  dElow=$(echo "$E1 - $pr" | bc)
  Elow=$(echo "scale=4; $E1 - 0.5 * $dElow" | bc)
  ne=${prev[$i1]}
  dEhigh=$(echo "$E2 - $ne" | bc)
  Ehigh=$(echo "scale=4; $E2 - 0.5 * $dEhigh" | bc)

#echo $E1 $E2 $dE $Esum $flux $pr $dElow $Elow $ne $dEhigh $Ehigh #diagnostic
  #generate TRANSCAR input file for this beam
  ThisPrecParam="$tstart\\n$Elow $flux\\n$Ehigh -1.0\\n$tfin\\n-1.0 -1.0"
  echo "writing $ThisPrecParam to $PrecFN" | teea $BMlog
  echo -e "$ThisPrecParam" > "$PrecFN" # the -e option was in the original script

  echo "transconvec_13.op is running for a differential number flux of " | teea $BMlog
  echo "$flux eV cm-2 s-1 sR-1" | teea $BMlog
  echo "in the energy range $E1 to $E2 eV" | teea $BMlog
  echo "Starting at t=$tstart hours until t=$tfin hours" | teea $BMlog
  echo "output to $CurrDir" | teea $BMlog

  #Run the sim
  ./run_transcar.sh "$CurrDir" | teea $BMlog

  #error trap
  LastErr=$?
  case "$LastErr" in
  0) echo "Energy ${energies[$i]} completed in $CurrDir" | teea $BMlog
     ;;
  98) echo "Skipped energy ${energies[$i]} in $CurrDir" | teea $BMlog
     ;;
  99) echo "Aborting run_beams of $CurrDir per user Ctrl+C" | teea $BMlog
     ;;
  *) echo "Transcar exit code $LastErr" |  teea $BMlog
     exit $LastErr
     ;;
  esac

done #for energy

echo "Simulations complete."
}

teea ()
{
tee --append "$1"
}

#input argument is root output directory
runBeams "$1"
