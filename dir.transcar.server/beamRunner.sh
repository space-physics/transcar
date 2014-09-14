#!/bin/bash

DebugMsg=0
#set +e #don't stop on any error, so that return codes from transcar will work
RODIR=$1
E1=$2
E2=$3
pr1=$4
pr2=$5

CurrDir="$RODIR/beam$E1"
# freshen simulation directory
[[ -d $CurrDir ]] && { \rm -rv $CurrDir; mkdir -v $CurrDir; }

BMlog=$RODIR/Beams.log
TCconfig=dir.input/DATCAR

flux0=10016000000.0 #FIXME where does this come from?

tstart=$(grep "precipitation start time (seconds)" $TCconfig | cut -f1)
tfin=$(grep "precipitation end time (seconds)" $TCconfig | cut -f1)

[[ -z $tstart ]] && { echo "error: could not find precip start in $TCconfig"; exit 1; }
[[ -z $tfin ]] && { echo "error: could not find precip end in $TCconfig"; exit 1; }  

# make a directory for this beam -- 
# everything relevant to sim will reside in this directory, including executable!
mkdir -pv $CurrDir/dir.input $CurrDir/dir.output \
          $CurrDir/dir.data/dir.linux/dir.geomag \
          $CurrDir/dir.data/dir.linux/dir.projection \
          $CurrDir/dir.data/dir.linux/dir.cine/dir.seff \
          $CurrDir/dir.data/dir.linux/dir.cine/dir.euvac \
          >>$BMlog 2>&1

PrecFN="$CurrDir/dir.input/precinput.dat"


  #calculate the beam characteristics
  dE=$(echo "$E2 - $E1" | bc)
  Esum=$(echo "$E2 + $E1" | bc)

  # bc will not give any decimal without scale argument in this case
  flux=$(echo "scale=4; $flux0 / 0.5 / $Esum / $dE" | bc) 

  #get 'padding' so that beams are truncated correctly in energy
  pr=$pr1
  dElow=$(echo "$E1 - $pr" | bc)
  Elow=$(echo "scale=4; $E1 - 0.5 * $dElow" | bc)
  ne=$pr2
  dEhigh=$(echo "$E2 - $ne" | bc)
  Ehigh=$(echo "scale=4; $E2 - 0.5 * $dEhigh" | bc)

[[ $DebugMsg -ne 0 ]] && echo "$E1 $E2 $dE $Esum $flux $pr $dElow $Elow $ne $dEhigh $Ehigh"

  #generate TRANSCAR input file for this beam
  ThisPrecParam="$tstart\\n$Elow $flux\\n$Ehigh -1.0\\n$tfin\\n-1.0 -1.0"
  #echo "writing $ThisPrecParam to $PrecFN"
  echo -e "$ThisPrecParam" > "$PrecFN" # the -e option was in the original script
if [[ $DebugMsg -ne 0 ]]; then  
  echo "transconvec_13.op is running for a differential number flux of "
  echo "$flux eV cm-2 s-1 sR-1"
  echo "in the energy range $E1 to $E2 eV"
  echo "Starting at t=$tstart sec. until t=$tfin sec."
  echo; echo "output to $CurrDir"
fi
  #Run the sim
  #(cd $exedir && exec ./run_transcar.sh $CurrDir)
  ./run_transcar.sh $CurrDir

  #error trap
  LastErr=$?
  case "$LastErr" in
  0) echo "run_beams: Energy $E1 completed in $CurrDir" | tee -a $BMlog;;
  98) echo "run_beams: Skipped energy $E1 in $CurrDir" | tee -a $BMlog;;
  99) echo "run_beams: Aborting run_beams of $CurrDir" | tee -a $BMlog; exit 99;;
  126) echo "run_beams: Transconvec was not given exec permissions" | tee -a $BMlog; exit 126;;
  *) echo "run_beams: Transcar exit code $LastErr" |  tee -a $BMlog; exit $LastErr;;
  esac
