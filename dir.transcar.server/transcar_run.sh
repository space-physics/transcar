#!/bin/bash
# brought up to BASH 4 conventions by Michael Hirsch mhirsch@bu.edu
# original by Matt Zettergren

# note: relative paths used here work because of "cd" step in beamRunner.sh

copyTranscarFiles()
{
[[ -z $1 ]] && { echo "error: must specify root directory" | tee -a $TClog; exit 1; }
ODIR=$1
[[ ! -d $ODIR ]] && { echo "Beam dir. $ODIR does not exist" | tee -a $TClog; exit 98; }

TClog=$ODIR/transcarBash.log
TCconfig=dir.input/DATCAR #this is hard-coded in transcar

initfn=$(grep "input file (initial ionospheric conditions)" $TCconfig | cut -d" " -f1)
if [[ -n $initfn ]]; then
 echo "using msis init file $initfn as specified in $TCconfig">>$TClog
else
 echo "initial conditions file not found, aborting"
 exit 1
fi

#copy transcar input
\cp -p -t "$ODIR/dir.input/" $TCconfig dir.input/$initfn 2>>$TClog
\cp -p -t "$ODIR/dir.data/" dir.data/type 2>>$TClog
\cp -p -t "$ODIR/dir.data/dir.linux/dir.geomag/" dir.data/dir.linux/dir.geomag/{data_geom.bin,igrf90.dat,igrf90s.dat} 2>>$TClog
\cp -p -t "$ODIR/dir.data/dir.linux/dir.projection/" dir.data/dir.linux/dir.projection/varpot.dat 2>>$TClog
#\cp -p -t "$ODIR/dir.data/dir.linux/dir.projection/" dir.data/dir.linux/dir.projection/varcourant.dat 2>>$TClog  #makes program crash
\cp -p -t "$ODIR/dir.data/dir.linux/dir.cine/" dir.data/dir.linux/dir.cine/{DAT{DEG,FEL,TRANS},FELTRANS,flux.flag} 2>>$TClog
\cp -p -t "$ODIR/dir.data/dir.linux/dir.cine/dir.euvac/" dir.data/dir.linux/dir.cine/dir.euvac/EUVAC.dat 2>>$TClog
\cp -p -t "$ODIR/dir.data/dir.linux/dir.cine/dir.seff/" dir.data/dir.linux/dir.cine/dir.seff/{crsb8,crsphot1.dat,rdtb8} 2>>$TClog
#copy transcar itself
\cp transconvec_13.op.out "$ODIR/"
}

runTranscar()
{
#time strace -e trace=open transconvec_13.op.out | tee $PWD/bid
#time transconvec_13.op.out | tee $PWD/bid | egrep -A3 -e "*******|Error|felin.f"
#time transconvec_13.op.out > $PWD/bid & tail -f $PWD/bid | egrep --line-buffered -A3 -e "*******|Error|felin.f"

#must have /usr/bin/time to use the system time command, as bash time command won't accept arguments
# egrep --line-buffered is necessary to avoid annoying random line truncations while waiting for stdout to fill
#/usr/bin/time --output="$ODIR/timing.txt" ./transconvec_13.op.out | \

# exec used to reduce memory usage http://stackoverflow.com/questions/786376
 (cd $ODIR && exec ./transconvec_13.op.out) 2>>$ODIR/transcarErrors.log >>$ODIR/transcar.log
    # egrep -i --line-buffered -e "Warning|Error|felin.f|trans.f|Input eV/cm2/s" | \

}

# this noMore is only used for non-parallel runs
#noMore ()
#{
#echo "run_transcar: aborting on user request Ctrl C" | tee -a $TClog
#return 99
#}
#trap noMore SIGINT

copyTranscarFiles $1
runTranscar

