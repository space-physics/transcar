#!/bin/bash
# brought up to BASH 4 conventions by Michael Hirsch mhirsch@bu.edu
# original by Matt Zettergren 

# note: relative paths used here work because of "cd" step in beamRunner.sh

runTranscar()
{
[[ -z $1 ]] && { echo "error: must specify root directory e.g. ./run_beams ../BT"; exit 1; } 

#set +e #don't stop on error unless I direct it to

#not necessary
#http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in 
#TDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 

ODIR=$1
TClog=$ODIR/transcarBash.log
TCconfig=dir.input/DATCAR
msisinitfn=$(grep "input file (initial ionospheric conditions)" $TCconfig | cut -d" " -f1)
echo "using msis init file $msisinitfn as specified in $TCconfig">>$TClog 2>&1

#check for unused output dir
[[ ! -d "$ODIR" ]] && { echo "Energy directory $ODIR does not exist" | tee -a $TClog; exit 98; }

#copy transcar input
cp -v -t "$ODIR/dir.input/" $TCconfig dir.input/$msisinitfn >>$TClog 2>&1
cp -v -t "$ODIR/dir.data/" dir.data/type >>$TClog 2>&1
cp -v -t "$ODIR/dir.data/dir.linux/dir.geomag/" dir.data/dir.linux/dir.geomag/{data_geom.bin,igrf90.dat,igrf90s.dat} >>$TClog 2>&1
cp -v -t "$ODIR/dir.data/dir.linux/dir.projection/" dir.data/dir.linux/dir.projection/varpot.dat >>$TClog 2>&1
#cp -v -t "$ODIR/dir.data/dir.linux/dir.projection/" dir.data/dir.linux/dir.projection/varcourant.dat >>$TClog 2>&1  #makes program crash
cp -v -t "$ODIR/dir.data/dir.linux/dir.cine/" dir.data/dir.linux/dir.cine/{DAT{DEG,FEL,TRANS},FELTRANS,flux.flag} >>$TClog 2>&1
cp -v -t "$ODIR/dir.data/dir.linux/dir.cine/dir.euvac/" dir.data/dir.linux/dir.cine/dir.euvac/EUVAC.dat >>$TClog 2>&1
cp -v -t "$ODIR/dir.data/dir.linux/dir.cine/dir.seff/" dir.data/dir.linux/dir.cine/dir.seff/{crsb8,crsphot1.dat,rdtb8} >>$TClog 2>&1

#copy transcar itself
cp -v transconvec_13.op.out "$ODIR/"


# run TRANSCAR
#echo "RUN_TRANSCAR.SH: transconvec_13.op is running"
#time strace -e trace=open transconvec_13.op.out | tee $PWD/bid
#time transconvec_13.op.out | tee $PWD/bid | egrep -A3 -e "*******|Error|felin.f"
#time transconvec_13.op.out > $PWD/bid & tail -f $PWD/bid | egrep --line-buffered -A3 -e "*******|Error|felin.f"

#must have /usr/bin/time to use the system time command, as bash time command won't accept arguments!
# egrep --line-buffered is necessary to avoid annoying random line truncations while waiting for stdout to fill!
#/usr/bin/time --output="$ODIR/timing.txt" ./transconvec_13.op.out | \


# exec used to reduce memory usage http://stackoverflow.com/questions/786376
 (cd $ODIR && exec ./transconvec_13.op.out) >>$ODIR/TranscarErrors.txt 2>&1
    # egrep -i --line-buffered -e "Warning|Error|felin.f|trans.f|Input eV/cm2/s" | \
      
}

# this noMore is only used for non-parallel runs
#noMore ()
#{
#echo "run_transcar: aborting on user request Ctrl C" | tee -a $TClog
#return 99
#}
#trap noMore SIGINT

runTranscar $1

