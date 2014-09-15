#!/bin/bash
# Michael Hirsch 2014

# Note: this was not possible using find/tail/grep without a for loop -- would need gawk
RODIR=$1
[[ -z $RODIR ]] && { echo "error you must specify directory to examine"; exit 1; }

abnormal=0
normal=0

for f in $(find $RODIR -mindepth 1 -maxdepth 2 -type f -name "TranscarErrors.txt"); do
  outcome=$(tail -n1 $f)
  [[ $outcome != *fin\ normale* ]] && { echo "abnormal completion in $f"; abnormal=$((abnormal+1)); }
  [[ $outcome == *fin\ normale* ]] && normal=$((normal+1))
done

echo "in directory $RODIR, we detected"
echo "$normal normal simulation completion"
echo "$abnormal ABNORMAL simulation completion"
