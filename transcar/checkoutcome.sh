#!/bin/bash
# Michael Hirsch 2014

# Note: this was not possible using find/tail/grep without a for loop -- would need gawk
RODIR=$1
[[ -z $RODIR ]] && { echo "error you must specify directory to examine"; exit 1; }

abnormal=0
normal=0

for f in $(find $RODIR -mindepth 1 -maxdepth 2 -type f -name "transcarError.log"); do
  outcome=$(tail -n1 $f)

#don't do this if as one-liner, gives impossible good and bad at same time for first file
  if [[ $outcome == *fin\ normale* ]]; then
     ((normal++)) 
  else 
     echo "abnormal completion in $f"
     ((abnormal++))
  fi
done

echo "in directory $RODIR, we detected"
echo "$normal normal simulation completion"
echo "$abnormal ABNORMAL simulation completion"
