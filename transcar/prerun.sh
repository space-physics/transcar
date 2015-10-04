#!/bin/bash
# Michael Hirsch
# tasks to do before running the sim

remotes=(labHST0 labHST1 irs4 irs3 swoboj)
exedir=code/transcar/transcar

freshenout()
{
[[ -d $1 ]] && echo $1
}

getcode()
{
# git pull requires a public repo
(cd $exedir && git pull && cd dir.source && make -s)

}

for remote in ${remotes[@]}; do
ssh $remote -t "(cd $exedir && git pull && cd dir.source && make -s && cd; cd code/transcar-utils && git pull && cd ../hist-utils && git pull)"
done

#freshenout $1/$2
#getcode $1
