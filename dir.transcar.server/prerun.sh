#!/bin/bash
# Michael Hirsch
# tasks to do before running the sim

freshenout()
{
[[ -d $1 ]] && echo $1
}

getcode()
{
# git pull requires a public repo
(cd $exedir && git pull)
}

freshenout $1/$2
getcode $1
