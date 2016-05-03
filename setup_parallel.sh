#!/bin/bash

(
wd=$(mktemp -d)
wget -nc -P $wd ftp://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2

cd $wd
tar -xf parallel-latest.tar.bz2
cd parallel-*
./configure && make && make install
)
