#!/bin/bash

wget -nc -P /tmp ftp://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2

(
cd /tmp
tar -xf parallel-latest.tar.bz2
cd parallel-*
./configure && make && make install
)
