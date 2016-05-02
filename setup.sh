#!/bin/bash
# Installer script for Python Transcar
# Michael Hirsch
#
# tested with gfortran, but should with ifort, pgf, et al. Let me know if it doesn't.

(
cd transcar/dir.source/dir.obj
cmake ..
make -j7 --quiet
)

conda install --file requirements.txt
python setup.py develop
