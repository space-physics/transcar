#!/bin/bash
# Installer script for Python Transcar
# Michael Hirsch
#
# tested with gfortran, but should with ifort, pgf, et al. Let me know if it doesn't.

(
cd transcar/dir.source
make
)

python setup.py develop
