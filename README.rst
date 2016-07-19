.. image:: https://travis-ci.com/scienceopen/transcar.svg?token=qSpy37WAwefTyjRcouQS&branch=master
    :target: https://travis-ci.com/scienceopen/transcar

========
transcar
========

:Fortran Authors: P.L. Blelly, J. Lilensten, M. Zettergren
:Python Author, Fortran cleanup: Michael Hirsch

:Runs on: Linux/Unix, Mac, Windows (Cygwin or Windows Subsytem for Linux)

TRANSCAR 1D flux tube ionospheric energy deposition flux transport model.
Considers solar input and background conditions via MSIS, HWM.
Models disturbance propagation in ionosphere via models including LCPFCT.

Note, despite a very substantial effort to clean up the code, numerous
deficiencies exist, that would require cleanup function by function.
This is probably too large an effort, given that alternative models are
available.

.. contents::

Install
=======
::

  apt-get install gfortran cmake make bc wget
  git clone https://github.com/scienceopen/transcar
  cd transcar/transcar/dir.source/dir.obj
  cmake ..
  make -j7 --quiet
  cd ../..

To Run
======
Simulations are configured in transcar/dir.input/DATCAR. Simulations are run by::
    
    cd transcar
    ./run_beams_parallel.sh /tmp/tc
    
where ``/tmp/tc`` is the output directory. Files are automatically erased there, so be careful!

Optional
========

Python installation
-------------------
If you wish to use the Python interfaces (alpha test)::

    python setup.py develop

Specify compiler
----------------
Cmake uses your system default Fortran compiler, but if you wish to use another compiler, set environment variable FC in the Cmake call. For example, Intel ``ifort``::

    FC=ifort cmake ..

Parallel remote execution
-------------------------
install the latest GNU Parallel into ~/bin by::

    wget -O - pi.dk/3 | bash

Code
====

`transconvec <https://github.com/scienceopen/transcar/blob/master/transcar/dir.source/transconvec_13.op.f>`_  Main Program
