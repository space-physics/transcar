.. image:: https://travis-ci.com/scienceopen/transcar.svg?token=qSpy37WAwefTyjRcouQS&branch=master
    :target: https://travis-ci.com/scienceopen/transcar

========
transcar
========

:Fortran Authors: P.L. Blelly, J. Lilensten, M. Zettergren
:Python Author, Fortran major cleanup: Michael Hirsch

TRANSCAR 1D flux tube ionospheric energy deposition flux transport model.
Considers solar input and background conditions via MSIS, HWM.
Models disturbance propagation in ionosphere via models including LCPFCT.

Note, despite a very substantial effort to clean up the code, numerous
deficiencies exist, that would require cleanup function by function.
This is probably too large an effort, given that alternative models are
available. Matt doesn't use Transcar anymore.

I suggest getting started with Bash-based (terminal shell) program first, then
if needed we can do to the networked poor-mans cluster version. There is extra time
needed for one-time setup of those nodes.

.. contents::

Prereqs
=======
::

    sudo apt-get install gfortran cmake make parallel bc wget
    
If you don't have GNU Parallel, you can install it via ``setup_parallel.sh`` on any system including Cygwin.
If it's not compiling for you let me know, it really should be platform-independent.

Fortran Installation
====================
The basic program uses Bash and Fortran code, so it runs anywhere (Linux/BSD/Mac/Windows)::

  git clone https://github.com/scienceopen/transcar
  cd transcar/transcar/dir.source/dir.obj
  cmake ..
  make -j7 --quiet
  cd ../..

If you then modify the Fortran source code, you just have to in dir.obj type::

    make
  
To Run
======
Simulations are configured in transcar/dir.input/DATCAR. Simulations are run by::
    
    cd transcar
    ./run_beams_parallel.sh /tmp/tc
    
where /tmp/tc is the output directory. Files are automatically erased there, so be careful!

Optional
========

Python installation
-------------------
If you wish to use the user-friendly Python interfaces (alpha test)::

    python setup.py develop

Specify compiler
----------------
Cmake uses your system default Fortran compiler, but if you wish to use another compiler, set environment variable FC in the Cmake call. For example, Intel ``ifort``::

    FC=ifort cmake ..


