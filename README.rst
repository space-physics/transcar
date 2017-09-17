.. image:: https://travis-ci.org/scivision/transcar.svg
    :target: https://travis-ci.org/scivision/transcar

========
transcar
========

:Fortran Authors: P.L. Blelly, J. Lilensten, M. Zettergren
:Python Author, Fortran cleanup: Michael Hirsch

:Runs on: Linux/Unix, Mac, Windows

TRANSCAR 1D flux tube ionospheric energy deposition flux transport model.
Considers solar input and background conditions via MSIS, HWM.
Models disturbance propagation in ionosphere via models including LCPFCT.

.. contents::

Prereqs
=======
Because Transcar is Python & Fortran based, it runs on any PC/Mac with Linux, OS X, Windows, etc.

You can use your preferred Python install. 
I use `Anaconda Python <http://continuum.io/downloads>`_.

Linux
-----
::

    apt install gfortran cmake make

Mac
---
::

    brew install gcc cmake make

Windows
-------
I often use `Windows Subsystem for Linux <https://www.scivision.co/install-windows-subsystem-for-linux/>`_, which is a full Ubuntu 16.04 system made by Microsoft in the Windows App Store.

However, you can also `install gfortran, cmake and make natively on Windows <https://www.scivision.co/windows-gcc-gfortran-cmake-make-install/>`_


Install
=======
::

  git clone https://github.com/scivision/transcar

  cd transcar/transcar/dir.source/dir.obj
  cmake ..
  make -j7

  cd ../..
  python setup.py develop

Run
======
Simulations are configured in `transcar/dir.input/DATCAR <transcar/dir.input/DATCAR>`_. Simulations are run from the `transcar/` directory.

The legacy method is using Bash (Linux/Mac specific).
The new method is running via Python (runs on any operating system).

Python
------
from the `transcar/` subdirectory::

    python transcar_run.py /tmp/tc 

Bash
----
This is the Mac/Linux only legacy way of running Transcar::

    ./run_beams_parallel.sh /tmp/tc
    
where ``/tmp/tc`` is the output directory. Files are automatically erased there, so be careful!

Optional
========

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

`transconvec <https://github.com/scivision/transcar/blob/master/transcar/dir.source/transconvec_13.op.f>`_  Main Program
