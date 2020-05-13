# transcar

[![image](https://travis-ci.org/scivision/transcar.svg)](https://travis-ci.org/scivision/transcar)

* Fortran Authors:   P.L. Blelly, J. Lilensten, M. Zettergren
* Python Author, Fortran cleanup:  Michael Hirsch

Runs on: Linux/Unix, Mac, Windows

TRANSCAR 1D flux tube ionospheric energy deposition flux transport
model. Considers solar input and background conditions via MSIS, HWM.
Models disturbance propagation in ionosphere via models including
LCPFCT.

## Prereqs

Because Transcar is Python & Fortran based, it runs on any PC/Mac with
Linux, OS X, Windows, etc.

You can use your preferred Python install. I use [Anaconda
Python](http://continuum.io/downloads).

### Linux

    apt install gfortran cmake make

### Mac

    brew install gcc cmake make

### Windows

Use MSYS2 or Windows Subsystem for Linux to get Gfortran compiler

## Install

CMake, Make/Ninja, Gfortran and Python 3:

```sh
python setup.py develop

cmake -S dir.source -B dir.source/dir.obj

cmake --build dir.source/dir.obj
```

## Run

Simulations are configured in
[transcar/dir.input/DATCAR](transcar/dir.input/DATCAR). Simulations are
run from the [transcar/]{.title-ref} directory.

The legacy method is using Bash (Linux/Mac specific). The new method is
running via Python (runs on any operating system).

### Python

from the [transcar/]{.title-ref} subdirectory:

    python transcar_run.py /tmp/tc

### Bash

This is the Mac/Linux only legacy way of running Transcar:

    ./run_beams_parallel.sh /tmp/tc

where `/tmp/tc` is the output directory. Files are automatically erased
there, so be careful!

install the latest GNU Parallel into \~/bin by:

    wget -O - pi.dk/3 | bash
