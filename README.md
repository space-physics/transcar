[![image](https://travis-ci.org/scivision/transcar.svg)](https://travis-ci.org/scivision/transcar)
[![image](https://coveralls.io/repos/github/scivision/transcar/badge.svg?branch=next)](https://coveralls.io/github/scivision/transcar?branch=next)
[![image](https://ci.appveyor.com/api/projects/status/d4y6eqqjjq4uq2sw?svg=true)](https://ci.appveyor.com/project/scivision/transcar)
[![Maintainability](https://api.codeclimate.com/v1/badges/7c237d2870d0611e5df6/maintainability)](https://codeclimate.com/github/scivision/transcar/maintainability)

# Transcar

Fortran Authors: P.L. Blelly, J. Lilensten, M. Zettergren

Python Author, Fortran cleanup:  Michael Hirsch

TRANSCAR 1D flux tube ionospheric energy deposition flux transport model. 
Considers solar input and background conditions via MSIS, HWM.
Models disturbance propagation in ionosphere via models including LCPFCT.

## Prereqs

Because Transcar is Python & Fortran based, it runs on any PC/Mac with Linux, OS X, Windows, etc.

You can use your preferred Python >= 3.6. 
I use [Anaconda Python](http://continuum.io/downloads).

-   Linux: `apt install gfortran cmake make`
-   Mac: `brew install gcc cmake make`
-   Windows: use either:
    -   [Windows Subsystem for Linux](https://www.scivision.co/install-windows-subsystem-for-linux/),
        a full Ubuntu system made by Microsoft in the Windows App Store.
    -   install gfortran, cmake and make [natively on Windows](https://www.scivision.co/brew-install-scoop-for-windows/)

## Install

You'll need Cmake, Gfortran and Python 3:

    python -m pip install -e .

## Usage

Simulations are configured in
[transcar/dir.input/DATCAR](transcar/dir.input/DATCAR). 
Simulations are run from the top directory. 
Python runs Transcar in parallel execution using `concurrent.futures`, dynamically adapting to the number
of CPU cores available:
```sh
python RunTranscar.py /tmp/tc
```

## Reference

[transconvec](https://github.com/scivision/transcar/blob/master/transcar/dir.source/transconvec_13.op.f)
Main Program

### Manual compile

This is not normally needed, just for reference:
```sh
cd dir.source/dir.obj
cmake ..
cmake --build .

cd ..
```

### Bash

This is the Mac/Linux only legacy way of running Transcar:

    ./run_beams_parallel.sh /tmp/tc

where `/tmp/tc` is the output directory. Files are automatically erased
there, so be careful!

### Specify compiler

Cmake uses your system default Fortran compiler, but if you wish to use
another compiler, set environment variable FC in the Cmake call. For
example, Intel `ifort`:

    FC=ifort cmake ..

### Parallel remote execution

install the latest GNU Parallel into ~/bin by:

    wget -O - pi.dk/3 | bash
