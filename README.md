[![Build Status](https://travis-ci.org/space-physics/transcar.svg?branch=next)](https://travis-ci.org/space-physics/transcar)
[![image](https://coveralls.io/repos/github/scivision/transcar/badge.svg?branch=next)](https://coveralls.io/github/scivision/transcar?branch=next)
[![Build status](https://ci.appveyor.com/api/projects/status/ij3jskpvfvprm185?svg=true)](https://ci.appveyor.com/project/scivision/transcar-yhafk)
[![Maintainability](https://api.codeclimate.com/v1/badges/7c237d2870d0611e5df6/maintainability)](https://codeclimate.com/github/scivision/transcar/maintainability)

# Transcar

Fortran Authors: P.L. Blelly, J. Lilensten, M. Zettergren

Python front-end, and Fortran interfacing:  Michael Hirsch

TRANSCAR 1D flux tube ionospheric energy deposition flux transport model.
Considers solar input and background conditions via MSIS, HWM.
Models disturbance propagation in ionosphere via models including LCPFCT.

## Prereqs

Because Transcar is Python & Fortran based, it runs on any PC/Mac with Linux, MacOS, Windows, etc.
Most Fortran compilers can be used, including Gfortran and Intel.

* Linux / Windows Subsystem for Linux: `apt install gfortran cmake make`
* Mac: `brew install gcc cmake make`

### Windows

From native Windows
[install CMake](https://cmake.org/download)
and either of:

* [Gfortran](https://www.scivision.dev/install-msys2-windows/)
* [Intel Parallel Studio](https://www.scivision.dev/install-intel-compiler-icc-icpc-ifort/)


## Install

from Terminal / Command Prompt

```sh
git clone https://github.com/scivision/transcar

cd transcar

python -m pip install -e .

python build.py
```

## Usage

Simulations are configured in
[dir.input/DATCAR](./dir.input/DATCAR).
Simulations are run from the top directory.
Python runs Transcar in parallel using
[concurrent.futures.ThreadPoolExecutor](https://docs.python.org/3/library/concurrent.futures.html),
dynamically adapting to the number of CPU cores available:

```sh
python MonoenergeticBeams.py /tmp/tc
```

## Notes

[transconvec](https://github.com/scivision/transcar/blob/master/transcar/dir.source/transconvec_13.op.f)

### Manual compile

This is not normally needed, just for reference:
```sh
cd dir.source/dir.obj
cmake ..
cmake --build . -j

cd ..
```


### Specify compiler

Cmake uses your system default Fortran compiler, but if you wish to use
another compiler, set environment variable FC in the Cmake call. For
example, Intel `ifort`:

    FC=ifort cmake ..
