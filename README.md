# Transcar

[![Build Status](https://travis-ci.org/space-physics/transcar.svg?branch=next)](https://travis-ci.org/space-physics/transcar)

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

cmake -B build

cmake --build build
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

### Plotting

The simulation results are loaded and plotted by the [transcarread](https://github.com/space-physics/transcarread) Python package.
