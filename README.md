# Transcar 1-D time-dependent ionosphere flux tube model

![Actions Status](https://github.com/space-physics/transcar/workflows/ci_unix/badge.svg)

Fortran Authors: P.L. Blelly, J. Lilensten, M. Zettergren

Python front-end, and Fortran interfacing:  Michael Hirsch

TRANSCAR 1D flux tube ionospheric energy deposition flux transport model.
Considers solar input and background conditions via MSIS, HWM.
Models disturbance propagation in ionosphere via models including LCPFCT.

## Prereqs

Because Transcar is Python & Fortran based, it runs on any PC/Mac with Linux, MacOS, Windows, etc.
Fortran compilers can be used, including Gfortran and Intel.
However, there are limitations because of the non-standard code used in Transcar, not all compiler vendors or versions work.

* Linux / Windows Subsystem for Linux: `apt install gfortran cmake ninja-build`
* MacOS / Homebrew: `brew install gcc cmake ninja`

### Known working

macOS, Linux, Windows :

* GCC 16.2.0 with -O1 or -O0

### Known not working

Apple:

* Flang 23.1 with -O3, -O2, -O1, -O0 (SIGABRT -6)
* GCC 16.2 with -O3 or -O2 (SIGBUS -10)
* GCC 15.3 with -O1 or -O0 (error 2); SIGBUS -10 with -O3 or -O2

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

python -m pip install -e ./transcar

cmake -B transcar/build -S transcar

cmake --build transcar/build
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
