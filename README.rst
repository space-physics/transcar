========
transcar
========

:Fortran Authors: P.L. Blelly, J. Lilensten, M. Zettergren
:Python Author, Fortran cleanup: Michael Hirsch

TRANSCAR 1D flux tube ionospheric energy deposition flux transport model.
Considers solar input and background conditions via MSIS, HWM.
Models disturbance propagation in ionosphere via models including LCPFCT.

Note, despite a very substantial effort to clean up the code, numerous
deficiencies exist, that would require cleanup function by function.
This is probably too large an effort, given that alternative models are
available. Matt doesn't use Transcar anymore.

Installation
============

.. code:: bash

  git clone https://github.com/scienceopen/transcar
  cd transcar
  ./setup.sh
  
  
To Run
======
Simulations are configured in transcar/dir.input/DATCAR. Simulations are run by::
    
    cd transcar
    ./run_beams_parallel.sh /tmp/tc
    
where /tmp/tc is the output directory. Files are automatically erased there, so be careful!
