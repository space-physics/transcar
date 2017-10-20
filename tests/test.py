#!/usr/bin/env python
import subprocess
from pathlib import Path
from numpy.testing import run_module_suite, assert_allclose
#
from transcarread import ExcitationRates

rdir = Path(__file__).parents[1]
print(rdir)
rdir = Path('/tmp/test/')
refdir = Path('tests/beam947.2')
odir = rdir/refdir.name
kinfn = 'dir.output/emissions.dat'

def test_transcar():
    subprocess.run(['python', 'RunTranscar.py', 'newdata', '-infn', 'tests/test_E1E2prev.csv'],#  cwd=rdir,
                          timeout=90)

    refexc, tref = ExcitationRates(refdir/kinfn)[:2]

    exc, t = ExcitationRates(odir/kinfn)[:2]

    ind=[[1,12,5],[0,62,8]]

    for i in ind:
        assert_allclose(refexc[i[0],i[1],i[2]],
                           exc[i[0],i[1],i[2]], rtol=1e-3)

    assert tref[0]==t[0]

if __name__ == '__main__':
    run_module_suite()
