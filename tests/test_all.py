#!/usr/bin/env python
import pandas
from pathlib import Path
from numpy.testing import assert_allclose
import tempfile
#
import transcar
from transcarread import ExcitationRates

rdir = Path(__file__).parents[1]
print(rdir)
rdir = Path(tempfile.gettempdir()) / 'newdata'
refdir = Path('tests/beam947.2')
odir = rdir/refdir.name
kinfn = Path('dir.output') / 'emissions.dat'

def test_transcar():
    odir.mkdir(parents=True, exist_ok=True)

    params = {'rodir': rdir,
              'Q0': 70114000000.0,
              'msgfn': 'transcar.log',
              'errfn': 'transcarError.log'
              }

    beams = pandas.read_csv('tests/test_E1E2prev.csv', header=None, names=['E1','E2','pr1','pr2']).squeeze()

    transcar.iterbeams(beams, params)

    refexc, tref = ExcitationRates(refdir/kinfn)[:2]

    exc, t = ExcitationRates(odir/kinfn)[:2]

    ind=[[1,12,5],[0,62,8]]

    for i in ind:
        assert_allclose(refexc[i[0],i[1],i[2]],
                           exc[i[0],i[1],i[2]], rtol=1e-3)

    assert tref[0]==t[0]

if __name__ == '__main__':
    test_transcar()
