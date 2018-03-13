#!/usr/bin/env python
import pandas
from pathlib import Path
from numpy.testing import assert_allclose
import tempfile
#
import transcar
import transcarread as tr

root = Path(__file__).parents[1]
beam =  'beam947.2'
odir = Path(tempfile.gettempdir()) / 'newdata'
refdir = root / 'tests'/beam
kinfn = 'dir.output/emissions.dat'

def test_transcar():
    odir.mkdir(parents=True, exist_ok=True)

    params = {'rodir': odir,
              'Q0': 70114000000.0,
              'msgfn': 'transcar.log',
              'errfn': 'transcarError.log'
              }

    beams = pandas.read_csv(root / 'tests/test_E1E2prev.csv', header=None, names=['E1','E2','pr1','pr2']).squeeze()

    transcar.iterbeams(beams, params)

    refexc = tr.ExcitationRates(refdir/kinfn)

    exc = tr.ExcitationRates(odir/beam/kinfn)

    ind=[[1,12,5],[0,62,8]]

    for i in ind:
        assert_allclose(refexc[i[0],i[1],i[2]],
                           exc[i[0],i[1],i[2]], rtol=1e-3)

    assert refexc.time == exc.time


if __name__ == '__main__':
    test_transcar()
