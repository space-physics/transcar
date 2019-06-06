#!/usr/bin/env python
import pandas
from pathlib import Path
import pytest
import tempfile
from pytest import approx
import transcar.base as transcar
import transcarread as tr

root = Path(__file__).parents[1]
beam = 'beam947.2'
refdir = root / 'tests'/beam
kinfn = 'dir.output/emissions.dat'


def test_transcar():

    with tempfile.TemporaryDirectory() as odir:

        odir = Path(odir).expanduser()

        params = {'rodir': odir,
                  'Q0': 70114000000.0,
                  'msgfn': 'transcar.log',
                  'errfn': 'transcarError.log'
                  }

        beams = pandas.read_csv(root / 'tests/test_E1E2prev.csv', header=None,
                                names=['E1', 'E2', 'pr1', 'pr2']).squeeze()

        transcar.mono_beam_arbiter(beams, params)

        refexc = tr.ExcitationRates(refdir/kinfn)

        exc = tr.ExcitationRates(odir/beam/kinfn)

    ind = [[1, 12, 5], [0, 62, 8]]

    for i in ind:
        assert refexc[i[0], i[1], i[2]].values == approx(exc[i[0], i[1], i[2]].values, rel=1e-3)

    assert refexc.time.shape == refexc.time.shape, 'did you rerun the test without clearing the output directory first?'
    assert (refexc.time == exc.time).all(), 'simultation time of current run did not match reference run'


if __name__ == '__main__':
    pytest.main(['-x', __file__])
