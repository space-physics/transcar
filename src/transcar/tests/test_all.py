from pathlib import Path

import pandas
from pytest import approx
import transcar.base as transcar
import transcarread as tr

root = Path(__file__).parent
beam = "beam947.2"
refdir = root / beam
kinfn = "dir.output/emissions.dat"


def test_run_transcar(tmp_path):

    odir = tmp_path

    params = {
        "rodir": odir,
        "Q0": 70114000000.0,
        "logfn": "transcar.log",
        "errfn": "transcarError.log",
    }

    beams = pandas.read_csv(root / "test_E1E2prev.csv", header=None, names=["E1", "E2", "pr1", "pr2"]).squeeze()
    beams = beams.to_dict()

    if not transcar.mono_beam_arbiter(beams, params):
        raise RuntimeError(f"Transcar run failed. See {tmp_path / 'transcarError.log'} for clues")

    refexc = tr.ExcitationRates(refdir / kinfn)

    exc = tr.ExcitationRates(odir / beam / kinfn)

    ind = [[1, 12, 5], [0, 62, 8]]

    for i in ind:
        assert refexc[i[0], i[1], i[2]].values == approx(exc[i[0], i[1], i[2]].values, rel=1e-3)

    assert refexc.time.shape == refexc.time.shape, "did you rerun the test without clearing the output directory first?"
    assert (refexc.time == exc.time).all(), "simulation time of current run did not match reference run"
