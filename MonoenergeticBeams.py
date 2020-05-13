#!/usr/bin/env python
"""
Executes Transcar to output "monoenergetic" electron beams.
Optionally, in parallel.
"""
import concurrent.futures
import logging
from pathlib import Path
from pandas import read_csv
from argparse import ArgumentParser
import time

#
import transcar
import signal

try:
    import psutil

    Ncpu = psutil.cpu_count(logical=False)
except ImportError:
    import os

    Ncpu = max(1, os.cpu_count() // 2)


def main():
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    p = ArgumentParser(description="parallel Transcar runner")
    p.add_argument("rodir", help="root of beam directory to output")
    p.add_argument("-Q0", help="Assumed particle flux", type=float, default=70114000000.0)
    p.add_argument("-infn", help="energy bin CSV file", default="BT_E1E2prev.csv")
    p.add_argument("--msgfn", help="file to write transcar messages to", default="transcar.log")
    p.add_argument("--errfn", help="file to write transcar Errors to", default="transcarError.log")
    p.add_argument("-np", help="number of concurrent processes", type=int, default=Ncpu)
    p = p.parse_args()

    rodir = Path(p.rodir).expanduser().resolve()
    infn = Path(p.infn).expanduser().resolve()

    params = {"rodir": rodir, "Q0": p.Q0, "msgfn": p.msgfn, "errfn": p.errfn}

    rodir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=rodir / "Beams.log",
        filemode="a",
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        level=logging.DEBUG,
    )

    tic = time.time()
    beams = read_csv(infn, header=None, names=["E1", "E2", "pr1", "pr2"])

    print("using", p.np, "concurrent Transcar runs")

    # %% do run
    if p.np == 1:
        for _, beam in beams.iterrows():
            transcar.mono_beam_arbiter(beam, params)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=p.np) as executor:
            future_beam = (executor.submit(transcar.mono_beam_arbiter, beam, params) for _, beam in beams.iterrows())
            for future in concurrent.futures.as_completed(future_beam):
                future.result()

    print(f"DONE in {time.time() - tic:.1f} seconds")


if __name__ == "__main__":
    main()
