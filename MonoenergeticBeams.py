#!/usr/bin/env python
"""
Executes Transcar to output "monoenergetic" electron beams.
Optionally, in parallel.
"""
import concurrent.futures
import logging
from pathlib import Path
from pandas import read_csv
import os
from argparse import ArgumentParser
#
import transcar.base as transcar
import signal


def main():
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    p = ArgumentParser(description='parallel instance transcar runner')
    p.add_argument('rodir', help='root of beam directory')
    p.add_argument('-Q0', help='Assumed particle flux', type=float, default=70114000000.0)
    p.add_argument('-infn', help='energy bin CSV file', default='BT_E1E2prev.csv')
    p.add_argument('--msgfn', help='file to write transcar messages to', default='transcar.log')
    p.add_argument('--errfn', help='file to write transcar Errors to', default='transcarError.log')
    p.add_argument('-np', help='number of concurrent processes', type=int,
                   default=max(1, os.cpu_count()-1))
    p = p.parse_args()

    rodir = Path(p.rodir).expanduser()
    infn = Path(p.infn).expanduser()

    params = {'rodir': rodir,
              'Q0': p.Q0,
              'msgfn': p.msgfn,
              'errfn': p.errfn
              }

    rodir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=rodir/'Beams.log',
                        filemode='a',
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.DEBUG)

    beams = read_csv(infn, header=None, names=['E1', 'E2', 'pr1', 'pr2'])

    print('using', p.np, 'concurrent Transcar runs')
# %%
    with concurrent.futures.ThreadPoolExecutor(max_workers=p.np) as executor:
        future_beam = (executor.submit(transcar.mono_beam_arbiter, beam, params) for _, beam in beams.iterrows())
        for future in concurrent.futures.as_completed(future_beam):
            future.result()


if __name__ == '__main__':
    main()
