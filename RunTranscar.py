#!/usr/bin/env python
"""
Executes Transcar to output "monoenergetic" electron beams.
Optionally, in parallel.
"""
import concurrent.futures
import logging
import itertools
from pathlib import Path
from pandas import read_csv
#
import transcar


def main():
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    from argparse import ArgumentParser
    p = ArgumentParser(description='parallel instance transcar runner')
    p.add_argument('rodir', help='root of beam directory')
    p.add_argument('-Q0', help='Assumed particle flux', type=float, default=70114000000.0)
    p.add_argument('-infn', help='energy bin CSV file', default='BT_E1E2prev.csv')
    p.add_argument('--msgfn', help='file to write transcar messages to', default='transcar.log')
    p.add_argument('--errfn', help='file to write transcar Errors to', default='transcarError.log')
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
# %%
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(transcar.iterbeams,
                     beams.iterrows(), itertools.repeat(params),
                     timeout=600)


if __name__ == '__main__':
    main()
