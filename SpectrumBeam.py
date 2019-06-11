#!/usr/bin/env python
"""
Executes Transcar with user-defined input particle flux spectrum

python SpectrumBeam.py flux.csv /tmp/spectrumout

python SpectrumBeam.py flux.csv c:/temp/spectrumout
"""
import logging
from pathlib import Path
import numpy as np
from pandas import read_csv
from argparse import ArgumentParser
import signal
#
import transcar.base as transcar


def main():
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    p = ArgumentParser(description='parallel instance transcar runner')
    p.add_argument('fluxfn',  help='particle flux csv filename')
    p.add_argument('rodir', help='output directory')
    p.add_argument('-Q0', help='Assumed particle flux', type=float, default=70114000000.0)
    p.add_argument('-infn', help='energy bin CSV file', default='BT_E1E2prev.csv')
    p.add_argument('--msgfn', help='file to write transcar messages to', default='transcar.log')
    p.add_argument('--errfn', help='file to write transcar Errors to', default='transcarError.log')
    p = p.parse_args()

    rodir = Path(p.rodir).expanduser()
    infn = Path(p.infn).expanduser()
    fluxfn = Path(p.fluxfn).expanduser()

    params = {'rodir': rodir,
              'Q0': p.Q0,
              'msgfn': p.msgfn,
              'errfn': p.errfn
              }

    logfn = rodir/'Beams.log'

    print('logging to', logfn)
    rodir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=logfn,
                        filemode='a',
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.DEBUG)

    beams = read_csv(infn, header=None, names=['E1', 'E2', 'pr1', 'pr2'])
    beams['flux'] = np.loadtxt(fluxfn, delimiter=',')

    transcar.beam_spectrum_arbiter(beams, params)


if __name__ == '__main__':
    main()
