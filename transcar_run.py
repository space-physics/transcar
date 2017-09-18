#!/usr/bin/env python
"""
to run in parallel on multiple PCs, use GNU Parallel to take the place of this serial script.
"""
from pathlib import Path
from pandas import read_csv
#
from transcar import setuptranscario,setupPrecipitation,runTranscar,transcaroutcheck
# %%
def runbeam(rodir:Path, Q0:float, beam, logfn:Path, errfn:Path):
# %% copy the Fortran static init files to this directory (simple but robust)
    datinp,odir = setuptranscario(rodir, beam['E1'])
    setupPrecipitation(odir, datinp, beam, Q0)
# %% run the compiled executable
    runTranscar(odir, errfn, logfn)
#%% check output trivially
    transcaroutcheck(odir, errfn)


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='parallel instance transcar runner')
    p.add_argument('rodir',help='root of beam directory')
    p.add_argument('-Q0',help='Assumed particle flux',type=float,default=70114000000.0)
    p.add_argument('-infn',help='energy bin CSV file',default='transcar/BT_E1E2prev.csv')
    p.add_argument('--msgfn',help='file to write transcar messages to',default='transcar.log')
    p.add_argument('--errfn',help='file to write transcar Errors to',default='transcarError.log')
    p = p.parse_args()

    rodir = Path(p.rodir).expanduser()
    infn = Path(p.infn).expanduser()

    beams = read_csv(infn, header=None, names=['E1','E2','pr1','pr2'])

    for i,beam in beams.iterrows():
        runbeam(rodir, p.Q0, beam, p.msgfn, p.errfn)