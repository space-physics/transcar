#!/usr/bin/env python
"""
Executes Transcar to output "monoenergetic" electron beams.
Optionally, in parallel.
"""
import multiprocessing
import concurrent.futures
import logging
from pathlib import Path
from pandas import read_csv
#
from transcar import setuptranscario,setupPrecipitation,runTranscar,transcaroutcheck
#
Ncpu=multiprocessing.cpu_count() // 2  # //2 makes one thread per CPU for 2 thread Intel Hyperthreading

# %%
def runbeam(rodir:Path, Q0:float, beam, logfn:Path, errfn:Path):
# %% copy the Fortran static init files to this directory (simple but robust)
    datinp,odir = setuptranscario(rodir, beam['E1'])
    setupPrecipitation(odir, datinp, beam, Q0)
# %% run the compiled executable
    runTranscar(odir, errfn, logfn)
#%% check output trivially
    isok = transcaroutcheck(odir, errfn)

    return isok


def okmain(beam):
    beam = beam[1]
    #rint(beam)
    isok = runbeam(rodir, p.Q0, beam, p.msgfn, p.errfn)


    if not isok:
        logging.warning(f'retrying beam{beam["E1"]}')
        isok = runbeam(rodir, p.Q0, beam, p.msgfn, p.errfn)
        if not isok:
            logging.error(f'failed on beam{beam["E1"]} on 2nd try, aborting')


if __name__ == '__main__':
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    from argparse import ArgumentParser
    p = ArgumentParser(description='parallel instance transcar runner')
    p.add_argument('rodir',help='root of beam directory')
    p.add_argument('-Q0',help='Assumed particle flux',type=float,default=70114000000.0)
    p.add_argument('-infn',help='energy bin CSV file',default='BT_E1E2prev.csv')
    p.add_argument('--msgfn',help='file to write transcar messages to', default='transcar.log')
    p.add_argument('--errfn',help='file to write transcar Errors to', default='transcarError.log')
    p = p.parse_args()

    rodir = Path(p.rodir).expanduser()
    infn = Path(p.infn).expanduser()

    rodir.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(filename=rodir/'Beams.log',
                            filemode='a',
                            format='%(asctime)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)

    beams = read_csv(infn, header=None, names=['E1','E2','pr1','pr2'])
    abeam = beams.values

    print(Ncpu,'workers')

    with concurrent.futures.ProcessPoolExecutor(max_workers=Ncpu) as executor:
        executor.map(okmain, beams.iterrows())

