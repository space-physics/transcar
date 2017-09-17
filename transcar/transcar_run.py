#!/usr/bin/env python
"""
This file is normally called from GNU Parallel
"""
from pathlib import Path
import logging
from collections import deque
from subprocess import Popen
#
from transcarread import readTranscarInput
from transcar import cp_parents
#
transcarexe = 'transconvec_13.op.out'
fileok = 'finish.status'
# hard-coded in Fortran
datcar = 'dir.input/DATCAR'
precfn = 'dir.input/precinput.dat'

def runTranscar(odir,errfn,msgfn):
    odir = Path(odir).expanduser()
    with (odir/errfn).open('w')  as ferr, (odir/msgfn).open('w') as fout:
    #we need cwd feature of Popen that call() doesn't have
        exe = (odir/transcarexe).resolve()

        proc = Popen(args=exe, cwd=odir, stdout=fout, stderr=ferr, shell=False)
        out,err = proc.communicate() #this makes popen block (we want this)
        #print(out,err) #will be none since we didn't use PIPE


def transcaroutcheck(odir,errfn):
    fok = odir/fileok
    try:
      with (odir/errfn).open('r') as ferr:
        last = deque(ferr,1)[0].rstrip('\n')
        with open(fok,'w') as f:
            if last == 'STOP fin normale':
                f.write('true')
            else:
                f.write('false')
                logging.warn(f'transcaroutcheck: {odir} got unexpected return value, transcar may not have finished the sim')
                raise AttributeError(last)
    except (IOError) as e:
        logging.error(f'transcaroutcheck: problem reading transcar output.  {e}' )
    except IndexError as e:
        with open(fok,'w') as f:
            f.write('false')
            logging.warn(f'transcaroutcheck: {odir} got unexpected return value, transcar may not have finished the sim')


def setuptranscario(rodir,beamEnergy):
    inp = readTranscarInput(datcar)

    odir = rodir/'beam{:.1f}'.format(beamEnergy)

    try:
        (odir/'dir.output').mkdir(parents=True,exist_ok=True)
    except OSError:
        raise ValueError(f'{beamEnergy} directory already existed, aborting')
# %% move files where needed for this instantiation
    flist = [datcar, 'dir.input'/inp['precfile'], 'dir.data/type',transcarexe]
    flist.extend(['dir.data/dir.linux/dir.geomag/' +s  for s in ['data_geom.bin','igrf90.dat','igrf90s.dat']])
    flist.append('dir.data/dir.linux/dir.projection/varpot.dat')
    #transcar sigsegv on val_fit_ if FELTRANS is blank!
    flist.extend(['dir.data/dir.linux/dir.cine/' +s  for s in ['DATDEG','DATFEL','DATTRANS','flux.flag','FELTRANS']])
    flist.append('dir.data/dir.linux/dir.cine/dir.euvac/EUVAC.dat')
    flist.extend(['dir.data/dir.linux/dir.cine/dir.seff/' +s  for s in ['crsb8','crsphot1.dat','rdtb8']])

    cp_parents(flist,odir)


    return inp,odir

def setupPrecipitation(odir,inp,E1,E2,pr1,pr2, flux0):
    dE =   E2-E1
    Esum = E2+E1
    flux = flux0 / 0.5 / Esum / dE
    Elow = E1 - 0.5*(E1 - pr1)
    Ehigh= E2 - 0.5*(E2 - pr2)

    precout = '\n'.join(('{}','{} {}','{} -1.0','{}','-1.0 -1.0')).format(
              inp['precipstartsec'], Elow, flux, Ehigh, inp['precipendsec'])

    with (odir/precfn).open('w') as f:
        f.write(precout)
#%%
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='parallel instance transcar runner')
    p.add_argument('rodir',help='root of beam directory')
    p.add_argument('flux0',help='flux0',type=float)
    p.add_argument('E1',help='energy in eV',type=float)
    p.add_argument('E2',help='energy in eV',type=float)
    p.add_argument('pr1',help='pr1',type=float)
    p.add_argument('pr2',help='pr2',type=float)
    p.add_argument('--msgfn',help='file to write transcar messages to',default='transcar.log')
    p.add_argument('--errfn',help='file to write transcar Errors to',default='transcarError.log')
    p = p.parse_args()

    datinp,odir = setuptranscario(p.rodir,p.E1)
    setupPrecipitation(odir,datinp, p.E1, p.E2, p.pr1, p.pr2, p.flux0)

    runTranscar(odir, p.errfn, p.msgfn)
#%% check output trivially
    transcaroutcheck(odir,p.errfn)
