from __future__ import division
import logging
from collections import deque
from os.path import join,abspath
from os import makedirs
from subprocess import Popen
#
import sys
sys.path.extend(['../../transcar-utils','../../hist-utils'])
from parseTranscar import readTranscarInput
from cp_parents import cp_parents
from empty_file import empty_file

transcarexe = 'transconvec_13.op.out'
fileok = 'finish.status'
# hard-coded in Fortran
datcar = 'dir.input/DATCAR'
precfn = 'dir.input/precinput.dat'

def runTranscar(odir,errfn,msgfn):
    with open(join(odir,errfn),'w')  as ferr, open(join(odir,msgfn),'w') as fout:
    #we need cwd feature of Popen that call() doesn't have
        exe = abspath(join(odir,transcarexe))
        try:
            proc = Popen(args=exe, cwd=odir, stdout=fout, stderr=ferr, shell = False)
            out,err = proc.communicate() #this makes popen block (we want this)
            #print(out,err) #will be none since we didn't use PIPE
        except IOError as e:
            exit(e)   
 
def transcaroutcheck(odir,errfn):
    fok = join(odir,fileok)
    try:
      with open(join(odir,errfn),'r') as ferr:
        last = deque(ferr,1)[0].rstrip('\n')
        with open(fok,'w') as f:
            if last == 'STOP fin normale':
                f.write('true')
            else:
                print('transcaroutcheck: {} got unexpected return value, transcar may not have finished the sim'.format(odir))
                print(last)
                f.write('false')
    except IOError as e:
        print('transcaroutcheck: problem reading transcar output.  {}'.format(e) )
                
def setuptranscario(rodir,beamEnergy):
    inp = readTranscarInput(datcar) 
    
    odir = join(rodir,'beam{:.1f}'.format(beamEnergy))
    
    try: 
        makedirs(join(odir,'dir.output'))
    except OSError: 
        exit('{} directory already existed, aborting'.format(beamEnergy))
    
    flist = [datcar, join('dir.input',inp['precfile']), 'dir.data/type',transcarexe]
    flist.extend(['dir.data/dir.linux/dir.geomag/' +s  for s in ['data_geom.bin','igrf90.dat','igrf90s.dat']])
    flist.append('dir.data/dir.linux/dir.projection/varpot.dat')
    flist.extend(['dir.data/dir.linux/dir.cine/' +s  for s in ['DATDEG','DATFEL','DATTRANS','flux.flag']])
    flist.append('dir.data/dir.linux/dir.cine/dir.euvac/EUVAC.dat')
    flist.extend(['dir.data/dir.linux/dir.cine/dir.seff/' +s  for s in ['crsb8','crsphot1.dat','rdtb8']])
    cp_parents(flist,odir)
    
    empty_file(join(odir,'dir.data/dir.linux/dir.cine/FELTRANS'))
    
    return inp,odir
        
def setupPrecipitation(odir,inp,E1,E2,pr1,pr2, flux0):
    dE = E2-E1
    Esum = E2+E1
    flux = flux0 / 0.5 / Esum / dE
    Elow = E1 - 0.5*(E1 - pr1)
    Ehigh =E2 - 0.5*(E2 - pr2)
    
    precout = '\n'.join(('{}','{} {}','{} -1.0','{}','-1.0 -1.0')).format(
              inp['precipstartsec'], Elow, flux, Ehigh, inp['precipendsec'])
    
    with open(join(odir,precfn),'w') as f:
        f.write(precout)
#%%    
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='parallel instance transcar runner')
    p.add_argument('rodir',help='root of beam directory',type=str)
    p.add_argument('flux0',help='flux0',type=float)
    p.add_argument('E1',help='energy in eV',type=float)
    p.add_argument('E2',help='energy in eV',type=float)
    p.add_argument('pr1',help='pr1',type=float)
    p.add_argument('pr2',help='pr2',type=float)
    p.add_argument('--msgfn',help='file to write transcar messages to',type=str,default='transcar.log')
    p.add_argument('--errfn',help='file to write transcar Errors to',type=str,default='transcarError.log')
    p = p.parse_args()

    datinp,odir = setuptranscario(p.rodir,p.E1)
    setupPrecipitation(odir,datinp, p.E1, p.E2, p.pr1, p.pr2, p.flux0)
    
    runTranscar(odir, p.errfn, p.msgfn)
#%% check output trivially
    transcaroutcheck(odir,p.errfn) 
    