import subprocess
from pathlib import Path
import logging
import pandas
from typing import Dict, Any
# %% constants dictacted by legacy Fortran code
from .io import setuptranscario, setupPrecipitation, transcaroutcheck, TRANSCAREXE


def iterbeams(beam: Dict[str, float], P: Dict[str, Any]):
    if isinstance(beam, pandas.Series):
        beam = beam.to_dict()

    isok = runbeam(beam, P)

    if isok:
        print(f'OK {beam["E1"]:.0f} eV')
    else:
        logging.warning(f'retrying beam{beam["E1"]}')
        isok = runbeam(beam, P)
        if not isok:
            logging.error(f'failed on beam{beam["E1"]} on 2nd try, aborting')


def runbeam(beam: Dict[str, float], P: Dict[str, Any]) -> bool:
    """Run a particular beam energy vs. time"""
# %% copy the Fortran static init files to this directory (simple but robust)
    datinp, odir = setuptranscario(P['rodir'], beam['E1'])
    setupPrecipitation(odir, datinp, beam, P['Q0'])
# %% run the compiled executable
    runTranscar(odir, P['errfn'], P['msgfn'])
# %% check output trivially
    isok = transcaroutcheck(odir, P['errfn'])

    return isok


def runTranscar(odir: Path, errfn: Path, msgfn: Path):
    """actually run Transcar exe"""
    odir = Path(odir).expanduser().resolve()  # MUST have resolve()!!

    with (odir/errfn).open('w') as ferr, (odir/msgfn).open('w') as fout:
        subprocess.run(str(TRANSCAREXE), cwd=odir, stdout=fout, stderr=ferr)
