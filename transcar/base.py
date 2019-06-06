import subprocess
from pathlib import Path
import logging
import pandas
from typing import Dict, Any
# %% constants dictacted by legacy Fortran code
from .io import setup_dirs, setup_monoprec, setup_spectrum_prec, transcaroutcheck, TRANSCAREXE


def beam_spectrum_arbiter(beam: pandas.DataFrame,
                          P: Dict[str, Any]):
    """
    run beam with user-defined flux spectrum
    """

    isok = run_spectrum(beam, P)

    if isok:
        print(f'OK {beam["E1"]:.1f} eV')
    else:
        logging.warning(f'retrying beam{beam["E1"]:.1f}')
        isok = run_spectrum(beam, P)
        if not isok:
            logging.error(f'failed on beam{beam["E1"]:.1f} on 2nd try, aborting')


def run_spectrum(beam: pandas.DataFrame, P: Dict[str, Any]) -> bool:
    """
    Run beam spectrum
    """
# %% copy the Fortran static init files to this directory (simple but robust)
    datinp, odir = setup_dirs(P['rodir'])
    setup_spectrum_prec(odir, datinp, beam)
# %% run the compiled executable
    runTranscar(odir, P['errfn'], P['msgfn'])
# %% check output trivially
    isok = transcaroutcheck(odir, P['errfn'])

    return isok


def mono_beam_arbiter(beam: Dict[str, float], P: Dict[str, Any]):
    """
    run monoenergetic beam
    """
    if isinstance(beam, pandas.Series):
        beam = beam.to_dict()

    isok = run_monobeam(beam, P)

    if isok:
        print(f'OK {beam["E1"]:.1f} eV')
    else:
        logging.warning(f'retrying beam{beam["E1"]}')
        isok = run_monobeam(beam, P)
        if not isok:
            logging.error(f'failed on beam{beam["E1"]} on 2nd try, aborting')


def run_monobeam(beam: Dict[str, float], P: Dict[str, Any]) -> bool:
    """Run a particular beam energy vs. time"""
# %% copy the Fortran static init files to this directory (simple but robust)
    datinp, odir = setup_dirs(P['rodir'] / f'beam{beam["E1"]:.1f}')
    setup_monoprec(odir, datinp, beam, P['Q0'])
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
