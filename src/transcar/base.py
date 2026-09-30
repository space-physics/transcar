import subprocess
from pathlib import Path
import logging
import typing as T
import os

import pandas

from .io import setup_dirs, setup_monoprec, setup_spectrum_prec, transcaroutcheck, transcar_paths


def beam_spectrum_arbiter(beam: pandas.DataFrame, P: dict[str, T.Any]):
    """
    run beam with user-defined flux spectrum
    """
    odir = P["rodir"]

    print("Running Transcar single-threaded.")
    print(
        odir / P["msgfn"],
        "logs the simulation output text, watch this file to see simulation progress.",
    )

    if run_spectrum(beam, P):
        print("OK: Transcar complete")
        return

    raise RuntimeError(f"Transcar run failed. See {odir / P['errfn']} for clues")


def run_spectrum(beam: pandas.DataFrame, P: dict[str, T.Any]) -> bool:
    """
    Run beam spectrum
    """
    # %% copy the Fortran static init files to this directory (simple but robust)
    datinp, odir = setup_dirs(P["rodir"], P)
    setup_spectrum_prec(odir, datinp, beam)
    # %% run the compiled executable
    isok = runTranscar(odir, P["errfn"], P["msgfn"])
    # %% check output trivially
    return isok and transcaroutcheck(odir, P["errfn"])


def mono_beam_arbiter(beam: dict[str, float], P: dict[str, T.Any]):
    """
    run monoenergetic beam
    """
    if isinstance(beam, pandas.Series):
        beam = beam.to_dict()

    isok = run_monobeam(beam, P)

    if isok:
        print(f"OK {beam['E1']:.1f} eV")
    else:
        logging.warning(f"retrying beam{beam['E1']:.1f}")
        isok = run_monobeam(beam, P)
        if not isok:
            logging.error(f"failed on beam{beam['E1']:.1f} on 2nd try, aborting")


def run_monobeam(beam: dict[str, float], P: dict[str, T.Any]) -> bool:
    """Run a particular beam energy vs. time"""
    # %% copy the Fortran static init files to this directory (simple but robust)
    datinp, odir = setup_dirs(P["rodir"] / f"beam{beam['E1']:.1f}", P)
    setup_monoprec(odir, datinp, beam, P["Q0"])
    # %% run the compiled executable
    isok = runTranscar(odir, P["errfn"], P["msgfn"])
    # %% check output trivially
    return isok and transcaroutcheck(odir, P["errfn"])


def runTranscar(odir: Path, errfn: Path, msgfn: Path) -> bool:
    """actually run Transcar exe"""
    odir = Path(odir).expanduser().resolve()  # MUST have resolve()!!

    exe = transcar_paths()["transconvec"]

    err_file = odir / errfn
    out_file = odir / msgfn

    with err_file.open("w") as ferr, out_file.open("w") as fout:
        ret = subprocess.run(exe, cwd=odir, stdout=fout, stderr=ferr)

    if ret.returncode != 0:
        logging.error(f"{odir.name} error code {ret.returncode} see {err_file}")
        if os.name == "nt":
            match ret.returncode:
                case 3221225725:
                    logging.error(f"{exe} stack overflow indicated")

    return ret.returncode == 0
