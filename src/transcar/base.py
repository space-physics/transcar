import subprocess
from pathlib import Path
import logging
import typing as T
import signal
import os

import pandas

from .io import setup_dirs, setup_monoprec, setup_spectrum_prec, transcaroutcheck, transcar_paths


def describe_returncode(rc: int) -> str:
    """a short reason for subprocess failure."""

    _win_code = {
        0xC0000135: "DLL not found - missing compiler runtime libraries?",
        0xC00000FD: "stack overflow",
        0xC0000005: "access violation (segfault)",
        0xC0000409: "stack buffer overrun",
        0xC0000374: "heap corruption",
        3: "abort (CRT abort, exit code 3)",
    }

    if os.name == "nt":
        return _win_code.get(rc, f"error code {rc}")
    if rc < 0:
        try:
            sig = signal.Signals(-rc)
        except ValueError:
            return f"killed by signal {-rc}"

        return f"killed by {sig.name}"

    return f"exited with status {rc}"


def beam_spectrum_arbiter(beam: pandas.DataFrame, P: dict[str, T.Any]) -> bool:
    """
    run beam with user-defined flux spectrum
    """
    odir = P["rodir"]

    print("Running Transcar single-threaded.")
    print(
        odir / P["msgfn"],
        "logs the simulation output text, watch this file to see simulation progress.",
    )

    if ok := run_spectrum(beam, P):
        print("OK: Transcar complete")
    else:
        logging.error(f"Transcar run failed. See {odir / P['errfn']} for clues")

    return ok


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


def mono_beam_arbiter(beam: dict[str, float], P: dict[str, T.Any]) -> bool:
    """
    run monoenergetic beam
    """
    if isinstance(beam, pandas.Series):
        beam = beam.to_dict()

    if isok := run_monobeam(beam, P):
        print(f"OK {beam['E1']:.1f} eV")
    else:
        logging.warning(f"retrying beam{beam['E1']:.1f}")
        isok = run_monobeam(beam, P)
        if not isok:
            logging.error(f"failed on beam{beam['E1']:.1f} on 2nd try, aborting")

    return isok


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
        logging.error(f"{odir.name}: error code {ret.returncode}\nError logfile: {err_file}\n{describe_returncode(ret.returncode)}")

    return ret.returncode == 0
