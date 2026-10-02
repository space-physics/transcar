from pathlib import Path
from datetime import datetime, timedelta
import logging
import shutil
import pandas
import os
import collections
import typing as T
import functools


@functools.cache
def transcar_paths() -> dict[str, Path]:

    root = Path(__file__).resolve().parents[2]
    tcpaths: list[str | Path | None] = [os.getenv("TRANSCAR_ROOT"), root, root / "build"]
    for p in tcpaths:
        if exe := shutil.which("transconvec", path=p):
            break

    if exe is None:
        raise FileNotFoundError(f"could not find transconvec executable under {root}")

    paths = {
        "transconvec": Path(exe),
        "root": root,
        "datain": root / "dir.input",
        "datadir": root / "dir.data",
        "datcar": root / "dir.input/DATCAR",
        "finish_status": root / "finish.status",
        "precip": Path("dir.input/precinput.asc"),
    }
    return paths


def cp_parents(files: list[Path], target_dir: Path, origin: Path | None = None) -> None:
    """
    inputs
    ------
    files: list of files to copy (source)
    target_dir: directory to copy into
    origin: non-relative path of source, that needs to be trimmed.

    This function requires Python >= 3.6.

    This acts like bash cp --parents in Python
    inspiration from
    http://stackoverflow.com/questions/15329223/copy-a-file-into-a-directory-with-its-original-leading-directories-appended

    example
    source: /tmp/e/f
    dest: /tmp/a/b/c/d/
    result: /tmp/a/b/c/d/tmp/e/f

    cp_parents('/tmp/a/b/c/d/boo','/tmp/e/f')
    cp_parents('x/hi','/tmp/e/f/g')  --> copies ./x/hi to /tmp/e/f/g/x/hi
    """
    # %% make list if it's a string
    if isinstance(files, (str, Path)):
        files = [files]
    # %% cleanup user
    files = [Path(f).expanduser() for f in files]  # relative path or absolute path is fine
    target_dir = Path(target_dir).expanduser()
    # %% work
    for f in files:
        if origin:
            fsource = f.relative_to(origin).parent
        else:
            fsource = f.parent

        newpath = target_dir / fsource  # to make it work like cp --parents, copying absolute paths if specified
        newpath.mkdir(parents=True, exist_ok=True)
        shutil.copy2(f, newpath)


def transcaroutcheck(odir: Path, logfile: str | Path = "transcar.log", ok: str = "fin normale") -> bool:
    """
    checks output for the normal completion message "fin normale"

    Compilers format the normal completion message differently and may print
    IEEE flag notes after, e.g. "STOP fin normale" or "Fortran STOP: fin normale".
    """
    tpaths = transcar_paths()

    isok = False
    fok = odir / tpaths["finish_status"]
    try:
        with (odir / logfile).open("r", errors="replace") as ferr:
            last = collections.deque(ferr, 1)[0]

        if ok in last:
            isok = True
            fok.write_text("true")
            logging.info(f"{odir} completed successfully.")
        else:
            fok.write_text("false")
            logging.warning(f"{odir} ended sim early:  {last}")
    except IOError as e:
        logging.error(f"problem reading transcar output.  {e}")
    except IndexError as e:  # empty file
        fok.write_text("false")
        logging.warning(f"{odir} Transcar may not have finished the sim   {e}")

    return isok


def setup_dirs(odir: Path, params: dict[str, T.Any]) -> tuple[dict[str, T.Any], Path]:
    """
    prepare output directory for a beam
    """

    tpaths = transcar_paths()

    datcar = params["datcar"] if "datcar" in params else tpaths["datcar"]

    inp = readTranscarInput(datcar)
    # %% cleanup bad runs
    out = odir / "dir.output"
    out.mkdir(parents=True, exist_ok=True)
    for fn in ("ediffnumflux.dat", "emissions.dat", "flux.output", "transcar_output"):
        if (out / fn).is_file():
            (out / fn).unlink()
    # %% move files where needed for this instantiation
    flist = []
    # precfn is NOT included here!

    # decided to just use the original executable rather than make so many copies
    # flist += [tpaths["transconvec"]]

    # does not operate consistently (segfault) if not in same directory, verified by hand June 2019.
    flist += [tpaths["datain"] / inp["precfile"], tpaths["datadir"] / "type"]
    flist += [tpaths["datadir"] / "dir.linux/dir.geomag" / s for s in ["data_geom.bin", "igrf90.dat", "igrf90s.dat"]]
    flist += [tpaths["datadir"] / "dir.linux/dir.projection/varpot.dat"]
    # transcar sigsegv on val_fit_ if FELTRANS is blank!
    flist += [tpaths["datadir"] / "dir.linux/dir.cine" / s for s in ["DATDEG", "DATFEL", "DATTRANS", "flux.flag", "FELTRANS"]]
    flist += [tpaths["datadir"] / "dir.linux/dir.cine/dir.euvac/EUVAC.dat"]
    flist += [tpaths["datadir"] / "dir.linux/dir.cine/dir.seff" / s for s in ["crsb8", "crsphot1.dat", "rdtb8"]]

    cp_parents(flist, odir, tpaths["root"])
    # may have uniquely named input DATCAR, that always needs to be in output dir as
    # dir.input/DATCAR due to legacy Fortran hardcoding
    shutil.copy2(datcar, odir / "dir.input/DATCAR")

    return inp, odir


def setup_monoprec(odir: Path, inp: dict[str, T.Any], beam: dict[str, float], flux0: float) -> None:
    """
    write dir.input/precinput.asc for monoenergetic beam case
    """
    tpaths = transcar_paths()

    ofn = Path(odir).expanduser() / tpaths["precip"]

    E1 = beam["E1"]
    E2 = beam["E2"]
    pr1 = beam["pr1"]
    pr2 = beam["pr2"]

    dE = E2 - E1
    Esum = E2 + E1
    flux = flux0 / 0.5 / Esum / dE
    Elow = E1 - 0.5 * (E1 - pr1)
    Ehigh = E2 - 0.5 * (E2 - pr2)

    precout = "\n".join(
        (
            f"{inp['precipstartsec']}",
            f"{Elow} {flux}",
            f"{Ehigh} -1.0",
            f"{inp['precipendsec']}",
            "-1.0 -1.0",
        )
    )

    print("writing", ofn)

    ofn.write_text(precout)


def setup_spectrum_prec(odir: Path, inp: dict[str, T.Any], beam: pandas.DataFrame) -> None:
    """
    write dir.input/precinput.asc for beam with shaped differential number flux
    """

    tpaths = transcar_paths()

    ofn = Path(odir).expanduser() / tpaths["precip"]

    dat = str(inp["precipstartsec"])

    for _, ebin in beam.iterrows():
        Elow, Ehigh, flux = compute_Ebin(ebin)

        dat += f"\n{Elow:.3f} {flux:.3f}"

    dat += f"\n{Ehigh:.3f} -1.0"

    dat += f"\n{inp['precipendsec']}\n-1.0 -1.0"

    print("writing", ofn)

    ofn.write_text(dat)


def compute_Ebin(ebin: pandas.Series) -> tuple[float, float, float]:
    E1 = ebin["E1"]
    E2 = ebin["E2"]
    pr1 = ebin["pr1"]
    pr2 = ebin["pr2"]

    dE = E2 - E1
    Esum = E2 + E1
    flux = ebin["flux"] / 0.5 / Esum / dE
    Elow = E1 - 0.5 * (E1 - pr1)
    Ehigh = E2 - 0.5 * (E2 - pr2)

    return Elow, Ehigh, flux


def readTranscarInput(infn: Path) -> dict[str, T.Any]:
    """
    The transcar input file is indexed by line number --this is what the Fortran
      #  code of transcar does, and it's what we do here as well.
    """
    infn = Path(infn).expanduser()
    hd: dict[str, T.Any] = {}

    with infn.open("r") as f:
        hd["kiappel"] = int(f.readline().split()[0])
        hd["precfile"] = f.readline().split()[0]
        hd["dtsim"] = float(f.readline().split()[0])  # "dto"
        hd["dtfluid"] = float(f.readline().split()[0])  # "sortie"
        hd["iyd_ini"] = int(f.readline().split()[0])
        hd["dayofsim"] = datetime.strptime(str(hd["iyd_ini"]), "%Y%j")
        hd["simstartUTCsec"] = float(f.readline().split()[0])  # "tempsini"
        hd["simlengthsec"] = float(f.readline().split()[0])  # "tempslim"
        hd["jpreci"] = int(f.readline().split()[0])
        # transconvec calls the next two latgeo_ini, longeo_ini
        hd["latgeo_ini"], hd["longeo_ini"] = [float(a) for a in f.readline().split(None)[0].split(",")]
        hd["tempsconv_1"] = float(f.readline().split()[0])  # from transconvec, time before precip
        hd["tempsconv"] = float(f.readline().split()[0])  # from transconvec, time after precip
        hd["step"] = float(f.readline().split()[0])
        hd["dtkinetic"] = float(f.readline().split()[0])  # transconvec calls this "postinto"
        hd["vparaB"] = float(f.readline().split()[0])
        hd["f107ind"] = float(f.readline().split()[0])
        hd["f107avg"] = float(f.readline().split()[0])
        hd["apind"] = float(f.readline().split()[0])
        hd["convecEfieldmVm"] = float(f.readline().split()[0])
        hd["cofo"] = float(f.readline().split()[0])
        hd["cofn2"] = float(f.readline().split()[0])
        hd["cofo2"] = float(f.readline().split()[0])
        hd["cofn"] = float(f.readline().split()[0])
        hd["cofh"] = float(f.readline().split()[0])
        hd["etopflux"] = float(f.readline().split()[0])
        hd["precinfn"] = f.readline().split()[0]
        hd["precint"] = int(f.readline().split()[0])
        hd["precext"] = int(f.readline().split()[0])
        hd["precipstartsec"] = float(f.readline().split()[0])
        hd["precipendsec"] = float(f.readline().split()[0])

        # %% derived parameters not in datcar file
        hd["tstartSim"] = hd["dayofsim"] + timedelta(seconds=hd["simstartUTCsec"])
        hd["tendSim"] = hd["dayofsim"] + timedelta(seconds=hd["simlengthsec"])
        # TODO verify this isn't added to start
        hd["tstartPrecip"] = hd["dayofsim"] + timedelta(seconds=hd["precipstartsec"])
        hd["tendPrecip"] = hd["dayofsim"] + timedelta(seconds=hd["precipendsec"])

    return hd
