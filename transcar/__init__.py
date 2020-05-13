from pathlib import Path
import shutil
import logging
from collections import deque
import subprocess

#
from transcarread import readTranscarInput

# %% constants dictacted by legacy Fortran code
R = Path(__file__).resolve().parents[1]
exe = shutil.which("transconvec_13.op.out", path=str(R))
fileok = "finish.status"
# hard-coded in Fortran
din = Path("dir.input")
ddat = Path("dir.data")
DATCAR = din / "DATCAR"
precfn = din / "precinput.dat"


def cp_parents(files, target_dir: Path):
    """
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
    files = (Path(f).expanduser() for f in files)  # relative path or absolute path is fine
    target_dir = Path(target_dir).expanduser()
    # %% work
    for f in files:
        newpath = target_dir / f.parent  # to make it work like cp --parents, copying absolute paths if specified
        newpath.mkdir(parents=True, exist_ok=True)
        shutil.copy2(f, newpath)


def runTranscar(odir: Path, errfn: Path, msgfn: Path):
    odir = Path(odir).expanduser()

    with (odir / errfn).open("w") as ferr, (odir / msgfn).open("w") as fout:
        subprocess.run(args=[exe], cwd=odir, stdout=fout, stderr=ferr, shell=False, universal_newlines=True)

        # current code error checks rely on serial operation.
        # subprocess.Popen(args=exe, cwd=odir, stdout=fout, stderr=ferr, shell=False)


def transcaroutcheck(odir, errfn):
    fok = odir / fileok
    try:
        with (odir / errfn).open("r") as ferr:
            last = deque(ferr, 1)[0].rstrip("\n")

        if last == "STOP fin normale":
            fok.write_text("true")
        else:
            fok.write_text("false")
            logging.warn(f"transcaroutcheck: {odir} got unexpected return value, transcar may not have finished the sim")
            raise AttributeError((odir / errfn).read_text())
    except (IOError) as e:
        logging.error(f"transcaroutcheck: problem reading transcar output.  {e}")
    except IndexError:
        with open(fok, "w") as f:
            f.write("false")
            logging.warn(f"transcaroutcheck: {odir} got unexpected return value, transcar may not have finished the sim")


def setuptranscario(rodir: Path, beamEnergy: float):
    inp = readTranscarInput(DATCAR)

    odir = rodir / f"beam{beamEnergy:.1f}"

    (odir / "dir.output").mkdir(parents=True, exist_ok=True)
    # %% move files where needed for this instantiation
    flist = [DATCAR, din / inp["precfile"], ddat / "type"]
    flist += [ddat / "dir.linux/dir.geomag" / s for s in ["data_geom.bin", "igrf90.dat", "igrf90s.dat"]]
    flist += [ddat / "dir.linux/dir.projection/varpot.dat"]
    # transcar sigsegv on val_fit_ if FELTRANS is blank!
    flist += [ddat / "dir.linux/dir.cine" / s for s in ["DATDEG", "DATFEL", "DATTRANS", "flux.flag", "FELTRANS"]]
    flist += [ddat / "dir.linux/dir.cine/dir.euvac/EUVAC.dat"]
    flist += [ddat / "dir.linux/dir.cine/dir.seff" / s for s in ["crsb8", "crsphot1.dat", "rdtb8"]]

    cp_parents(flist, odir)

    return inp, odir


def setupPrecipitation(odir, inp, beam, flux0):
    ofn = odir / precfn

    E1 = beam["E1"]
    E2 = beam["E2"]
    pr1 = beam["pr1"]
    pr2 = beam["pr2"]

    dE = E2 - E1
    Esum = E2 + E1
    flux = flux0 / 0.5 / Esum / dE
    Elow = E1 - 0.5 * (E1 - pr1)
    Ehigh = E2 - 0.5 * (E2 - pr2)

    precout = "\n".join((f"{inp['precipstartsec']}", f"{Elow} {flux}", f"{Ehigh} -1.0", f"{inp['precipendsec']}", "-1.0 -1.0"))

    print("writing", ofn)

    ofn.write_text(precout)
