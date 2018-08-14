#!/usr/bin/env python
"""
Checks Transcar output "monoenergetic" electron beams.
"""
from pathlib import Path
#
from transcar.io import transcaroutcheck
from argparse import ArgumentParser


def main():

    p = ArgumentParser(description='Check Transcar output')
    p.add_argument('rodir', help='root of beam directory')
    p.add_argument('--errfn', help='file to write transcar Errors to', default='transcarError.log')
    p = p.parse_args()
    # %%
    rodir = Path(p.rodir).expanduser()

    dlist = [d for d in rodir.iterdir() if d.is_dir() and d.stem.startswith('beam')]

    for d in dlist:
        transcaroutcheck(d, p.errfn, 'STOP fin normale')


if __name__ == '__main__':
    main()
