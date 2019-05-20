#!/usr/bin/env python
import os
import subprocess

OBJ = 'dir.source/dir.obj'

if os.name == 'nt':
    subprocess.check_call(['cmake', '-G', 'MinGW Makefiles',
                           '-DCMAKE_SH="CMAKE_SH-NOTFOUND', '..'], cwd=OBJ)
else:
    subprocess.check_call(['cmake', '..'], cwd=OBJ)

subprocess.check_call(['cmake', '--build', str(OBJ), '-j'])
