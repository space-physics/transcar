#!/usr/bin/env python
import os
from setuptools import setup
import subprocess

OBJ = 'dir.source/dir.obj'


setup()

if os.name == 'nt':
    subprocess.check_call(['cmake', '-G', 'MinGW Makefiles',
                           '-DCMAKE_SH="CMAKE_SH-NOTFOUND', '..'], cwd=OBJ)
else:
    subprocess.check_call(['cmake', '..'], cwd=OBJ)

subprocess.check_call(['cmake', '--build', '.'], cwd=OBJ)
