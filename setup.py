#!/usr/bin/env python
install_requires=['python-dateutil','pytz','pandas']
tests_require=['pytest','nose','coveralls','transcarread']
OBJ = 'dir.source/dir.obj'
# %%
import os
from setuptools import setup,find_packages
import subprocess

setup(name='transcar',
      packages=find_packages(),
      version='0.2.0',
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scivision/transcar',
      long_description=open('README.rst').read(),
      description='Parallel-executing Transcar 1-D particle precipitation ionospheric model.',
      classifiers=[
      'Development Status :: 4 - Beta',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 3',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      ],
      python_requires='>=3.6',
      install_requires=install_requires,
      tests_require=tests_require,
      extras_require={'tests':tests_require,},
      scripts=['RunTranscar.py'],
	  )

if os.name =='nt':
    subprocess.check_call(['cmake','-G','MinGW Makefiles','..'],cwd=OBJ)
    subprocess.check_call(['mingw32-make'],cwd=OBJ)
else:
    subprocess.check_call(['cmake','..'],cwd=OBJ)
    subprocess.check_call(['make'],cwd=OBJ)
