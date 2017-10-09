#!/usr/bin/env python

req=['nose','python-dateutil','pandas']

import pip
try:
    import conda.cli
    conda.cli.main('install',*req)
except Exception as e:
    pip.main(['install'] + req)

# %%
from setuptools import setup
import subprocess

setup(name='transcar',
      packages=['transcar'],
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scivision/transcar',
      dependency_links = [],
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3',
      ],
     install_requires=req,
	  )

subprocess.check_call(['cmake','..'],cwd='dir.source/dir.obj')
subprocess.check_call(['make'],cwd='dir.source/dir.obj')
