#!/usr/bin/env python

req=['nose']
pipreq=['transcarread']

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
      dependency_links = ['https://github.com/scivision/transcarread/tarball/master#egg=transcarread',],
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3.6',
      ],
     install_requires=req + pipreq,
	  )

subprocess.run(['cmake','..'],cwd='dir.source/dir.obj')
subprocess.run(['make','-j2'],cwd='dir.source/dir.obj')