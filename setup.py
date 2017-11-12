#!/usr/bin/env python

req=['nose','python-dateutil','pandas']
pipreq=['transcarread']

try:
    import conda.cli
    conda.cli.main('install',*req)
except Exception:
    pass
# %%
from setuptools import setup
import subprocess

setup(name='transcar',
      packages=['transcar'],
      version='0.1.0',
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scivision/transcar',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 3 - Alpha',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3',
      ],
      python_requires='>=3.6',
      install_requires=req+pipreq,
	  )

subprocess.check_call(['cmake','..'],cwd='dir.source/dir.obj')
subprocess.check_call(['make'],cwd='dir.source/dir.obj')
