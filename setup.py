#!/usr/bin/env python
from setuptools import setup
import subprocess

try:
    subprocess.call(['conda','install','--quiet','--file','requirements.txt'])
except Exception:
    pass


setup(name='transcar',
	  description='Transcar 1-D flux tube model',
	  author='Michael Hirsch',
	  url='https://github.com/scienceopen/transcar',
      dependency_links = ['https://github.com/scienceopen/histutils/tarball/master#egg=histutils',
                          'https://github.com/scienceopen/transcarread/tarball/master#egg=transcarread'],
     install_requires=['histutils','transcarread',
     		       'pathlib2'],
      packages=['transcar'],
	  )
	  
