#!/usr/bin/env python
from setuptools import setup
try:
    import conda.cli
    conda.cli.main('install','--file','requirements.txt')
except Exception as e:
    print(e)


setup(name='transcar',
      dependency_links = ['https://github.com/scienceopen/histutils/tarball/master#egg=histutils',
                          'https://github.com/scienceopen/transcarread/tarball/master#egg=transcarread'],
     install_requires=['histutils','transcarread'],
      packages=['transcar'],
	  )
	  
