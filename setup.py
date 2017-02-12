#!/usr/bin/env python
from setuptools import setup

req=['h5py','scipy','pandas','pytz','numpy','nose','astropy','scikit-image','python-dateutil','matplotlib',
    'histutils','transcarread']

setup(name='transcar',
      packages=['transcar'],
      url='https://github.com/scienceopen/transcar',
      dependency_links = ['https://github.com/scienceopen/transcarread/tarball/master#egg=transcarread'],
     install_requires=req,
	  )
	  
