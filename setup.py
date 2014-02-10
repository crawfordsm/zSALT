#!/usr/bin/env python

from distutils.core import setup
DISTUTILS_DEBUG=True

# Set affiliated package-specific settings
PACKAGENAME = 'zsalt'
DESCRIPTION = 'Measuring redshifts from SALT longslit observations'
LONG_DESCRIPTION = """
The zsalt package is set up to provide a nearly
automatic interface to measure redshifts from spectra 
observed with SALT using the Long Slit mode on RSS 
"""
AUTHOR = 'Steve Crawford'
AUTHOR_EMAIL = 'crawford@saao.ac.za'
LICENSE = 'BSD'
URL = 'http://pysalt.salt.ac.za'

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = '0.1.dev'


setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      url=URL,
      license=LICENSE,
      requires=['PySpectrograph', 'pyraf.iraf'],
      packages=['zsalt']
     )

