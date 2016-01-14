
"""
AGNEXTRACT

extract galaxy spectra from image or images and combine them

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter


from extract import extract, write_extract
from calfunc import calfunc

import spectools as st

from PySpectrograph.Spectra import findobj, Spectrum

from salt_extract import salt_extract

if __name__=='__main__':
   if len(sys.argv)>=5:
      calfile=sys.argv[4]
   else: 
      calfile=None
   print calfile
   specformat = 'lcogt'
   salt_extract(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), specformat=specformat, convert=True, calfile=calfile, normalize=False)
