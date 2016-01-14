
"""
GALEXTRACT

extract galaxy spectra from image or images and combine them

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter


from pyraf import iraf
from iraf import pysalt

from specsky import skysubtract
from specextract import extract, write_extract
from specslitnormalize import specslitnormalize
from specsens import specsens
from speccal import speccal

from PySpectrograph.Spectra import findobj

from agnextract import salt_extract

if __name__=='__main__':
    import argparse 
    parser = argparse.ArgumentParser(description='Extract spectra from SALT 2D rectified image')
    parser.add_argument('objfile', help='SALT 2D rectified image')
    parser.add_argument('yc', help='Central row of object', type=int)
    parser.add_argument('dy', help='Half width of object', type=int)
    parser.add_argument('--spst', dest='cal_file', default='',
                   help='SPST calibration file')
    parser.add_argument('--f', dest='format', default='ascii', choices=['ascii','lcogt'],
                   help='Format for output file')
    args = parser.parse_args()

   
    salt_extract(args.objfile, args.yc, args.dy, specformat=args.format, calfile=args.cal_file, convert=True, normalize=False)
