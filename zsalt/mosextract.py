
"""
MOSEXTRACT

extract galaxy spectra from MOS image

"""

import os, sys, glob, shutil

import numpy as np
from astropy.io import fits
from scipy.ndimage.filters import median_filter
from PySpectrograph.Spectra import findobj, Spectrum


import pylab as pl

from extract import extract, write_extract
from calfunc import calfunc
from salt_extract import write_ascii, clean_spectra

import spectools as st


def mos_extract(img, y1, y2, sy1, sy2, ext=None, normalize=True, calfile=None, specformat='ascii', convert=True):
    hdu = fits.open(img)
   
    data = hdu[ext].data

    ap = extract(hdu, method='normal', section=[(y1, y2)], ext=ext, convert=convert)
    sk = extract(hdu, method='normal', section=[(sy1, sy2)], ext=ext, convert=convert)
      
    flux = ap[0].ldata/(y2-y1) - sk[0].ldata/(sy2-sy1)
    flux[flux < 1] = np.median(flux)
    flux_spec = Spectrum.Spectrum(ap[0].wave, flux, abs(ap[0].lvar)**0.5, stype='continuum')
    outfile = 'CHILES_%s_%s.txt' % (hdu[ext].header['SLITNAME'].strip(), img.split('.')[0][-2:])
    write_ascii(outfile, flux_spec, clobber=True)

    pl.plot(ap[0].wave, np.convolve(ap[0].ldata/(y2-y1) - sk[0].ldata/(sy2-sy1), np.ones(5), mode='same'))
    pl.figure()
    pl.plot(ap[0].wave, np.convolve(ap[0].ldata, np.ones(5), mode='same'))
    pl.plot(ap[0].wave, np.convolve(sk[0].ldata, np.ones(5), mode='same'))
    pl.show()

    
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Extract spectra from SALT 2D rectified image')
    parser.add_argument('objfile', help='SALT 2D rectified image')
    parser.add_argument('y1', help='y1', type=int)
    parser.add_argument('y2', help='y1', type=int)
    parser.add_argument('sy1', help='y1', type=int)
    parser.add_argument('sy2', help='y1', type=int)
    parser.add_argument('--spst', dest='cal_file', default='',
                   help='SPST calibration file')
    parser.add_argument('--f', dest='format', default='ascii', choices=['ascii','lcogt'],
                   help='Format for output file')
    parser.add_argument('--e', dest='ext', default=1, type=int,
                   help='Extension to extract')
    args = parser.parse_args()


    mos_extract(args.objfile, args.y1, args.y2, args.sy1, args.sy2, specformat=args.format, ext=args.ext, calfile=args.cal_file, convert=True, normalize=False)

