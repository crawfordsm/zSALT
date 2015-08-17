
"""
AGNEXTRACT

extract galaxy spectra from image or images and combine them

"""

import os, sys, glob, shutil

import numpy as np
from astropy.io import fits
from scipy.ndimage.filters import median_filter


from pyraf import iraf
from iraf import pysalt

from specsky import skysubtract
from specextract import extract, write_extract
from specslitnormalize import specslitnormalize
from specsens import specsens
from speccal import speccal, calfunc
import spectools as st

from PySpectrograph.Spectra import findobj, Spectrum

from agnextract import write_lcogt

def agncalibrate(img, outfile, calfile, specformat='lcogt'):

    
    #set up some files that will be needed
    logfile='specext.log'

    hdu = fits.open(img)
    w1 = hdu[0].header['CRVAL1']
    p1 = hdu[0].header['CRPIX1']
    dw = hdu[0].header['CD1_1']
    f = hdu[0].data[0][0]
    e = hdu[0].data[3][0]
    xarr = np.arange(len(f))
    w = (xarr)*dw+w1


    cal_spectra=st.readspectrum(calfile, error=False, ftype='ascii')
    airmass=hdu[0].header['AIRMASS']
    exptime=hdu[0].header['EXPTIME']
    extfile=iraf.osfn("pysalt$data/site/suth_extinct.dat")
    ext_spectra=st.readspectrum(extfile, error=False, ftype='ascii')

    flux_spec=Spectrum.Spectrum(w, f, e, stype='continuum')
    flux_spec=calfunc(flux_spec, cal_spectra, ext_spectra, airmass, exptime, True)
    hdu[0].data[0][0] = flux_spec.flux
    hdu[0].data[3][0] = flux_spec.var
    hdu.writeto(outfile, clobber=True)


if __name__=='__main__':
   specformat = 'lcogt'
   agncalibrate(sys.argv[1], sys.argv[2], sys.argv[3], specformat=specformat)
