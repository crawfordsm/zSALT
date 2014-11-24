
#Script for making BPM mask 

import os, sys
from astropy.io import fits

from pyraf import iraf
from iraf import pysalt

from saltprepare import saltprepare
from saltbias import saltbias

outfile = 'tmp.fits'

saltprepare(sys.argv[1], outfile, '', createvar=False, badpixelimage='', clobber=True, logfile='tmp.log', verbose=True)
saltbias(outfile, outfile, '', subover=True, trim=True, subbias=False, masterbias='', median=False, function='polynomial', order=5,rej_lo=3.0,rej_hi=5.0, niter=10, plotover=False, turbo=False, clobber=True, logfile='tmp.log', verbose=True)


hdu = fits.open(outfile)
limit = float(sys.argv[2])
for i in range(1,len(hdu)): 
    mask = (hdu[i].data < limit)
    hdu[i].data = hdu[i].data * 0.0
    hdu[i].data[mask] = 1


if os.path.isfile('bpm.fits'): os.remove('bpm.fits')
hdu.writeto('bpm.fits')

os.remove('tmp.fits')
os.remove('tmp.log')
