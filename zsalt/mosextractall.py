
import sys
import numpy as np
from astropy.io import fits

from pyraf import iraf
from iraf import pysalt 

import mostools as mt


def mos_extract_all(img):
    hdu = fits.open(img)
    
    for i in range(1, len(hdu)):
        if hdu[i].name == 'SCI':
           print(i)
           try:
              flux, sky = mt.extract_fit_flux(hdu[i].data, padding=1)
           except:
              continue
                

           xarr = np.arange(len(flux))
           # convert using the WCS information
           try:
               w0 = hdu[i].header['CRVAL1']
               dw = hdu[i].header['CD1_1']
           except Exception as e:
               msg = 'Error on Ext %i: %s' % (i, e)
               raise Exception(msg)
           warr = w0 + dw * xarr

           #clean up spectra
           err = abs(flux)**0.5
 
           #write out spectra
           outfile = '{}_{}_{}.txt'.format(hdu[0].header['OBJECT'], hdu[i].header['SLITNAME'].strip(), img.split('.')[0][-2:])
           fout=open(outfile, 'w')
           for j in range(len(flux)):
               fout.write('{} {} {}\n'.format(warr[j], flux[j], err[j]))
           fout.close()
               


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Extract spectra from SALT 2D rectified image')
    parser.add_argument('objfile', help='SALT 2D rectified image')
    args = parser.parse_args()
    mos_extract_all(args.objfile)
