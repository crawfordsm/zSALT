
import sys
import numpy as np
from astropy.io import fits

from pyraf import iraf
from iraf import pysalt 

import mostools as mt


def source_extract(img, yc, dy, ext=1, trim=0):
    hdu = fits.open(img)
    
    data = hdu[ext].data[yc-dy:yc+dy,:]
    if 'VAREXT' in hdu[ext].header:
       varext = hdu[ext].header['VAREXT']
       err = hdu[varext].data[yc-dy:yc+dy,:]
    if 'BPMEXT' in hdu[ext].header:
       bpmext = hdu[ext].header['BPMEXT']
       mask = hdu[bpmext].data[yc-dy:yc+dy,:]

    flux, sky = mt.extract_fit_flux(data, mask=mask, err=err, padding=0, bounds=False)

    xarr = np.arange(len(flux))
    # convert using the WCS information
    try:
        w0 = hdu[ext].header['CRVAL1']
        dw = hdu[ext].header['CD1_1']
    except Exception as e:
        msg = 'Error on Ext {}: {:s}'.format(ext, str(e))
        raise Exception(msg)
    warr = w0 + dw * xarr

    #clean up spectra
    err = abs(flux)**0.5
 
    #write out spectra
    obsdate = img.split('P')[-1][0:8]
    outfile = '{}_{}_{}.txt'.format(hdu[0].header['OBJECT'], obsdate, yc)
    fout=open(outfile, 'w')
    for j in range(0+trim, len(flux)-trim):
        if not np.isnan(flux[j]):
           fout.write('{} {} {}\n'.format(warr[j], flux[j], err[j]))
    fout.close()
               
    return outfile


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Extract spectra from SALT 2D rectified image')
    parser.add_argument('objfile', help='SALT 2D rectified image')
    parser.add_argument('yc', help='Y-position of source')
    parser.add_argument('dy', help='Y-window')
    parser.add_argument('--trim', dest='trim', default=0, type=int, help='Amount to trim spectra')
    args = parser.parse_args()
    source_extract(args.objfile, int(args.yc), int(args.dy), trim=args.trim)
