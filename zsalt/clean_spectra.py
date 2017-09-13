
import sys
import numpy as np
from astropy.io import fits

from pyraf import iraf
from iraf import pysalt 

import mostools as mt


def clean_spectra(img, region1=None, region2=None):
    w, f, e = np.loadtxt(img, usecols=(0,1,2), unpack=True)

    mask = np.ones(len(w))
    if region1 is not None:
       x1, x2 = region1
       mask[x1:x2]=0
    if region2 is not None:
       x1, x2 = region2
       mask[x1:x2]=0

    mask[0:10] = 0 
    mask[-10:] = 0
    mask =  (mask==1)
    warr = w[mask]
    flux = f[mask]
    err = e[mask]

    #clean bad spectra lines
    mask = (warr > 6295) * (warr < 6305) + (warr > 5570) * (warr < 5585) + (warr>6360)*(warr<6370)
    mask = (mask==0)
    warr = warr[mask]
    flux = flux[mask]
    err = err[mask]
    
    return warr, flux, err
    

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Clean an extracted spectra from SALT 2D rectified image')
    parser.add_argument('objfile', help='SALT 2D rectified image')
    parser.add_argument('--r1',  dest='r1', help='Region1: should be given as x1,x2', default=None)
    parser.add_argument('--r2',  dest='r2', help='Region2: should be given as x1,x2', default=None)
    parser.add_argument('--outfile',  dest='outfile', help='Output file to write to', default=None)
    args = parser.parse_args()
    if args.r1 is not None: region1=[int(x) for x in args.r1.split(',')]
    if args.r2 is not None: region2=[int(x) for x in args.r2.split(',')]
    warr, flux, err = clean_spectra(args.objfile, region1=region1, region2=region2)

    if args.outfile is None:
       outfile = args.objfile
    else:
       outfile = args.outfile
    fout=open(outfile, 'w')
    for j in range(len(flux)):
        fout.write('{} {} {}\n'.format(warr[j], flux[j], err[j]))
    fout.close()
