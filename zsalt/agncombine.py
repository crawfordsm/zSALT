#!/usr/local/bin/python
import os
import sys
import numpy as np
from astropy.io import fits
from pylab import *

print os.path.dirname(__file__)

figure()
n=0
hdu_list = []
flux_list = []
raw_list = []
sky_list = []
var_list = []
for infile in sys.argv[1:]:
   hdu = fits.open(infile)
   hdu_list.append(hdu)
   flux_list.append(hdu[0].data[0][0])
   raw_list.append(hdu[0].data[1][0])
   sky_list.append(hdu[0].data[2][0])
   var_list.append(hdu[0].data[3][0])


   
   w1 = hdu[0].header['CRVAL1']
   p1 = hdu[0].header['CRPIX1']
   dw = hdu[0].header['CD1_1']
   f = hdu[0].data[0][0]
   r = hdu[0].data[1][0]
   s = hdu[0].data[2][0]
   e = hdu[0].data[3][0]
   xarr = np.arange(len(f))
   w = (xarr)*dw+w1

   f=np.convolve(f, np.ones(5), mode='same')/5
   #if n == 0: n = f.mean()
   #plot(w,f * n/f.mean())
   plot(w,e)
   #plot(w,e)

var_list = np.array(var_list)
var_list[var_list==0] = var_list[var_list>0].min()
flux = np.average(np.array(flux_list), axis=0, weights=1.0/var_list)
raw  = np.average(np.array(raw_list), axis=0, weights=1.0/var_list)
sky  = np.average(np.array(sky_list), axis=0, weights=1.0/var_list)
var = np.std(np.array(flux_list), axis=0)

hdu = hdu_list[0]
hdu[0].data[0][0] = flux
hdu[0].data[1][0] = raw 
hdu[0].data[2][0] = sky 
hdu[0].data[3][0] = var

name_list = sys.argv[1].split('.')
obsdate =  hdu[0].header['DATE-OBS'].replace('-', '')
name = name_list[0]+'.' + obsdate + '.fits'

hdu.writeto(name, clobber=True)

print name #len(w), len(flux)
figure()
plot(w, var)

show()
