import sys, os, string, argparse
import pyfits
import numpy as np
from PySpectrograph import Spectrum
from PySpectrograph.Utilities.fit import interfit

from redshift import xcor_redshift, loadtext, loadiraf, loadsdss

import pylab as pl


if __name__=='__main__':

   parser = argparse.ArgumentParser(description='Match redshift to template')
   parser.add_argument('spectra', help='Spectra to measure redshift')
   parser.add_argument('-n', dest='noplot', default=True, action='store_false',
                    help='do not plot the data')
   parser.add_argument('--z1', dest='z1', default=0.0001, help='default lower redshift', type=float)
   parser.add_argument('--z2', dest='z2', default=1.2001, help='default lower redshift', type=float)

   args = parser.parse_args()

   infile = args.spectra

   if infile.count('fits'):
      hdu=pyfits.open(infile)
      spec=loadiraf(hdu)
   else:
      spec=loadtext(infile)
 
   dirpath = os.path.dirname(__file__)
 
   best_cc = 0
   z1 = args.z1
   z2 = args.z2
   for i in range(23, 33):
       template_name = 'spDR2-0{}.fit'.format(string.zfill(i, 2))
       thdu=pyfits.open(dirpath+'/template/'+template_name)
       template=loadsdss(thdu)
   

   #z_arr, cc_arr=xcor_redshift(spec, template, z1=0.0, z2=5.5, zstep=0.01)
       z_arr, cc_arr=xcor_redshift(spec, template, z1=z1, z2=z2, zstep=0.0001)
       z=z_arr[cc_arr.argmax()]
       print i, z, cc_arr.max(), cc_arr.argmax()
       if best_cc < cc_arr.max():
          best_cc = cc_arr.max()
          best_z = z
          best_i = i
          best_arr = cc_arr
   #z_arr, cc_arr=xcor_redshift(spec, template, z1=z-0.05, z2=z+0.05, zstep=0.0001)
   #z=z_arr[cc_arr.argmax()]
   #print z
   z = best_z
   i = best_i
   cc_arr = best_arr
   template_name = 'spDR2-0{}.fit'.format(string.zfill(i, 2))
   thdu=pyfits.open(dirpath+'/template/'+template_name)
   template=loadsdss(thdu)
   print best_i, best_z, best_arr.max()
   if args.noplot: 
      pl.figure()
      pl.plot(z_arr, cc_arr)
      pl.figure()
      cflux=np.convolve(spec.flux, np.ones(10), mode='same')
      pl.plot(spec.wavelength, cflux)
      nflux=np.interp(spec.wavelength, (1+z)*template.wavelength, template.flux)
      #pl.plot((1+z)*template.wavelength, template.flux*spec.flux.mean()/template.flux.mean())
      pl.plot(spec.wavelength, nflux*cflux.mean()/nflux.mean())
      pl.show()

   
   
    
