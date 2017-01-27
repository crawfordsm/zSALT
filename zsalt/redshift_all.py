import sys, os, string
import pyfits
import numpy as np
from PySpectrograph import Spectrum
from PySpectrograph.Utilities.fit import interfit

from redshift import xcor_redshift, loadtext, loadiraf, loadsdss

import pylab as pl


if __name__=='__main__':
 
   if sys.argv[1].count('fits'):
      hdu=pyfits.open(sys.argv[1])
      spec=loadiraf(hdu)
   else:
      spec=loadtext(sys.argv[1])
 
   dirpath = os.path.dirname(__file__)
 
   best_cc = 0
   z1 = 0.0
   z2 = 1.2
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
   print best_z
   pl.figure()
   pl.plot(z_arr, cc_arr)
   pl.figure()
   cflux=np.convolve(spec.flux, np.ones(10), mode='same')
   pl.plot(spec.wavelength, cflux)
   nflux=np.interp(spec.wavelength, (1+z)*template.wavelength, template.flux)
   #pl.plot((1+z)*template.wavelength, template.flux*spec.flux.mean()/template.flux.mean())
   pl.plot(spec.wavelength, nflux*cflux.mean()/nflux.mean())
   pl.show()

   
   
    
