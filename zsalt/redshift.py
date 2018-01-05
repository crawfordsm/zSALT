import sys
import pyfits
import numpy as np
from PySpectrograph import Spectrum
from PySpectrograph.Utilities.fit import interfit

import pylab as pl

def ncor(x, y):
    """Calculate the normalized correlation of two arrays"""
    d=np.correlate(x,x)*np.correlate(y,y)
    if d<=0: return 0
    return np.correlate(x,y)/d**0.5


def xcor_redshift(spectra, template, sub=False, z1=0, z2=1, zstep=0.001):
    """Meaure the redshift of a spectra by cross correlating it 
       with a template

       returns an array of correlation values
    """
    zvalue=np.arange(z1,z2,zstep)
    cc_arr=np.zeros(len(zvalue))
    sflux=continuum_subtract(spectra)
    tflux=continuum_subtract(template)
    for i,z in enumerate(zvalue):
        nflux=np.interp(spectra.wavelength, template.wavelength*(1+z), tflux)
        cc_arr[i]=ncor(sflux, nflux)
  
    return  zvalue, cc_arr

def continuum_subtract(spec, function='polynomial', order=7):
    """Fit a function to a spectra and subtract the continuum"""
    wc=interfit(spec.wavelength, spec.flux, function=function, order=order)
    wc.interfit()
 
    return spec.flux-wc(spec.wavelength)
    
def loadtext(infile):
   warr, farr=np.loadtxt(infile, usecols=(0,1), unpack=True)
   spec=Spectrum.Spectrum(warr, farr, stype='continuum')
   return spec
    

def loadiraf(hdu):
   farr=hdu[0].data
   xarr=np.arange(len(farr))
   warr=hdu[0].header['CRVAL1']+hdu[0].header['CDELT1']*(xarr+hdu[0].header['CRPIX1'])
   mask=(farr>10)
   spec=Spectrum.Spectrum(warr[mask], farr[mask], stype='continuum')
   return spec

def loadsdss(hdu):
   farr=hdu[0].data[0]
   xarr=np.arange(len(farr))
   warr=10**(hdu[0].header['CRVAL1']+hdu[0].header['CD1_1']*(xarr+1))
   spec=Spectrum.Spectrum(warr, farr, stype='continuum')
   return spec

if __name__=='__main__':
 
   if sys.argv[1].count('fits'):
      hdu=pyfits.open(sys.argv[1])
      spec=loadiraf(hdu)
   else:
      spec=loadtext(sys.argv[1])
 
 
   thdu=pyfits.open(sys.argv[2])
   template=loadsdss(thdu)


   z_arr, cc_arr=xcor_redshift(spec, template, z1=0.0001, z2=1.20, zstep=0.0001)
   z=z_arr[cc_arr.argmax()]
   #z_arr, cc_arr=xcor_redshift(spec, template, z1=z-0.05, z2=z+0.05, zstep=0.0001)
   #z=z_arr[cc_arr.argmax()]
   print z
   pl.figure()
   pl.plot(z_arr, cc_arr)
   pl.figure()
   cflux=np.convolve(spec.flux, np.ones(10), mode='same')
   pl.plot(spec.wavelength, cflux)
   nflux=np.interp(spec.wavelength, (1+z)*template.wavelength, template.flux)
   #pl.plot((1+z)*template.wavelength, template.flux*spec.flux.mean()/template.flux.mean())
   pl.plot(spec.wavelength, nflux*cflux.mean()/nflux.mean())
   pl.show()

   
   
    
