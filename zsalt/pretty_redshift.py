import os, sys
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
   print len(farr)
   xarr=np.arange(len(farr))
   warr=10**(hdu[0].header['CRVAL1']+hdu[0].header['CD1_1']*(xarr+1))
   spec=Spectrum.Spectrum(warr, farr, stype='continuum')
   return spec

def readlinelist(infile):
    line_wave=[]
    line_name=[]
    for lines in open(infile).readlines():
        l=lines.split()
        line_wave.append(l[0])
        line_name.append(l[1])
    return line_wave, line_name

if __name__=='__main__':
 
   if sys.argv[1].count('fits'):
      hdu=pyfits.open(sys.argv[1])
      spec=loadiraf(hdu)
   else:
      spec=loadtext(sys.argv[1])
 

   thdu=pyfits.open(sys.argv[2])
   zc=float(sys.argv[3])
   template=loadsdss(thdu)
   z1=max(0,zc-0.15)
   z2=max(0,zc+0.15)

   z_arr, cc_arr=xcor_redshift(spec, template, z1=z1, z2=z2, zstep=0.0001)
   z=z_arr[cc_arr.argmax()]
   z = zc
   print z
   #z_arr, cc_arr=xcor_redshift(spec, template, z1=z-0.05, z2=z+0.05, zstep=0.0001)
   #z=z_arr[cc_arr.argmax()]
   #print z
   pl.figure()
   sp=pl.axes([0.15,0.15,0.8,0.8])
   cflux=np.convolve(spec.flux, np.ones(10), mode='same')
   #cflux*=1000e16
   sp.plot(spec.wavelength, cflux, color='#000000')
   coef=np.polyfit(spec.wavelength, cflux, 3)
   #sp.plot(spec.wavelength, np.polyval(coef, spec.wavelength))

   nflux=np.interp(spec.wavelength, (1+z)*template.wavelength, template.flux)

   tcoef=np.polyfit(spec.wavelength, nflux*cflux.mean()/nflux.mean(), 2)
   #ratio=cflux.mean()/nflux.mean()#*np.polyval(coef, spec.wavelength)/np.polyval(tcoef, spec.wavelength)
   ratio=cflux.mean()/nflux.mean()*np.polyval(coef, spec.wavelength)/np.polyval(tcoef, spec.wavelength)
   #sp.plot(spec.wavelength, nflux*ratio*0.5-0.4e-16, color='#FF0000')
   sp.plot(spec.wavelength, nflux*ratio, color='#FF0000')
   #sp.plot(spec.wavelength, np.polyval(tcoef, spec.wavelength))
   #pl.plot((1+z)*template.wavelength, template.flux*spec.flux.mean()/template.flux.mean())
   spname=sys.argv[1].split('_')[0]
   #sp.set_ylim([0,2000])
   x1,x2=sp.get_xlim()
   y1,y2=sp.get_ylim()
   print y1,y2, x1,x2
   line_wave, line_name=readlinelist(os.path.dirname(__file__)+'/sdss.linelist')
   dx=10
   for w,n in zip(line_wave, line_name):
     w=float(w)*(1+z)
     if w>x1 and w< x2:
       sp.plot([w,w],[y1,y2],ls='--', color='#AAAAAA')
       sp.text(w, y2-dx, '$%s$' % n.replace('_', '\\'), color='#AAAAAA', fontsize=8)
       #if dx<300:  
       #   dx+=100
       #else:
       #   dx=100
   spname=sys.argv[4]
   sp.text(x1+0.1*(x2-x1),0.8*y2,spname, fontsize=24)
   sp.text(x1+0.1*(x2-x1),0.70*y2,'z=%5.4f' % zc, fontsize=24)
   sp.set_ylabel('Counts')
   sp.set_xlabel('$\lambda \ (\AA)$')
   if len(sys.argv)>5:
      sy1=float(sys.argv[5])
      sy2=float(sys.argv[6])
   else:
      sy1=0.7
      sy2=0.7

   if False:
     cc=pl.axes([sy1, sy2,0.2,0.2])
     cc.plot(z_arr, cc_arr, color='#777777')
     xticks=np.arange(100*z1,100*z2+1,10, dtype=int)/100.0
     print xticks
     cc.set_xticks(xticks)
     cc.set_yticklabels([])
     cc.set_xlabel('z')
     cc.set_title('X-corr Function')

   pl.savefig(spname+'.png')
   pl.show()
   
   
    
