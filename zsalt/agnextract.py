
"""
AGNEXTRACT

extract galaxy spectra from image or images and combine them

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter


from extract import extract, write_extract
from calfunc import calfunc

import spectools as st

from PySpectrograph.Spectra import findobj, Spectrum

def agnextract(img, yc=None, dy=None, normalize=True, calfile=None, convert=True, specformat='ascii'):

    
    #set up some files that will be needed
    logfile='specext.log'


    #create the spectra text files for all of our objects
    spec_list=[]
    #skynormalize the data
    if normalize:
       from specslitnormalize import specslitnormalize
       specslitnormalize(img, 'n'+img, '', response=None, response_output=None, order=3, conv=1e-2, niter=20,
                     startext=0, clobber=True,logfile='salt.log',verbose=True)
       img = 'n'+img
    hdu=pyfits.open(img)
    target=hdu[0].header['OBJECT']
    ofile='%s.%s_%i_%i.ltxt' % (target, extract_date(img), extract_number(img), yc)
    if specformat=='lcogt': ofile=ofile.replace('ltxt', 'fits')
    #ofile = img.replace('fits', 'txt')

    extract_spectra(hdu, yc, dy, ofile, smooth=False, grow=10, clobber=True, specformat=specformat, convert=convert, calfile=calfile)

 

def speccombine(spec_list, sfile):
   """Combine N spectra"""

   w1,f1, e1=np.loadtxt(spec_list[0], usecols=(0,1,2), unpack=True)

   w=w1
   f=1.0*f1
   e=e1**2

   for sfile in spec_list[1:]:
      w2,f2, e2=np.loadtxt(sfile, usecols=(0,1,2), unpack=True)
      if2=np.interp(w1, w2, f2)
      ie2=np.interp(w1, w2, e2)
      f2=f2*np.median(f1/if2)
      f+=if2
      e+=ie2**2

   f=f/len(spec_list)
   e=e**0.5/len(spec_list)

   fout=open(sfile, 'w')
   for i in range(len(w)):
           fout.write('%f %e %e\n' % (w[i], f[i], e[i]))
   fout.close()


def extract_spectra(hdu, yc, dy, outfile, minsize=5, thresh=3, grow=0, smooth=False, maskzeros=False, 
                    convert=True,  cleanspectra=True, calfile=None, clobber=True, specformat='ascii'):
    """From an image, extract a spectra.   

    """
    data=hdu[1].data

    #replace the zeros with the average from the frame
    if maskzeros:
       mean,std=iterstat(data[data>0])
       #rdata=mean  np.random.normal(mean, std, size=data.shape)
       data[data<=0]=mean #rdata[data<=0]

    y1=yc-dy
    y2=yc+dy
    ap_list=extract(hdu, method='normal', section=[(y1,y2)], minsize=minsize, thresh=thresh, convert=convert)
    sy1a=y2
    sy2a=sy1a+2.0*dy
    ska_list=extract(hdu, method='normal', section=[(sy1a,sy2a)], minsize=minsize, thresh=thresh, convert=convert)
    sy2b=y1-dy
    sy1b=sy2b-2.0*dy
    skb_list=extract(hdu, method='normal', section=[(sy1b,sy2b)], minsize=minsize, thresh=thresh, convert=convert)
    print sy1b, sy2b

    sdata = 0.5*(ska_list[0].ldata/(sy2a-sy1a) + skb_list[0].ldata/(sy2b-sy1b))
    #sdata = ska_list[0].ldata/(sy2a-sy1a)
    #sdata = skb_list[0].ldata/(sy2b-sy1b)
    raw = 1.0 * ap_list[0].ldata
    print 'extract:', ap_list[0].ldata[1124]
    ap_list[0].ldata=ap_list[0].ldata-float(y2-y1) * sdata
    print 'sky:', ap_list[0].ldata[1124]
 
    print ap_list[0].wave[10], ap_list[0].ldata[10], ap_list[0].lvar[10]
    flux_spec=Spectrum.Spectrum(ap_list[0].wave, ap_list[0].ldata, abs(ap_list[0].lvar)**0.5, stype='continuum')
    print flux_spec.wavelength[10], flux_spec.flux[10], flux_spec.var[10]

    if cleanspectra:
       clean_spectra(ap_list[0], grow=grow)
    print 'clean:', ap_list[0].ldata[1124]

    if calfile:
           cal_spectra=st.readspectrum(calfile, error=False, ftype='ascii')
           airmass=hdu[0].header['AIRMASS']
           exptime=hdu[0].header['EXPTIME']
           extfile=os.path.dirname(st.__file__)+"/suth_extinct.dat"
           print extfile
           ext_spectra=st.readspectrum(extfile, error=False, ftype='ascii')

           flux_spec=Spectrum.Spectrum(ap_list[0].wave, ap_list[0].ldata, abs(ap_list[0].lvar)**0.5, stype='continuum')
           print flux_spec.flux[10], flux_spec.var[10]
           flux_spec=calfunc(flux_spec, cal_spectra, ext_spectra, airmass, exptime, True)
           print flux_spec.flux[10], flux_spec.var[10]
    else:
        flux_spec = Spectrum.Spectrum(ap_list[0].wave, ap_list[0].ldata, abs(ap_list[0].lvar)**0.5, stype='continuum')
    
    if specformat == 'ascii':
        write_ascii(outfile, flux_spec, clobber=clobber)
    elif specformat == 'lcogt':
        write_lcogt(outfile, flux_spec, hdu, sky=float(y2-y1) * sdata, raw = raw, clobber=clobber)

sky_lines = [] #[4358.34, 5577.338, 6300.304,6363.8000,7715.0116]

def write_ascii(outfile, ap_spec, clobber=True):
    if os.path.isfile(outfile): os.remove(outfile)
    fout = open(outfile, 'w')
    for i in range(len(ap_spec.wavelength)):
        fout.write('%7.3f %e %e\n' % (ap_spec.wavelength[i], ap_spec.flux[i], ap_spec.var[i]))
    fout.close()

def write_lcogt(outfile, ap, hdu, sky, raw, clobber=True):
    """Write data in the format of an LCOGT output file

    """

    header = hdu[0].header
    data = np.array([[ap.flux], [raw], [sky], [ap.var]])
    header['CTYPE1'] = 'LINEAR'
    header['CTYPE2'] = 'LINEAR'
    header['CRVAL1'] = hdu[1].header['CRVAL1']
    header['CRPIX1'] = 1
    header['CD1_1'] = hdu[1].header['CD1_1']
    header['CD2_2'] = 1.0
    header['WAT0_001']= 'system=equispec'
    header['WAT1_001']= 'wtype=linear label=Wavelength units=angstroms'
    header['WAT2_001']= 'wtype=linear'
    header['APNUM1'] = '1 1 1 1'
    header['NAXIS']   = 3 
    header['NAXIS1']  = len(ap.flux)
    header['NAXIS2']  = 1
    header['NAXIS3']  = 4 
    header['UTSTART'] = hdu[0].header['TIME-OBS']

   
 
    hdulist = pyfits.PrimaryHDU(data, header = header)
    hdulist.writeto(outfile, clobber=clobber)

    return 

def clean_spectra(ap_spec, dres=2, grow=0, mask=False, replace=True):
    """Clean a spectrum and make it more presenatable"""
    xmax=len(ap_spec.ldata)
    for w in sky_lines:
        mask=(abs(ap_spec.wave-w)<dres)
        ap_spec.ldata[mask] = 0

    ap_spec.ldata[0] = 0
    ap_spec.ldata[-1] = 0

    ap_spec.ldata[ap_spec.ldata<0] = 0 


    #grow the spectra
    if grow:
      nid = np.where(ap_spec.ldata==0)[0]
      for i in nid:
          x1=max(0,int(i-grow))
          x2=min(int(i+grow),xmax-1)
          ap_spec.ldata[x1:x2]=0

    #replace 
    if replace:
       nid = np.where(ap_spec.ldata==0)[0]
       mask = (ap_spec.ldata>0)
       for i in nid:
           ap_spec.ldata[i] = np.interp(ap_spec.wave[i], ap_spec.wave[mask], ap_spec.ldata[mask])
 

    return ap_spec

def extract_number(img):
    """Get the image number only"""
    img=img.split('.fits')
    nimg=int(img[0][-4:])
    return nimg

def extract_date(img):
    """Get the date"""
    img=img.split('.fits')
    obsdate=int(img[0][-8:-4])
    print obsdate
    return obsdate

def iterstat(data, thresh=3, niter=5):
    mean=data.mean()
    std=data.std()
    for i in range(niter):
        mask=(abs(data-mean)<std*thresh)
        mean=data[mask].mean()
        std=data[mask].std()
    return mean, std

    

def findskysection(section, skysection=[800,900], skylimit=100):
    """Look through the section lists and determine a section to measure the sky in

       It should be as close as possible to the center and about 200 pixels wide
    """
    #check to make sure it doesn't overlap any existing spectra
    #and adjust if necessary
    for y1, y2 in section:
        if -30< (skysection[1]-y1)<0:
           skysection[1]=y1-30
        if 0< (skysection[0]-y2)<30:
           skysection[0]=y2+30
    if skysection[1]-skysection[0] < skylimit: print "WARNING SMALL SKY SECTION"
    return skysection
    


if __name__=='__main__':
   if len(sys.argv)>=5:
      calfile=sys.argv[4]
   else: 
      calfile=None
   print calfile
   specformat = 'lcogt'
   agnextract(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), specformat=specformat, convert=True, calfile=calfile, normalize=False)
