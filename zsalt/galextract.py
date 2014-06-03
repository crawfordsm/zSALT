
"""
GALEXTRACT

extract galaxy spectra from image or images and combine them

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter


from pyraf import iraf
from iraf import pysalt

from specsky import skysubtract
from specextract import extract, write_extract
from specslitnormalize import specslitnormalize
from specsens import specsens
from speccal import speccal
from quickspec import clean_spectra

from PySpectrograph.Spectra import findobj

def galextract(img, yc=None, dy=None, normalize=True, calfile=None):

    
    #set up some files that will be needed
    logfile='specext.log'


    #create the spectra text files for all of our objects
    spec_list=[]
    #skynormalize the data
    if normalize:
       specslitnormalize(img, 'n'+img, '', response=None, response_output=None, order=3, conv=1e-2, niter=20,
                     startext=0, clobber=True,logfile='salt.log',verbose=True)

    hdu=pyfits.open('n'+img)
    target=hdu[0].header['OBJECT']
    ofile='%s.%s_%i_%i.ltxt' % (target, extract_date(img), extract_number(img), yc)

    extract_spectra(hdu, yc, dy, ofile, smooth=False, grow=10, clobber=True)

    if calfile is not None: 
           airmass=hdu[0].header['AIRMASS']
           exptime=hdu[0].header['EXPTIME']
           extfile=iraf.osfn("pysalt$data/site/suth_extinct.dat")
           speccal(ofile, ofile.replace("txt", "spec"), calfile, extfile, airmass, exptime, clobber=True, logfile='salt.log', verbose=True)

 

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


def extract_spectra(hdu, yc, dy, outfile, minsize=5, thresh=3, grow=0, smooth=False, maskzeros=False, cleanspectra=True, clobber=True):
    """From an image, extract a spectra.   

    """
    if smooth:
       data=smooth_data(hdu[1].data)
    else:
       data=hdu[1].data

    #replace the zeros with the average from the frame
    if maskzeros:
       mean,std=iterstat(data[data>0])
       #rdata=mean  np.random.normal(mean, std, size=data.shape)
       data[data<=0]=mean #rdata[data<=0]

    y1=yc-dy
    y2=yc+dy

    sy1=y2+5.0*dy
    sy2=sy1+10.0*dy

    ap_list=extract(hdu, method='normal', section=[(y1,y2)], minsize=minsize, thresh=thresh, convert=True)
    sk_list=extract(hdu, method='normal', section=[(sy1,sy2)], minsize=minsize, thresh=thresh, convert=True)
    
    ap_list[0].ldata=ap_list[0].ldata-float(y2-y1)/(sy2-sy1)*sk_list[0].ldata

    if cleanspectra:
       clean_spectra(ap_list[0], grow=grow)
    
    write_extract(outfile, [ap_list[0]], outformat='ascii', clobber=clobber)

def smooth_data(data, mbox=25):
    mdata=median_filter(data, size=(mbox, mbox))
    return data-mdata

def find_section(section, y):
    """Find the section closest to y"""
    best_i=-1
    dist=1e5
    for i in range(len(section)):
        d=min(abs(section[i][0]-y), abs(section[i][1]-y))
        if d < dist:
           best_i=i
           dist=d
    return best_i

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
   galextract(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
