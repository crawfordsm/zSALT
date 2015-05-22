
"""
MOSEXTRACT

extract galaxy spectra from MOS image

"""

import os, sys, glob, shutil

import numpy as np
from astropy.io import fits
from scipy.ndimage.filters import median_filter

from astropy import modeling as md

import pylab as pl

#from pyraf import iraf
#from iraf import pysalt

#from specsky import skysubtract
#from specextract import extract, write_extract
#from specslitnormalize import specslitnormalize
#from specsens import specsens
#from speccal import speccal
#from quickspec import clean_spectra

from PySpectrograph.Spectra import findobj

def mosextract(img, yc=None, dy=None, normalize=True, calfile=None):

    
    #set up some files that will be needed
    logfile='specext.log'


    #create the spectra text files for all of our objects
    spec_list=[]
    #skynormalize the data
    #if normalize:
    #    specslitnormalize(img, 'n'+img, '', response=None, response_output=None, order=3, conv=1e-2, niter=20,
    #                 startext=0, clobber=True,logfile='salt.log',verbose=True)

    hdu=fits.open(img)
    ofile = img.replace('fits', 'txt')

    #pl.imshow(hdu[1].data[:,1340:1640], aspect='auto', origin='lower')
    #pl.show()

    extract_spectra(hdu)

    #extract_spectra(hdu, yc, dy, ofile, smooth=False, grow=10, clobber=True)

@md.models.custom_model_1d
def model_of_skys(x, amplitude1=1., mean1=-1., sigma1=1.,
                        amplitude2=1., mean2=1., sigma2=1.):
    return (amplitude1 * np.exp(-0.5 * ((x - mean1) / sigma1)**2) +
            amplitude2 * np.exp(-0.5 * ((x - mean2) / sigma2)**2))



def extract_spectra(hdu, dy=7):
    """Extract a spectrum for a source 

    """

    #first define the good region around the spectra

    max_list=[]
    x_arr = np.arange(len(hdu[1].data[0]))
    max_arr = x_arr*0.0 - 1
    for xc in x_arr: 
       f = 1.0 * hdu[1].data[dy:-dy,xc]
       m = (hdu[3].data[dy:-dy,xc]==0)
       if m.sum():
          max_arr[xc] = f.argmax() 

    m_init = md.models.Legendre1D(5) 
    fit = md.fitting.LevMarLSQFitter()
    fmax = fit(m_init, x_arr[max_arr>0], max_arr[max_arr>0])
    #pl.plot(x_arr, max_arr)
    #pl.plot(m(x_arr))
    #pl.show()
       

    #xc = 1402
    #xc = 1802
    for xc in range(len(hdu[1].data[0])):
        f = 1.0 * hdu[1].data[dy:-dy,xc]
        e = hdu[varext].data[dy:-dy,xc]
        m = (hdu[3].data[dy:-dy,xc]>0)
        y = np.arange(len(f))
        i = fmax(xc)
        mask = (f>0) * (m==0)
        mask[i-dy:i+dy] = 0
        if mask.sum():
           f[i-dy:i+dy] = np.interp(y[i-dy:i+dy], y[mask], f[mask]) 
        
        try:
           hdu[1].data[dy:-dy,xc] -=f #- m(y)
           hdu[1].data[:dy,xc] = 0
           hdu[1].data[-dy:,xc] = 0
        except:
           pass
    hdu.writeto('tmp.fits', clobber=True)
    f = hdu[1].data.sum(axis=0)
    m = hdu[3].data.sum(axis=0)
    warr = x_arr * hdu[1].header['CD1_1']  + hdu[1].header['CRVAL1'] 
    fout = open('tmp.txt', 'w')
    for i in range(len(f)):
        if f[i]>0: 
           fout.write('%6.2f %e\n' % (warr[i], f[i]))
    exit()
 
    #pl.plot(m(y))

    niter=5
    #for i in range(niter):
    #for i in range(niter):
    #chi = abs(f - m(y))
    #    mask = mask * (chi < 3*chi[mask].mean())
    #    m = fit(m_init, y[mask], f[mask])

        
    print m

    pl.plot(f[mask])
    pl.plot(m(y[mask]))
    #pl.plot(f-m(y))
    #pl.plot(m(y[mask]))
    print chi[mask].mean()
    pl.figure()
    pl.plot(chi[mask])
    pl.show()
    


def old(hdu, yc, dy, outfile, minsize=5, thresh=3, grow=0, smooth=False, maskzeros=False, cleanspectra=True, clobber=True):
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
   mosextract(sys.argv[1]) #, int(sys.argv[2]), int(sys.argv[3]))
