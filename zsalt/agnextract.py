
"""
SNEXTRACT

extract spectra from image or images and combine them

"""

import os, sys, glob, shutil, string

import numpy as np
import pyfits
from scipy import ndimage as nd


from pyraf import iraf
from iraf import pysalt

from specsky import skysubtract
from specextract import extract, write_extract
from specslitnormalize import specslitnormalize
from specsens import specsens
from speccal import speccal

from PySpectrograph.Spectra import findobj

def specextract(img, calfile=None):

    #set up some files that will be needed
    logfile='specext.log'

    #create the spectra text files for all of our objects
    spec_list=[]
    spec_list.extend(extract_spectra(img, calfile=calfile, smooth=False, clobber=True))
    print spec_list

 

def speccombine(spec_list, obsdate):
   """Combine N spectra"""
   print spec_list[0]
   w1,f1, e1=np.loadtxt(spec_list[0][0], usecols=(0,1,2), unpack=True)

   w=w1
   f=1.0*f1
   e=e1**2

   for sfile in spec_list:
      w2,f2, e2=np.loadtxt(sfile[0], usecols=(0,1,2), unpack=True)
      if2=np.interp(w1, w2, f2)
      ie2=np.interp(w1, w2, e2)
      f2=f2*np.median(f1/if2)
      f+=if2
      e+=ie2**2

   f=f/len(spec_list)
   e=e**0.5/len(spec_list)
 
   return w,f,e


def cleanspectra(w, f, e=None, dw=2,  grow=6, neg=False):
    """Remove possible bad pixels"""
    if e is None: e = f*0.0+f.std()
 
    m = (w>0)
    #set a few unreliable sky lines to zero
    for l in [5577, 6300, 6364]:
        m[abs(w-l)<dw]=0

    if neg: m = m * (f>0)

    #remove and grow the bad areas
    m = nd.minimum_filter(m, size=grow)

    return w[m],f[m],e[m]

  
def write_spectra(sfile, w, f, e):
    fout=open(sfile, 'w')
    for i in range(len(w)):
        if f[i]!=np.nan:
           fout.write('%f %e %e\n' % (w[i], f[i], e[i]))
    fout.close()
 
def normalizespectra(sfile, compfile):
    """Normalize spectra by the comparison object"""

    #read in the spectra
    w,f,e=np.loadtxt(sfile, usecols=(0,1,2), unpack=True)
   
    #read in the comparison spectra
    cfile=sfile.replace('MCG-6-30-15', 'COMP')
    print cfile
    wc,fc,ec=np.loadtxt(cfile, usecols=(0,1,2), unpack=True)

    #read in the base star
    ws,fs,es=np.loadtxt(compfile, usecols=(0,1,2), unpack=True)
 
    #calcualte the normalization
    ifc=np.interp(ws, wc, fc) 
    norm=np.median(fs/ifc)
    print norm
    f=norm*f
    e=norm*e

    #write out the result
    fout=open(sfile, 'w')
    for i in range(len(w)):
        fout.write('%f %e %e\n' % (w[i], f[i], e[i]))
    fout.close()

    #copy

    
 

def extract_spectra(img, yc=None, oy=10, dy=50, minsize=5, thresh=3, findobject=False, 
                    niter=5, calfile=None, smooth=False, maskzeros=False, clobber=True):
    """Create a list of spectra for each of the objects in the images"""
    #okay we need to identify the objects for extraction and identify the regions for sky extraction
    #first find the objects in the image

    #skynormalize the data
    #specslitnormalize(img, 'n'+img, '', response=None, response_output=None, order=3, conv=1e-2, niter=20,
    #                 startext=0, clobber=False,logfile='salt.log',verbose=True)

    print 'Extract Spectra from ', img
    hdu=pyfits.open(img)
    target=hdu[0].header['OBJECT'].replace(' ', '')
    propcode=hdu[0].header['PROPID']
    airmass=hdu[0].header['AIRMASS']
    exptime=hdu[0].header['EXPTIME']

    if smooth:
       data=smooth_data(hdu[1].data)
    else:
       data=hdu[1].data

    #replace the zeros with the average from the frame
    if maskzeros:
       mean,std=iterstat(data[data>0])
       #rdata=mean  np.random.normal(mean, std, size=data.shape)
       print mean, std
       data[data<=0]=mean #rdata[data<=0]

    #use manual intervention to get section
    if findobject:
       section=findobj.findObjects(data, method='median', specaxis=1, minsize=minsize, thresh=thresh, niter=niter)
         
       if yc is None and len(section)>0:
          yc = np.mean(section[0])
       elif len(section)>0:
          diff = 1e6
          for i in range(len(section)):
              y = np.mean(section[i])
              if abs(y-yc) < diff: 
                 bestyc = y
                 diff = abs(y-yc)
          yc = bestyc
              

    if yc is None:
        os.system('ds9 %s &' % img)
    print len(hdu)
    if len(hdu)==2: 
        print 'Using basic extraction'
        if yc is None:
           y1=int(raw_input('y1:'))
           y2=int(raw_input('y2:'))
           sy1=int(raw_input('sky y1:'))
           sy2=int(raw_input('sky y2:'))
        ap_list=extract(hdu, method='normal', section=[(y1,y2)], minsize=minsize, thresh=thresh, convert=True)
        sk_list=extract(hdu, method='normal', section=[(sy1,sy2)], minsize=minsize, thresh=thresh, convert=True)
        
        ap_list[0].ldata=ap_list[0].ldata-float(y2-y1)/(sy2-sy1)*sk_list[0].ldata
    
        ofile='%s.%s_%i.txt' % (target, extract_date(img), extract_number(img))
        write_extract(ofile, [ap_list[0]], outformat='ascii', clobber=clobber)

        w, f, e = np.loadtxt(ofile, usecols=(0,1,2), unpack=True)
        w, f, e=cleanspectra(w, f, e, neg=True) 
        m = (w>3900)*(w<8100)
        write_spectra(ofile, w[m], f[m], e[m])
  
    else: 
        print 'Using advanced extraction'

        if yc is None: yc=int(raw_input('yc:'))
        

        w0=hdu[1].header['CRVAL1']
        dw=hdu[1].header['CD1_1']
        xarr = np.arange(len(hdu[1].data[0]))
        warr=w0+dw*xarr
       
        print hdu[1].data[yc, 1462], hdu[2].data[yc,1462]
        warr, madata, var = skysub_region(warr, hdu[1].data, hdu[2].data, hdu[3].data, yc, oy, dy)
        w, f, e = masked_extract(warr, madata[yc-oy:yc+oy, :], var[yc-oy:yc+oy, :])
        print yc
        ofile='%s.%s_%i_%i.txt' % (target, extract_date(img), extract_number(img), yc)
        write_spectra(ofile, w, f, e)
        

    if calfile is not None: 
       extfile=iraf.osfn("pysalt$data/site/suth_extinct.dat")
       speccal(ofile, ofile.replace("txt", "spec"), calfile, extfile, airmass, exptime, clobber=True, logfile='salt.log', verbose=True)
    
    spec_list=[ofile, airmass, exptime, propcode]

    return spec_list

def masked_extract(w, madata, var=None, grow=10):
    """Extraction of spectra from an array.  The extration
       returns a weighted average of an array where the weighting
       is based on the distance from the center line and the
       variance of the frame
    """
    print w.min(), w.max()
    print madata[10, 1062], var[10,1062]
    ylen = len(madata)
    ywei = abs((np.arange(ylen) - 0.5*ylen)) + 1
    ywei = 1.0 / ywei
    ywei.shape = (ylen, 1)
    if var is None:
       ewei = abs(f)
    else:
       ewei = abs(madata)/var
    weights =  ewei * ywei
    f = np.ma.average(madata, axis=0, weights=weights)
    e = np.ma.average(var/abs(madata), axis=0, weights=weights)
    e = e * abs(f)
    e = e**0.5   / ywei.sum()**0.5
    w, f, e = cleanspectra(w, f, e, grow=grow)

    return w,f,e

def skysub_region(warr, data, var, mask,  yc, oy, dy, order=1, x1=670, x2=2100, grow=10):
    """Within the region of interested specified by yc-dy to yc+dy, exclude  
       any possible object and fit a polynomial to the remaining sky 
       background pixels.   

    """
    madata = np.ma.array(data, mask=mask)

    #fit the background along each column
    for i in range(x1,x2):
        y = data[yc-dy:yc+dy, i]
        x = np.arange(len(y))
        m = (mask[yc-dy:yc+dy, i]==0) * (abs (x - 0.5*len(x)) > oy)
        if y[m].any():
           #g = fitdata(x[m], y[m]) 
           #madata[yc-dy:yc+dy,i] = madata[yc-dy:yc+dy,i] - g(x)
           coef = np.polyfit(x[m],y[m],1)
           madata[yc-dy:yc+dy,i] = madata[yc-dy:yc+dy,i] - np.polyval(coef, x)
    return warr[x1:x2], madata[:, x1:x2], var[:,x1:x2]


def smooth_data(data, mbox=25):
    mdata=nd.filters.median_filter(data, size=(mbox, mbox))
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
    obsdate=int(img[0][-12:-4])
    return string.zfill(obsdate,8)

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
   snextract(sys.argv[1], calfile=sys.argv[2])
