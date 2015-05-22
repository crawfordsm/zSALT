
"""
IMRED 

Reduction script for SALT data -- this is 
for science level reductions with variance frames

This includes step that are not yet included in the pipeline 
and can be used for extended reductions of SALT data. 

It does require the pysalt package to be installed 
and up to date.

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter

 
from pyraf import iraf
from iraf import pysalt

from saltobslog import obslog
from saltprepare import  *
from saltbias import saltbias
from saltgain import saltgain
from saltxtalk import saltxtalk
from saltcrclean import saltcrclean
from saltcombine import saltcombine
from saltflat import saltflat
from saltmosaic import saltmosaic
from saltillum import saltillum

def imred(infile_list, prodir, bpmfile=None, gaindb = None, cleanup=True):

    #get the name of the files
    infiles=','.join(['%s' % x for x in infile_list])
    

    #get the current date for the files
    obsdate=os.path.basename(infile_list[0])[1:9]
    print obsdate

    #set up some files that will be needed
    logfile='spec'+obsdate+'.log'
    flatimage='FLAT%s.fits' % (obsdate)
    dbfile='spec%s.db' % obsdate

    #create the observation log
    obs_dict=obslog(infile_list)

 
    #prepare the data
    saltprepare(infiles, '', 'p', createvar=False, badpixelimage='', clobber=True, logfile=logfile, verbose=True)

    #bias subtract the data
    saltbias('pP*fits', '', 'b', subover=True, trim=True, subbias=False, masterbias='',  
              median=False, function='polynomial', order=5, rej_lo=3.0, rej_hi=5.0, 
              niter=10, plotover=False, turbo=False, 
              clobber=True, logfile=logfile, verbose=True)

    add_variance('bpP*fits', bpmfile)

    #gain correct the data 
    usedb = False
    if gaindb: usedb = True
    saltgain('bpP*fits', '', 'g', gaindb=gaindb, usedb=usedb, mult=True, clobber=True, logfile=logfile, verbose=True)

    #cross talk correct the data
    saltxtalk('gbpP*fits', '', 'x', xtalkfile = "", usedb=False, clobber=True, logfile=logfile, verbose=True)

 
    #flat field correct the data
    flat_imgs=''
    for i in range(len(infile_list)):
        if obs_dict['CCDTYPE'][i].count('FLAT'):
           if flat_imgs: flat_imgs += ','
           flat_imgs += 'xgbp'+os.path.basename(infile_list[i])

    if 0: #len(flat_imgs)!=0:
         saltcombine(flat_imgs,flatimage, method='median', reject=None, mask=False,    \
                weight=False, blank=0, scale=None, statsec='[200:300, 600:800]', lthresh=3,    \
                hthresh=3, clobber=True, logfile=logfile, verbose=True)
         saltillum(flatimage, flatimage, '', mbox=11, clobber=True, logfile=logfile, verbose=True)

         saltflat('xgbpP*fits', '', 'f', flatimage, minflat=0.8, allext=False, clobber=True, logfile=logfile, verbose=True)
    else:
         flats=None
         imfiles=glob.glob('xgbpP*fits')
         for f in imfiles:
             shutil.copy(f, 'f'+f)

    #cosmic ray clean the data
    #only clean the object data
    for i in range(len(infile_list)):
        if obs_dict['CCDTYPE'][i].count('OBJECT') and obs_dict['INSTRUME'][i].count('RSS'):
          img='fxgbp'+os.path.basename(infile_list[i])
          saltcrclean(img, img, '', crtype='edge', thresh=5, mbox=11, bthresh=5.0,
                flux_ratio=0.2, bbox=25, gain=1.0, rdnoise=5.0, fthresh=5.0, bfactor=2,
                gbox=3, maxiter=5, multithread=True,  clobber=True, logfile=logfile, verbose=True)

    #mosaic the data
    geomfile=iraf.osfn("pysalt$data/rss/RSSgeom.dat")
    try:
       saltmosaic('fxgbpP*fits', '', 'm', geomfile, interp='linear', cleanup=True, geotran=True, clobber=True, logfile=logfile, verbose=True)
    except:
       saltmosaic('fxgbpP*fits', '', 'm', geomfile, interp='linear', cleanup=True, geotran=True, clobber=True, logfile=logfile, verbose=True)

    #clean up the images
    if cleanup:
           for f in glob.glob('p*fits'): os.remove(f)
           for f in glob.glob('bp*fits'): os.remove(f)
           for f in glob.glob('gbp*fits'): os.remove(f)
           for f in glob.glob('xgbp*fits'): os.remove(f)
           for f in glob.glob('fxgbp*fits'): os.remove(f)


def add_variance(filenames, bpmfile):
    file_list=glob.glob(filenames)
    badpixelstruct = saltio.openfits(bpmfile)
    for f in file_list:
        struct = pyfits.open(f)
        nsciext=len(struct)-1
        nextend=nsciext
        for i in range(1, nsciext+1):
            hdu=CreateVariance(struct[i], i, nextend+i)
            struct[i].header.set('VAREXT',nextend+i, comment='Extension for Variance Frame')
            struct.append(hdu)
        nextend+=nsciext
        for i in range(1, nsciext+1):
            hdu=createbadpixel(struct, badpixelstruct, i, nextend+i)
            struct[i].header.set('BPMEXT',nextend+i, comment='Extension for Bad Pixel Mask')
            struct.append(hdu)
        nextend+=nsciext
        struct[0].header.set('NEXTEND', nextend)
        if os.path.isfile(f): os.remove(f)
        struct.writeto(f)

if __name__=='__main__':
   raw_files=glob.glob(sys.argv[1]+'P*fits')
   bpmfile = sys.argv[2]
   prodir=os.path.curdir+'/'
   imred(raw_files, prodir, cleanup=True, bpmfile=bpmfile)
