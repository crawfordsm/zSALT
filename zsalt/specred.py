
"""
SPECREDUCE

General data reduction script for SALT long slit data.

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
from saltprepare import saltprepare
from saltbias import saltbias
from saltgain import saltgain
from saltxtalk import saltxtalk
from saltcrclean import saltcrclean
from saltcombine import saltcombine
from saltflat import saltflat
from saltmosaic import saltmosaic
from saltillum import saltillum

from specidentify import specidentify
from specrectify import specrectify
from specsky import skysubtract
from specextract import extract, write_extract
from specsens import specsens
from speccal import speccal

from PySpectrograph.Spectra import findobj

def specred(infiles, logprefix=None, imreduce=True, specreduce=True, flatfield=False, 
            fringe=False, automethod='Matchlines', inter=True, badpixelimage=None, lampdir="pysalt$data/linelists/",
            cleanup=True):
    """specred is a script to produce 2D-wavelength calibrated data 
       from  SALT longslit observations.  It includes steps not
       currently in the pipeline as well as additional steps to wavelength 
       calibrate the data.   

       Parameters
       ----------
       infiles: list
          input list of files to reduce.  This should be the list of raw data files

       logprefix: string
          Prefix for logifles.   Suggested is to set it to the obsdate
 
       imreduce: boolean
          Run the basic ccd reductions on the data

       specreduce: boolean
          Wavelength calibrate the data

       flatfield: boolean
          Apply a flatfield correction based on data included in infiles. This applies
          an illumination correction and so only corrects for pixel to pixel variance

       fringe: boolean
          Apply a fringe correction to the data.   
          THIS IS CURRENTLY NOT IMPLIMENTED

       specreduce: boolean
          Wavelength calibrate the data

       automethod: string	 
          Method for automatic line identifications.  Options are 'Matchlines'
          and 'Zeropoint'
    
       inter: boolean
          Interactively identify arc lines

       cleanup: boolean
          Remove any temporary files made during the process due to 
          intermediate steps
  
       logfile: string or None
       

       TODO 
       ----
       * Include fringe correction

    """

    obsdate=logprefix
    if logprefix is None: obsdate=''
    logfile='spec'+obsdate+'.log'
    flatimage='FLAT%s.fits' % (obsdate)
    dbfile='spec%s.db' % obsdate
    infile_list=infiles.split(',')

    #create the observation log
    obs_dict=obslog(infile_list)

 
    if imreduce:   
      #prepare the data
      createvar=True
      if badpixelimage is None: createvar=False
      saltprepare(infiles, '', 'p', createvar=createvar, badpixelimage=badpixelimage, clobber=True, logfile=logfile, verbose=True)

      #bias subtract the data
      saltbias('pP*fits', '', 'b', subover=True, trim=True, subbias=False, masterbias='',  
              median=False, function='polynomial', order=5, rej_lo=3.0, rej_hi=5.0, 
              niter=10, plotover=False, turbo=False, 
              clobber=True, logfile=logfile, verbose=True)

      #gain correct the data
      saltgain('bpP*fits', '', 'g', usedb=False, mult=True, clobber=True, logfile=logfile, verbose=True)

      #cross talk correct the data
      saltxtalk('gbpP*fits', '', 'x', xtalkfile = "", usedb=False, clobber=True, logfile=logfile, verbose=True)

 
      #flat field correct the data
      flat_imgs=''
      for i in range(len(infile_list)):
        if obs_dict['CCDTYPE'][i].count('FLAT'):
           if flat_imgs: flat_imgs += ','
           flat_imgs += 'xgbp'+os.path.basename(infile_list[i])

      if flatfield:
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
      saltmosaic('fxgbpP*fits', '', 'm', geomfile, interp='linear', cleanup=True, geotran=True, clobber=True, logfile=logfile, verbose=True)

      #clean up the images
      if cleanup:
           for f in glob.glob('p*fits'): os.remove(f)
           for f in glob.glob('bp*fits'): os.remove(f)
           for f in glob.glob('gbp*fits'): os.remove(f)
           for f in glob.glob('xgbp*fits'): os.remove(f)
           for f in glob.glob('fxgbp*fits'): os.remove(f)


    #set up the name of the images
    if specreduce:
       for i in range(len(infile_list)):
           if obs_dict['OBJECT'][i].upper().strip()=='ARC':
               lamp=obs_dict['LAMPID'][i].strip().replace(' ', '')
               arcimage='mfxgbp'+os.path.basename(infile_list[i])
               lampfile=iraf.osfn("%s%s.txt" % (lampdir, lamp))
               specidentify(arcimage, lampfile, dbfile, guesstype='rss', 
                  guessfile='', automethod=automethod,  function='legendre',  order=3, 
                  rstep=100, rstart='middlerow', mdiff=10, thresh=2, niter=5, smooth=3,
                  inter=inter, clobber=True, logfile=logfile, verbose=True)

               specrectify(arcimage, outimages='', outpref='x', solfile=dbfile, caltype='line', 
                   function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
                   blank=0.0, clobber=True, logfile=logfile, verbose=True)
     

    objimages=''
    for i in range(len(infile_list)):
       if obs_dict['CCDTYPE'][i].count('OBJECT') and obs_dict['INSTRUME'][i].count('RSS'):
          if objimages: objimages += ','
          objimages+='mfxgbp'+os.path.basename(infile_list[i])

    if specreduce:
      #run specidentify on the arc files
      specrectify(objimages, outimages='', outpref='x', solfile=dbfile, caltype='line', 
           function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
           blank=0.0, clobber=True, logfile=logfile, verbose=True)


    return



if __name__=='__main__':

   rawdir=sys.argv[1]
   prodir=os.path.curdir+'/'

   #get the name of the files
   infile_list=glob.glob(rawdir+'*.fits')
   infiles=','.join(['%s' % x for x in infile_list])
    

   #get the current date for the files
   obsdate=os.path.basename(infile_list[0])[1:9]
   print obsdate

   #set up some files that will be needed
   specred(infiles, logprefix=obsdate, imreduce=True, specreduce=True, flatfield=False, cleanup=True)
