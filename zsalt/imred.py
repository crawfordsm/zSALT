
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
from astropy.io import fits
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

def imred(infile_list, prodir, bpmfile=None, gaindb = None, geomfile = None, cleanup=True):

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
    if geomfile is None: geomfile=iraf.osfn("pysalt$data/rss/RSSgeom.dat")
    try:
       saltmosaic('fxgbpP*fits', '', 'm', geomfile, interp='linear', cleanup=True, geotran=True, clobber=True, logfile=logfile, verbose=True)
    except:
       saltmosaic('fxgbpP*fits', '', 'm', geomfile, interp='linear', cleanup=True, geotran=True, clobber=True, logfile=logfile, verbose=True)

    #add needed header columns and make a better mask
    for img in infile_list:
        filename = 'mfxgbp'+os.path.basename(img)
        hdu = fits.open(filename, 'update')
        hdu[2].header.update('EXTNAME','VAR')
        hdu[3].header.update('EXTNAME','BPM')
        bpm_rc = (hdu[3].data>0).astype('uint8')
        zeroscicol = hdu['SCI'].data.sum(axis=0) == 0
        bpmgapcol = bpm_rc.mean(axis=0) == 1
        addbpmcol = zeroscicol & ~bpmgapcol
        addbpmcol[np.argmax(addbpmcol)-4:np.argmax(addbpmcol)] = True    # allow for chip tilt
        bpm_rc[:,addbpmcol] = 1
        hdu[3].data = bpm_rc
        hdu.writeto(filename,clobber=True)


    #clean up the images
    if cleanup:
           for f in glob.glob('p*fits'): os.remove(f)
           for f in glob.glob('bp*fits'): os.remove(f)
           for f in glob.glob('gbp*fits'): os.remove(f)
           for f in glob.glob('xgbp*fits'): os.remove(f)
           for f in glob.glob('fxgbp*fits'): os.remove(f)


def add_variance(filenames, bpmfile):
    file_list=glob.glob(filenames)
    badpixelstruct = fits.open(bpmfile)
    for f in file_list:
        struct = fits.open(f)
        nsciext=len(struct)-1
        nextend=nsciext
        for i in range(1, nsciext+1):
            hdu=CreateVariance(struct[i], i, nextend+i)
            hdu.header.update('EXTNAME','VAR')
            struct[i].header['VAREXT'] = (nextend+i, 'Extension for Variance Frame')
            struct.append(hdu)
        nextend+=nsciext
        for i in range(1, nsciext+1):
            if os.path.basename(bpmfile)=='bpm_rss_11.fits':
                hdu=masterbadpixel(struct, badpixelstruct, i, nextend+i)
            else:
                hdu=createbadpixel(struct, badpixelstruct, i, nextend+i)
            struct[i].header['BPMEXT']=(nextend+i, 'Extension for Bad Pixel Mask')
            hdu.header.update('EXTNAME','BPM')
            struct.append(hdu)
        nextend+=nsciext
        struct[0].header['NEXTEND'] = nextend
        if os.path.isfile(f): os.remove(f)
        struct.writeto(f)

def masterbadpixel(inhdu, bphdu, sci_ext, bp_ext):
    """Create the bad pixel hdu bp_ext for inhdu[sci_ext] from a master, bphdu
    """

    if bphdu is None:
        data=np.zeros_like(inhdu[sci_ext].data).astype("uint8")
    else:
        infile=inhdu.fileinfo(0)['filename']
        bpfile=bphdu.fileinfo(0)['filename']
        masternext = len(bphdu)-1
        masterext = (sci_ext-1) % masternext + 1        # allow for windows

        if not saltkey.compare(inhdu[0], bphdu[0], 'INSTRUME', infile, bpfile):
            message = '%s and %s are not the same %s' % (infile,bpfile, 'INSTRUME')
            raise SaltError(message)
        else:
            rows,cols = inhdu[sci_ext].data.shape
            cbin,rbin = np.array(inhdu[sci_ext].header["CCDSUM"].split(" ")).astype(int)
            masterrows,mastercols = bphdu[masterext].data.shape
            master_rc = np.ones((masterrows+(masterrows % rbin),mastercols+(mastercols % cbin)))
            master_rc[:masterrows,:mastercols] = bphdu[masterext].data
            masterrows,mastercols=(masterrows+(masterrows % rbin),mastercols+(mastercols % cbin))
            ampsec = inhdu[sci_ext].header["AMPSEC"].strip("[]").split(",")
            r1,r2 = (np.array(ampsec[1].split(":")).astype(int) / rbin) * rbin
            c1,c2 = (np.array(ampsec[0].split(":")).astype(int) / cbin) * cbin
            if c1 > c2: c1,c2 = c2,c1
            bin_rc = (master_rc.reshape(masterrows/rbin,rbin,mastercols/cbin,cbin).sum(axis=3).sum(axis=1) > 0)
            data = bin_rc[ r1:r2, c1:c2 ].astype('uint8')

    header=inhdu[sci_ext].header.copy()
    header.update('EXTVER',bp_ext)
    header.update('SCIEXT',sci_ext,comment='Extension of science frame')

    return fits.ImageHDU(data=data, header=header, name='BPM')


if __name__=='__main__':
   raw_files=glob.glob(sys.argv[1]+'P*fits')
   bpmfile = sys.argv[2]
   prodir=os.path.curdir+'/'
   imred(raw_files, prodir, cleanup=True, bpmfile=bpmfile)
