
"""
mosred 

Process MOS spectral reductions of the data and produce
the output spectra for each object

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter

from pyraf import iraf
from iraf import pysalt

from saltobslog import obslog

from specslit import specslit
from specidentify import specidentify
from specrectify import specrectify

def mosred(infile_list, slitmask,propcode=None, dy=0, inter=True, automethod='Matchlines'):

    #set up the files
    infiles=','.join(['%s' % x for x in infile_list])
    obsdate=os.path.basename(infile_list[0])[7:15]

    #set up some files that will be needed
    logfile='spec'+obsdate+'.log'
    dbfile='spec%s.db' % obsdate

    #create the observation log
    obs_dict=obslog(infile_list)


    #apply the mask to the data sets
    for i in range(len(infile_list)):
        specslit(image=infile_list[i], outimage='', outpref='s', exttype='rsmt', slitfile='../../P001423N01.xml', 
                 outputslitfile='', regprefix='ds_', sections=3, width=25.0, sigma=2.2, thres=6.0, order=1, padding=5, yoffset=dy, 
                 inter=False, clobber=True, logfile=logfile, verbose=True)

        
    exit()

    for i in range(len(infile_list)):
           if obs_dict['OBJECT'][i].upper().strip()=='ARC' and obs_dict['PROPID'][i].upper().strip()==propcode:
               lamp=obs_dict['LAMPID'][i].strip().replace(' ', '')
               arcimage=os.path.basename(infile_list[i])
               if lamp == 'NONE': lamp='CuAr'
               lampfile=iraf.osfn("pysalt$data/linelists/%s.salt" % lamp)

               specidentify(arcimage, lampfile, dbfile, guesstype='rss',
                  guessfile='', automethod=automethod,  function='legendre',  order=3,
                  rstep=100, rstart='middlerow', mdiff=20, thresh=3, niter=5, smooth=3,
                  inter=inter, clobber=True, logfile=logfile, verbose=True)

               #specrectify(arcimage, outimages='', outpref='x', solfile=dbfile, caltype='line',
               #    function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
               #    blank=0.0, clobber=True, logfile=logfile, verbose=True)

    objimages=''
    spec_list=[]
    for i in range(len(infile_list)):
       if obs_dict['CCDTYPE'][i].count('OBJECT') and obs_dict['INSTRUME'][i].count('RSS')  and obs_dict['PROPID'][i].upper().strip()==propcode:
          img = infile_list[i]
          ##rectify it
          specrectify(img, outimages='', outpref='x', solfile=dbfile, caltype='line',
            function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
            blank=0.0, clobber=True, logfile=logfile, verbose=True)


  


    


