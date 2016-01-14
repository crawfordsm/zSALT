
"""
specred

Process spectral reductions of the data and produce
the output spectra for each object

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter

from pyraf import iraf
from iraf import pysalt

from saltobslog import obslog

from specidentify import specidentify
from specrectify import specrectify
from specarcstraighten import specarcstraighten


def specred(infile_list, propcode=None, inter=True, guessfile = None, automethod='Matchlines'):

    #set up the files
    infiles=','.join(['%s' % x for x in infile_list])
    obsdate=os.path.basename(infile_list[0])[7:15]

    #set up some files that will be needed
    logfile='spec'+obsdate+'.log'
    dbfile='spec%s.db' % obsdate
    sdbfile='spec_straight_%s.db' % obsdate

    #create the observation log
    obs_dict=obslog(infile_list)
    for i in range(len(infile_list)):
           print infile_list[i], obs_dict['OBJECT'][i].upper().strip(), obs_dict['PROPID'][i].upper().strip()
           if obs_dict['OBJECT'][i].upper().strip()=='ARC' and obs_dict['PROPID'][i].upper().strip()==propcode:
               lamp=obs_dict['LAMPID'][i].strip().replace(' ', '')
               arcimage=os.path.basename(infile_list[i])
               if lamp == 'NONE': lamp='CuAr'
               lampfile=iraf.osfn("pysalt$data/linelists/%s.salt" % lamp)
               #lampfile='/Users/crawford/research/kepler/Xe.salt' 

               #straighten the arc
               specarcstraighten(arcimage, sdbfile , function='legendre', order=2, rstep=10,
                      rstart='middlerow', nrows=1, dcoef=[0.5, 1.0, 0.0], ndstep=10,
                      startext=0, clobber=False, logfile='salt.log', verbose=True)

               #rectify it
               specrectify(arcimage, outimages='', outpref='s', solfile=sdbfile, caltype='line',
                   function='legendre',  order=2, inttype='interp', w1=None, w2=None, dw=None, nw=None,
                   blank=0.0, clobber=True, logfile=logfile, verbose=True)

               #idnetify the line
               if guessfile is None:
                   guesstype='rss'
                   guessfile=''
               else:
                   guesstype='file'
                   
               
               specidentify('s'+arcimage, lampfile, dbfile, guesstype=guesstype,
                  guessfile=guessfile, automethod=automethod,  function='legendre',  order=3,
                  rstep=100, rstart='middlerow', mdiff=20, thresh=5, niter=5, smooth=3,
                  inter=inter, clobber=True, logfile=logfile, verbose=True)

               #apply the final rectification
               specrectify('s'+arcimage, outimages='', outpref='x', solfile=dbfile, caltype='line',
                   function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
                   blank=0.0, clobber=True, logfile=logfile, verbose=True)

    objimages=''
    spec_list=[]
    for i in range(len(infile_list)):
       if obs_dict['CCDTYPE'][i].count('OBJECT') and obs_dict['INSTRUME'][i].count('RSS')  and obs_dict['PROPID'][i].upper().strip()==propcode:
          img = infile_list[i]
          # rectify it
          specrectify(arcimage, outimages='', outpref='s', solfile=sdbfile, caltype='line',
              function='legendre',  order=2, inttype='interp', w1=None, w2=None, dw=None, nw=None,
              blank=0.0, clobber=True, logfile=logfile, verbose=True)

          specrectify('s'+img, outimages='', outpref='x', solfile=dbfile, caltype='line',
            function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
            blank=0.0, nearest=True, clobber=True, logfile=logfile, verbose=True)


  


    


