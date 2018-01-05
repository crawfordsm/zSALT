import os, sys, glob
import argparse
import string 

import numpy as np
from astropy.io import fits
from scipy import signal

from pyraf import iraf
from iraf import pysalt
from specrectify import specrectify
from saltobslog import obslog


from imred import imred
from autoarclens import auto_arc_lens
from combine_dithers import combine_dithers
from source_extract import source_extract
from redshift_all import redshift_all
from pretty_redshift import pretty_redshift

#bpmfile = os.path.dirname(os.path.abspath(__file__))+'/bpm_2x4.fits'
bpmfile = os.path.dirname(__file__)+'/data/bpm_rss_11.fits'
lampfile = os.path.dirname(__file__)+'/data/Xe.lens'
guessfile = os.path.dirname(__file__)+'/data/lens.db'

parser = argparse.ArgumentParser(description='Reduce SALT Lens Data')
parser.add_argument('ddir', help='Top level directory with SALT data')
parser.add_argument('-s', dest='basic_red', default=True, action='store_false',
                    help='Skip basic reduction')
args = parser.parse_args()


ddir = args.ddir
os.chdir(ddir)
if not os.path.isdir('sci'): os.mkdir('sci')
os.chdir('sci')

#basic image reuctions
infile_list = glob.glob('../raw/P*fits')
if args.basic_red: imred(infile_list, './', bpmfile, cleanup=True)


#proceed with the spectroscopic reductions
logfile = 'spec{}.log'.format(ddir)
dbfile = 'spec{}.db'.format(ddir)

infile_list = glob.glob('m*fits')
infiles=','.join(['%s' % x for x in infile_list])
obsdate=os.path.basename(infile_list[0])[7:15]
obs_dict=obslog(infile_list)
for i in range(len(infile_list)):
    if obs_dict['CCDTYPE'][i].upper().strip()=='ARC':
       arc_image = infile_list[i]
       auto_arc_lens(arc_image, dbfile=dbfile, ndstep=20, logfile=logfile)

obj_dict = {}
for i in range(len(infile_list)):
    if obs_dict['CCDTYPE'][i].count('OBJECT') and obs_dict['INSTRUME'][i].count('RSS'):
          img = infile_list[i]
          specrectify(img, outimages='', outpref='x', solfile=dbfile, caltype='line',
            function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
            blank=0.0, nearest=True, clobber=True, logfile=logfile, verbose=True)
      
          
          obj = obs_dict['OBJECT'][i].strip()
          if obj in obj_dict.keys():
             obj_dict[obj].append('x'+img)
          else:
             obj_dict[obj] = ['x'+img]


# combine the data if necessary
red_list = []
for obj in obj_dict.keys():
    if len(obj_dict[obj])>1:
       print(obj_dict[obj])
       output = combine_dithers(obj_dict[obj])
       red_list.append(output)
    else:
       red_list.append(obj_dict[obj][0])


# extract all of the sources
spec_list = []
for img in red_list:
    hdu = fits.open(img)
    farr = hdu[1].data.sum(axis=1)
    farr = farr - np.median(farr)
    farr[farr<0] = 0
    xp = signal.find_peaks_cwt(farr, np.array([5]))
    for x in xp:
        print('Extracting position {} from {}'.format(x, img))
        outfile = source_extract(img, int(x), 10)
        spec_list.append(outfile)

    

# run the redshift on each source

fout = open('results_{}.txt'.format(obsdate), 'w')
for spec in  spec_list:
    t, z, cc = redshift_all(spec, z1=0.0001, z2=1.200, plot_result=False, show=False)
    template = os.path.dirname(__file__) +  '/template/spDR2-0{}.fit'.format(string.zfill(t, 2))
    name = spec.split('_')[0]
    xp = spec.split('_')[2].strip('.txt')
    pretty_redshift(spec, template, z, name, False)
    out_str = '{} {} {} {} {}\n'.format(name, obsdate, xp, z, t)
    fout.write(out_str)
    print(out_str)
fout.close()   
