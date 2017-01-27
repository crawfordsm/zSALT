import os, sys, glob

from astropy.io import fits

import imred
from mosred import mosred
import argparse

bpmfile = os.path.dirname(imred.__file__)+'/data/bpm_rss_11.fits'


parser = argparse.ArgumentParser(description='Reduce SALT Lens Data')
parser.add_argument('ddir', help='Top level directory with SALT data')
parser.add_argument('-s', dest='basic_red', default=True, action='store_false',
                    help='Skip basic reduction')
parser.add_argument('--dy', dest='dy', default=-12, type=float, 
                    help='Mask offset')
parser.add_argument('slitmask', help='Slitmask for observations')


args = parser.parse_args()
args.slitmask = os.path.abspath(args.slitmask)
ddir = args.ddir


os.chdir(ddir)
if not os.path.isdir('sci'): os.mkdir('sci')
os.chdir('sci')

#basic image reuctions
infile_list = glob.glob('../raw/P*fits')
# remove alignment images
spec_list = []
for infile in infile_list:
    hdu = fits.open(infile)
    if hdu[0].header['OBSMODE'].strip() != 'IMAGING':
       spec_list.append(infile)
infile_list = spec_list
print spec_list

if args.basic_red: imred.imred(infile_list, './', bpmfile, cleanup=True)

#spectroscopic reductions
propid=None
infile_list = glob.glob('m*fits')
mosred(infile_list, args.slitmask, dy=args.dy)

#extract the data



