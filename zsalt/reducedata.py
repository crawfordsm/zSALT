import os, sys, glob

from imred import imred
from specred import specred
import argparse

#bpmfile = os.path.dirname(os.path.abspath(__file__))+'/bpm_2x4.fits'
bpmfile = os.path.dirname(__file__)+'/data/bpm_rss_11.fits'


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

#spectroscopic reductions
propid=None
infile_list = glob.glob('m*fits')
specred(infile_list, propid, inter=True, guessfile=None)

#extract the data



