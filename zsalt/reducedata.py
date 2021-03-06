import os, sys, glob

from imred import imred
from specred import specred
import argparse

#bpmfile = os.path.dirname(os.path.abspath(__file__))+'/bpm_2x4.fits'
bpmfile = os.path.dirname(__file__)+'/data/bpm_rss_11.fits'


parser = argparse.ArgumentParser(description='Reduce SALT Lens Data')
parser.add_argument('ddir', help='Top level directory with SALT data')
parser.add_argument('-p', dest='preprocess', default=False, action='store_true',
                    help='prepocess the line identification')
parser.add_argument('-a', dest='auto', default=False, action='store_true',
                    help='Autoidentify the lines')
parser.add_argument('-s', dest='basic_red', default=True, action='store_false',
                    help='Skip basic reduction')
parser.add_argument('--g', dest='guess_file', default=None, type=str, 
                    help='Initial guess file for reductions')

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
specred(infile_list, propid, inter=not args.auto, guessfile=args.guess_file, preprocess=args.preprocess)

#extract the data


