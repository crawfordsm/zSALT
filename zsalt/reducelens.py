import os, sys, glob
reddir = '/Users/crawford/programs/zSALT/zsalt/'
sys.path.append(reddir)

from imred import imred
from specred import specred

bpmfile = '/Users/crawford/research/kepler/bpm_agn.fits'

ddir = sys.argv[1]

os.chdir(ddir)
if not os.path.isdir('sci'): os.mkdir('sci')
os.chdir('sci')

#basic image reuctions
infile_list = glob.glob('../raw/P*fits')
imred(infile_list, './', bpmfile, cleanup=True)

#spectroscopic reductions
#propid = '2014-1-RSA_OTH-012'
propid = '2015-1-MLT-006'
infile_list = glob.glob('m*fits')
#specred(infile_list, propid, inter=True)

#extract the data



