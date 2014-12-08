import os, sys, glob
reddir = '/Users/crawford/programs/zSALT/zsalt/'
sys.path.append(reddir)

from imred import imred
from specred import specred

bpmfile = os.getcwd() + '/bpm.fits'
print bpmfile

ddir = sys.argv[1]
propid = sys.argv[2]

os.chdir(ddir)
if not os.path.isdir('sci'): os.mkdir('sci')
os.chdir('sci')

#basic image reuctions
infile_list = glob.glob('../raw/P*fits')
#imred(infile_list, './', bpmfile, cleanup=True)

#spectroscopic reductions
infile_list = glob.glob('m*fits')
specred(infile_list, propid, inter=True)

#extract the data



