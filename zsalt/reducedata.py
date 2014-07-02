import os, sys, glob
reddir = '[PATHTOZSALT]/zSALT/zsalt/'
sys.path.append(reddir)
bpmfile = '[PATHTOBPM]/bpm.fits'
propid = '[PROPID]'

from imred import imred
from specred import specred


ddir = sys.argv[1]

os.chdir(ddir)
if not os.path.isdir('sci'): os.mkdir('sci')
os.chdir('sci')

#basic image reuctions
infile_list = glob.glob('../raw/P*fits')
imred(infile_list, './', bpmfile, cleanup=True)

#spectroscopic reductions
infile_list = glob.glob('m*fits')
specred(infile_list, propid, inter=True)

#UNEDIT IF YOU WANT TO REDUCE MOS DATA
#mosxml = [MOSXMLFILE]
#mosred = infile_list, propid, mosxml, inter=True)





