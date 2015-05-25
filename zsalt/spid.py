import sys

from pyraf import iraf
from iraf import pysalt
from specidentify import specidentify

specidentify(images=sys.argv[1], linelist=sys.argv[2], outfile='wav.db', guesstype='rss', guessfile='', automethod='Matchlines', function='poly', order=3, rstep=100, rstart='middlerow', mdiff=20, thresh=3, niter=5, smooth=3, inter=True, startext=0 ,clobber=True, textcolor='black', logfile='salt.log',
verbose=True)
