import sys

from pyraf import iraf
from iraf import pysalt
from specidentify import specidentify

from astropy.io import fits

header = fits.open(sys.argv[1])[0].header
lampid = header['LAMPID'].strip().replace(' ', '')
print lampid
specidentify(images=sys.argv[1], linelist='/Users/crawford/programs/pysalt/data/linelists/%s.salt' % lampid, outfile='wav.db', guesstype='rss', guessfile='', automethod='Matchlines', function='poly', order=3, rstep=100, rstart='middlerow', mdiff=20, thresh=3, niter=5, smooth=5, inter=True, startext=0 ,clobber=True, textcolor='black', logfile='salt.log',
verbose=True)
