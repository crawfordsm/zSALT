#! /usr/bin/env python

import sys, os
from pyraf import iraf

badpixel = os.path.dirname(os.path.abspath(__file__))+'/box.txt'

obsdate = sys.argv[1]

def masking():
    os.chdir('%s/sci' % (obsdate))

    # create temporary directory to keep masked images

    if not os.path.isdir('temp'):
       os.mkdir('temp')

    # copy images to avoid overwriting the original files

    os.system('cp x*fits temp')
    os.chdir('temp')

    # apply masking

    iraf.obsolete(_doprint=0)
    iraf.ofixpix(images='x*fits[1]', badpixel=badpixel)

if __name__ == '__main__':
   masking()
