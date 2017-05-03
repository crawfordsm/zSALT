#! /usr/bin/env python

import sys, os, time, pyfits

str_line = os.path.dirname(os.path.abspath(__file__))+'/str_line.reg'

obsdate = sys.argv[1]

os.chdir('%s/sci' % (obsdate))

images = [img for img in os.listdir('.') if img.startswith('x')]

for img in images:
    hdulist = pyfits.open(img)
    prihdr = hdulist[0].header
    if prihdr['OBJECT'] == 'ARC':
       print img
       os.system('ds9 %s -regions %s' % (img, str_line))
       time.sleep(2)
