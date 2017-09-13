import os 
import sys 

import numpy as np
from scipy.ndimage import interpolation
from scipy.signal import correlate
from astropy.io import fits

import pylab as pl

def calc_offset(hdu_list):
 
    y0 = None
    pos = []
    n_images = len(hdu_list)
    for i in range(n_images):
        yarr = hdu_list[i][1].data.sum(axis=1) 
        yarr = yarr - np.median(yarr)
        yarr[yarr<0] = 0
        if y0 is None: 
           y0 = 1.0 * yarr
           pos = [0.00]
        else:
           z = correlate(y0, yarr)
           pos.append(len(y0) - z.argmax())
   
    return pos

def combine_dithers(image_list, pos=None):

    # open the set of hdus
    hdu_list = [fits.open(img) for img in image_list]
    n_images = len(hdu_list)

    # determine the position of each image
    if pos is None: pos = calc_offset(hdu_list)

    # shift each image
    # assumes there is a variance and bpm, but just does a simple shift
    for i in range(n_images):
        hdu_list[i][1].data = interpolation.shift(hdu_list[i][1].data, np.array([pos[0] - pos[i], 0]))
        try:
           hdu_list[i][2].data = interpolation.shift(hdu_list[i][2].data, np.array([pos[0] - pos[i], 0]))**2
           hdu_list[i][3].data = interpolation.shift(hdu_list[i][3].data, np.array([pos[0] - pos[i], 0])) 
        except:
           pass
  
        # sum the results
        if i>0:
           hdu_list[0][1].data += hdu_list[i][1].data
           hdu_list[0][2].data += hdu_list[i][2].data
           hdu_list[0][3].data *= hdu_list[i][3].data
 
    hdu_list[0][1].data /= n_images
    hdu_list[0][2].data /= hdu_list[0][2].data**0.5 / n_images
    obsdate = hdu_list[0][0].header['DATE-OBS'].replace('-','')
    output = '{}_P{}.fits'.format(hdu_list[0][0].header['OBJECT'].replace(' ', '_'), obsdate)
    hdu_list[0].writeto(output, overwrite=True)
    return output


if __name__ == "__main__":

    combine_dithers(sys.argv[1:])
    exit()
