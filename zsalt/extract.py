#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved. See LICENSE file for more detail.                   #
"""
SPECEXTRACT extracts a 1-D spectrum from a 2-D data file.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       15 Nov 2010

TODO
----
1. The task still needs to be written

LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import time
import numpy as np
from astropy.io import fits


from PySpectrograph.Spectra import Spectrum

from PySpectrograph.Spectra import apext
from PySpectrograph.Spectra import findobj

debug = True


def extract(hdu, ext=1, method='normal', section=[],
            minsize=3.0, thresh=3.0, convert=True):
    """For a given image, extract a 1D spectra from the image
       and write the spectra to the output file

    """

    ap_list = []
    i = ext
    if hdu[i].name == 'SCI':
        # set up the data, variance, and bad pixel frames
        # first step is to find the region to extract
        data_arr = hdu[i].data
        try:
            var_arr = hdu[hdu[i].header['VAREXT']].data
        except:
            var_arr = None
        try:
            bpm_arr = hdu[hdu[i].header['BPMEXT']].data
        except:
            bpm_arr = None
        var_arr = None
        bpm_arr = None

        xarr = np.arange(len(data_arr[0]))

        # convert using the WCS information
        try:
            w0 = hdu[i].header['CRVAL1']
            dw = hdu[i].header['CD1_1']
        except Exception as e:
            msg = 'Error on Ext %i: %s' % (i, e)
            raise Exception(msg)
        warr = w0 + dw * xarr

        # convert from air to vacuum
        if convert:
            warr = Spectrum.air2vac(warr)

        # set up the sections in case of findobj
        if section is None:
            section = findobj.findObjects(
                data_arr,
                method='median',
                specaxis=1,
                minsize=minsize,
                thresh=thresh,
                niter=5)

        # extract all of the  regions
        for sec in section:
            ap = apext.apext(warr, data_arr, ivar=var_arr)
            y1, y2 = sec
            ap.flatten(y1, y2)
            ap_list.append(ap)

    return ap_list


def write_extract(ofile, ap_list, outformat='ascii', fvar=None, clobber=False):
    """Write out to either a txt file or fits file depending on the extension
       of ofile

    """
    if outformat == 'FITS':
        write_extract_fits(ofile, ap_list, clobber)
    elif outformat == 'ascii':
        write_extract_text(ofile, ap_list, clobber)
    else:
        msg = '%s is not a supported output format' % outformat
        raise Exception(msg)
    return


def write_extract_text(ofile, ap_list, clobber=False):
    """Write out the extracted spectrum to a text file.  If the file already
       exists, this will not overwrite it.  The first

       For each spectrum in ap_list, it will add a columns onto the output file
       so that the first column is always wavelength, the second column is
       flux, and the third column is sigma, and then repeat the flux and sigma
       columns

       ofile: Output file to write

       ap_list:  List of extracted spectrum

       clobber: delete ofile if it already exists


    """
    if os.path.isfile(ofile) and not clobber:
        return

    # open the file
    dout = saltio.openascii(ofile, 'w')

    # first extract warr, assume it is the same for all frames
    warr = ap_list[0].wave

    # write out the spectrum
    for i in range(len(warr)):
        outstr = '%7.3f ' % warr[i]
        for ap in ap_list:
            flux = ap.ldata[i]
            try:
                fvar = abs(ap.lvar[i]) ** 0.5
            except:
                fvar = 1
            outstr += "%7.3f %7.3f " % (flux, fvar)
        outstr += '\n'
        dout.write(outstr)
    dout.close()
    return


def write_extract_fits(ofile, ap_list, clobber=False):
    """Write out the extracted spectrum to a FITS table.  If the file already
        exists, this will not overwrite it.

        For each spectrum in ap_list, it will add another extension to the
        fits file.  Each extension will have the first column as wavelength,
        the second column as counts, and the third column as sigma on the
        counts.

        ofile: Output file to write

        ap_list:  List of extracted spectrum

        clobber: delete ofile if it already exists


    """
    # delete the file
    if os.path.isfile(ofile) and clobber:
        saltio.delete(ofile)

    # create the primary array
    hdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([hdu])

    # create the columns and the
    for ap in ap_list:
        fvar = abs(ap.lvar) ** 0.5
        # create the columns
        col1 = pyfits.Column(
            name='wavelength',
            format='D',
            unit='Angstroms',
            array=ap.wave)
        col2 = pyfits.Column(
            name='counts',
            format='D',
            unit='Counts',
            array=ap.ldata)
        col3 = pyfits.Column(name='counts_err', format='D', array=fvar)

        # add to the table
        tbhdu = pyfits.new_table([col1, col2, col3])
        hdulist.append(tbhdu)

    # write it out
    hdulist.writeto(ofile)
    return
