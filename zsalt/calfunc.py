#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
"""
SPECCAL corrections a given observation by a calibration curve
and the extinction curve for the site.  The task assumes a 1-D spectrum that
has already been caled from the original observations.


Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       21 Mar 2011

TODO
----

LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import time
import numpy as np
import pyfits

from PySpectrograph.Spectra import Spectrum
from PySpectrograph.Utilities.fit import interfit


debug = True


def calfunc(obs_spectra, std_spectra, ext_spectra,
            airmass, exptime, error=False):
    """Given an observe spectra, calculate the calibration curve for the
       spectra.  All data is interpolated to the
       binning of the obs_spectra.  The calibrated spectra is then calculated
       from:
       C =  F_obs/ F_std / 10**(-0.4*A*E)/T/dW
       where F_obs is the observed flux from the source,  F_std  is the
       standard spectra, A is the airmass, E is the
       extinction in mags, T is the exposure time and dW is the bandpass

    Parameters
    -----------
    obs_spectra--spectrum of the observed star (counts/A)
    std_spectra--know spectrum of the standard star (ergs/s/cm2/A)
    ext_spectra--spectrum of the extinction curve (in mags)
    airmass--airmass of the observations
    exptime--exposure time of the observations
    function
    """

    # re-interpt the std_spectra over the same wavelength
    std_spectra.interp(obs_spectra.wavelength)

    # re-interp the ext_spetra over the sam ewavelength
    ext_spectra.interp(obs_spectra.wavelength)

    # create the calibration spectra
    cal_spectra = Spectrum.Spectrum(
        obs_spectra.wavelength,
        obs_spectra.flux.copy(),
        stype='continuum')

    # set up the bandpass
    bandpass = np.diff(obs_spectra.wavelength).mean()

    # correct for extinction
    cal_spectra.flux = obs_spectra.flux / \
        10 ** (-0.4 * airmass * ext_spectra.flux)

    # correct for the exposure time and calculation the calitivity curve
    cal_spectra.flux = cal_spectra.flux / exptime / bandpass / std_spectra.flux

    # correct the error calc
    if error:
        cal_spectra.var = obs_spectra.var * cal_spectra.flux / obs_spectra.flux

    return cal_spectra

