import os
import numpy as np

from astropy.io import fits

from pyraf import iraf
from iraf import pysalt

from PySpectrograph.Spectra import Spectrum


import WavelengthSolution
import spectools as st
import AutoIdentify as ai
from specidentify import specidentify, writeIS
from specrectify import readsolascii

from saltsafelog import logging
import saltsafekey as saltkey


def auto_arc_lens(arc_image, dbfile='wav.db', ndstep=20, logfile='salt.log'):
    """Automatically process an arc image for the SALT Lens project

    """
    hdu = fits.open(arc_image)
    hdr = hdu[0].header
    if hdr['LAMPID'] == 'Xe' and hdr['GRATING'] == 'PG0900' and hdr['GRTILT'] == 15.875:
        print('Automatically processing arc image')
        data = hdu[1].data
        ystart = int(0.5 * len(data))
        xarr = np.arange(len(data[ystart]), dtype='int64')
        farr = data[ystart]

        lampfile = os.path.dirname(__file__)+'/data/Xe.lens'
        guessfile = os.path.dirname(__file__)+'/data/lens.db'

        slines, sfluxes = st.readlinelist(lampfile)
        spectrum = Spectrum.Spectrum(
            slines, sfluxes, dw=0.1, stype='line', sigma=6)
        swarr = spectrum.wavelength
        sfarr = spectrum.flux * farr.max() / spectrum.flux.max()

        soldict = readsolascii(guessfile, {})
        soldict = (soldict[soldict.keys()[0]])
        ws = WavelengthSolution.WavelengthSolution(
             xarr, xarr, function=soldict[7], order=soldict[8])
        ws.func.func.domain = soldict[11]
        ws.set_coef(soldict[10][2])
 
        # start the pre processing
        dcoef = ws.coef * 0.0
        dcoef[0] = 0.5 * 6 *  ndstep
        dcoef[1] = 0.1 * ws.coef[1]
        ws = st.findxcor(xarr, farr, swarr, sfarr, ws,
                         dcoef=dcoef, ndstep=ndstep, best=False, inttype='interp')
        xp, wp = st.crosslinematch(xarr, farr, slines, sfluxes, ws,
                                   res=6, mdiff=20, wdiff=10,
                                   sections=3, sigma=5, niter=5)
        ws = st.findfit(np.array(xp), np.array(wp), ws=ws, thresh=ws.thresh)
        print(ws)
        with logging(logfile, True) as log:
            iws = ai.AutoIdentify(xarr, data, slines, sfluxes, ws, farr=farr,
                      method='Matchlines', rstep=100, istart=ystart, nrows=1,
                      res=6, dres=0.25, mdiff=20, sigma=5,
                      smooth=3, niter=5, dc=5, ndstep=ndstep, 
                      oneline=False, log=log, verbose=True)

            # get the basic information about the spectrograph
            dateobs = saltkey.get('DATE-OBS', hdu[0])
            try:
                utctime = saltkey.get('UTC-OBS', hdu[0])
            except SaltError:
                utctime = saltkey.get('TIME-OBS', hdu[0])

            instrume = saltkey.get('INSTRUME', hdu[0]).strip()
            grating = saltkey.get('GRATING', hdu[0]).strip()
            grang = saltkey.get('GR-ANGLE', hdu[0])
            grasteps = saltkey.get('GRTILT', hdu[0])
            arang = saltkey.get('AR-ANGLE', hdu[0])
            arsteps = saltkey.get('CAMANG', hdu[0])
            rssfilter = saltkey.get('FILTER', hdu[0])
            specmode = saltkey.get('OBSMODE', hdu[0])
            masktype = saltkey.get('MASKTYP', hdu[0]).strip().upper()
            slitname = saltkey.get('MASKID', hdu[0])
            slit = st.getslitsize(slitname)
            xbin, ybin = saltkey.ccdbin(hdu[0], arc_image)
            writeIS(iws, dbfile, dateobs=dateobs, utctime=utctime, instrume=instrume,
                    grating=grating, grang=grang, grasteps=grasteps, arsteps=arsteps,
                    arang=arang, rfilter=rssfilter, slit=slit, xbin=xbin,
                    ybin=ybin, objid=None, filename=arc_image, log=log, verbose=True)
        print(iws)

        #self.findfit()

 
    else:
        lamp = hdr['LAMPID']
        lampfile=iraf.osfn("pysalt$data/linelists/%s.salt" % lamp)
        lampfile = 'Xe.lens'
        specidentify(arc_image, lampfile, dbfile, guesstype='rss', 
                  guessfile=None, automethod='Matchlines',  function='legendre',  order=3,
                  rstep=100, rstart='middlerow', mdiff=20, thresh=5, niter=5, smooth=3,
                  inter=True, clobber=True,  preprocess=True, logfile=logfile, verbose=True)
        print("Running specidenity in interactive mode")


if __name__=='__main__':
   import sys
   auto_arc_lens(sys.argv[1])
