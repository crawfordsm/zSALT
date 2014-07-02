zSALT
=====

Measuring redshifts from SALT longslit and MOS observations.  


Instructions
============

1. Clone the git repository onto your device

   git clone https://github.com/crawfordsm/zSALT.git

2. Run the script 'gettemplates.sh' in the directory to download
   the SDSS templates. 

3. Copy the script 'reducedata.py' to the directory that has the
   data that you want to reduce.   Edit this file to point to the 
   zsalt directory, the proposal code, and the bad pixel file map.
   If you are using MOS data, edit the file to point to the xml
   file describing the MOS mask and to use mosred instead of
   specred.

4. The directory where you will run the data should have the format
   such that in the top level directory there should be the
   obsdate (e.g. 20131227) and then the raw directory inside that.
   The raw directory should only include the files that you want
   reduced so any extra files which are sometimes included in your
   directory should be removed.

5. Run reducedata.py with

   python reducedata.py 20131227

   Replace the obsdate with the appropriate observations date

6. Reducedata should step through all the tasks needed to produce
   wavelength calibrated images.   There will be one interactive
   step to identify arc lines.

7. Once complete, inspect the final product in the sci directory
   that was created.  Skylines should be straight and cosmic
   rays should be cleaned. 

8. To extract spectra, use extractobject.py.  This can be called just
   by giving the frame name in which to extract the object.  It will
   then ask for the central row of the object to be extract.  This 
   can be repeated for different rows. 

   python extractobject.py xmfxgbpP201302010093.fits

   You will either need the full path to extractobject.py or copy it 
   into the directory.  This will perform a sky subtraction and 
   extraction of the spectrum.

9. To calculate the redshift, use redshift.py. You will have to give
   it a template to match and different templates are available 
   for matching.  The format of it is:
 
   python redshift.py [SPECTRUMFILE] [TEMPLATEFILE]

   This will display the spectrum, measure a redshift, and overplot
   the template at the best fit redshift.  This may need to be 
   repeated for different redshifts


