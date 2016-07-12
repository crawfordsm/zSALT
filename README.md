zSALT
=====

Measuring redshifts from SALT longslit and MOS observations.  


Instructions
============

1. Clone the git repository onto your device

        git clone https://github.com/crawfordsm/zSALT.git

2. Run the script 'gettemplates.sh' in the directory to download
   the SDSS templates. 

3. Run the script 'reducedata.py' from the zsalt directory on the directory that has the data that you want to reduce.  


4. The directory where you will run the data should have the format
   such that in the top level directory there should be the
   obsdate (e.g. 20131227) and then the raw directory inside that.
   The raw directory should only include the files that you want
   reduced so any extra files which are sometimes included in your
   directory should be removed.

5. Run reducedata.py with

        python reducedata.py 20131227

   Replace the obsdate with the appropriate observations date and location of that director

6. Reducedata should step through all the tasks needed to produce
   wavelength calibrated images.   There will be one interactive
   step to identify arc lines.

7. Once complete, inspect the final product in the sci directory
   that was created.  Skylines should be straight and cosmic
   rays should be cleaned. 

To extract spectra and measure the redshift
-------------------------------------------

   A. use extractobject.py.  This can be called just
   by giving the frame name in which to extract the object.  It will
   then ask for the central row of the object to be extract.  This 
   can be repeated for different rows. 

        python galextract xmfxgbpP201302010093.fits 513 5

   You will either need the full path to extractobject.py. This will perform a sky subtraction and 
   extraction of the spectrum.  Provide the filename, the y-center of the spectrum to extract, and the 
   half-width of the aperture to extract. 

   B. To calculate the redshift, use redshift.py. You will have to give
   it a template to match and different templates are available 
   for matching.  The format of it is:
 
        python redshift.py [SPECTRUMFILE] [TEMPLATEFILE]

   This will display the spectrum, measure a redshift, and overplot
   the template at the best fit redshift.  This may need to be 
   repeated for different redshifts

To extract the AGN for the LCOGT project
----------------------------------------

Run agnextract.py. Give it a filename, y-center (central row of the object ) and dy (half width of the object) like so:

         python agnextract.py xsmfxgbpP201505200016.fits 552 5
    
This will extract galaxy spectra from the image.
