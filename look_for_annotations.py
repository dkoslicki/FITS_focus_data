# This is a simple script that, given a directory,
# will look for any FITS files that have the ANNOTATE keyword set to 'True' and will print the filename
# and the value of the ANNOTATION keyword.
import astropy
import astropy.io.fits as fits
import os
import sys
import photutils
import shutil

folder = "D:\\Astrophotography subs depreciated\\ghost\\lights"
#folder = 'D:\\Dropbox\\Astrophotography\\Sadr region\\mosaic\\ForAsteroids\\SubFolder5'
#folder = "D:\\Dropbox\\Astrophotography\\Ephemeris\\1663 van den Bos\\subs"
print(folder)
# iterate over all files in the folder
for filename in os.listdir(folder):
    filename = os.path.join(folder, filename)
    # if the file is a FITS file
    if filename.endswith('.fit'):
        # import the fits file
        hdulist = fits.open(filename)
        # get the header
        header = hdulist[0].header
        # get the value of the ANNOTATE keyword, if it exists
        annotation = header.get('ANNOTATE')
        # if the ANNOTATION keyword exists, print the filename and the value of the ANNOTATION keyword
        if annotation is not None:
            print(filename, annotation)
            # get the ephemeris name
            ephemeris = annotation.split(' ')[-1].split(';')[0]
            # create a subfolder for the files that have the ANNOTATE keyword set to 'True'
            new_folder = os.path.join(folder, ephemeris)
            if not os.path.exists(new_folder):
                os.makedirs(new_folder)
            # copy the file to the new folder
            shutil.copy(filename, new_folder)
        # close the fits file
        hdulist.close()
