import astropy
import astropy.io.fits as fits
import os
import sys
import photutils
import shutil
from collections import Counter
import numpy as np


file = "D:\Dropbox\Astrophotography\Calibration\SX-46\Old BPM\BPM-sxASCOM_driver_Copyright_c_Dad_Dog_Development_Ltd.-4540x3640-1.fits"
# load the fits file
hdulist = fits.open(file)
# get the data
data = hdulist[0].data
vals = Counter(np.ndarray.flatten(data))
good_val = 127
bad_val = 255
print(f"BPM: {os.path.dirname(file)}")
print("Bad pixels: ", vals[bad_val])
print("Good pixels: ", vals[good_val])
print("Ratio: ", vals[bad_val]/vals[good_val])


file = "D:\Dropbox\Astrophotography\Calibration\SX-46\Mid June 2022\BPM-sxASCOM_driver_Copyright_c_Dad_Dog_Development_Ltd.-4540x3640-1.fits"
# load the fits file
hdulist = fits.open(file)
# get the data
data = hdulist[0].data
vals = Counter(np.ndarray.flatten(data))
good_val = 127
bad_val = 255
print(f"BPM: {os.path.dirname(file)}")
print("Bad pixels: ", vals[bad_val])
print("Good pixels: ", vals[good_val])
print("Ratio: ", vals[bad_val]/vals[good_val])

# and a more recent one
file = "D:\Dropbox\Astrophotography\Calibration\SX-46\August 2022 EL panel flats\BPM-sxASCOM_driver_Copyright_c_Dad_Dog_Development_Ltd.-4540x3640.fits"
# load the fits file
hdulist = fits.open(file)
# get the data
data = hdulist[0].data
vals = Counter(np.ndarray.flatten(data))
good_val = 127
bad_val = 255
print(f"BPM: {os.path.dirname(file)}")
print("Bad pixels: ", vals[bad_val])
print("Good pixels: ", vals[good_val])
print("Ratio: ", vals[bad_val]/vals[good_val])

# And the most recent one
file = "D:\Dropbox\Astrophotography\Calibration\SX-46\September 2022\Darks default\BPM-sxASCOM_driver_Copyright_c_Dad_Dog_Development_Ltd.-4540x3640.fits"
# load the fits file
hdulist = fits.open(file)
# get the data
data = hdulist[0].data
vals = Counter(np.ndarray.flatten(data))
good_val = 127
bad_val = 255
print(f"BPM: {os.path.dirname(file)}")
print("Bad pixels: ", vals[bad_val])
print("Good pixels: ", vals[good_val])
print("Ratio: ", vals[bad_val]/vals[good_val])


folder = "D:\\Dropbox\\Astrophotography\\Calibration\\SX-46\\September 2022\\frames"
temps = []
# iterate over all files in the folder
for filename in os.listdir(folder):
    if filename.startswith("D_"):
        hdulist = fits.open(os.path.join(folder, filename))
        header = hdulist[0].header
        temp = header.get('CCD-TEMP')
        temps.append(temp)
        hdulist.close()
print("Mean temp: ", np.mean(temps))

folder = "D:\\Dropbox\\Astrophotography\\Calibration\\SX-46\\Old BPM\\August darks"
temps = []
# iterate over all files in the folder
for filename in os.listdir(folder):
    if filename.startswith("D_"):
        hdulist = fits.open(os.path.join(folder, filename))
        header = hdulist[0].header
        temp = header.get('CCD-TEMP')
        temps.append(temp)
        hdulist.close()
print("Mean temp: ", np.mean(temps))
