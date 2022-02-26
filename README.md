# FITS_focus_data
Extract FWHM, temp, focus positions, etc. from a directory of astronomical FITS files

# Methods

## Astropy default
The keys that start with `their_` are computed using Astropy/Photutils: first an image is segmented to extract the stars, then `source_catalog.fwhm` and similar methods are used to get the stats. I haven't looked into how Astropy/Photutils does this, hence rolling my own solution.

## My own solution
First, I segment the image as above to extract the stars, but then I fit a 2D Gaussian distribution with `opt.curve_fit`. I then take the conservative/pessimistic approach of taking the maximum of the two axes' 1D FWHM.

# Notes
No optimization has been done save for parallelizing over the FITS files
