import photutils
from photutils.segmentation import detect_threshold
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils.segmentation import detect_sources
from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use('module://backend_interagg')
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as opt
import os
import glob
import pandas as pd
import multiprocessing
from multiprocessing import get_context


def segment_image(fits_file_name, nsigma=3., FWHM=5, kernel_size=5, npixels=7, padding=20, contrast=0.001, take_first=False, plot=False):
    """
    This function will find and cutout the stars from a FITS file
    :param fits_file_name: the file name for a fits file
    :param nsigma: how far above the backround a star needs to be in order to be detected
    :param FWHM: a rough estimate of the FWHM of the stars in the image, in order to smooth things
    :param kernel_size: how big the smoothing kernel should be
    :param npixels: how many contiguous pixels are required in order to detect a star
    :param padding: only relevant if plotting in order to show cutout and surrounding area
    :param contrast: something to do with the contrast between the star and the background
    :param take_first: int, if you want to take the first NxN pixels of an image for speed
    :param plot: bool, if you want to plot the star cutouts
    :return: a list of arrays of ints that represent the star cutouts, the astropy source_catalog that accompanies those
    """
    # open the file
    hdul = fits.open(fits_file_name)
    data = hdul[0].data
    # subset for speed if requested
    if take_first:
        data = data[0:take_first, 0:take_first]
    # set detection above background
    threshold = detect_threshold(data, nsigma=nsigma)
    # segment the image using a smoothing kernel
    sigma = FWHM * gaussian_fwhm_to_sigma
    # smoothing kernel
    kernel = Gaussian2DKernel(sigma, x_size=kernel_size, y_size=kernel_size)
    kernel.normalize()
    # do the segmentation
    segm = detect_sources(data, threshold, npixels=npixels, kernel=kernel)
    # deblend sources
    segm_deblend = photutils.segmentation.deblend_sources(data, segm, npixels=npixels, kernel=kernel, nlevels=32,
                                                          contrast=contrast)
    # for the fun of it, even though I don't trust the FWHM calculations
    source_catalog = photutils.segmentation.SourceCatalog(data, segm_deblend)
    num_stars = len(segm_deblend.segments)
    star_nums = range(num_stars)
    cutouts = list()
    for star_it in range(num_stars):
        star_ind = star_nums[star_it]
        star_slice = segm_deblend.segments[star_ind].slices
        # plot if that was asked for
        if plot:
            fig, axes = plt.subplots(2, 1, figsize=(10 * num_stars, 12.5))
            sns.heatmap(data[star_slice], ax=axes[0])
            axes[0, star_it].set_title(f'Extracted star\n measured FWHM={source_catalog[star_ind].fwhm.value}')
            x_start = np.max([0, star_slice[0].start - padding])
            x_stop = np.min([data.shape[0], star_slice[0].start + padding])
            y_start = np.max([0, star_slice[1].start - padding])
            y_stop = np.min([data.shape[1], star_slice[1].start + padding])
            # ax2.imshow(data[x_start:x_stop,y_start:y_stop], origin='lower', cmap='Greys_r', norm=norm)
            sns.heatmap(data[x_start:x_stop, y_start:y_stop], ax=axes[1])
            axes[1, star_it].set_title('Surrounding area')
        # Return the data itself
        cutouts.append(data[star_slice])
    return cutouts, source_catalog


def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """
    A two dimensional normal distribution
    :param xy: [1D array, 1D array] of x-y coordinates
    :param amplitude: amplitude of the distribution
    :param xo: mean x
    :param yo: mean y
    :param sigma_x: x std dev
    :param sigma_y: y std dev
    :param theta: rotation
    :param offset: height offset
    :return: 1D unravelled plot of the gaussian
    """
    x = xy[0]
    y = xy[1]
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()


def get_stats_cutout(cutout, plot=False):
    """
    Fits a 2D Gaussian to a cutout and returns stats about them
    :param cutout: a 2D array that represents a cutout of a star
    :param plot: if you want to plot the fit to check if it makes sense
    :return: max_FWHM (taking each dimension idependently), roundness (ratio of sigma_x/sigma_y), fit (all parameters that can be passed to the twoD_Gaussian function)
    """
    xp, yp = cutout.shape
    x, y = np.mgrid[:xp, :yp]
    # Initial guess
    amplitude = np.max(cutout)
    xo = np.floor(cutout.shape[0] / 2)
    yo = np.floor(cutout.shape[1] / 2)
    sigma_x = 2
    sigma_y = 2
    theta = 0
    offset = 0
    initial_guess = (amplitude, xo, yo, sigma_x, sigma_y, theta, offset)
    # Do the fitting
    try:
        popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), cutout.ravel(), p0=initial_guess, maxfev=1000)
        data_fitted = twoD_Gaussian((x, y), *popt)
    except RuntimeError:
        popt = (np.NAN, np.NAN, np.NAN, np.NAN, np.NAN, np.NAN, np.NAN)
        data_fitted = np.NAN
    if plot:
        plt.figure()
        fig, axes = plt.subplots(1, 3, figsize=(30, 12.5))
        sns.heatmap(cutout, ax=axes[0])
        plt.title("Data")
        f = twoD_Gaussian((x, y), *popt).reshape(xp, yp)
        sns.heatmap(f, ax=axes[1])
        plt.title("Model")
        sns.heatmap(np.abs(cutout - f), ax=axes[2])
        plt.title("Residual")
        plt.show()
    fit_amp = popt[0]
    fit_xo = popt[1]
    fit_yo = popt[2]
    fit_sigma_x = popt[3]
    fit_sigma_y = popt[4]
    fit_theta = popt[5]
    fit_offset = popt[6]
    fit = (fit_amp, fit_xo, fit_yo, fit_sigma_x, fit_sigma_y, fit_theta, fit_offset)
    # max FWHM
    max_FWHM = np.max(np.abs([2 * np.sqrt(2 * np.log(2)) * fit_sigma_x, 2 * np.sqrt(2 * np.log(2)) * fit_sigma_y]))
    # roundness
    try:
        roundness = fit_sigma_x / fit_sigma_y
    except:
        roundness = 0
    return max_FWHM, roundness, fit


def get_star_stats_fits(fits_file_name, nsigma=3., FWHM=5, kernel_size=5, npixels=7, contrast=0.001, take_first=False,
                        plot_stats=False):
    """
    Get FWHM and other data from a FITS file
    :param nsigma: how far above the backround a star needs to be in order to be detected
    :param FWHM: a rough estimate of the FWHM of the stars in the image, in order to smooth things
    :param kernel_size: how big the smoothing kernel should be
    :param npixels: how many contiguous pixels are required in order to detect a star
    :param padding: only relevant if plotting in order to show cutout and surrounding area
    :param contrast: something to do with the contrast between the star and the background
    :param take_first: int, if you want to take the first NxN pixels of an image for speed
    :param plot_stats: bool, if you want to plot the statistics (histogram, over all stars)
    :return: dictionary with keys: file_name, my_FWHM, my_roundness, their_FWHM,
    their_elongation, num_stars, filter, temp, focus_pos
    """
    # get the cutouts
    cutouts, source_catalog = segment_image(fits_file_name, nsigma=nsigma, FWHM=FWHM, kernel_size=kernel_size,
                                            npixels=npixels, contrast=contrast, take_first=take_first)
    # gather the stats for each cutout
    my_stats = list()
    for cutout in cutouts:
        my_stats.append(get_stats_cutout(cutout))
    if plot_stats:
        fig, axes = plt.subplots(1, 4, figsize=(40, 12.5))
        my_fwhm = list(map(lambda x: x[0] if ~np.isnan(x[0]) else 0, my_stats))
        #my_fwhm = my_fwhm[~np.isnan(my_fwhm)]
        sns.histplot(my_fwhm, ax=axes[0])
        axes[0].set_xlim(0, 5)
        axes[0].set_title('My FWHM')
        my_roundnes = list(map(lambda x: x[1] if ~np.isnan(x[1]) else 0, my_stats))
        #my_roundnes = my_roundnes[~np.isnan(my_roundnes)]
        sns.histplot(my_roundnes, ax=axes[1])
        axes[1].set_xlim(0, 5)
        axes[1].set_title('My roundness')
        sns.histplot(source_catalog.fwhm, ax=axes[2])
        axes[2].set_xlim(0, 5)
        axes[2].set_title('Their FWHM')
        sns.histplot(source_catalog.elongation, ax=axes[3])
        axes[3].set_title('Their elongation')
        plt.show()
    to_return = dict()
    to_return['my_FWHM'] = np.median(list(map(lambda x: x[0] if ~np.isnan(x[0]) else 0, my_stats)))
    to_return['my_roundness'] = np.median(list(map(lambda x: x[1] if ~np.isnan(x[1]) else 0, my_stats)))
    to_return['their_FWHM'] = np.median(source_catalog.fwhm).value
    to_return['their_elongation'] = np.median(source_catalog.elongation).value
    to_return['num_stars'] = len(cutouts)
    to_return['file_name'] = fits_file_name
    to_return['filter'] = fits.getval(fits_file_name, 'filter')
    to_return['temp'] = fits.getval(fits_file_name, 'amb-temp')
    to_return['focus_pos'] = fits.getval(fits_file_name, 'focuspos')
    to_return['object'] = fits.getval(fits_file_name, 'object')
    return to_return#, my_stats


def add_stats_to_fits(fits_file_name, take_first=False):
    """
    Writes the stats to the FITS header
    :param fits_file_name: FITS file name
    :param take_first: int N/bool, if you want to only analyze the first NXN pixels
    :return: none
    """
    res = get_star_stats_fits(fits_file_name, take_first=take_first, plot_stats=False)
    fits.setval(fits_file_name, 'my_FWHM', value=res["my_FWHM"])
    fits.setval(fits_file_name, 'my_roundness', value=res["my_roundness"])
    fits.setval(fits_file_name, 'their_FWHM', value=res["their_FWHM"])
    fits.setval(fits_file_name, 'their_elongation', value=res["their_elongation"])


def map_func(file_name, take_first=False):
    """
    Helper function for parallel processing
    :param file_name: FITS file name
    :param take_first: int/bool if you want to analyze only the first NXN pixels of the image
    :return: same as get_star_stats_fits
    """
    return get_star_stats_fits(file_name, take_first=take_first)


def get_stats_dir(dir_name, recursive=False, num_threads=10):
    """
    Compute a data frame with all the focus data from FITS files in a directory
    :param dir_name: a directory containing FITS files
    :param recursive: if you want to get all sub-folders too
    :param num_threads: number of threads to use, default=10
    :return: pandas data frame
    """
    #file_names = glob.glob("*.fit", root_dir=dir_name, recursive=recursive)
    # glob somehow got downgraded?
    file_names = glob.glob(os.path.join(dir_name, "*.fit"))
    file_names = list(map(lambda x: os.path.join(dir_name, x), file_names))
    # M1 mac requires forking
    pool = get_context("fork").Pool(num_threads)
    res = list(pool.map(map_func, file_names))
    pool.close()
    # turn it into a data-frame
    df = pd.DataFrame.from_records(res, index='file_name', columns=res[0].keys())
    return df

##################################################################
# tests
def test_segment_image():
    fits_file_name = "/Users/dmk333/Dropbox/Astrophotography/FocusTest/L_Ha_2509_300s.fit"
    cutouts, source_catalog = segment_image(fits_file_name)
    sns.heatmap(cutouts[0])
    plt.show()


def test_get_stats():
    fits_file_name = "/Users/dmk333/Dropbox/Astrophotography/FocusTest/L_Ha_2509_300s.fit"
    cutouts, source_catalog = segment_image(fits_file_name)
    my_stats = get_stats_cutout(cutouts[0])
    #print(my_stats)
    get_stats_cutout(cutouts[0], plot=True)


def test_calculate_star_stats():
    fits_file_name = "/Users/dmk333/Dropbox/Astrophotography/FocusTest/L_Ha_2509_300s.fit"
    res = get_star_stats_fits(fits_file_name, take_first=500, plot_stats=False)
    print(res)
    res = get_star_stats_fits(fits_file_name, take_first=500, plot_stats=True)


def test_add_stats_to_fits():
    fits_file_name = "/Users/dmk333/Dropbox/Astrophotography/FocusTest/L_Ha_2509_300s.fit"
    add_stats_to_fits(fits_file_name, take_first=500)

#test_get_stats()
#test_calculate_star_stats()
#fits_file_name = "/Users/dmk333/Dropbox/Astrophotography/FocusTest/L_Ha_2509_300s.fit"
#res = calculate_star_stats(fits_file_name, plot_stats=True)
#test_add_stats_to_fits()

