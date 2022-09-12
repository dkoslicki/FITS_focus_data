# Investigate if the number of warm/hot pixels is increaseing

from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use('module://backend_interagg')
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# import for q q plot
import scipy.stats as stats
import numbers


def qqplot(x, y, quantiles=None, interpolation='nearest', ax=None, rug=False,
           rug_length=0.05, rug_kwargs=None, **kwargs):
    """Draw a quantile-quantile plot for `x` versus `y`.

    Parameters
    ----------
    x, y : array-like
        One-dimensional numeric arrays.

    ax : matplotlib.axes.Axes, optional
        Axes on which to plot. If not provided, the current axes will be used.

    quantiles : int or array-like, optional
        Quantiles to include in the plot. This can be an array of quantiles, in
        which case only the specified quantiles of `x` and `y` will be plotted.
        If this is an int `n`, then the quantiles will be `n` evenly spaced
        points between 0 and 1. If this is None, then `min(len(x), len(y))`
        evenly spaced quantiles between 0 and 1 will be computed.

    interpolation : {‘linear’, ‘lower’, ‘higher’, ‘midpoint’, ‘nearest’}
        Specify the interpolation method used to find quantiles when `quantiles`
        is an int or None. See the documentation for numpy.quantile().

    rug : bool, optional
        If True, draw a rug plot representing both samples on the horizontal and
        vertical axes. If False, no rug plot is drawn.

    rug_length : float in [0, 1], optional
        Specifies the length of the rug plot lines as a fraction of the total
        vertical or horizontal length.

    rug_kwargs : dict of keyword arguments
        Keyword arguments to pass to matplotlib.axes.Axes.axvline() and
        matplotlib.axes.Axes.axhline() when drawing rug plots.

    kwargs : dict of keyword arguments
        Keyword arguments to pass to matplotlib.axes.Axes.scatter() when drawing
        the q-q plot.
    """
    # Get current axes if none are provided
    if ax is None:
        ax = plt.gca()

    if quantiles is None:
        quantiles = min(len(x), len(y))

    # Compute quantiles of the two samples
    if isinstance(quantiles, numbers.Integral):
        quantiles = np.linspace(start=0, stop=1, num=int(quantiles))
    else:
        quantiles = np.atleast_1d(np.sort(quantiles))
    x_quantiles = np.quantile(x, quantiles, interpolation=interpolation)
    y_quantiles = np.quantile(y, quantiles, interpolation=interpolation)

    # Draw the rug plots if requested
    if rug:
        # Default rug plot settings
        rug_x_params = dict(ymin=0, ymax=rug_length, c='gray', alpha=0.5)
        rug_y_params = dict(xmin=0, xmax=rug_length, c='gray', alpha=0.5)

        # Override default setting by any user-specified settings
        if rug_kwargs is not None:
            rug_x_params.update(rug_kwargs)
            rug_y_params.update(rug_kwargs)

        # Draw the rug plots
        for point in x:
            ax.axvline(point, **rug_x_params)
        for point in y:
            ax.axhline(point, **rug_y_params)

    # Draw the q-q plot
    ax.scatter(x_quantiles, y_quantiles, **kwargs)



# two random darks
#files = ["/Users/dmk333/Dropbox/Astrophotography/Calibration/SX-46/Old BPM/June darks/D_11274_600s.fit", "/Users/dmk333/Dropbox/Astrophotography/Calibration/SX-46/September 2022/frames/D_17929_300s.fit"]
# two master darks
files = ["/Users/dmk333/Dropbox/Astrophotography/Calibration/SX-46/Mid June 2022/MD-IG_-1.0-E600.0s-sxASCOM_driver_Copyright_c_Dad_Dog_Development_Ltd.-4540x3640-all_channels-session_1.fits", "/Users/dmk333/Dropbox/Astrophotography/Calibration/SX-46/September 2022/Darks default/MD-IG_-1.0-E600.0s-sxASCOM_driver_Copyright_c_Dad_Dog_Development_Ltd.-4540x3640-all_channels-session_1.fits"]
dark_vals = []
for file in files:
    # Import the fits file
    hdul = fits.open(file)
    # Get the data
    data = hdul[0].data
    # flatten the array
    data = data.flatten()
    # plot seaborn histogram
    #data_large = data[data > 2*np.median(data)]
    dark_vals.append(data)
    data_large = data
    sns.distplot(data_large, hist=True, kde=False,
                    bins=int(180/5), color = 'darkblue',
                    hist_kws={'edgecolor':'black'},
                    kde_kws={'linewidth': 4})
    # take log of y axis
    plt.yscale('log')
    # Add labels
    plt.title('Histogram of Pixel Values')
    plt.xlabel('Pixel Value')
    plt.ylabel('Count')
    # Show the plot
    plt.show()
    print("Number of pixels above 2*median: ", len(data_large))
    print("median: ", np.median(data))
    print("Skew: ", pd.Series(data_large).skew())
    # number of pixels above threshold
    print("Number of pixels above threshold: ", len(data[data > 0.5]))

plt.figure()
# randomly sample the data
data1_subsample = dark_vals[0][np.random.randint(0, len(dark_vals[0]), 10000)]
data2_subsample = dark_vals[1][np.random.randint(0, len(dark_vals[1]), 10000)]
qqplot(data1_subsample, data2_subsample, quantiles=[.10,.20,.30,.40,.50,.60,.70,.80,.90,1], interpolation='nearest', ax=None, rug=True,
           rug_length=0.05)
plt.xlabel('June')
plt.ylabel('September')
plt.show()


# scatter plot of two dark frames
plt.figure()
plt.scatter(dark_vals[0], dark_vals[1], s=0.1)
plt.xlabel('June')
plt.ylabel('September')
plt.show()

# plot of difference between two dark frames
plt.figure()
plt.scatter(dark_vals[0], dark_vals[1] - dark_vals[0], s=0.1)
# add an y=0 line
plt.plot([0, 1], [0, 0], 'r-', lw=2)
plt.xlabel('June')
plt.ylabel('September - June')
plt.title('Master Dark Frame Comparison')
# save the plot
plt.savefig('dark_diff.png', dpi=300)
plt.show()





# Try adding histograms to the axes
def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    binwidth = 0.01
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    #bins = np.arange(-lim, lim + binwidth, binwidth)
    bins = np.arange(0, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    # scale x axis to log
    ax_histx.set_yscale('log')
    ax_histy.hist(y, bins=bins, orientation='horizontal')
    ax_histy.set_xscale('log')

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]
# start with a square Figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_axes(rect_scatter)
plt.plot([0, 1], [0, 0], 'r-', lw=2)
plt.xlabel('June')
plt.ylabel('September - June')
plt.title('Master Dark Frame Comparison')
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)
# use the previously defined function
scatter_hist(dark_vals[0], dark_vals[1] - dark_vals[0], ax, ax_histx, ax_histy)
plt.savefig('dark_diff_with_hist.png', dpi=300)
plt.show()

