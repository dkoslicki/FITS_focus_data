import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import os
import seaborn as sns
plt.rcParams['figure.figsize'] = (10, 8)

# data based on single pixel
dir_name = "/Users/dmk333/Dropbox/Astrophotography/Headphones nebula/alignedNormalizedForKalman/"
file_names = glob.glob(os.path.join(dir_name, "*.fits"))
file_names = list(map(lambda x: os.path.join(dir_name, x), file_names))
data_list = []
for file_name in file_names:
    hdul = fits.open(file_name)
    data = np.array(hdul[0].data)
    data_list.append(data)

z = []
#i = 1739
i = 8
j = 343
for data in data_list:
    z.append(data[i, j])

x = np.median(z)
print(f"median value={x}")
# intial parameters
n_iter = len(z)
sz = (n_iter,) # size of array
#x = -0.37727 # truth value (typo in example at top of p. 13 calls this z)
#z = np.random.normal(x,0.1,size=sz) # observations (normal about x, sigma=0.1)

Q = 1e-5  # process variance

# allocate space for arrays
xhat = np.zeros(sz)       # a posteri estimate of x
P = np.zeros(sz)          # a posteri error estimate
xhatminus = np.zeros(sz)  # a priori estimate of x
Pminus = np.zeros(sz)     # a priori error estimate
K = np.zeros(sz)          # gain or blending factor

R = 0.1**2  # estimate of measurement variance, change to see effect
#R = 100*np.std(z)

# initial guesses
xhat[0] = 0.0
#xhat[0] = np.median(z)
#P[0] = 1.0
P[0] = .01

for k in range(1, n_iter):
    # time update
    xhatminus[k] = xhat[k-1]
    Pminus[k] = P[k-1]+Q

    # measurement update
    K[k] = Pminus[k]/(Pminus[k]+R)
    xhat[k] = xhatminus[k]+K[k]*(z[k]-xhatminus[k])
    P[k] = (1-K[k])*Pminus[k]

print(f"Median pixel={np.median(z)}")
print(f"Median kalman={np.median(xhat)}")
print(f"Last Kalman={xhat[-1]}")

plt.figure()
plt.plot(z, 'k+', label='noisy measurements')
plt.plot(xhat, 'b-', label='a posteri estimate')
plt.axhline(x, color='g', label='median value')
plt.legend()
plt.title('Estimate vs. iteration step', fontweight='bold')
plt.xlabel('Iteration')
plt.ylabel('Pixel value')
plt.show()


plt.figure()
valid_iter = range(1, n_iter)  # Pminus not valid at step 0
plt.plot(valid_iter, Pminus[valid_iter], label='a priori error estimate')
plt.title('Estimated $\it{\mathbf{a \ priori}}$ error vs. iteration step', fontweight='bold')
plt.xlabel('Iteration')
plt.ylabel('Pixel value(^2?)')
plt.setp(plt.gca(), 'ylim', [0, .01])
plt.show()

#############################################################################
# Now let's see if we can do it over the entire image

dir_name = "/Users/dmk333/Dropbox/Astrophotography/Headphones nebula/alignedNormalizedForKalman/"
file_names = glob.glob(os.path.join(dir_name, "*.fits"))
file_names = list(map(lambda x: os.path.join(dir_name, x), file_names))
data_list = []
ms = []
ns = []
# for subsampling
subsample = False
sub_m = 1000
sub_n = 1000
for file_name in file_names:
    hdul = fits.open(file_name)
    data = np.array(hdul[0].data)
    m, n = data.shape
    if subsample:
        m = sub_m
        n = sub_n
        data_list.append(data[:m, :n])
    else:
        data_list.append(data)
    ms.append(m)
    ns.append(n)


if len(np.unique(ms)) != 1 or len(np.unique(ms)) != 1:
    raise Exception("Dimensions of images don't match")
else:
    max_m = ms[0]
    max_n = ns[0]

z = np.array(data_list)
num_images = len(z)

####
# Let's see if I can remove outliers (i.e. hot pixels)
# Idea: identify hot pixel locations through those locations whose variances are higher than expected
z_no_hot = z
med = np.median(z, axis=0)
var = np.var(z, axis=0)
sigma = 3
bad_i, bad_j = np.where(var > sigma*np.mean(var))
# now for each of the locations that have a hot pixel show up, let's replace it with the median
for i, j in zip(bad_i, bad_j):
    series = z[:, i, j]
    series_med = np.median(series)
    series_var = np.var(series)
    upper_lim = series_med + sigma*series_var
    lower_lim = series_med - sigma*series_var
    spike_locs = np.where(np.logical_or(series > upper_lim, series < lower_lim))
    series[spike_locs] = series_med
    z_no_hot[:, i, j] = series

z = z_no_hot

# median for each pixel across the series of images
x = np.median(z, axis=1)
# intial parameters
n_iter = num_images
sz = z.shape  # size of array

#Q = 1e-5  # process variance
Q = 1e-3

# allocate space for arrays
xhat = np.zeros(sz)       # a posteri estimate of x
P = np.zeros(sz)          # a posteri error estimate
xhatminus = np.zeros(sz)  # a priori estimate of x
Pminus = np.zeros(sz)     # a priori error estimate
K = np.zeros(sz)          # gain or blending factor

R = 0.1**2  # estimate of measurement variance, change to see effect

# initial guesses
xhat[0] = np.median(z, axis=0)
#P[0] = 1.0
P[0] = .01

for k in range(1, n_iter):
    # time update
    xhatminus[k] = xhat[k-1]
    Pminus[k] = P[k-1]+Q

    # measurement update
    K[k] = Pminus[k]/(Pminus[k]+R)
    xhat[k] = xhatminus[k]+K[k]*(z[k]-xhatminus[k])
    P[k] = (1-K[k])*Pminus[k]

if False:
    fig, axes = plt.subplots(1, 3)
    k_median_estimate = np.median(xhat, axis=0)
    k_last_estimate = xhat[-1]
    median_estimate = np.median(z, axis=0)
    sns.heatmap(median_estimate, ax=axes[0])
    sns.heatmap(k_median_estimate, ax=axes[1])
    sns.heatmap(k_last_estimate, ax=axes[2])
    plt.show()

# Then save it as a new fits file
out_path = "/Users/dmk333/Dropbox/Astrophotography/Headphones nebula/KalmanOut/kalman_median.fit"
xhat_median = np.median(xhat, axis=0).astype('float32')
fits.writeto(out_path, xhat_median, overwrite=True)

out_path = "/Users/dmk333/Dropbox/Astrophotography/Headphones nebula/KalmanOut/kalman_last.fit"
xhat_last = xhat[-1].astype('float32')
fits.writeto(out_path, xhat_last, overwrite=True)

out_path = "/Users/dmk333/Dropbox/Astrophotography/Headphones nebula/KalmanOut/median.fit"
just_median = np.median(z, axis=0)
fits.writeto(out_path, just_median, overwrite=True)

print("finished")


########################################################
# Trying the approach of https://github.com/SRC877/Comparing-and-combining-Kalman-and-Wiener-Filter-for-video-denoising
# which used the code from
# % Purpose
# % Implements a predictive Kalman-like filter in the time domain of the image
# % stack. Algorithm taken from Java code by C.P. Mauer.
# % http://rsb.info.nih.gov/ij/plugins/kalman.html
percentvar = .005
gain = 0.9
imageStack = z
imageStack = np.concatenate((imageStack, [z[-1]]))
width = imageStack.shape[1]
height = imageStack.shape[2]
stacksize = imageStack.shape[0]
tmp = np.ones((width, height))
predicted = imageStack[1]
predictedvar = tmp * percentvar
noisevar = predictedvar

for i in range(1, stacksize-1):
    stackslice = imageStack[i + 1]
    observed = stackslice
    Kalman = predictedvar / (predictedvar + noisevar)
    corrected = gain * predicted + (1.0 - gain) * observed + Kalman * (observed - predicted)
    correctedvar = predictedvar * (tmp - Kalman)
    predictedvar = correctedvar
    predicted = corrected
    imageStack[i] = corrected

res = imageStack[1:-2]

out_path = "/Users/dmk333/Dropbox/Astrophotography/Headphones nebula/KalmanOut/kalman2.fit"
kalman2 = res[-1].astype('float32')
fits.writeto(out_path, kalman2, overwrite=True)
print('finished')



#########################################
# outlier hunting
x_start = 1743
y_start = 1303
pad = 10
for i in range(x_start-pad,x_start+pad):
    for j in range(y_start-pad, y_start+pad):
        sns.histplot(z[:,i,j])
        plt.title(f"({i,j})")
        plt.show()
        plt.pause(.1)

plt.imshow(np.max(z[:,x_start-pad:x_start+pad,y_start-pad:y_start+pad], axis=0)); plt.show()
plt.imshow(z[5,x_start-pad:x_start+pad,y_start-pad:y_start+pad]); plt.show()

i = 1733 #1750  # 1733
j = 1306 #1298  1294 1293 # 1306
plt.plot(z[:,i,j]); plt.show()