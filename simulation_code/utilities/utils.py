"""
Utility functions for the optical simulation package
"""
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize as opt
from scipy import ndimage
from scipy.optimize import curve_fit


def cart2pol(y, x):
    """
    utility function for converting from cartesian to polar
    Parameters: y, x;
    Returns: rho, theta; theta = [-pi, pi]; unit: radian;
    """
    theta = np.arctan2(y, x)
    rho = np.hypot(y, x)
    return rho, theta


def pol2cart(rho, theta):
    """
    utility function for converting from polar to cartesian
    Parameters: rho, theta; theta = [-pi, pi]; unit: radian;
    Returns: y, x
    """
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return y, x


def sph2cart(rho, theta, phi):
    """
    utility function for converting from spherical to cartesian
    Parameters: rho, theta, phi; theta = [0, pi]; phi = [-pi, pi]
    Returns: z, y, x
    """
    x = rho * np.sin(theta) * np.cos(phi)
    y = rho * np.sin(theta) * np.sin(phi)
    z = rho * np.cos(theta)
    return z, y, x


def cart2sph(z, y, x):
    """
    utility function for converting from cartesian to spherical
    Parameters: z, y, x;
    Returns: rho, theta, phi; theta = [0, pi]; phi = [-pi, pi]
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    phi = np.arctan2(y, x)
    theta = np.arccos(z/r)
    return r, theta, phi


def rotation(x, y, z, m):
    """3D rotation of xyz by matrix m, used by mapping_pupil"""
    x1 = x * m[0][0] + y * m[0][1] + z * m[0][2]
    y1 = x * m[1][0] + y * m[1][1] + z * m[1][2]
    z1 = x * m[2][0] + y * m[2][1] + z * m[2][2]
    return x1, y1, z1


def psqrt(data):
    """Take the positive square root, negative values will be set to zero."""
    # make zero array
    sdata = np.zeros_like(data, float)
    # fill only sqrt of positive values
    sdata[data > 0] = np.sqrt(data[data > 0])
    return sdata


def cross_section_3dpsf_2dfit(data, path=None, size=None, res=None, zsize=None, zres=None, plt_profiles=False, title=None, iso_voxel=True, figsize=(6, 4)):
    """show the cross section of a 3D psf"""
    if data.ndim != 3:
        raise ValueError("PSF is not 3D!")
    if res is None:
        res = 1
    if zres is None:
        zres = res

    # interpolate the data to have isotropic voxel
    if iso_voxel:
        data = ndimage.zoom(data, (zres/res, 1, 1))
        zres = res

    data = center_data(data)
    z_len, y_len, x_len = data.shape
    z_center, y_center, x_center = z_len // 2, y_len // 2, x_len // 2

    if size is None:
        size = min(y_len, x_len)
    if zsize is None:
        zsize = z_len

    z_start, z_end = z_center - zsize // 2, z_center + zsize // 2
    y_start, y_end = y_center - size // 2, y_center + size // 2
    x_start, x_end = x_center - size // 2, x_center + size // 2

    figure, axes = plt.subplots(nrows=2, ncols=2, figsize=figsize)
    if title is not None:
        plt.suptitle(title)
    axes[0, 0].imshow(data[z_center, y_start:y_end, x_start:x_end])
    axes[0, 1].imshow(np.rot90(data[z_start:z_end, y_start:y_end, x_center], 1, axes=(1, 0)))  # rotation needed for yz section
    axes[1, 0].imshow(data[z_start:z_end, y_center, x_start:x_end])

    # display a white image in axes[1,1] in the corner
    x_size = data.shape[2]
    z_size = data.shape[0]
    data_empty = np.zeros((z_size, x_size), dtype=np.uint8)
    data_empty.fill(255)
    axes[1, 1].imshow(data_empty, vmin=0, vmax=255, cmap='gray')
    axes[1, 1].axis('off')

    if plt_profiles:
        # get the FWHMs and angles, and display
        p1 = fitgaussian(data[z_center, :, :])  # p[3] = y_width, p[4] = x_width
        print("xy section", p1)
        p2 = fitgaussian(data[:, :, x_center])  # p[3] = z_width, p[4] = y_width
        print("yz section", p2)
        p3 = fitgaussian(data[:, y_center, :])  # p[3] = z_width, p[4] = x_width
        print("xz section", p3)

        # axes[1, 1].axis('off')
        FWHM_x = min(p1[4], p3[4]) * 2 * np.sqrt(2 * np.log(2)) * res
        FWHM_y = min(p1[3], p2[4]) * 2 * np.sqrt(2 * np.log(2)) * res
        FWHM_z = max(p2[3], p3[3]) * 2 * np.sqrt(2 * np.log(2)) * zres
        angle = p3[5]
        axes[1, 1].text(0, 0, "FWHM_x' = %.1f nm" % FWHM_x)
        axes[1, 1].text(0, z_size // 4, "FWHM_y' = %.1f nm" % FWHM_y)
        axes[1, 1].text(0, 2 * z_size // 4, "FWHM_z' = %.1f nm" % FWHM_z)
        axes[1, 1].text(0, 3 * z_size // 4, "Angle = %.1f degree" % angle)
        figure.tight_layout()

    if path is not None:
        plt.savefig(path + '/cross_sections.pdf', dpi=300)


def fit_profiles_3d(data, pixel_xy, pixel_z, path):
    ind = np.unravel_index(np.argmax(data, axis=None), data.shape)
    z = ind[0]
    y = ind[1]
    x = ind[2]

    figure, axes = plt.subplots(nrows=2, ncols=2)
    profile_x = data[z, y, :]
    profile_y = data[z, :, x]
    # get the z profile by choosing the max value in each slice; works well for simulated data
    profile_z = np.amax(data, axis=0)

    axes[0, 0].imshow(data[z, :, :])
    axes[0, 1].imshow(data[:, :, x])
    axes[1, 0].imshow(data[:, y, :])
    axes[1, 1].axis('off')
    figure.tight_layout()
    plt.savefig(path + '/fitted_profiles.pdf', dpi=300)
    plt.savefig(path + '/fitted_profiles.png', dpi=300)


def gaussian(height, center_x, center_y, width_x, width_y, rot_by_angle):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)

    rot_by_angle = np.deg2rad(rot_by_angle)
    center_x = center_x * np.cos(rot_by_angle) - center_y * np.sin(rot_by_angle)
    center_y = center_x * np.sin(rot_by_angle) + center_y * np.cos(rot_by_angle)

    def rotgauss(x, y):
        xp = x * np.cos(rot_by_angle) - y * np.sin(rot_by_angle)
        yp = x * np.sin(rot_by_angle) + y * np.cos(rot_by_angle)
        g = height * np.exp(
            -(((center_x - xp) / width_x) ** 2 +
              ((center_y - yp) / width_y) ** 2) / 2.)
        return g

    return rotgauss


def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X * data).sum() / total
    y = (Y * data).sum() / total
    col = data[:, int(y)]
    width_x = np.sqrt(abs((np.arange(col.size) - y) ** 2 * col).sum() / col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(abs((np.arange(row.size) - x) ** 2 * row).sum() / row.sum())
    height = data.max()
    return height, x, y, width_x, width_y, 0.0


def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)  *  x and y seem to be reversed for images
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    error_function = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
    p, success = opt.leastsq(error_function, params)
    return p


def center_data(data):
    """Utility to center the data

    Parameters
    ----------
    data : ndarray
        Array of data points

    Returns
    -------
    centered_data : ndarray same shape as data
        data with max value at the central location of the array
    """
    # copy data
    centered_data = data.copy()
    # extract shape and max location
    data_shape = data.shape
    max_loc = np.unravel_index(data.argmax(), data_shape)
    # iterate through dimensions and roll data to the right place
    for i, (x0, nx) in enumerate(zip(max_loc, data_shape)):
        centered_data = np.roll(centered_data, nx // 2 - x0, i)
    return centered_data


# def cross_section_xz(data, path=None, size=None, res=None, zsize=None, zres=None, plt_profiles=False):
#     """show the cross section of a 3D psf"""
#     if data.ndim != 3:
#         raise ValueError("PSF is not 3D!")
#     if res is None:
#         res = 1
#     if zres is None:
#         zres = res
#
#     data = center_data(data)
#     z_len, y_len, x_len = data.shape
#     z_center, y_center, x_center = z_len // 2, y_len // 2, x_len // 2
#
#     if size is None:
#         size = min(y_len, x_len)
#     if zsize is None:
#         zsize = z_len
#
#     z_start, z_end = z_center - zsize // 2, z_center + zsize // 2
#     y_start, y_end = y_center - size // 2, y_center + size // 2
#     x_start, x_end = x_center - size // 2, x_center + size // 2
#
#     plt.figure()
#     plt.imshow(data[z_start:z_end, y_center, x_start:x_end])
#
#     # todo, check for correctness
#     if plt_profiles is True:
#         # get the FWHMs and angles, and display
#         p1 = fitgaussian(data[z_center, :, :])  # p[3] = y_width, p[4] = x_width
#         print("xy section", p1)
#         p2 = fitgaussian(data[:, :, x_center])  # p[3] = z_width, p[4] = y_width
#         print("yz section", p2)
#         p3 = fitgaussian(data[:, y_center, :])  # p[3] = z_width, p[4] = x_width
#         print("xz section", p3)
#
#         # axes[1, 1].axis('off')
#         FWHM_x = min(p1[4], p3[4]) * 2 * np.sqrt(2 * np.log(2)) * res * 1000
#         FWHM_y = min(p1[3], p2[4]) * 2 * np.sqrt(2 * np.log(2)) * res * 1000
#         FWHM_z = max(p2[3], p3[3]) * 2 * np.sqrt(2 * np.log(2)) * zres * 1000
#         angle = p3[5]
#         axes[1, 1].text(0, 0, "FWHM_x' = %.1f nm" % FWHM_x)
#         axes[1, 1].text(0, z_size // 4, "FWHM_y' = %.1f nm" % FWHM_y)
#         axes[1, 1].text(0, 2 * z_size // 4, "FWHM_z' = %.1f nm" % FWHM_z)
#         axes[1, 1].text(0, 3 * z_size // 4, "Angle = %.1f degree" % angle)
#         figure.tight_layout()
#
#     if path is not None:
#         plt.savefig(path + '/cross_section_xz.pdf', dpi=300)


def normalize(data):
    datamax = np.amax(data.flatten())
    return data / datamax


# def cross_section_3d(data, pixelsize_xy, pixelsize_z, title=None):
#     """get cross sections xy, yz, xz from a 3D stack"""
#     # for cross section, show scaled data
#     data_scaled = axial_rescale(data, pixelsize_z, pixelsize_xy)
#     ind = np.unravel_index(np.argmax(data_scaled, axis=None), data_scaled.shape)
#     z = ind[0]
#     y = ind[1]
#     x = ind[2]
#     figure, axes = plt.subplots(nrows=2, ncols=2, figsize=(6, 6))
#
#     axes[0, 0].imshow(data_scaled[z, :, :])
#     axes[0, 1].imshow(np.rot90(data_scaled[:, :, x], -1, axes=(1, 0)))  # rotation needed for yz section
#     axes[1, 0].imshow(data_scaled[:, y, :])
#
#     # for fitting, use original data
#     z_size = data.shape[0]
#     xy_size = data.shape[1]
#     plot_and_fit_gaussian_1d(data, pixelsize_xy, pixelsize_z, xy_size - 10, z_size - 10)
#
#     if title is not None:
#         plt.title(title)
#
#     figure.tight_layout()


def cross_section_3d_1dfit(data, path=None, size=None, res=None, zsize=None, zres=None, plt_profiles=False, title=None, iso_voxel=True, figsize=(6, 4)):
    """show the cross section of a 3D psf and fit the plots with 1d fit"""
    if data.ndim != 3:
        raise ValueError("PSF is not 3D!")
    if res is None:
        res = 1
    if zres is None:
        zres = res

    # interpolate the data to have isotropic voxel
    if iso_voxel:
        data = ndimage.zoom(data, (zres/res, 1, 1))
        zres = res

    data = center_data(data)
    z_len, y_len, x_len = data.shape
    z_center, y_center, x_center = z_len // 2, y_len // 2, x_len // 2

    if size is None:
        size = min(y_len, x_len)
    if zsize is None:
        zsize = z_len

    z_start, z_end = z_center - zsize // 2, z_center + zsize // 2
    y_start, y_end = y_center - size // 2, y_center + size // 2
    x_start, x_end = x_center - size // 2, x_center + size // 2

    figure, axes = plt.subplots(nrows=2, ncols=2, figsize=figsize)
    if title is not None:
        plt.suptitle(title)
    axes[0, 0].imshow(data[z_center, y_start:y_end, x_start:x_end])
    axes[0, 1].imshow(np.rot90(data[z_start:z_end, y_start:y_end, x_center], 1, axes=(1, 0)))  # rotation needed for yz section
    axes[1, 0].imshow(data[z_start:z_end, y_center, x_start:x_end])

    if plt_profiles:
        plot_and_fit_gaussian_1d(data, res, zres, size, zsize)

    if path is not None:
        plt.savefig(path + '/cross_sections.pdf', dpi=300)



def axial_rescale(imArray, oldz, newz):
    """scale the 3D data along z to have isotropic voxel"""
    imArrayRescaled = ndimage.zoom(imArray, (oldz/newz, 1, 1))
    # imsave(inpath.replace('.tif','') + '_rescaled.tif', imArrayRescaled)
    return imArrayRescaled


def plot_and_fit_gaussian_1d(imArray, pixelsize_xy, pixelsize_z, xy_length, z_length):
    """fit 1d gaussian profile of a 3D stack"""
    # find out the coordinates of the bead
    ind = np.unravel_index(np.argmax(imArray, axis=None), imArray.shape)
    zCenter = ind[0]
    yCenter = ind[1]
    xCenter = ind[2]

    # get the normalized profiles
    colors = ['crimson', 'blue', 'green', 'lightblue', 'magenta', 'red', 'black']
    xLine = imArray[zCenter, yCenter, :] / imArray[zCenter, yCenter, :].max()
    yLine = imArray[zCenter, :, xCenter] / imArray[zCenter, :, xCenter].max()
    zLine = imArray[:, yCenter, xCenter] / imArray[:, yCenter, xCenter].max()

    # recenter the profiles
    xLine = recenter(xLine)
    yLine = recenter(yLine)
    zLine = recenter(zLine)

    # fit x profile
    xCoord = np.linspace(-xLine.shape[0] / 2, xLine.shape[0] / 2, xLine.shape[0], endpoint=True) * pixelsize_xy
    lengthToShow = range(imArray.shape[2] // 2 - xy_length // 2, imArray.shape[2] // 2 + xy_length // 2)
    ampl = xLine.max() - xLine.min()
    mean = xCoord[xLine.argmax()]  # note this correction
    background = xLine.min()
    sigma = 200  # note this correction
    popt, pcov = curve_fit(gaus, xCoord, xLine, p0=[ampl, mean, sigma, background])
    FWHM = 2 * np.sqrt(2 * np.log(2)) * popt[2]

    plt.plot(xCoord[lengthToShow], xLine[lengthToShow], linestyle='None', marker='.', markersize=5, color=colors[0])
    plt.plot(xCoord[lengthToShow], gaus(xCoord, popt[0], popt[1], popt[2], popt[3])[lengthToShow],
             label="FWHM_x = %.1f nm" % FWHM, color=colors[0])

    # fit y profile
    yCoord = np.linspace(-yLine.shape[0] / 2, yLine.shape[0] / 2, yLine.shape[0], endpoint=True) * pixelsize_xy
    lengthToShow = range(imArray.shape[1] // 2 - xy_length // 2, imArray.shape[1] // 2 + xy_length // 2)
    ampl = yLine.max() - yLine.min()
    mean = yCoord[yLine.argmax()]  # note this correction
    background = yLine.min()
    sigma = 200  # note this correction
    popt, pcov = curve_fit(gaus, yCoord, yLine, p0=[ampl, mean, sigma, background])
    FWHM = 2 * np.sqrt(2 * np.log(2)) * popt[2]

    plt.plot(yCoord[lengthToShow], yLine[lengthToShow], linestyle='None', marker='.', markersize=5, color=colors[1])
    plt.plot(yCoord[lengthToShow], gaus(yCoord, popt[0], popt[1], popt[2], popt[3])[lengthToShow],
             label="FWHM_y = %.1f nm" % FWHM, color=colors[1])

    # fit z profile
    lengthToShow = range(imArray.shape[0] // 2 - z_length // 2, imArray.shape[0] // 2 + z_length // 2)
    zCoord = np.linspace(-zLine.shape[0] / 2, zLine.shape[0] / 2, zLine.shape[0], endpoint=True) * pixelsize_z
    ampl = zLine.max() - zLine.min()
    mean = zCoord[zLine.argmax()]                   #note this correction
    background = zLine.min()
    sigma = 200        #note this correction
    popt, pcov = curve_fit(gaus, zCoord, zLine, p0=[ampl, mean, sigma, background])
    FWHM = 2*np.sqrt(2*np.log(2))*popt[2]
    plt.plot(zCoord[lengthToShow], zLine[lengthToShow], linestyle = 'None', marker='.', markersize = 5, color = colors[2])
    plt.plot(zCoord[lengthToShow], gaus(zCoord, popt[0], popt[1], popt[2], popt[3])[lengthToShow],
             label="FWHM_z = %.1f nm" % FWHM, color = colors[2])

    plt.legend(loc='upper left')
    plt.xlabel('Distance (nm)')
    plt.ylabel('Intensity (a.u.)')
    # plt.show()


# recenter a line to its peak intensity
def recenter(data):
    return np.roll(data, len(data) // 2 - np.argmax(data))


# 1D-gauss function
def gaus(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b
