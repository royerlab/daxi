import numpy as np
from matplotlib import pyplot as plt
from utilities.utils import cart2pol, pol2cart, sph2cart, cart2sph, rotation
from skimage.transform import downscale_local_mean
from pyotf.zernike import zernike, noll2degrees, noll2name


class RemoteFocusing(object):
    """A class defining a remote focusing module with two objective at the remote space of the
    primary objectives O1. To avoid confusion, the two remote objectives are called O2 and O3,
    to be consistent with that usually used in the papers.
    """

    def __init__(self, wavelength=0.525, n2=1, na2=0.75, n3=1.33, na3=1.0, alpha=0, res=0.1,
                 size=1024, zres=None, zsize=80, binning=4, numZern=15, zern_mode=0, zern_coeff=0, is_zern_o2=True):
        """Generate an object for the remote focusing module based on the given parameters

        Input Parameters
        ----------
        wavelength : numeric
            wavelength of the emission light, unit: um
        n2 : numeric
            refractive index of O2, usually is 1 (air objective)
        na2 : numeric
            numerical aperture of O2
        n3 : numeric
            refractive index of O3, usually is 1.33 (water objective)
        na3 : numeric
            numerical aperture of O3
        alpha : numeric
            tilting angle between O2 and O3, unit: degree
        res : numeric
            pixel size in xy plane;
        size : numeric
            number of pixels in the xy plane, need to be large if angle is not 0 for sampling reason
        binning : numeric
            number of bin to scale down the pupil function to generate the psf, to avoid discontinuity
            in the pupul function due to mapping
        zres : numeric
            pixel size along the z axis;
        zsize : numeric
            number of pixels along the z axis
        numZern : numeric
            up to this number of zernike modes can be included
        zern_mode : numeric
            which zernike mode exist in the system
        zern_coeff : numeric
            coefficient of that zernike mode

        Auto-generated parameters
        ----------
        delta_k : numeric
            unit frequency of the pupil function, defined by 1 / size / res
        phase_o2 : ndarray (n, n)
            phase component of the pupil function before O2, unit: wave; multiple by 2pi to get radian
        mag_o2 : ndarray (n, n)
            magnitude component of the pupil function before O2
        phase_o3 : ndarray (n, n)
            phase component of the pupil function after O3, unit: wave; multiple by 2pi to get radian
        mag_o3 : ndarray (n, n)
            magnitude component of the pupil function after O3
        psf_o2 : ndarray (n, n)
            psf at focal space of o2
        psf_o3 : ndarray (n, n)
            "equivalent" psf at focal space of o3
        """
        self.n2 = n2
        self.na2 = na2
        self.n3 = n3
        self.na3 = na3
        self.alpha = np.deg2rad(alpha)
        self.wavelength = wavelength
        self.res = res
        if zres is None:
            zres = res
        self.zres = zres
        self.size = size
        self.delta_k = 1 / (size * res)
        self.zsize = zsize
        self.binning = binning
        self.numZern = numZern
        self.zern_mode = zern_mode
        self.zern_coeff = zern_coeff
        self.is_zern_o2 = is_zern_o2
        self.zerns = self._generate_zernikes()
        self.phase_o2 = self._calculate_o2_phase()
        self.mag_o2 = self._calculate_o2_mag()
        if self.is_zern_o2:
            if self.zern_coeff != 0:
                self.phase_o2 = self.zern_coeff * self.zerns[self.zern_mode-1]
            self.phase_o3 = self._calculate_o3_phase()
        else:
            self.phase_o3 = self._calculate_o3_phase()
            self.phase_o3 += self.zern_coeff * downscale_local_mean(self.zerns[self.zern_mode-1], (binning, binning))
        self.mag_o3 = self._calculate_o3_mag()
        self.psf3d = self._calculate_3d_psf_o3()
        # self.phase_o2_bin = downscale_local_mean(self.phase_o2, (self.binning, self.binning))
        # self.phase_o3_bin = downscale_local_mean(self.phase_o3, (self.binning, self.binning))
        # self.mag_o2_bin = downscale_local_mean(self.mag_o2, (self.binning, self.binning))
        # self.mag_o3_bin = downscale_local_mean(self.mag_o3, (self.binning, self.binning))

    def _calculate_o2_phase(self):
        """calculate the phase component of the pupil function before O2, without aberrations"""
        delta_k = self.delta_k  # unit frequency
        x = np.linspace(-delta_k * self.size / 2, delta_k * self.size / 2, self.size, endpoint=False)
        (X, Y) = np.meshgrid(x, x)
        rho, theta = cart2pol(Y, X)  # rho map of the pupil function
        k_max = self.na2 / self.wavelength
        return (rho < k_max).astype(float)  # convert from boolean to float

    def _calculate_o2_mag(self):
        """calculate the magnitude component of the pupil function before O2, without apodizations"""
        delta_k = self.delta_k  # unit frequency
        x = np.linspace(-delta_k * self.size / 2, delta_k * self.size / 2, self.size, endpoint=False)
        (X, Y) = np.meshgrid(x, x)
        rho, theta = cart2pol(Y, X)  # rho map of the pupil function
        k_max = self.na2 / self.wavelength
        return (rho < k_max).astype(float)  # convert from boolean to float

    def _calculate_o3_phase(self):
        """calculate the phase component of the pupil function after O3, without apodizations"""
        return self.mapping(self.phase_o2, is_mag=False)

    def _calculate_o3_mag(self):
        """calculate the magnitude component of the pupil function after O3, without apodizations"""
        return self.mapping(self.mag_o2, is_mag=True)

    def mapping(self, pupil, is_mag):
        """mapping the pupil function from O2 to O3, with or without an angle
         Note: only consider tilt between O2 and O2 along the x axis

        Parameters
        ----------
        pupil : ndarray (n, n)
            pupil function of O2, phase or magnitude
        is_mag : boolean
            is this to calculate magnitude or phase

        Returns
        ___________
        pupil2 : ndarray (n, n)
            the pupil function of O3
        """
        #######
        # convert pupil function from k space to spherical coordinates, rotate by alpha,
        # pass through interface n2-n3 and get back to k space
        delta_k = self.delta_k  # unit frequency
        k_min = -delta_k * self.size / 2
        k_max = delta_k * self.size / 2
        k = np.linspace(k_min, k_max, self.size, endpoint=False)
        (kx, ky) = np.meshgrid(k, k)
        kr, phi = cart2pol(ky, kx)   # from k space to spherical coordinates
        # get the theta angle from kr; min(1, VAL) for arcsin; using k = n / lambda * sin(theta)
        theta = np.arcsin(np.minimum(kr * self.wavelength / self.n2, 1))
        z, y, x = sph2cart(1, theta, phi)   # convert the spherical coords to cartesian coords
        # define the rotation matrix along y by an angle alpha
        rot = np.array([[np.cos(self.alpha), 0, np.sin(self.alpha)], [0, 1, 0], (-np.sin(self.alpha), 0, np.cos(self.alpha))])
        # rotate the incoming light cone from O2 by the angle defined by the rotation matrix
        x1, y1, z1 = rotation(x, y, z, rot)
        # transform back to spherical coordinate of the light cone from O2, theta < [0, pi], could be > pi/2
        r1, theta1, phi1 = cart2sph(z1, y1, x1)
        theta1 = (theta1 <= np.pi / 2) * theta1  # convert theta1 > pi/2 to zeros, later remove all the zeros values

        # light cone after interface n2-n3
        phi2 = -phi1
        theta2 = np.arcsin(np.sin(theta1) * self.n2 / self.n3)
        kr2 = self.n3 / self.wavelength * np.sin(theta2)
        ky2, kx2 = pol2cart(kr2, phi2)  # get back to k space

        k_dl = self.na2 / self.wavelength  # bandwidth supported by o2
        ind_i, ind_j = np.where(kr <= k_dl)  # calculate iterators
        ####
        # not binning inside
        # reconstruct the effective pupil plane
        # pupil2 = np.zeros((self.size, self.size), pupil.dtype)
        # counters = np.zeros((self.size, self.size))
        # for i, j in zip(ind_i, ind_j):
        #     pupil2[int(ky2[i][j] / delta_k + self.size / 2)][int(kx2[i][j] / delta_k + self.size / 2)] \
        #             += pupil[i][j]
        #     counters[int(ky2[i][j] / delta_k + self.size / 2)][int(kx2[i][j] / delta_k + self.size / 2)] += 1
        # # reset the center in pupil2 as it it overestimated
        # pupil2[self.size // 2][self.size // 2] = pupil2[self.size // 2 + 1][self.size // 2 + 1]
        # if is_mag is True:
        #     return pupil2          # return the cumulated magnitude
        # return pupil2 / np.maximum(counters, 1)  # get the average phase
        ####
        # with internal binning
        # reconstruct the effective pupil plane
        new_size = self.size // self.binning
        pupil2 = np.zeros((new_size, new_size), pupil.dtype)
        counters = np.zeros((new_size, new_size))
        for i, j in zip(ind_i, ind_j):
            ind_ky2 = np.rint(ky2[i][j] / delta_k / self.binning).astype(int)
            ind_kx2 = np.rint(kx2[i][j] / delta_k / self.binning).astype(int)
            pupil2[ind_ky2 + new_size // 2][ind_kx2 + new_size // 2] += pupil[i][j]
            counters[ind_ky2 + new_size // 2][ind_kx2 + new_size // 2] += 1
        if is_mag is True:
            # reset the center in magnitude as it it overestimated
            pupil2[new_size // 2][new_size // 2] = pupil2[new_size // 2 + 1][new_size // 2 + 1]
            return pupil2 / self.binning ** 2          # return the cumulated magnitude
        return pupil2 / np.maximum(counters, 1)        # get the average phase

    def _calculate_3d_psf_o3(self):
        """calculate the 3d psf after o3"""
        delta_k = self.delta_k * self.binning  # unit frequency for O3 pupil function
        size = self.phase_o3.shape[0]

        k_min = -delta_k * size / 2
        k_max = delta_k * size / 2
        k = np.linspace(k_min, k_max, size, endpoint=False)
        (kx, ky) = np.meshgrid(k, k)
        kr, phi = cart2pol(ky, kx)

        kmag = self.n3 / self.wavelength
        kz = np.sqrt(np.maximum(kmag ** 2 - kr ** 2, 0))  # set the negative values to 0
        z_range = np.linspace(-self.zres * self.zsize / 2, self.zres * self.zsize / 2, self.zsize)
        defocus = np.exp(2 * np.pi * 1j * kz * z_range[:, np.newaxis, np.newaxis])

        mag = self.mag_o3[np.newaxis]
        phase = self.phase_o3[np.newaxis]
        pupils = mag * np.exp(1j * 2 * np.pi * phase) * defocus
        # psf_a = np.fft.fftshift(np.fft.ifftn(np.fft.fftshift(pupils, axes=(1, 2)), axes=(1, 2)), axes=(1, 2))
        psf_a = np.fft.fftshift(np.fft.ifftn(pupils, axes=(1, 2)), axes=(1, 2))
        return np.abs(psf_a) ** 2  # return PSF

    def _generate_zernikes(self):
        """calculate the zernike modes"""
        delta_k = self.delta_k  # unit frequency
        k_min = -delta_k * self.size / 2
        k_max = delta_k * self.size / 2
        k = np.linspace(k_min, k_max, self.size, endpoint=False)
        (kx, ky) = np.meshgrid(k, k)
        kr, theta = cart2pol(ky, kx)   # from k space to spherical coordinates

        if self.is_zern_o2:
            na = self.na2
        else:
            na = self.na3
        k_dl = na / self.wavelength  # bandwidth supported by o2
        kr = kr / k_dl  # normalize kr
        return zernike(kr, theta, np.arange(1, self.numZern + 1))


if __name__ == "__main__":
    """example client"""
    from tifffile import imwrite
    from utilities.utils import cross_section_3dpsf_2dfit
    mysetup = RemoteFocusing(n2=1.0, na2=0.90, n3=1.525, alpha=30, size=1024, zern_mode=8, zern_coeff=0.3)

    plt.figure(0)
    plt.subplot(2, 2, 1)
    plt.imshow(mysetup.mag_o2)
    imwrite("temp" + '/mag_o2.tiff', np.single(mysetup.mag_o2))
    plt.title("mag_o2")
    plt.subplot(2, 2, 2)
    plt.imshow(mysetup.mag_o3)
    plt.title("mag_o3")
    imwrite("temp" + '/mag_o3.tiff', np.single(mysetup.mag_o3))
    plt.subplot(2, 2, 3)
    plt.imshow(mysetup.phase_o2)
    imwrite("temp" + '/phase_o2.tiff', np.single(mysetup.phase_o2))
    plt.title("phase_o2")
    plt.subplot(2, 2, 4)
    plt.imshow(mysetup.phase_o3)
    imwrite("temp" + '/phase_o3.tiff', np.single(mysetup.phase_o3))
    plt.title("phase_o3")
    plt.savefig("before_bin")

    imwrite("temp" + '/psf3d.tiff', np.single(mysetup.psf3d))

    plt.figure(1)
    cross_section_3dpsf_2dfit(mysetup.psf3d, size=64)

    # plt.figure(2)
    # fig, axs = plt.subplots(3, 5, figsize=(20, 12))
    # # fill out plot
    # print("zern shape : ", mysetup.zerns.shape)
    # for ax, (k, v) in zip(axs.ravel(), noll2name.items()):
    #     # ax.matshow(np.ma.array(zern, mask=r > 1), vmin=-1, vmax=1, cmap="coolwarm")
    #     ax.matshow(mysetup.zerns[k-1], vmin=-1, vmax=1, cmap="coolwarm")
    #     ax.set_title(v + r", $Z_{{{}}}^{{{}}}$".format(*noll2degrees(k)))
    #     ax.axis("off")
    # fig.tight_layout()
    plt.show()


