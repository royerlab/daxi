import numpy as np
from matplotlib import pyplot as plt
from pyotf.zernike import zernike, noll2name
from utilities.utils import cart2pol, cross_section_3dpsf
from numpy.linalg import lstsq


class PupilFunction(object):
    """A class defining the pupil function of an optical system"""

    def __init__(self, mag=None, phase=None, numZern=50, unit_circle=True, na=1, wl=500, res=200):
        """Generate a pupil function object based on the given parameters
        most of the operations are for the phase component, magnitude component in pupil function plays a less
        important part

        Input Parameters
        ----------
        mag : ndarray (n, n)
            magnitude component of the pupil function
        phase : ndarray (n, n)
            phase component of the pupil function, unit: wave; multiple by 2pi to get radian
            note that pyotf will calculate phase in the unit of radian
        numZern : numeric
            number of zernike modes to use


        Auto-generated parameters
        ----------
        size : numeric
            number of pixels from the input pupil function
        zerns: ndarray (size, size)
            unit zernike modes on a circular plane
        zern_coef : ndarray (numZern,)
            zernike coefficients
        strehl : numeric
            strehl ratio based on the phase error
        strehl_subtracted : numeric
            strehl ratio based on the phase error with the first four zernike terms removed
        phase_subtracted : ndarray (n, n)
            phase component of the pupil function, after removing the first four zernike terms
        kr : numeric
            maximally supported bandwidth
        """

        if phase is None:
            raise ValueError("No input for phase!")
        self.phase = phase
        if phase.shape[0] != phase.shape[1]:
            raise ValueError("Pupil plane have different dimensions along x and y")
        self.size = phase.shape[0]
        # calculate the normalized k space max value
        self.k_dl = na / wl  # bandwidth limited by diffraction
        if unit_circle:
            self.k_max_nor = 1
        else:
            kmax = 1 / res / 2  # max k determined by the input real space parameters
            self.k_max_nor = kmax / self.k_dl   # normalized max k
        self.kr, self.theta = self._calculate_kr()
        self.valid_points = self.kr <= 1.0
        if mag is None:
            mag = np.ones(phase.shape)
            mag[self.kr > 1.0] = 0
        self.mag = mag
        self.rms = np.std(self.phase[self.valid_points])   # std same as rms when mean is zero
        # self.rms = np.sqrt(np.mean(self.phase[self.valid_points] ** 2))
        # self.rms = np.sqrt(np.mean(self.phase[self.valid_points] ** 2) - np.mean(self.phase[self.valid_points]) ** 2)
        self.strehl = np.exp(- 4 * np.pi ** 2 * self.rms ** 2)
        self.numZern = numZern
        self.zerns = self._generate_zernikes()
        self.zern_coef = self._fit_to_zerns()
        self.phase_subtracted = self._calculate_phase_withoutfirstfour()
        self.rms_substracted = np.std(self.phase_subtracted)
        self.strehl_substracted = np.exp(- 4 * np.pi ** 2 * self.rms_substracted ** 2)

    def _crop_pupil(self):
        """crop the pupil to where the supporting bandwidth ends"""
        raise NotImplementedError

    def _calculate_kr(self):
        """calculate the kr coordinates of the pupil function"""
        x = np.linspace(-self.k_max_nor, self.k_max_nor, self.size, endpoint=False)
        xx, yy = np.meshgrid(x, x)  # xy indexing is default
        kr, theta = cart2pol(yy, xx)
        return kr, theta

    def _generate_zernikes(self):
        """calculate the zernike modes"""
        # x = np.linspace(-1, 1, self.size, endpoint=False)  # generate the base according to pupil size
        # xx, yy = np.meshgrid(x, x)  # xy indexing is default
        # r, theta = cart2pol(yy, xx)
        return zernike(self.kr, self.theta, np.arange(1, self.numZern))

    def _fit_to_zerns(self):
        """get the zernike coeffs by decomposing the pupil function into zernike modes"""
        # x = np.linspace(-1, 1, self.size, endpoint=False)  # generate the base according to pupil size
        # xx, yy = np.meshgrid(x, x)  # xy indexing is default
        # r, theta = cart2pol(yy, xx)
        valid_points = self.valid_points
        data2fit = self.phase[valid_points]
        zerns2fit = self.zerns[:, valid_points].T
        # fit the points
        coefs, _, _, _ = lstsq(zerns2fit, data2fit, rcond=None)
        # return the coefficients
        return coefs

    def plot_named_coefs(self):
        """Plot the first 15 zernike mode coefficients
        These coefficients correspond to the classical abberations
        """
        # set up the subplot
        fig, ax = plt.subplots(1, 1, sharex=True, figsize=(6, 6))
        # get the ordered names
        ordered_names = [noll2name[i + 1] for i in range(len(noll2name))]
        # make an x range for the bar plot
        x = np.arange(len(ordered_names)) + 1
        # pull the data
        data = self.zern_coef[: len(ordered_names)]
        # make the bar plot
        ax.bar(x, data, align="center", tick_label=ordered_names)
        # set up axes
        ax.axis("tight")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_ylabel("Phase Coefficient")
        fig.tight_layout()
        # return figure handles
        return fig, ax

    def plot_coefs(self):
        """Same as `plot_named_coefs` but for all coefs"""
        fig, ax = plt.subplots(1, 1, sharex=True)
        ax.bar(np.arange(self.zern_coef.size) + 1, self.zern_coef)
        ax.axis("tight")
        ax.set_xlabel("Noll's Number")
        fig.tight_layout()
        return fig, ax

    def _calculate_phase_withoutfirstfour(self):
        """calculate the phase without the first four terms"""
        # reconstruct the phase by adding zern5 to zernN
        # copy_coeffs = np.copy(self.zern_coef)
        # copy_coeffs[0:4] = 0  # set the first four to zero
        # return self._recon_from_zerns(copy_coeffs)
        # reconstruct the phase by substracting the first four
        copy_coeffs = self.zern_coef[0:4]
        return self.phase - (copy_coeffs[:, np.newaxis, np.newaxis] * self.zerns[0:4]).sum(0)

    def _recon_from_zerns(self, coeffs):
        """reconstruct pupil from zernike coefs"""
        return (coeffs[:, np.newaxis, np.newaxis] * self.zerns).sum(0)

    def calculate_3d_psf(self, zres, zsize):
        """calculate the 3d psf after o3"""
        "todo, check for correctness with unit circle phase"
        kr, phi = self.kr, self.theta
        kmag = self.k_dl

        kz = np.sqrt(np.maximum(kmag ** 2 - kr ** 2, 0))  # set the negative values to 0
        z_range = np.linspace(-zres * zsize / 2, zres * zsize / 2, zsize)
        defocus = np.exp(2 * np.pi * 1j * kz * z_range[:, np.newaxis, np.newaxis])

        mag = self.mag[np.newaxis]
        phase = self.phase[np.newaxis]
        pupils = mag * np.exp(1j * 2 * np.pi * phase) * defocus
        # psf_a = np.fft.fftshift(np.fft.ifftn(np.fft.fftshift(pupils, axes=(1, 2)), axes=(1, 2)), axes=(1, 2))
        psf_a = np.fft.fftshift(np.fft.ifftn(pupils, axes=(1, 2)), axes=(1, 2))
        return np.abs(psf_a) ** 2  # return PSF


if __name__ == "__main__":
    from opd.opd_by_coverslip import OPDbyCoverslip
    #######
    # generate the two modes using air and water objective
    wavelength = 0.5  # wavelength of the light, unit: um
    opd_water = OPDbyCoverslip(na=1.33 * 0.75, n_media=1.33)
    opd_air = OPDbyCoverslip(na=0.75, n_media=1.0)

    phaseError = (opd_air.opd2d - opd_water.opd2d) / wavelength  # convert the unit from um to wave
    plt.figure()
    plt.imshow(phaseError)
    plt.colorbar()
    # plt.show()

    pf = PupilFunction(phase=phaseError, numZern=20)
    print("raw rms : ", pf.rms)
    print("raw strehl ratio : ", pf.strehl)
    print("the zernike coeffients: ", pf.zern_coef)
    print("rms after removing 1-4th zern : ", pf.rms_substracted)
    print("strehl ratio after removing 1-4th zern : ", pf.strehl_substracted)

    # plt.figure()
    fig, ax = pf.plot_named_coefs()
    # # pf.plot_coefs()

    # plt.figure()
    fig, ax = pf.plot_coefs()

    plt.figure()
    plt.subplot(1, 2, 1)
    plt.imshow(pf.phase)
    plt.title("Original phase")
    plt.colorbar()
    plt.subplot(1, 2, 2)
    plt.imshow(pf.phase_subtracted)
    plt.title("Phase after removal")
    plt.colorbar()
    plt.show()

    # data = pf.calculate_3d_psf(zres=20000, zsize=20)
    # cross_section_3dpsf(data)
    plt.show()