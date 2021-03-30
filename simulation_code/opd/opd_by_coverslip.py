import numpy as np
from matplotlib import pyplot as plt
from utilities.utils import cart2pol


class OPDbyCoverslip(object):
    """A class defining the optical path difference when putting a coverslip
    along the path from the object to the lens
    """

    def __init__(self, n_media=1.0, n_cs=1.525, d=170, na=0.75, size=128, wavelength=0.5):
        """Generate a opd object based on the given parameters

        Input Parameters
        ----------
        n_media : numeric
            refractive index of the immersion media, e.g. water, air, oil, etc.
        n_cs : numeric
            refractive index of the coverslip, mostly glass
        d : numeric
            thickness of the coverslip, unit: um
        na : numeric
            numerical aperture of the objective lens
        size : int
            the sampling size
        wavelength : numeric
            wavelength of the light, unit: um


        Auto generated parameters
        ----------
        k_na  : ndarray (size, )
            sampling of the NA
        opd1d: ndarray (1, size)
            the calculated opd along 1D
        opd2d: ndarray (size, size)
            the calculated opd along 2D
        """
        self.n_media = n_media
        self.n_cs = n_cs
        self.d = d
        self.na = na
        self.size = size
        self.wavelength = wavelength
        self.k_na = np.linspace(0, na, size)
        self.k_angle = np.degrees(np.arcsin(self.k_na / self.n_media))
        self.opd1d = self._calculate_opd1d()
        self.opd2d = self._calculate_opd2d()

    def _calculate_opd1d(self):
        """calculate the opd along 1d"""
        opd1d = np.zeros(self.size)
        angle_media = np.arcsin(self.k_na[1:] / self.n_media)
        angle_cs = np.arcsin(self.k_na[1:] / self.n_cs)
        opd1d[1:] = self.d * self.n_cs * (1 / np.cos(angle_cs) - 1) + self.d * self.n_media * np.tan(angle_cs) \
            / np.sin(angle_media) * (np.cos(angle_media) - 1)
        # opd1d = np.insert(_array, 0, 0)
        return opd1d

    def _calculate_opd2d(self):
        """calculate the 2D opd by assuming rotational symmetry of opd1d"""
        k1d_half = np.arange(self.size)
        k1d = np.concatenate(([self.size], np.flip(k1d_half[1:]), k1d_half))
        kxx, kyy = np.meshgrid(k1d, k1d)
        kr, phi = cart2pol(kyy, kxx)
        opd2d = np.zeros((2*self.size, 2*self.size))
        for i in range(2 * self.size):
            for j in range(2 * self.size):
                ind = int(round(kr[i, j]))
                if ind < self.size:
                    opd2d[i, j] = self.opd1d[ind]
        return opd2d

    def plot1d(self, label="none"):
        """plot the opd1d as a function of the angle, opd is calculated as a function as NA, but plotting over
        angle makes it easier to compare over objectives"""
        plt.plot(self.k_angle, self.opd1d, label=label)
        plt.xlabel("angle (degree)")
        plt.ylabel("Relative optical path length (um)")


if __name__ == "__main__":
    # na = 0.9975, it is to assure the objectives have the same max collection angle
    opd_water = OPDbyCoverslip(na=1.33 * 0.75, n_media=1.33)
    opd_air = OPDbyCoverslip(na=0.75, n_media=1.0)
    # opdModel.plot1d()
    # opdModel.plot2d()
    # plt.plot(opdModel.k_na, opdModel.opd1d)
    plt.figure(0)
    plt.subplot(2, 1, 1)
    opd_water.plot1d("water")
    opd_air.plot1d("air")
    plt.legend()
    # plt.show()

    # plt.figure(1)
    plt.subplot(2, 1, 2)
    plt.plot(opd_air.k_angle, opd_air.opd1d - opd_water.opd1d)
    plt.xlabel("angle (degree)")
    plt.ylabel("Relative OPD (um)")

    plt.figure(1)
    plt.suptitle("OPD when putting a coverslip along the path")
    plt.subplot(1, 3, 1)
    plt.imshow(opd_air.opd2d)
    plt.colorbar()
    plt.title("OPD air")

    plt.subplot(1, 3, 2)
    plt.imshow(opd_water.opd2d)
    plt.colorbar()
    plt.title("OPD water")

    plt.subplot(1, 3, 3)
    plt.imshow(opd_air.opd2d - opd_water.opd2d)
    plt.colorbar()
    plt.title("relative")

    plt.show()
    # print("max angle in air :", opd_air.k_angle[-1])
    # print("marginal ray OPD for air objective : ", opd_air.opd1d[-1])
    # print("max angle in water :", opd_water.k_angle[-1])
    # print("marginal ray OPD for water objective : ", opd_water.opd1d[-1])
