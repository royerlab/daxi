from opd.opd_by_coverslip import OPDbyCoverslip
from opd.pupil_function import PupilFunction
import matplotlib.pyplot as plt
import numpy as np


wavelength = 0.5
opd_water = OPDbyCoverslip(na=1.33 * 0.75, n_media=1.33, wavelength=wavelength)
opd_air = OPDbyCoverslip(na=0.75, n_media=1.0, wavelength=wavelength)

# plot 1D optical path length
figure, axes = plt.subplots(figsize=(6, 4))
axes.plot(opd_water.k_angle, opd_water.opd1d / wavelength, label="water")
axes.plot(opd_air.k_angle, opd_air.opd1d / wavelength, label="air")
axes.set_xlabel("Angle (degree)")
axes.set_ylabel("Relative optical path length (lambda)")
axes.legend()
plt.savefig('opl.pdf', dpi=300)
plt.show()

# plot 2D optical path difference
figure, axes2 = plt.subplots(figsize=(6, 4))
im = axes2.imshow((opd_air.opd2d - opd_water.opd2d) / wavelength)
plt.colorbar(im)
plt.savefig('opd.pdf', dpi=300)
plt.show()


# fit zernike modes
phaseError = (opd_air.opd2d - opd_water.opd2d) / wavelength # convert the unit from um to wave

pf = PupilFunction(phase=phaseError, numZern=20)
print("raw rms : ", pf.rms)
print("raw strehl ratio : ", pf.strehl)
print("the zernike coeffients: ", pf.zern_coef)
print("rms after removing 1-4th zern : ", pf.rms_substracted)
print("strehl ratio after removing 1-4th zern : ", pf.strehl_substracted)

# plt.figure()
fig, ax = pf.plot_named_coefs()
plt.savefig('zernikes_withName.pdf', dpi=300)
plt.show()

fig, ax = pf.plot_coefs()
plt.savefig('zernikes_noName.pdf', dpi=300)
plt.show()