from tifffile import imwrite
from utilities.utils import cross_section_3dpsf_2dfit, cross_section_3d_1dfit
from psf.remote_focusing import RemoteFocusing
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage


res = 0.1  # in um
mag = 1.33

######################
# High NA straight
mysetup = RemoteFocusing(wavelength=0.525, n2=1, na2=0.95, n3=1.525, na3=1.0, alpha=0, size=1024, zsize=128, res=res)

figure, axes = plt.subplots(figsize=(3, 2))

im = axes.imshow((mysetup.mag_o3 > 0).astype(float), cmap='gray')
plt.savefig('high_na_straight_pupil.pdf', dpi=300)

cross_section_3d_1dfit(mysetup.psf3d, size=128, res=res * 1000 / 1.33, plt_profiles=True)
# cross_section_3dpsf_2dfit(mysetup.psf3d, size=128, res=res * 1000 / 1.33, plt_profiles=True)
plt.savefig('high_na_straight_psf.pdf', dpi=300)


# High NA oblique
mysetup = RemoteFocusing(wavelength=0.525, n2=1, na2=0.95, n3=1.525, na3=1.0, alpha=30, size=1024, zsize=128, res=res)

figure, axes = plt.subplots(figsize=(3, 2))
im = axes.imshow(mysetup.mag_o3 / np.amax(mysetup.mag_o3) * 1.0, cmap='gray')
plt.savefig('high_na_oblique_pupil.pdf', dpi=300)

cross_section_3d_1dfit(mysetup.psf3d, size=128, res=res * 1000 / 1.33, plt_profiles=True)
plt.savefig('high_na_oblique_psf.pdf', dpi=300)

# rotate psf to get the plots
angle = -10  # degree
rotated_psf = ndimage.rotate(mysetup.psf3d, angle=angle, axes=(2, 0))
cross_section_3d_1dfit(rotated_psf, size=128, res=res * 1000 / 1.33, plt_profiles=True)
plt.savefig('high_na_oblique_psf_rotated.pdf', dpi=300)
