# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 16:08:57 2019

@author: dingxu
"""

from photutils.datasets import make_100gaussians_image
from photutils import CircularAperture, CircularAnnulus
import matplotlib.pyplot as plt
from photutils import aperture_photometry

data = make_100gaussians_image()
positions = [(145.1, 168.3), (84.5, 224.1), (48.3, 200.3)]
aperture = CircularAperture(positions, r=5)
annulus_aperture = CircularAnnulus(positions, r_in=10, r_out=15)

apers = [aperture, annulus_aperture]
#norm = simple_norm(data, 'sqrt', percent=99)
plt.imshow(data, cmap = 'gray')
aperture.plot(color='white', lw=2)
annulus_aperture.plot(color='red', lw=2)

phot_table = aperture_photometry(data, apers)
print(phot_table)

bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
bkg_sum = bkg_mean * aperture.area
final_sum = phot_table['aperture_sum_0'] - bkg_sum

phot_table['residual_aperture_sum'] = final_sum

print(phot_table['residual_aperture_sum'])
 