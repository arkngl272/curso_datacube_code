#!/usr/bin/env python
# -*- coding: utf8 -*-
print 'Importing modules...'
import numpy as np
from astropy.io import fits
print 'Done!'


def Av_func(ha,hb):
	av = 1./(Av_b - Av_a)*2.5*(np.log10(ha/hb)-np.log10(2.86))
	if av > 6.:
		av = 0.
	return av

# define constants used to correct for dust extinction
av1 = 0.8136172	# galactic extinction curve at 6570 Å
av2 = 0.8167424 # galactic extinction curve at 6550 Å
av3 = 1.271163  # galactic extinction curve at 4870 Å
av4 = 1.278603  # galactic extinction curve at 4850 Å
ha_l = 6562.8	# Halpha central wavelength
hb_l = 4861.3	# Hbeta central wavelength
Av_a = (av1-av2)/(6570.-6550.)*(ha_l-6550.)+av2
Av_b = (av3-av4)/(4870.-4850.)*(hb_l-4850.)+av4
kfile = 'JO36_DATACUBE_FINAL_v1_ec_eo_res_gau_spax_noflag.fits'


hdukin = fits.open(kfile)
k_head = hdukin[0].header
k_cube = hdukin[0].data
nx = k_head["NAXIS1"]
ny = k_head["NAXIS2"]
Ha_mat = k_cube[3,:,:]
Hb_mat = k_cube[6,:,:]
av_mat = np.zeros([ny,nx])
for _x in range(nx):
	for _y in range(ny):
		if (Ha_mat[_y,_x] > 0.0 and Hb_mat[_y,_x] > 0.0):
			av_mat[_y,_x] = Av_func(Ha_mat[_y,_x],Hb_mat[_y,_x])
		if (Hb_mat[_y,_x] <= 0.0 or av_mat[_y,_x] < 0.):
			av_mat[_y,_x] = 0.0
newheader = k_head.copy()
del newheader[5:19]
del newheader[19:]
newheader["NAXIS"] = 2
newheader["COMMENT"]= "Extinction map for the galaxy JO36"
newheader["COMMENT"]= "Made with the routine av.py"
hdu_ext = fits.PrimaryHDU(av_mat,header=newheader)
hdu_ext.writeto('ext_map.fits')
			
