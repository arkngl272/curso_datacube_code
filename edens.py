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

def sfr_func(ha):
	sfr = 4.*np.pi*(dist*mpc2cm)**2*ha*4.6e-42*1e-20
	return sfr

def edens(R):
	if R < 0.4:
		R = 0.4
	if R > 1.435:
		R = 1.435
	rho = 0.0543*np.tan(-3.0553*R+2.8506)+6.98-10.6905*R+9.9186*R**2-3.5422*R**3		
	return rho

## define constants
dist = 180.
mpc2cm = 3.0856775815E+24

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
Sb_mat = k_cube[9,:,:]
Sr_mat = k_cube[10,:,:]
av_mat = np.zeros([ny,nx])
sfr_map = np.zeros([ny,nx])
eld_map = np.zeros([ny,nx])
for _x in range(nx):
	for _y in range(ny):
		if (Ha_mat[_y,_x] > 0.0 and Hb_mat[_y,_x] > 0.0):
			av_mat[_y,_x] = Av_func(Ha_mat[_y,_x],Hb_mat[_y,_x])
		if (Hb_mat[_y,_x] <= 0.0 or av_mat[_y,_x] < 0.):
			av_mat[_y,_x] = 0.0
		if Ha_mat[_y,_x] > 0.0 :
			Ha_corr = Ha_mat[_y,_x]*10**(0.4*av_mat[_y,_x]*Av_a)
			sfr_map[_y,_x] = sfr_func(Ha_corr)
		else:
			sfr_map[_y,_x] = -999.
		if  (Sb_mat[_y,_x] > 0. and Sr_mat[_y,_x] > 0.):
			r = Sb_mat[_y,_x]/Sr_mat[_y,_x]
			eld_map[_y,_x] = edens(r)
		else:
			eld_map[_y,_x] = -999.
newheader1 = k_head.copy()
del newheader1[5:19]
del newheader1[19:]
newheader1["NAXIS"] = 2
newheader2 = newheader1.copy()
newheader3 = newheader1.copy()
newheader1["COMMENT"]= "Extinction map for the galaxy JO36"
newheader1["COMMENT"]= "Made with the routine edens.py"
newheader2["COMMENT"]= "SFR map for the galaxy JO36"
newheader2["COMMENT"]= "Made with the routine edens.py"
newheader3["COMMENT"]= "Electron density map for the galaxy JO36"
newheader3["COMMENT"]= "Made with the routine edens.py"
hdu_ext = fits.PrimaryHDU(av_mat,header=newheader1)
hdu_ext.writeto('ext_map.fits')
hdu_s = fits.PrimaryHDU(sfr_map,header=newheader2)
hdu_s.writeto('sfr_map.fits')
hdu_s = fits.PrimaryHDU(sfr_map,header=newheader3)
hdu_s.writeto('edens_map.fits')
			
