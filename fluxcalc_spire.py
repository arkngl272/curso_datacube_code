#!/usr/bin/env python
# -*- coding: utf8 -*-
#import ds9   #	doesn't bloody work
import numpy as np
from astropy.io import fits
import pyregion
import pyregion._region_filter as filter

# define some SPIRE-typical constants
beam = [465, 822, 1769]						# beam area, in arcsec^2
ebar = [1.0432705,1.0521442,1.1025926] 		# effective beam area ratio (beta=2 SPIRE UM, table 5.8, page 96; see also M. Smith email, 16/04/14)
cnoise = [0.0058,0.0063,0.0068]				# confusion noise, in Jy/beam, as from Nguyen et al. 2010
kk = 2.*np.pi/(360.*3600.)					# 1 arcsec in rad
img = raw_input('Name of the SPIRE data: ')
reg_name = raw_input('Name of the DS9 file containing the regions for flux and sky measurements: ')
hdulist = fits.open(img)
wavel = input('Wavelength of the image: ')
if wavel == 250:
	wind = 0
if wavel == 350:
	wind = 1
if wavel == 500:
	wind = 2
header = hdulist[1].header
fluxmap = hdulist[1].data
errmap = hdulist[2].data
nx = header["NAXIS1"]
ny = header["NAXIS2"]
size = (ny, nx)		#	this is the size of the image
pixsz = np.abs(header["CDELT1"])*3600.
## converts from MJy/sr to Jy/pix
fluxmap = fluxmap*1e6*kk**2*pixsz**2
##	first of all measure the on-source flux
source =  pyregion.open(reg_name)
nsky = len(source) - 1	#	total apertures to measure the sky level
del source[1:]
smask = source.get_mask(shape=size)
flux = 0.0
errsq = 0.0
on_s_pix = 0
for i in range(nx):
	for n in range(ny):
		if smask[n,i] == True:
			if str(fluxmap[n,i]) != 'nan':
				flux = flux + fluxmap[n,i]
				errsq = errsq + errmap[n,i]**2
				on_s_pix = on_s_pix + 1  	#on-source pixels
sky_ap = np.zeros(nsky)
pix_ap = np.zeros(nsky)
sky_pix = []		# 	This is the array which will contain all the pixels' values in the sky apertures
for k in range(nsky):
	sky = pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	sflx = 0.0
	pix = 0
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				if str(fluxmap[n,i]) != 'nan':
					sflx = sflx + fluxmap[n,i]
					pix = pix + 1  
					sky_pix.append(fluxmap[n,i])
	sky_ap[k] = sflx
	pix_ap[k] = pix
sk_p_pix = sum(sky_ap)/sum(pix_ap)	#	this is the average sky per pixel
skyave = sky_ap/pix_ap				#	average sky per pixel
tflux = flux-on_s_pix*sk_p_pix		#	this is the background subtracted flux 
#	Now calculating the various uncertainties components
#  ######################################################
#	this is the calibration uncertainty
err1 = 0.07*tflux	
#	this is the instrumental uncertainty
err2 = np.sqrt(errsq)
#	this is the confusion noise
err3 = cnoise[wind]*np.sqrt(on_s_pix*pixsz**2/beam[wind])
#	this is the background uncertainty
#		instrumental background
eflx = 0.0
pix_b = 0
for k in range(nsky):
	sky = pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				if str(fluxmap[n,i]) != 'nan':
					eflx = eflx + errmap[n,i]**2	# sums pixels in the errormap, in the "sky"regions
					pix_b = pix_b + 1
err4_1 = np.sqrt(eflx)*on_s_pix/pix_b
print "Instrumental background:",err4_1," Jy"
#		confusion background
err4_2 = cnoise[wind]*np.sqrt(pixsz**2/beam[wind])*on_s_pix/np.sqrt(pix_b)
print "Confusion background:",err4_2," Jy"
#		large-scale background
err4_3 = np.std(sky_ap)/np.sqrt(nsky)
print "Large-scale background:",err4_3," Jy"
err4 = np.sqrt(err4_1**2+err4_2**2+err4_3**2)
print "Calibration error:", err1," Jy"
print "Instrumental error:", err2," Jy"
print "Confusion noise:", err3," Jy"
print "Background error:", err4," Jy"
tot_err = np.sqrt(err1**2 + err2**2 + err3**2 + err4**2)
rel_err = tot_err/tflux*100
print 'Background-subtracted flux:', tflux,'±',tot_err,' Jy  (',rel_err,' %)'
print 'Average sky per pixel value:', sk_p_pix,' Jy'
