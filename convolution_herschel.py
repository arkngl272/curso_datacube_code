#!/usr/bin/env python
# -*- coding: utf8 -*-
print 'Importing modules...'
import numpy as np
from astropy.io import fits
from FITS_tools.hcongrid import hcongrid
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve

print 'Done!'

wimg = input("Which image are you convolving? (1-4) ")

print ''
print 'Now reading data and setting constants...'
psf = [5.5,6.7,18.,24.,36.]
kk = 2.*np.pi/(360.*3600.)					# one arcsec in rad units
fw2sig = 1./(2.*np.sqrt(2.*np.log(2.)))
# pacs @ 70
if wimg == 1:
	img = 'M51_pacs70.fits'
	outfile = 'M51_pacs70_convolved.fits'
	wind = 0
	print ''
	print 'Processing pacs data at 70 microns'
	print ''
# pacs @ 160
if wimg == 2:
	img = 'M51_pacs160.fits'
	outfile = 'M51_pacs160_convolved.fits'
	wind = 0
	print ''
	print 'Processing pacs data at 160 microns'
	print ''
# spire @ 250
if wimg == 3:
	img = 'M51_spire250.fits'
	outfile = 'M51_spire250_convolved.fits'
	wind = 2
	print ''
	print 'Processing spire data at 250 microns'
	print ''
# spire @ 350
if wimg == 4:
	img = 'M51_spire350.fits'
	outfile = 'M51_spire350_convolved.fits'
	wind = 3
	print ''
	print 'Processing spire data at 350 microns'
	print ''
# read only the header of spire @ 500 to be used as WCS reference
img500 = 'M51_spire500.fits'
hdulist500 = fits.open(img500)
header500 = hdulist500[1].header
pixsz500 = np.abs(header500["CDELT1"])*3600.

hdulist = fits.open(img)
header = hdulist[1].header
flux = hdulist[1].data
nx = header["NAXIS1"]
ny = header["NAXIS2"]
flx_zero = flux.copy()
flx_zero[np.isnan(flux)] = 0.0	# converts NaNs into zeroes
pixsz = np.abs(header["CDELT1"])*3600.
print 'Done...!'
print ''

# Note the regridding process must be performed in surface brightness units!
# We first define the convolution kernel.
# We will use a gaussian kernel for simplicity, but in reality a proper
# PSF kernel should be used (See, e.g., Aniano+ 2011).
##
# First of all, all maps must be in surface brightness units. 
# PACS maps need hence to be converted. We choose the My/sr conversion
# to keep consistency with SPIRE maps.
if wind < 2:
	print "Now converting the flux into surface brightness units..."
	flx_sb = flx_zero/(1e6*kk**2*pixsz**2)	# converted in MJy/sr from Jy/pix 
else:
	flx_sb = flx_zero
print "Building the gaussian simulating the PSF"
print ''
# the following is the sigma of the gaussian, in pixels
sigma = np.sqrt(psf[4]**2-psf[wind]**2)/pixsz*fw2sig
kernel = Gaussian2DKernel(stddev=sigma,x_size=20,y_size=20) 
print ''
print "Convolving..."
convolved = convolve(flx_sb, kernel, mode='same')  #, method='direct')
print "Done!"
# now regrids to the new header
print ''
print 'Regridding...'
flxmap = hcongrid(convolved,header,header500)
# Now converting from MJy/sr to Jy/px
flx_reg = flxmap*1e6*kk**2*pixsz500**2
print ''
print 'Writing the file with the convolved inmage...'
hdu_s = fits.PrimaryHDU(flx_reg,header=header500)
hdu_s.writeto(outfile)
