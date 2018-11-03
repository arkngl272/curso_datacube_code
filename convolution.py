#!/usr/bin/env python
# -*- coding: utf8 -*-
print 'Importing modules...'
import numpy as np
from astropy.io import fits
from FITS_tools.hcongrid import hcongrid
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
print 'Done!'

print 'Reading data and setting constant...'

psf1 = # put here the value of the psf of the hi-res image in (for example) arcminute
psf2 = # put here the value of the psf of the low-res image in (for example) arcminute
kk = 2.*np.pi/(360.*3600.)					# one arcsec in rad units
img1 = 'image1'
img2 = 'image2'
hdulist1 = fits.open(img1)
hdulist2 = fits.open(img2)
header1 = hdulist1[0].header
header2 = hdulist2[0].header
cube = hdulist1[0].data		#here goes the data at the highest resolution

print 'Done...!'

nl = header1["NAXIS1"]
map = cube[0,m,:,:]
# Note the a convolution must be performed in surface brightness units
# We first define the convolution kernel.
# We will use a gaussian kernel for simplicity, but in reality a proper
# PSF kernel should be used (See, e.g., Aniano+ 2011)
# the sigma of the gaussian (PSF) should be defined in pixels
pixsz1 = np.abs(header1["CDELT1"])*3600.
pixsz2 = np.abs(header2["CDELT1"])*3600.
# first of all, all maps must be in surface brightness units. 
print "I am converting the flux into surface brightness"
flx_sb = map/(kk**2*pixsz1**2)	# converted in K/sr from K 
print "Building the gaussian simulating the PSF"
# the following is the gaussian FWHM in pixels
sigma = np.sqrt(psf2**2-psf1**2)/pixsz1
#kernel = Gaussian2DKernel(stddev=sigma,x_size=xs,y_size=ys) 
kernel = Gaussian2DKernel(stddev=sigma,x_size=40,y_size=40) 
print "Convolving..."
convolved = scipy_convolve(flx_sb, kernel, mode='same')  #, method='direct')
print "Done!"
# now regrids to the new header
flxmap = hcongrid(convolved,header1,header2)
# Now converting from K/sr to K
flx_reg = flxmap*kk**2*pixsz2**2
outfile = 'convolved.fits'

hdu_s = fits.PrimaryHDU(flx_reg,header=header2)
hdu_s.writeto(outfile)
