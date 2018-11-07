#!/usr/bin/env python
# -*- coding: utf8 -*-
import numpy as np
from astropy.io import fits
import pyregion
from scipy.optimize import minimize

Mpc2cm = 3.08568e24			#conversion factor from Mpc to cm
h_p = 6.62607e-34			#planck constant in m^2 kg s^-1
k_b = 1.3806485e-23     	#boltzmann constant in m^2 kg s^-2 K^-1
c_l = 2.99792458e8      	#speed of light in m s^-1
msol = 1.988e30				#solar mass value in kg
dist = 8. 					#distance in Mpc
knu0 = 0.192   				#dust emissivity constant at 350 micron in M^2 kg^-1
nu0 = c_l*1e6/350.			#frequency corresponding to 350 micron
const = msol/(dist*Mpc2cm)**2
fit_thresh = 0.001


def B_nu(tt):
	return (2.*h_p*nu**3)/c_l**2*1./(np.exp((h_p*nu)/(k_b*tt))-1)	#in W sr^-1 m^-2 Hz^-1

def ffunc_nu(md,tt,bb):
	mbb_f = md*const*knu0*(nu/nu0)**bb*B_nu(tt)		#in W cm^-2 Hz^-1
	mbb_erg = 1e7*mbb_f  	#in erg s^-1 cm^-2 Hz^-1
	return mbb_erg*1e23  	#in Jy

def chi2(params):
	md, tt, bb = params
	return np.sum(((ffunc_nu(md,tt,bb)-oflux)*check/error)**2)/(sum(check)-3.)

print "Reading images..."
wavel = [70., 160., 250., 350., 500.]
nu = np.zeros(5)
nu = c_l*1e6*np.divide(1.,wavel)		#frequency array in Hz
## PACS 70
img1 = 'M51_pacs70_convolved.fits'
hdulist1 = fits.open(img1)
flux1 = (hdulist1[0].data[:,:])
## PACS 160
img2 = 'M51_pacs160_convolved.fits'
hdulist2 = fits.open(img2)
flux2 = (hdulist2[0].data[:,:])
## SPIRE250
img3 = 'M51_spire250_convolved.fits'
hdulist3 = fits.open(img3)
flux3 = (hdulist3[0].data[:,:])
## SPIRE 350
img4 = 'M51_spire350_convolved.fits'
hdulist4 = fits.open(img4)
flux4 = (hdulist4[0].data[:,:])
## SPIRE 500
img5 = 'M51_spire500_convolved.fits'
hdulist5 = fits.open(img5)
header = (hdulist5[0].header)
flux5 = (hdulist5[0].data[:,:])
nx = header['NAXIS1']
ny = header['NAXIS2']
print 'Done...!'
print ''
print 'Now calculating the value of the background per pixel...'
## now calculating the average value of the sky in one (now big) pixel
bg = np.zeros(5)
sig = np.zeros(5)
size = (ny, nx)
# PACS 70
reg_name = 'pacs_lowres.reg'
source =  pyregion.open(reg_name)
nsky = len(source) - 1	#	total apertures to measure the sky level
sky_ap = np.zeros(nsky)
sky_pix = []
pix_ap = np.zeros(nsky)
for k in range(nsky):
	sky =  pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	sflx = 0.0
	pix = 0
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				sflx = sflx + flux1[n,i]
				sky_pix.append(flux1[n,i])
				pix = pix + 1  
	sky_ap[k] = sflx
	pix_ap[k] = pix
sig[0] = np.std(sky_pix)
bg[0] =  sum(sky_ap)/sum(pix_ap)	#	this is the average sky per pixel
# PACS 160
reg_name = 'pacs_lowres.reg'
source =  pyregion.open(reg_name)
nsky = len(source) - 1	#	total apertures to measure the sky level
sky_ap = np.zeros(nsky)
sky_pix = []
pix_ap = np.zeros(nsky)
for k in range(nsky):
	sky =  pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	sflx = 0.0
	pix = 0
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				sflx = sflx + flux2[n,i]
				sky_pix.append(flux2[n,i])
				pix = pix + 1  
	sky_ap[k] = sflx
	pix_ap[k] = pix
sig[1] = np.std(sky_pix)
bg[1] =  sum(sky_ap)/sum(pix_ap)	#	this is the average sky per pixel
# SPIRE 250
reg_name = 'spire_lowres.reg'
source =  pyregion.open(reg_name)
nsky = len(source) - 1	#	total apertures to measure the sky level
sky_ap = np.zeros(nsky)
sky_pix = []
pix_ap = np.zeros(nsky)
for k in range(nsky):
	sky =  pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	sflx = 0.0
	pix = 0
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				sflx = sflx + flux3[n,i]
				sky_pix.append(flux3[n,i])
				pix = pix + 1  
	sky_ap[k] = sflx
	pix_ap[k] = pix
sig[2] = np.std(sky_pix)
bg[2] =  sum(sky_ap)/sum(pix_ap)	#	this is the average sky per pixel
# SPIRE 350
reg_name = 'spire_lowres.reg'
source =  pyregion.open(reg_name)
nsky = len(source) - 1	#	total apertures to measure the sky level
sky_ap = np.zeros(nsky)
sky_pix = []
pix_ap = np.zeros(nsky)
for k in range(nsky):
	sky =  pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	sflx = 0.0
	pix = 0
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				sflx = sflx + flux4[n,i]
				sky_pix.append(flux4[n,i])
				pix = pix + 1  
	sky_ap[k] = sflx
	pix_ap[k] = pix
sig[3] = np.std(sky_pix)
bg[3] =  sum(sky_ap)/sum(pix_ap)	#	this is the average sky per pixel
# SPIRE 500
reg_name = 'spire_lowres.reg'
source =  pyregion.open(reg_name)
nsky = len(source) - 1	#	total apertures to measure the sky level
sky_ap = np.zeros(nsky)
sky_pix = []
pix_ap = np.zeros(nsky)
for k in range(nsky):
	sky =  pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	sflx = 0.0
	pix = 0
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				sflx = sflx + flux5[n,i]
				sky_pix.append(flux5[n,i])
				pix = pix + 1  
	sky_ap[k] = sflx
	pix_ap[k] = pix
sig[4] = np.std(sky_pix)
bg[4] =  sum(sky_ap)/sum(pix_ap)	#	this is the average sky per pixel
print 'Done...!'
print ''
print 'Now performing the fit for every pixel...'
bounds = [(1e1,1e6), (5.,60.), (0.5, 3.0)]
x0 = [1e2, 20.,2.0] 
results = np.zeros([4,ny,nx])
for yyy in range(ny):
	for xxx in range(nx):
		oflux = np.zeros(5)
		error = np.zeros(5)
		oflux[0] = flux1[yyy,xxx]-bg[0]
		oflux[1] = flux2[yyy,xxx]-bg[1]
		oflux[2] = flux3[yyy,xxx]-bg[2]
		oflux[3] = flux4[yyy,xxx]-bg[3]
		oflux[4] = flux5[yyy,xxx]-bg[4]
		error[0] = oflux[0]*0.05
		error[1] = oflux[1]*0.05
		error[2] = oflux[2]*0.07
		error[3] = oflux[3]*0.07
		error[4] = oflux[4]*0.07
		check = np.zeros(5)
		for nnn in range(5):			
			if oflux[nnn] >= 3.*sig[nnn]:
				check[nnn] = 1
		if sum(check) > 3:
			print 'Fitting pixel ',xxx,yyy
			chisq = minimize(chi2,x0,bounds=bounds,method='TNC')	
			results[0,yyy,xxx] = chisq.fun
			results[1,yyy,xxx] = chisq.x[0]
			results[2,yyy,xxx] = chisq.x[1]
			results[3,yyy,xxx] = chisq.x[2]
			print 'Reduced chi2=',chisq.fun
		else:
			print 'Skipping...'
			results[0,yyy,xxx] = -999.
			results[1,yyy,xxx] = -999.
			results[2,yyy,xxx] = -999.
			results[3,yyy,xxx] = -999.
print 'Done...!'
print ''
print 'Now writing files'
outfile = "results_map.fits"
newheader = header.copy() 
del newheader[5:24]
del newheader[14:]
newheader['NAXIS'] = 3
newheader.insert(5, ('NAXIS3', 4))
newheader['PLANE0'] = 'chi2'
newheader['PLANE1'] = 'Dust mass [Msol]'
newheader['PLANE2'] = 'Dust temperature [K]'
newheader['PLANE3'] = 'Emissivity index'
hdu_s = fits.PrimaryHDU(results,header=newheader)
hdu_s.writeto(outfile)

