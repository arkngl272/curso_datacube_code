#!/usr/bin/env python
# -*- coding: utf8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from astropy.io import fits

Mpc2cm = 3.08568e24			#conversion factor from kpc to cm
h_p = 6.62607e-34			#planck constant in m^2 kg s^-1
k_b = 1.3806485e-23     	#boltzmann constant in m^2 kg s^-2 K^-1
c_l = 2.99792458e8      	#speed of light in m s^-1
msol = 1.988e30				#solar mass value in kg
dist = 8. 					#distance in kpc
knu0 = 0.192   				#dust emissivity constant at 350 micron in M^2 kg^-1
nu0 = c_l*1e6/350.			#frequency corresponding to 350 micron
const = msol/(dist*Mpc2cm)**2

def B_nu(tt):
	return (2.*h_p*nu**3)/c_l**2*1./(np.exp((h_p*nu)/(k_b*tt))-1)	#in W sr^-1 m^-2 Hz^-1

def ffunc_nu(md,tt,bb):
	mbb_f = md*msol*knu0*(nu/nu0)**bb*B_nu(tt)/(dist*Mpc2cm)**2		#in W cm^-2 Hz^-1
	mbb_erg = 1e7*mbb_f  	#in erg s^-1 cm^-2 Hz^-1
	return mbb_erg*1e23  	#in Jy

def chi2_f(md,tt,bb):
	return np.sum(((ffunc_nu(md,tt,bb)-oflux)/error)**2)

wavel = [70., 160., 250., 350., 500.]
bg = [-0.00064733, -0.00640535,  0.00408735,  0.00482752,  0.00264053]  #this is the background per pixel in Jy
lwavel = np.log10(wavel)
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
results = 'results_map.fits'
hdulist = fits.open(results)
header = hdulist[0].header
chi2_map = hdulist[0].data[0,:,:]
mass_map = hdulist[0].data[1,:,:]
temp_map = hdulist[0].data[2,:,:]
beta_map = hdulist[0].data[3,:,:]
nx = header["NAXIS1"]
ny = header["NAXIS2"]
# plotting section
nwav = 221
llmin = 1.4
llmax = 3.0
dl = (llmax-llmin)/float(nwav-1)
lwav = np.zeros(nwav)
for _i in range(nwav):
	lwav[_i]=llmin+dl*float(_i)
wav = 10**lwav		# wavelength array in microns
nu = np.zeros(nwav)
nu = c_l*1e6*np.divide(1.,wav)
xma = max(wav)*1.1
xmi = min(wav)*0.9
inloop = True
while inloop:
	print ""
	print "    Choose a pixel coordinates with 1<=x<=",nx," and 1<=y<=",ny,"."
	print ""
	wx = input("x-coordinate pixel? ")
	wy = input("y-coordinate pixel? ")
	print ""
	if wx > nx or wy > ny or wx < 1 or wy < 1:
		print ""
		print ""
		print "                     WARNING!!!"
		print "  (",wx,",",wy,") is not a valid pixel coordinate!"
		print ""
		print "        Check that you are within the image!"
	else:
		if chi2_map[wy-1,wx-1] >= 0:
			print ""
			print "      Plotting..."
			print ""
			oflux = np.zeros(5)
			error = np.zeros(5)
			oflux[0] = flux1[wy,wx]-bg[0]
			oflux[1] = flux2[wy,wx]-bg[1]
			oflux[2] = flux3[wy,wx]-bg[2]
			oflux[3] = flux4[wy,wx]-bg[3]
			oflux[4] = flux5[wy,wx]-bg[4]
			error[0] = oflux[0]*0.05
			error[1] = oflux[1]*0.05
			error[2] = oflux[2]*0.07
			error[3] = oflux[3]*0.07
			error[4] = oflux[4]*0.07
			loflux = np.log10(oflux)
			lerror = np.log10(error)
			dmass = mass_map[wy-1,wx-1]
			dtemp = temp_map[wy-1,wx-1]
			dbeta = beta_map[wy-1,wx-1]
			chi2 = chi2_map[wy-1,wx-1]
			bb_emission = ffunc_nu(dmass,dtemp,dbeta)
			print "beta=",dbeta
			print "Mass=",dmass
			print "Temperature=",dtemp
			print "Chi2=",chi2
			print oflux
			plt.errorbar(wavel,oflux, yerr=error, fmt='o')
			yma = max(oflux)*(1.5)
			ymi = yma/7e2
			plt.loglog(wav,bb_emission,linestyle='-',color='black')
			plt.loglog(wavel,oflux, 'ro', color='red')
			plt.xlabel(r'$\lambda$  [$\mu$ m]',fontsize=16)
			plt.ylabel(r'F$_\nu$  [Jy]',fontsize=16)
			plt.tick_params(labelsize=16)
			plt.xlim(xmi,xma)
			plt.ylim(ymi,yma)
			fig = plt.gcf()
			plt.show()
		else:
			print ""
			print "This pixel was not fitted"
			print ""
		whatnow = input("Continue (1) or exit (0) ")
		if whatnow == 0:
			inloop = False
print "Done!"			
			
