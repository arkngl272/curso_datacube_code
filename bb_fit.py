#!/usr/bin/env python
# -*- coding: utf8 -*-
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

mpc2cm = 3.08568e24			#conversion factor from Mpc to cm
h_p = 6.62607e-34			#planck constant in m^2 kg s^-1
k_b = 1.3806485e-23     	#boltzmann constant in m^2 kg s^-2 K^-1
c_l = 2.99792458e8      	#speed of light in m s^-1
msol = 1.988e30				#solar mass value in kg
dist = 8.
knu0 = 0.192   				#dust emissivity constant at 350 micron in M^2 kg^-1
nu0 = c_l*1e6/350.			#frequency corresponding to 350 micron

def B_nu(tt):
	return (2.*h_p*nu**3)/c_l**2*1./(np.exp((h_p*nu)/(k_b*tt))-1)	#in W sr^-1 m^-2 Hz^-1

def ffunc_nu(md,tt,bb):
	mbb_f = md*msol*knu0*(nu/nu0)**bb*B_nu(tt)/(dist*mpc2cm)**2		#in W cm^-2 Hz^-1
	mbb_erg = 1e7*mbb_f  	#in erg s^-1 cm^-2 Hz^-1
	return mbb_erg*1e23  	#in Jy

def chi2(params):
	md, tt, bb = params
	return np.sum(((ffunc_nu(md,tt,bb)-fluxes)/errors)**2)

wavel = [70., 160., 250., 350., 500.]
nu = np.zeros(5)
nu = c_l*1e6*np.divide(1.,wavel)		#frequency array in Hz
ofluxes = [194.3, 422.0, 204.8, 84.3, 28.3]
errors =  [11.1, 21.6, 14.3, 5.9, 2.0]

lwavel = np.log10(wavel)
bounds = [(1e5,5e8), (5.,50.), (0.5,3.)]
x0 = [1e7, 25., 2.] 
fluxes = ofluxes
chisq = minimize(chi2,x0,bounds=bounds,method='TNC')	
chi2 = chisq.fun/2.    #this is the reduced chi2: 5 constraints - 3 degrees of freedom
mass = chisq.x[0]
temp  = chisq.x[1]
beta  = chisq.x[2]
print "Best-fit reduced chi2: ",chi2
print "Best-fit dust mass: ",mass," [Msol]"
print "Best-fit dust temperature: ",temp," [K]" 
print "Best-fit dust emissivity index: ",beta

yma = max(ofluxes)*2
ymi = yma/80.
xma = max(wavel)*1.5
xmi = min(wavel)*0.7
plt.xscale('log')
plt.yscale('log')
plt.plot(wavel, ofluxes, 'ro', color='red')
plt.errorbar(wavel, ofluxes, yerr=errors, fmt='o')
# for plotting the best fit model
nwav = 201
llmin = 1.5
llmax = 2.9
dl = (llmax-llmin)/float(nwav-1)
lwav = np.zeros(nwav)
for _i in range(nwav):
	lwav[_i]=llmin+dl*float(_i)
wavel = 10**lwav		# wavelength array in microns
del nu
nu = np.zeros(nwav)
nu = c_l*1e6*np.divide(1.,wavel)
bb_emission = ffunc_nu(mass,temp,beta)
plt.plot(wavel, bb_emission, linestyle='-', color='black')
plt.xlabel(r'$\lambda$  [$\mu$ m]',fontsize=16)
plt.ylabel(r'F$_\nu$  [Jy]',fontsize=16)
plt.tick_params(labelsize=16)
plt.xlim(xmi,xma)
plt.ylim(ymi,yma)
plt.show()

