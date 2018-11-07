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
bounds = [(1e4,1e9), (5.,50.), (0.5,3.)]
x0 = [1e7, 25., 1.8] 
nrun = 1000
chi2_arr = np.zeros(nrun)
dmass_arr = np.zeros(nrun)
dtemp_arr = np.zeros(nrun)
dbeta_arr = np.zeros(nrun)
print "Now finding the best fit and errorbars..."
for i in range(nrun):
	fluxes = np.random.normal(ofluxes, errors)
	chisq = minimize(chi2,x0,bounds=bounds,method='TNC')	
	chi2_arr[i] = chisq.fun/2.    #this is the reduced chi2: 5 constraints - 3 degrees of freedom
	dmass_arr[i] = chisq.x[0]
	dtemp_arr[i] = chisq.x[1]
	dbeta_arr[i] = chisq.x[2]
print "Done...!"
print ""
print "Now finding acceptable models..."
print ""
deltachi2 = 9.21 + np.min(chi2_arr)    #this is at 99% confidence level  
accept_chi2 = []
accept_mass = []
accept_temp = []
accept_beta = []
for i in range(nrun):
	if chi2_arr[i] <= deltachi2:
		accept_chi2.append(chi2_arr[i])
		accept_mass.append(dmass_arr[i])
		accept_temp.append(dtemp_arr[i])
		accept_beta.append(dbeta_arr[i])
mass_best = np.median(accept_mass)
mass_sigm = np.std(accept_mass)
temp_best = np.median(accept_temp)
temp_sigm = np.std(accept_temp)
beta_best = np.median(accept_beta)
beta_sigm = np.std(accept_beta)
nbins = nrun/50
print "Best fit dust mass: ",mass_best," ± ",mass_sigm
print "Best fit dust temperature: ",temp_best," ± ",temp_sigm
print "Best fit dust emissivity: ",beta_best," ± ",beta_sigm
print "Mass with mean value: ",np.mean(accept_mass)
print "Temp with mean value: ",np.mean(accept_temp)
print "Beta with mean value: ",np.mean(accept_beta)
print "Now plotting..."
####################
# plotting section
####################

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
#fig.subplots_adjust(wspace=0.2,hspace=-0.1) #,top=None,bottom=None)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
ax.set_axis_bgcolor('none')

# plots histogram of the dust mass values
fig1 = fig.add_subplot(2,2,1)
n_m, bins_m, patches = plt.hist(accept_mass, nbins, facecolor='green')
plt.xlabel(r'M$_D$  [M$_\odot$]',fontsize=16)
plt.ylabel(r'N',fontsize=16)

# plots histogram of the dust temperature values
fig1 = fig.add_subplot(2,2,2)
n_t, bins_t, patches = plt.hist(accept_temp, nbins, facecolor='green')
plt.xlabel(r'T$_D$  [K]',fontsize=16)
plt.ylabel(r'N',fontsize=16)

# plots histogram of the dust mass values
fig1 = fig.add_subplot(2,2,3)
n_b, bins_b, patches = plt.hist(accept_beta, nbins, facecolor='green')
plt.xlabel(r'$\beta$',fontsize=16)
plt.ylabel(r'N',fontsize=16)

# plots best model fit to the data
fig4 = fig.add_subplot(2,2,4)
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
bb_emission = ffunc_nu(mass_best,temp_best,beta_best)
plt.plot(wavel, bb_emission, linestyle='-', color='black')
plt.xlabel(r'$\lambda$  [$\mu$ m]',fontsize=16)
plt.ylabel(r'F$_\nu$  [Jy]',fontsize=16)
plt.tick_params(labelsize=16)
plt.xlim(xmi,xma)
plt.ylim(ymi,yma)
plt.show()

