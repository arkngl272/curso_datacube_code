#!/usr/bin/env python
# -*- coding: utf8 -*-
print 'Importing modules...'
import numpy as np
from astropy.io import fits
import pyregion
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
print 'Done!'

def f_k03(x): 
	return ((0.61/(x-0.05))+1.3)

def f_k01(x): 
	return ((0.61/(x-0.47))+1.19)

def f_s10(x):
	return (-0.4-1.5)/(-0.6-0.46)*(x-.46)+1.5

xk03=np.arange(-2,0.04,0.01)
xk01=np.arange(-2,0.46,0.01)
xs10=np.arange(-0.65,0.5,0.01)
k03=f_k03(xk03)
k01=f_k01(xk01)
s10=f_s10(xs10)

kfile = 'JO36_DATACUBE_FINAL_v1_ec_eo_res_gau_spax_noflag.fits'
hdukin = fits.open(kfile)
k_head = hdukin[0].header
k_cube = hdukin[0].data
nx = k_head["NAXIS1"]
ny = k_head["NAXIS2"]
size = (ny,nx)
Ha = k_cube[3,:,:]
Hb = k_cube[6,:,:]
SIIb = k_cube[9,:,:]
SIIr = k_cube[10,:,:]
NIIr = k_cube[5,:,:]
OIII = k_cube[8,:,:]
OIb = k_cube[11,:,:]
OIr = k_cube[12,:,:]
source =  pyregion.open("diag3.reg")
mask = source.get_mask(shape=size)
x_rat1 = []
y_rat1 = []
x_rat2 = []
y_rat2 = []
x_rat3 = []
y_rat3 = []
x_rat4 = []
y_rat4 = []
NII_dia = np.zeros([ny,nx])
for n in range(nx):
	for i in range(ny):
		if mask[i,n] == True:
	# [NII] diagnostic
			if (Ha[i,n] > 0. and Hb[i,n] > 0. and OIII[i,n] > 0. and NIIr[i,n] > 0.):
				xrat = 0.0
				yrat = 0.0
				xrat = np.log10(NIIr[i,n]/Ha[i,n])
				yrat = np.log10(OIII[i,n]/Hb[i,n])
				cond1=((xrat<0.05) & (yrat<f_k03(xrat)))
				cond3=((xrat<0.47) & (yrat>f_k01(xrat)) & (yrat>f_s10(xrat)))
				cond4=((xrat<0.47) & (yrat>f_k01(xrat)) & (yrat<f_s10(xrat)))
				cond2=((~cond1)&(~cond3)&(~cond4))
				if cond1 == True:
					x_rat1.append(xrat)
					y_rat1.append(yrat)
					NII_dia[i,n] = 1
				if cond3 == True:
					x_rat3.append(xrat)
					y_rat3.append(yrat)
					NII_dia[i,n] = 2
				if cond4 == True:
					x_rat4.append(xrat)
					y_rat4.append(yrat)
					NII_dia[i,n] = 0
				if cond2 == True:
					x_rat2.append(xrat)
					y_rat2.append(yrat)
					NII_dia[i,n] = 6
			
#fig = plt.figure(figsize=(14,6))
#fig1 = fig.add_subplot(1,2,1)
plt.plot(x_rat1, y_rat1, 'ro', color='red')
plt.plot(x_rat2, y_rat2, 'ro', color='orange')
plt.plot(x_rat3, y_rat3, 'ro', color='blue')
plt.plot(x_rat4, y_rat4, 'ro', color='green')
plt.plot(xk03,k03,linestyle='-', color='blue')
plt.plot(xk01,k01,linestyle=':', color='black',linewidth=4)
plt.plot(xs10,s10,linestyle='--', color='green',linewidth=3)
plt.xlabel(r'log([NII]/H$\alpha$)',fontsize=16)
plt.ylabel(r'log([OIII]/H$\beta$)',fontsize=16)
plt.xlim(-1,0.4)
plt.ylim(-1.15,1.)

#prova = np.zeros([ny,nx])+10
#fig1 = fig.add_subplot(1,2,2)
#plt.imshow(NII_dia)
#plt.xlabel(r'Pixel Coordinates',fontsize=16)
#plt.ylabel(r'Pixel Coordinates',fontsize=16)
#plt.xlim(0,nx)
#plt.ylim(0,ny)

plt.show()


#
#newheader1 = k_head.copy()
#del newheader1[5:19]
#del newheader1[19:]
#newheader1["NAXIS"] = 2
#newheader2 = newheader1.copy()
#newheader3 = newheader1.copy()
#newheader1["COMMENT"]= "Extinction map for the galaxy JO36"
#newheader1["COMMENT"]= "Made with the routine edens.py"
#newheader2["COMMENT"]= "SFR map for the galaxy JO36"
#newheader2["COMMENT"]= "Made with the routine edens.py"
#newheader3["COMMENT"]= "Electron density map for the galaxy JO36"
#newheader3["COMMENT"]= "Made with the routine edens.py"
#hdu_ext = fits.PrimaryHDU(av_mat,header=newheader1)
#hdu_ext.writeto('ext_map.fits')
#hdu_s = fits.PrimaryHDU(sfr_map,header=newheader2)
#hdu_s.writeto('sfr_map.fits')
#hdu_s = fits.PrimaryHDU(sfr_map,header=newheader3)
#hdu_s.writeto('edens_map.fits')
			
