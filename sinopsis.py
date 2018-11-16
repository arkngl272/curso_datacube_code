#!/usr/bin/env python
# -*- coding: utf8 -*-
print ''
print ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print ' %                                                                           % '
print ' %      This routine performs a series of checks on the fits to MUSE         % '
print ' %   data using SINOPSIS. It will create a total chi2 map, plus a chi2 map   % '
print ' %       for each continuum band. It will then provide some statistics       % '
print ' %    over the general goodness of the fitting run, which might be used to   % '
print ' %     more carefully and directly check the spectra (included routine).     % '
print ' %                                                                           % '
print ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print ''
print ''
print 'Importing modules'
#import pyfits
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import fnmatch
import os
rc('text',usetex=True)
print 'Done!'
print 'Data aquisition...'
fitmask = 'null'
z_mask = 'null'
z_mask_abs = 'null'
z_mask_em = 'null'
for file in os.listdir('.'):
	if fnmatch.fnmatch(file, '*DATACUBE_FINAL*'):
		obscube = file
	if fnmatch.fnmatch(file, '*z_mask.fits'):
		z_mask = file
	if fnmatch.fnmatch(file, '*z_abs_mask.fits'):
		z_mask_abs = file
	if fnmatch.fnmatch(file, '*z_em_mask.fits'):
		z_mask_em = file
	if fnmatch.fnmatch(file, '*modelcube.fits*'):
		model = file
	if fnmatch.fnmatch(file, '*nolines*'):
		model_nolines = file
	if fnmatch.fnmatch(file, '*fitmask*'):
		fitmask = file
	if fnmatch.fnmatch(file, '*.dat'):
		for i in range(len(file)-2):
			if file[i:i+4] == '.dat':
				break
		rootname = file[0:i]
try: obscube
except: 
	print ""
	print "                   WARNING!!!"
	print ""
	print "The observed datacube doesn't have a standard name"
	print ""
	obscube = raw_input("Give me the observed datacube filename: ")
z_g = True
if z_mask == 'null':
	z_g = False
if (z_mask_abs == 'null' and z_mask_em == 'null'):
	z_e = False
	z_a = False
if (z_mask_abs != 'null' and z_mask_em == 'null'):
	z_e = False
	z_a = True
if (z_mask_abs == 'null' and z_mask_em != 'null'):
	z_e = True
	z_a = False
if (z_mask_abs != 'null' and z_mask_em != 'null'):
	z_e = True
	z_a = True
if (z_g and (not z_a) and (not z_e)):	
	print ""
	print "                   WARNING!!!"
	print ""
	print "The redshift mask doesn't have a standard name"
	print "A common redshift mask is will be used for gas and stars:"
	print ""
	z_g = True
	z_mask = raw_input("Give me the redhisft mask filename: ")
hdu_data = fits.open(obscube)
hdu_model = fits.open(model)
hdu_model_nolines = fits.open(model_nolines)
obs = hdu_data[1].data
if (not z_a and not z_e):
	hdu_zmask = fits.open(z_mask)
	z_mat = hdu_zmask[0].data
if (not z_g):
	if z_e:
		hdu_zmask_e = fits.open(z_mask_em)
		z_mat_e = hdu_zmask_e[0].data
	hdu_zmask_a = fits.open(z_mask_abs)
	z_mat_a = hdu_zmask_a[0].data	
m_mat = hdu_model[0].data
mn_mat = hdu_model_nolines[0].data
obs_head = hdu_data[1].header
nx = obs_head["NAXIS1"]
ny = obs_head["NAXIS2"]
nl = obs_head["NAXIS3"]
dl = obs_head["CDELT3"]
lmin = obs_head["CRVAL3"]
wavel = np.zeros(nl)
for x in range(nl):
	wavel[x] = lmin+x*dl
print "Done!"
print ""
ibound = [4600.,4845.,4858.,4870.,5040.,5210.,5400.,5650.,5955.,6150.,6400.,6620.,6820.,7110.]
ubound = [4750.,4853.,4864.,4878.,5140.,5310.,5500.,5800.,6055.,6250.,6490.,6690.,6920.,7210.]
nb = len(ibound)
if fitmask != 'null':
	hdu_fitmask = fits.open(fitmask)
	fm_mat = hdu_fitmask[0].data
else:
	fm_mat = np.zeros((ny,nx))
	for _y in range(ny):
		for _x in range(nx):
			if m_mat[0,_y,_x] >-10:
				fm_mat[_y,_x] = 1
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
		if fm_mat[wy-1,wx-1] == 1:
			print ""
			print "      Plotting..."
			print ""
			if z_g:
				redshift = z_mat[wy-1,wx-1]
				redshift_e = redshift
			else:
				redshift = z_mat_a[wy-1,wx-1]
				if z_e:
					redshift_e = z_mat_e[wy-1,wx-1]
				else:
					redshift_e = z_mat_a[wy-1,wx-1]
				if redshift_e < 0.:
					redshift_e = z_mat_a[wy-1,wx-1]
				if redshift < 0.:
					redshift = z_mat_e[wy-1,wx-1]
			obs_spec = obs[:,wy-1,wx-1]
			bf = m_mat[:,wy-1,wx-1]
			bf_nl = mn_mat[:,wy-1,wx-1]
			yma1 = np.mean(obs_spec)*3.5
			alpha_c = 6562.80*(1.0+redshift_e)
			for _i in range(nl):
				if wavel[_i] >= alpha_c:
					break		
			yma2 = obs_spec[_i]*1.3
			yma = max(yma1,yma2)
			ymi = -0.05*yma
			xma = max(wavel)*1.01
			xmi = min(wavel)*0.97
			f,axarr = plt.subplots(2,sharex=True,figsize=(18,15))
# upper plot		
			axarr[0].set_xlim(xmi,xma)
			axarr[0].set_ylim(ymi,yma)
			axarr[0].set_ylabel(r'F$_\lambda$  [1e20 erg/s/cm$^2$/\AA]',fontsize=20)
			axarr[0].tick_params(labelsize=20)
			tstring = 'Pixel coordinates: x=' + str(wx)+' y='+str(wy)	
			axarr[0].set_title(tstring)
			axarr[0].plot(wavel,bf,label="Best Fit",linestyle='-',color='red')
			axarr[0].plot(wavel,obs_spec,label="Observed",linestyle=':',color='black')
			axarr[0].legend()
			for _i in range(nb):
				xlow = ibound[_i]*(1.0+redshift)
				xhig = ubound[_i]*(1.0+redshift)
				x_b = [xlow,xhig]
				avey = np.mean(obs_spec)*0.75
				y_b = [avey,avey]
				axarr[0].plot(x_b,y_b,linestyle='-',color='black')		
				axarr[1].plot(x_b,y_b,linestyle='-',color='black')		
			
# mid-panel plot
			axarr[1].set_xlim(xmi,xma)
			axarr[1].set_ylim(ymi,yma)
			axarr[1].set_ylabel(r'F$_\lambda$  [1e20 erg/s/cm$^2$/\AA]',fontsize=20)
			axarr[1].plot(wavel,bf_nl,label="Best Fit without emission lines",linestyle='-',color='green')
			axarr[1].plot(wavel,obs_spec,label="Observed",linestyle=':',color='black')
			axarr[1].tick_params(labelsize=20)
			axarr[1].legend()

			#ofname = 'spec_'+str(wx)+'_'+str(wy)+'.dat'
			#mfname = 'mod_'+str(wx)+'_'+str(wy)+'.dat'
			#spectrum = np.array([wavel,obs_spec])
			#spectrum = np.array([wavel,bf])
			#spectrum = spectrum.T
			#np.savetxt(mfname,spectrum,fmt=['%10.4e','%12.5e'])
# plot reference lines
	#[NII]I
			x_lin = 6583.5*(1.+redshift_e)
			xmi_ma = [x_lin,x_lin]
			ymi_ma = [ymi,yma]
			axarr[0].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
			axarr[1].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
	#Halpha
			x_lin = 6562.8*(1.+redshift_e)
			xmi_ma = [x_lin,x_lin]
			ymi_ma = [ymi,yma]
			axarr[0].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
			axarr[1].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
	#[NII]II
			x_lin = 6548.0*(1.+redshift_e)
			xmi_ma = [x_lin,x_lin]
			ymi_ma = [ymi,yma]
			axarr[0].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
			axarr[1].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
	#Na(D)I
			x_lin = 5895.92*(1.+redshift)
			xmi_ma = [x_lin,x_lin]
			ymi_ma = [ymi,yma]
			axarr[0].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
			axarr[1].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
	#Na(D)II
			x_lin = 5889.95*(1.+redshift)
			xmi_ma = [x_lin,x_lin]
			ymi_ma = [ymi,yma]
			axarr[0].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
			axarr[1].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
	#Mag
			x_lin = 5176.7*(1.+redshift)
			xmi_ma = [x_lin,x_lin]
			ymi_ma = [ymi,yma]
			axarr[0].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
			axarr[1].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
	#Hbeta
			x_lin = 4861.3*(1.+redshift_e)
			xmi_ma = [x_lin,x_lin]
			ymi_ma = [ymi,yma]
			axarr[0].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
			axarr[1].plot(xmi_ma,ymi_ma,linestyle='--',color='blue')		
			fig = plt.gcf()
			plt.show()
			whatnow = input("Continue (1) or exit (0) ")
			if whatnow == 0:
				inloop = False
		elif fm_mat[wy-1,wx-1] == 0:
			print "No valid redshift avalilable for this pixel!"
			print "Please chose another one"
		elif fm_mat[wy-1,wx-1] == -1:
			print "The average flux was too low for this pixel," 
			print "hence sinopsis did not attempt a fit."
			print "Please chose another one"
		else:
			print "Uh-Oh... something went wrong here...?"	
print "Done!"
