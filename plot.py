#!/usr/bin/env python
# -*- coding: utf8 -*-
print ''
print ''
print 'Importing modules...'
import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
rc('text', usetex=True)
print 'DONE!'

Dist = 10.

fitsimg = "TYPE1_MODELS/t5_p0_q0_oa30_R10_Mcl0.97_i0_total.fits"
hdulist = fits.open(fitsimg)
header = (hdulist[0].header)
pix = header["CDELT1"]

sedfile = "TYPE1_MODELS/t5_p0_q0_oa30_R10_Mcl0.97_i0_sed.dat"
tdata = ascii.read(sedfile)
lamb = tdata['col1']
tflx = tdata['col2']
agn = tdata['col3']
dflx = tflx - agn

conv = pix/(Dist*1e6/206265.)

regfile = "sed.dat"
regdata = ascii.read(regfile)
reg1 = regdata[0][:]*conv**2
reg2 = regdata[1][:]*conv**2
reg3 = regdata[2][:]*conv**2
reg4 = regdata[3][:]*conv**2
tot = regdata[4][:]*conv**2


yma = max(tflx)*1.2
ymi = yma*0.001
xma = max(lamb)*1.01
xmi = min(lamb)*0.97
plt.figure(figsize=(15,10))
plt.loglog(lamb,tflx,linestyle='-',color='black')
plt.loglog(lamb,agn,linestyle='-',color='blue')
plt.loglog(lamb,dflx,linestyle='-',color='red')
plt.loglog(lamb,tot,linestyle='--',color='black',linewidth=4)
plt.loglog(lamb,reg1,linestyle='-',color='green')
plt.loglog(lamb,reg2,linestyle='-',color='magenta')
plt.loglog(lamb,reg3,linestyle='-',color='purple')
plt.loglog(lamb,reg4,linestyle='-',color='gray')
plt.xlabel(r'$\lambda$    [$\mu$m]',fontsize=20)
plt.ylabel(r'$\lambda$F$_\lambda$   [W m$^{-2}$]',fontsize=20)
plt.tick_params(labelsize=20)
plt.xlim(xmi,xma)
plt.ylim(ymi,yma)
fig = plt.gcf()
plt.show()
