#!/usr/bin/env python
# -*- coding: utf8 -*-

import pyregion 
import numpy as np
import astropy.io.fits as fits
import pyregion._region_filter as filter

modelname="../data/oa30_i0-zip/t5_p0_q0_oa30_R10_Mcl0.97_i0_total.fits"
hdulist =fits.open(modelname)
header = hdulist[0].header
cube = hdulist[0].data[:,:,:]
nx = header ["NAXIS1"]
ny = header ["NAXIS2"]
nl = header ["NAXIS3"]

source =pyregion.open("regions.reg")
nreg = len(source)
size = (ny, nx)
tot_sed = np.zeros((nl,nreg))

for k in range(nreg):
    print "Region ",k
    source = pyregion.open("regions.reg")
    del source[0:k]
    del source[k+1:nreg]
    mask = source.get_mask(shape=size)
    sflx = 0.0
    for m in range(nl):
        for i in range(nx):
            for n in range(ny):
                if mask[n,i] == True:
                    tot_sed[m,k] = tot_sed[m,k] + cube[m,n,i]

"""
fmt = "%12.4e %12.4e %12.4e %12.4e %12.4e \n"
ff = open("sed.dat",'w')
ff.write("reg1 reg2 reg3 reg4 tot \n")
for v in zip(sed1, sed2, sed3, sed4, sedt):
    ff.write(fmt%v)
ff.close()
"""