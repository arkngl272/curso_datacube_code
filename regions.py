
import numpy as np
import astropy.io.fits as fits 

import pyregion as filter

modelname="../data/oa30_i0-zip/t5_p0_q0_oa30_R10_Mcl0.97_i0_total.fits"
hdulist =fits.open(modelname)
header = hdulist[0].header
cube = hdulist[0].data[:,:,:]
nx = header ["NAXIS1"]
ny = header ["NAXIS2"]
nl = header ["NAXIS3"]



