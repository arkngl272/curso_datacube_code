#!/usr/bin/env python
# -*- coding: utf8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

wavel = [70., 160., 250., 350., 500.]
ofluxes = [316.05, 374.66, 165.91, 76.98, 30.58]
errors =  [30.3, 61.66, 27.5, 12.0, 5.0]
yma = max(ofluxes)*2
ymi = yma/80.
xma = max(wavel)*1.5
xmi = min(wavel)*0.7
plt.xscale('log')
plt.yscale('log')
plt.plot(wavel, ofluxes, 'ro', color='red')
plt.errorbar(wavel, ofluxes, yerr=errors, fmt='o')
plt.xlabel(r'$\lambda$  [$\mu$ m]',fontsize=16)
plt.ylabel(r'F$_\nu$  [Jy]',fontsize=16)
plt.tick_params(labelsize=16)
plt.xlim(xmi,xma)
plt.ylim(ymi,yma)
plt.show()

