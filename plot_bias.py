#!/usr/bin/env python3

# © 2020 University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Bias Plotter

# Import Python Libraries
import ccdproc
from ccdproc import CCDData
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt

# Reads the bias file
filename = "Calibration/Bias-S008-R002-C001-B2.fts"
ccd = CCDData.read(filename, unit = u.adu)

# Calculates the mean value averaging along the y and x axis
bias_x = np.mean(ccd.data, axis=0)
bias_y = np.mean(ccd.data, axis=1)

# Save the results as text ascii files
np.savetxt("bias_x.txt", bias_x)
np.savetxt("bias_y.txt", bias_y)

# and as figures
plt.close()
plt.xlabel('x pixel')
plt.ylabel('value')
plt.scatter(np.arange(len(bias_x)), bias_x)
plt.savefig("bias_x.png")

plt.close()
plt.xlabel('y pixel')
plt.ylabel('value')
plt.scatter(np.arange(len(bias_y)), bias_y)
plt.savefig("bias_y.png")