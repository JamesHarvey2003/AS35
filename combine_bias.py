#!/usr/bin/env python3

# © 2020 University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Master Bias Generator

# Import Python Libraries
import ccdproc
from ccdproc import CCDData
from astropy import units as u
import numpy as np
import os
import sys

list_file = "bias_files.txt"

# Check that the bias_files exists
if os.path.isfile(list_file) != True:
    print ("ERROR: " + list_file + " does not exist")
    sys.exit()

# Create the master calibration directory if needed
if not os.path.exists("master"):
    os.makedirs("master")

# Read the bias files
bias_list = []
for bias in open(list_file, "r"):
    bias = bias.strip()
    if len(bias) == 0 or bias[0] == "#":
        continue

    if os.path.isfile(bias) != True:
        print ("ERROR: The " + bias + " file listed in " + list_file + " does not exist")
        sys.exit()

    # Read the frame
    ccd = CCDData.read(bias, unit = u.adu)

    # Check that it is a bias frame
    if ccd.header["IMAGETYP"] != "Bias Frame":
        print ("ERROR: The " + bias + " file does not seem to be a bias file")
        sys.exit()

    bias_list.append(ccd)

# Check that there is at least 1 bias to be combined
if len(bias_list) == 0:
    print ("ERROR: " + list_file + " does not contain any valid file")
    sys.exit()

# Combine the bias
master_bias = ccdproc.combine(bias_list, method='average', dtype="float32")

# Calculate the mean and the standard deviation in the area delimited by (x1,y1) (x2,y2)
# EDIT the values of x1, x2, y1, and y2
x1 = 400
x2 = 700
y1 = 141
y2 = 142
print ("mean = ", np.mean(master_bias.data[y1:y2, x1:x2]))
print ("std = ", np.std(master_bias.data[y1:y2, x1:x2]))

# Save the master bias
master_bias.write("master/master_bias.fits", overwrite=True)
print ("Created master_bias.fits")
