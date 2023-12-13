#!/usr/bin/env python3

# © 2020 University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Combine Dark Frames

# Import Python Libraries
import ccdproc
from ccdproc import CCDData
from astropy import units as u
from astropy.stats import sigma_clipped_stats
import numpy as np
import os.path
import sys

list_file = "dark_files.txt"

# Check that file exists
if os.path.isfile(list_file) != True:
	print ("ERROR: " + list_file + " does not exist")
	sys.exit()

# Read the dark files
dark_list = []
for dark in open(list_file, "r"):
	dark = dark.strip()
	if len(dark) == 0 or dark[0] == "#":
		continue

	if os.path.isfile(dark) != True:
		print ("ERROR: The " + dark + " file listed in " + list_file + " does not exist")
		sys.exit()

	# Read the frame
	ccd = CCDData.read(dark, unit = u.adu)

	# Check that it is a dark frame
	if ccd.header["IMAGETYP"] != "Dark Frame":
		print ("ERROR: The " + dark + " file does not seem to be a dark file")
		sys.exit()

	dark_list.append(ccd)

# Check that there is at least 1 dark to be combined
if len(dark_list) == 0:
	print ("ERROR: " + list_file + " does not contain any valid file")
	sys.exit()

# Combine the dark
master_dark = ccdproc.combine(dark_list, method='average', dtype="float32")

# Calculate the mean and the standard deviation in the area delimited by (x1,y1) (x2,y2)
# EDIT the values of x1, x2, y1, and y2
x1 = 400
x2 = 700
y1 = 141
y2 = 142
mean, median, std = sigma_clipped_stats(master_dark.data[y1:y2,x1:x2], sigma=3.0, maxiters=5)

print ("mean = ", mean)
print ("std = ", std)

# Read the exposure time from the FITS header
print ("EXPTIME = ", master_dark.header["EXPTIME"])
