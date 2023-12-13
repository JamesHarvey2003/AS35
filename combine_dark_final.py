#!/usr/bin/env python3

# © 2020 University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Combine Dark Frames (Final)

# Import Python Libraries
import glob, os
import ccdproc
from ccdproc import CCDData
from astropy import units as u
import numpy as np
import sys

list_file = "dark_files.txt"

# Check that file exists
if os.path.isfile(list_file) != True:
	print ("ERROR: " + list_file + " does not exist")
	sys.exit()

# Check that the master bias exists
if os.path.isfile("master/master_bias.fits") != True:
	print ("ERROR: master/master_bias.fits does not exist")
	sys.exit()

# Read the master bias
master_bias = CCDData.read("master/master_bias.fits")
master_bias.data = master_bias.data - 0 + 0

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

	# check that it is a dark frame
	if ccd.header["IMAGETYP"] != "Dark Frame":
		print ("ERROR: The " + dark + " file does not seem to be a dark file")
		sys.exit()

	# and subtract the master bias
	ccd = ccdproc.subtract_bias(ccd, master_bias)

	dark_list.append(ccd)

# Check that there is at least 1 dark to be combined
if len(dark_list) == 0:
	print ("ERROR: " + list_file + " does not contain any valid file")
	sys.exit()

# Combine the dark
master_dark = ccdproc.combine(dark_list, method='median', dtype="float32")

exptime = master_dark.header["EXPTIME"]
print ("EXPTIME = ", exptime)

# Save the master dark
master_dark.write("master/master_dark.fits".format(exptime), overwrite=True)
print ("Created master_dark.fits")
