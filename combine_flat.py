#!/usr/bin/env python3

# © 2020 University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Combine Flat Field Frames

# Import Python Libraries
import glob, os
import ccdproc
from ccdproc import CCDData
from astropy import units as u
import numpy as np
import sys
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

list_file = "flat_files.txt"

# Check that the flat list exists
if os.path.isfile(list_file) != True:
	print ("ERROR: " + list_file + " does not exist")
	sys.exit()

# Check that the master bias exists
if os.path.isfile("master/master_bias.fits") != True:
	print ("ERROR: master/master_bias.fits does not exist")
	sys.exit()
	
# Read the master bias
master_bias = CCDData.read("master/master_bias.fits")


# Check that the master dark exists
if os.path.isfile("master/master_dark.fits") != True:
	print ("ERROR: master/master_dark.fits does not exist")
	sys.exit()
	
# Read the master dark
master_dark = CCDData.read("master/master_dark.fits")

# Read the flat files
flat_list = []
filter_name = None
for flat in open(list_file, "r"):
	flat = flat.strip()
	if len(flat) == 0 or flat[0] == "#":
		continue
	
	# Read the frame
	ccd = CCDData.read(flat, unit = u.adu)
	
	# Check that it is a flat field frame
	if ccd.header["IMAGETYP"].strip() != "FLAT":
		print ("ERROR: The " + flat + " file does not seem to be a flat file")
		sys.exit()
	
	if filter_name == None:
		filter_name = ccd.header["FILTER"].strip()
	else:
		if filter_name != ccd.header["FILTER"].strip():
			print ("ERROR: Creating a flat for filter " + filter_name + ", but the flat " + flat + " is for filter " + ccd.header["FILTER"].strip())
			sys.exit()
	
	# Subtract the master bias
	ccd = ccdproc.subtract_bias(ccd, master_bias)
	# and subtract the dark current scaled by the exposure time
	ccd = ccdproc.subtract_dark(ccd, master_dark, dark_exposure=master_dark.header["EXPTIME"]*u.s, data_exposure=ccd.header["EXPTIME"]*u.s, scale=True)
	# Normalize the flat field
	ccd.data = ccd.data/np.median(ccd.data)
	
	flat_list.append(ccd)
	
# Check that there is at least 1 dark to be combined
if len(flat_list) == 0:
	print ("ERROR: " + list_file + " does not contain any valid file")
	sys.exit()
	
# Combine the flats
master_flat = ccdproc.combine(flat_list, method='median', dtype="float32")

# Save the master flat
filter_name = ccd.header["FILTER"].strip()
master_flat.write("master/master_flat_" + filter_name + ".fits", clobber=True)
print ("Created master_flat_" + filter_name + ".fits")