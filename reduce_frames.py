#!/usr/bin/env python3

# © 2020 University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Reduction of science frames

# Import Python Libraries
import glob, os
import sys
import ccdproc
from ccdproc import CCDData
from astropy import units as u
from astropy.stats import sigma_clipped_stats
import numpy as np
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore')

# EDIT the name of the cluster
target = "NGC6939"
##

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

# Create the output directory if needed
if not os.path.exists(target + "_frames"):
	os.makedirs(target + "_frames")

# Find raw science frames
sci_files = glob.glob(target + "/" + target + "*")
# and reduce
for i, sci in enumerate(sci_files):
	# Read the science frame
	ccd = CCDData.read(sci, unit = u.adu)

	# Mask saturated pixels
	mask_saturated = (ccd.data > 50000)
	ccd.data = np.array(ccd.data, dtype=np.float32)
	ccd.data[mask_saturated] = np.nan

	# Load the appropriate flat field for this frame
	filter_name = ccd.header["FILTER"].strip()

	# Check that the master dark exists
	if os.path.isfile("master/master_flat_" + filter_name + ".fits") != True:
		print ("ERROR: master/master_flat_" + filter_name + ".fits")
		sys.exit()

	master_flat =  CCDData.read("master/master_flat_" + filter_name + ".fits")

	# Subtract bias
	ccd = ccdproc.subtract_bias(ccd, master_bias)

	# Subtract dark current
	ccd = ccdproc.subtract_dark(ccd, master_dark, dark_exposure=master_dark.header["EXPTIME"]*u.s, data_exposure=ccd.header["EXPTIME"]*u.s, scale=True)

	# Divide by flat
	ccd = ccdproc.flat_correct(ccd, master_flat, min_value=0.5)

	# Subtract global sky background
	mean, background, std = sigma_clipped_stats(ccd.data, sigma=3.0, maxiters=5)
	ccd.data = ccd.data - background

	# Divide by the exposure time
	ccd.data = ccd.data/ccd.header["EXPTIME"]
	ccd.unit = u.adu/u.s

	# Add keywords to the header
	ccd.header['SKY'] = background
	ccd.header['RAWFILE'] = sci

	# Save the calibrated frame
	ccd.write(target + "_frames/" + target + "_" + filter_name + "_" + '{:04d}'.format(i)  + ".fits" , clobber=True)
	print ("Created " + target + "_frames/" + target + "_" + filter_name + "_" + '{:04d}'.format(i)  + ".fits")
