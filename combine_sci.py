#!/usr/bin/env python3

# © 2020 University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Combine Science Frames

# Import Python Libraries
import glob, os
import sys
import ccdproc
import copy
from ccdproc import CCDData
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from reproject import reproject_interp
import numpy as np
from scipy.signal import medfilt
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore')


# EDIT the name of the cluster and the filters you want to image
target = "NGC0663"
filters = ["B", "V"]

# Create the output directory if needed
if not os.path.exists(target + "_combined"):
	os.makedirs(target + "_combined")

for filter_name in filters:
	# Find science frames for this filter
	sci_files = glob.glob(target + "_frames/" + target + "_" + filter_name + "_*.fits")

	sci_list = []
	ra = []
	dec = []
	for sci in sci_files:
		# Read the file
		ccd = CCDData.read(sci)
		# Store the coordinates of the frame corners
		w = wcs.WCS(ccd.header)
		
		(ra1, dec1) = w.wcs_pix2world(1, 1, 1)
		(ra2, dec2) = w.wcs_pix2world(ccd.header["NAXIS1"], ccd.header["NAXIS2"], 1)
		
		ra.append(ra1)
		ra.append(ra2)
		dec.append(dec1)
		dec.append(dec2)

		sci_list.append(ccd)
	
	# Check that there is at least 1 file to be combined
	if len(sci_list) == 0:
		print ("ERROR: no frames found for target = " + target + " and filter = " + filter_name)
		sys.exit()

	# Determine the reference image for the combination
	# convert lists to numpy arrays
	ra = np.array(ra)
	dec = np.array(dec)
	# Calculate average RA and Dec of the frames
	mean_ra = 0.5*(max(ra) + min(ra))
	mean_dec = 0.5*(max(dec) + min(dec))
	
	# Create reference header
	ref_header = copy.copy(sci_list[0].header)
	dist_ra = (((min(ra) - max(ra))*np.cos(np.radians(mean_dec)))**2 + (mean_dec - mean_dec)**2)**0.5
	dist_dec = (((mean_ra - mean_ra)*np.cos(np.radians(mean_dec)))**2 + (max(dec) - min(dec))**2)**0.5

	pix_size = abs(ref_header["CD1_1"])

	npix_ra = int(dist_ra/pix_size)
	npix_dec = int(dist_dec/pix_size)

	ref_header["CRVAL1"] = mean_ra
	ref_header["CRPIX1"] = npix_ra*0.5
	ref_header["CRVAL2"] = mean_dec
	ref_header["CRPIX2"] = npix_dec*0.5
	
	ref_header["NAXIS1"] = npix_ra
	ref_header["NAXIS2"] = npix_dec

	#Reproject the frames
	exposure_map = np.zeros((npix_dec, npix_ra))
	for i, ccd in enumerate(sci_list):
		print ("projecting" + str(i+1) + "/" + str(len(sci_list)))
		ccd.data, _ = reproject_interp((ccd.data, ccd.header), ref_header)
		
		exptime, _ = reproject_interp((np.zeros_like(ccd.data) + ccd.header["EXPTIME"], ccd.header), ref_header)
		mask_exptime = (~np.isfinite(exptime)) + ~(np.isfinite(ccd.data))
		exptime[mask_exptime] = 0.
		exposure_map = exposure_map + exptime
		
		# and mask nan values
		ccd.mask = (~(np.isfinite(ccd.data)))

	# Combine all the frames
	combined_image = ccdproc.combine(sci_list, method='median', dtype="float32")
	# Save the combined frame
	hdu = combined_image.to_hdu()
	hdu[0].header["CRVAL1"] = mean_ra
	hdu[0].header["CRPIX1"] = npix_ra*0.5
	hdu[0].header["CRVAL2"] = mean_dec
	hdu[0].header["CRPIX2"] = npix_dec*0.5
	
	# Mask pixel with low integration times
	exposure_map = medfilt(exposure_map, (15, 15)) # but keep bright stars
	mask_lowexposure = (exposure_map < np.median(exposure_map)*0.5)
	
	hdu[0].data[mask_lowexposure] = np.nan
	
	hdu.writeto(target + "_combined/" + target + "_" + filter_name + "_combined.fits", clobber=True)
	
	#hdu[0].data = exposure_map
	#hdu.writeto(target + "_combined/" + target + "_" + filter_name + "_combined_expmap.fits", clobber=True)
	
	print ("Created " + target + "_combined/" + target + "_" + filter_name + "_combined.fits")
