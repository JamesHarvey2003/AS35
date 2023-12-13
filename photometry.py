#!/usr/bin/env python3

# © 2020 University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Aperture Photometry

# Import Python Libraries
import glob, os
import sys
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from photutils import DAOStarFinder
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from astropy import wcs
import numpy as np
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore')



# EDIT the name of the cluster, the file the positions of the stars, the zeropoints and the aperture radius
target = "NGC0663"
stars_file = "stars_NGC0663_B.txt"
zeropoint_b = 19.51
zeropoint_v = 19.75
aperture_radius = 7.2

def do_photometry(target, filter_name, stars):
	filename = target + "_combined/" + target + "_" + filter_name + "_combined.fits"
	
	if os.path.isfile(filename) != True:
		print ("ERROR: " + filename + " does not exist")
		sys.exit()
	
	# Read the image
	hdulist = fits.open(filename)
	data = hdulist[0].data
	header = hdulist[0].header
	hdulist.close()
	
	# Convert from sky coordinates to pixel
	w = wcs.WCS(header)
	xpix,ypix = w.all_world2pix(stars, 0).T
	pixel_coordinates = list(zip(xpix,ypix))

	# Generate the annulus aperture for the local background
	annulus_apertures = CircularAnnulus(pixel_coordinates, r_in=6., r_out=12.)
	annulus_mask = annulus_apertures.to_mask()
	# and calculate the background level per pixel
	background_median = np.zeros(len(annulus_mask))
	for i, mask in enumerate(annulus_mask):
		data_cutout_aper = mask.cutout(data)
		background_median[i] = np.nanmedian(data_cutout_aper[data_cutout_aper != 0])

	# Generate the circular apertures
	apertures = CircularAperture(pixel_coordinates, r=aperture_radius)
	# and do the photometry
	fluxes = aperture_photometry(data, apertures)
	
	# Ignore the fluxes of stars too close the image borders
	mask_size = 15
	mask = (xpix < mask_size) | (xpix > data.shape[1] - mask_size) | (ypix < mask_size) | (ypix > data.shape[0] - mask_size)
	fluxes["aperture_sum"][mask] = np.nan
	
	# Return the background subtracted fluxes
	return fluxes["aperture_sum"] - background_median*apertures.area


if os.path.isfile(stars_file) != True:
	print ("ERROR: " + stars_file + " does not exist")
	sys.exit()
	
# Load the location of the stars
stars = np.loadtxt(stars_file)

# Measure the fluxes in counts/s in the B and V images
flux_b = do_photometry(target, "B", stars)
flux_v = do_photometry(target, "V", stars)

# Convert from fluxes to magnitudes using the provided zeropoint
mag_b = zeropoint_b -2.5*np.log10(flux_b)
mag_v = zeropoint_v -2.5*np.log10(flux_v)

# Save the positions and fluxes
mask = ((np.isfinite(mag_b) & np.isfinite(mag_v)))
data = np.array([stars[mask,0], stars[mask,1], mag_b[mask], mag_v[mask]]).T
np.savetxt("mag_" + target + ".txt", data, fmt=['%le','%le','%7.3f','%7.3f'], header="RA Dec magB magV")
print ("Magnitudes saved in mag_" + target + ".txt")
