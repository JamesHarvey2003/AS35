#!/usr/bin/env python3

# © 2020 University of Oxford - Department of Physics
# Physics Practical Course AS35 - Colour-magnitude diagrams of open clusters
# Star Detection

# Import Python Libraries
import glob, os
import sys
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from photutils import DAOStarFinder
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from astropy import wcs
import numpy as np
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore')

# EDIT the name of the cluster and the filter name
target = "NGC0663"
ref_filter = "V" # B or V


def find_stars(target, filter_name):
	filename = target + "_combined/" + target + "_" + filter_name + "_combined.fits"

	if os.path.isfile(filename) != True:
		print ("ERROR: " + filename + " does not exist")
		sys.exit()
	
	# Read the image
	hdulist = fits.open(filename)
	data = hdulist[0].data
	header = hdulist[0].header
	hdulist.close()
	
	# Compute the noise level
	mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
	
	# Find stars above 10sigma
	daofind = DAOStarFinder(fwhm=3.0, threshold=10*std)
	sources = daofind(data - median)

	print ("Found " + str(len(sources)) + " sources")
	
	xpix = sources['xcentroid']
	ypix = sources['ycentroid']
	
	# Plot the stars in the image
	pixpos = list(zip(xpix,ypix))
	apertures = CircularAperture(pixpos, r=5)
	norm = ImageNormalize(vmin=-std, vmax=20.*std, stretch=SqrtStretch())
	plt.close()
	plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
	apertures.plot(color="red", lw=1.5, alpha=0.5)
	plt.savefig(target + "_filter_" + filter_name + ".png", dpi=250)
	print ("Sources plotted in " + target + "_filter_" + filter_name + ".png")
	# Convert from pixel to sky coordinates
	w = wcs.WCS(header)

	return w.all_pix2world(np.array([xpix, ypix]).T, 0)
	

# Look for the stars
stars = find_stars(target, ref_filter)

# Save the positions of the stars in a text file
np.savetxt("stars_" + target + "_" + ref_filter + ".txt", stars)
print ("Coordinates saved in stars_" + target + "_" + ref_filter + ".txt")
