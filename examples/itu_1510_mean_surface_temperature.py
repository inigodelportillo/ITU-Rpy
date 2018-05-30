# -*- coding: utf-8 -*-
""" This example shows how to compute the anuaml mean surface temperature
over a large region of the Earth.

This image is similar to the one plotted in page 3 of Recommendation
ITU-R P.1510-1.
"""
import itur
import matplotlib.pyplot as plt

# Set Recommendation ITU-R P.837 to version 7
itur.models.itu1510.change_version(1)

# Generate a regular grid of latitude and longitudes with 0.1 degree resolution
# for the region of interest.
lat, lon = itur.utils.regular_lat_lon_grid(resolution_lat=0.1,
                                           resolution_lon=0.1)

# Compute the surface mean temperature
T = itur.models.itu1510.surface_mean_temperature(lat, lon)

# Display the results in a map
fig = plt.figure(figsize=(16, 8))
ax = fig.add_subplot(1, 1, 1)
m = itur.utils.plot_in_map(
        T, lat, lon, cmap='jet', vmin=230, vmax=310, ax=ax,
        cbar_text='Annual mean surface temperature [K]')
