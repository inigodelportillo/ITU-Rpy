# -*- coding: utf-8 -*-
""" This example shows how to compute the rainfall-rate (mm/hr) exceeded
for 0.01 % of the time of the average year over a large region of the Earth.

This image is similar to the one plotted in page 5 of Recommendation
ITU-R P.837-7.
"""
import itur
import matplotlib.pyplot as plt

# Set Recommendation ITU-R P.837 to version 7
itur.models.itu837.change_version(7)

# Generate a regular grid of latitude and longitudes with 0.1 degree resolution
# for the region of interest.
lat, lon = itur.utils.regular_lat_lon_grid(resolution_lat=0.1,
                                           resolution_lon=0.1)

# Compute the rainfall rate exceeded for 0.01 % of the time.
p = 0.01
R001 = itur.models.itu837.rainfall_rate(lat, lon, p)

# Display the results in a map
fig = plt.figure(figsize=(16, 8))
ax = fig.add_subplot(1, 1, 1)
m = itur.utils.plot_in_map(
    R001, lat, lon, cmap='jet', vmin=0, vmax=90, ax=ax,
    cbar_text='Rainfall rate exceeded for 0.01% of an average year [mm/hr]')
