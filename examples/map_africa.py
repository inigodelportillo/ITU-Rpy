# -*- coding: utf-8 -*-
""" This example shows how to compute the atmospheric attenuation exceeded
for 0.1 % of the time of the average year over a large region of the Earth.

It is assumed that the satellite is located in geostationary orbit, at the
4 E slot, and the link operates at 22.5 GHz with receiver-dishes of 1.2 m
diameter.

The satellite covers Africa, Europe and the Middle East.

Finally, we also plot the surface mean temperature distribution to illustrate
that other variables can also be computed using vectorized operations.
"""
import itur

# Generate a regular grid of latitude and longitudes with 0.1 degree resolution
# for the region of interest.
lat, lon = itur.utils.regular_lat_lon_grid(lat_max=60,
                                           lat_min=-60,
                                           lon_max=65,
                                           lon_min=-35,
                                           resolution_lon=0.1,
                                           resolution_lat=0.1)

# Satellite coordinates (GEO, 4 E)
lat_sat = 0
lon_sat = 4
h_sat = 35786 * itur.u.km

# Compute the elevation angle between satellite and ground stations
el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat, lat, lon)

# Set the link parameters
f = 22.5 * itur.u.GHz    # Link frequency
D = 1.2 * itur.u.m       # Antenna diameters
p = 0.1                  # Unavailability (Values exceeded 0.1% of time)

# Compute the atmospheric attenuation
Att = itur.atmospheric_attenuation_slant_path(lat, lon, f, el, p, D)

# Plot the results
m = itur.utils.plot_in_map(Att, lat, lon,
                           cbar_text='Atmospheric attenuation [dB]',
                           cmap='magma')

# Plot the satellite location
m.scatter(lon_sat, lat_sat, c='white', s=20)

# Now we show the surface mean temperature distribution
T = itur.surface_mean_temperature(lat, lon)\
    .to(itur.u.Celsius, equivalencies=itur.u.temperature())
m = itur.utils.plot_in_map(T, lat, lon,
                           cbar_text='Surface mean temperature [C]',
                           cmap='RdBu_r')
