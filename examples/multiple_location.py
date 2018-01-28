# -*- coding: utf-8 -*-
""" This example shows how to compute the atmospheric attenuation exceeded
for 0.1 % of the time for multiple ground stations.

It is assumed that the satellite is located in geostationary orbit, at the
77 W slot, and the link operates at 22.5 GHz with receiver-dishes of 1.2 m
diameter.

Finally, we also plot the surface mean temperature distribution to illustrate
that other variables can also be computed using vectorized operations.
"""
import itur
import numpy as np
import matplotlib.pyplot as plt

# Obtain the coordinates of the different cities
cities = {'Boston': (42.36, -71.06),
          'New York': (40.71, -74.01),
          'Los Angeles': (34.05, -118.24),
          'Denver': (39.74, -104.99),
          'Las Vegas': (36.20, -115.14),
          'Seattle': (47.61, -122.33),
          'Washington DC': (38.91, -77.04)}

lat = [coords[0] for coords in cities.values()]
lon = [coords[1] for coords in cities.values()]

# Satellite coordinates (GEO, 4 E)
lat_sat = 0
lon_sat = -77
h_sat = 35786 * itur.u.km

# Compute the elevation angle between satellite and ground stations
el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat, lat, lon)

# Set the link parameters
f = 22.5 * itur.u.GHz    # Link frequency
D = 1.2 * itur.u.m       # Antenna diameters
p = 0.1                  # Unavailability (Values exceeded 0.1% of time)

# Compute the atmospheric attenuation
Ag, Ac, Ar, As, Att = itur.atmospheric_attenuation_slant_path(
    lat, lon, f, el, p, D, return_contributions=True)

# Plot the results
city_idx = np.arange(len(cities))
width = 0.15

fig, ax = plt.subplots(1, 1)
ax.bar(city_idx, Att.value, 0.6, label='Total atmospheric Attenuation')
ax.bar(city_idx - 1.5 * width, Ar.value, width, label='Rain attenuation')
ax.bar(city_idx - 0.5 * width, Ag.value, width, label='Gaseous attenuation')
ax.bar(city_idx + 0.5 * width, Ac.value, width, label='Clouds attenuation')
ax.bar(city_idx + 1.5 * width, As.value, width,
       label='Scintillation attenuation')

# Set the labels
ticks = ax.set_xticklabels([''] + list(cities.keys()))
for t in ticks:
    t.set_rotation(45)
ax.set_ylabel('Atmospheric attenuation exceeded for 0.1% [dB]')

# Format image
ax.yaxis.grid(which='both', linestyle=':')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=2)
plt.tight_layout(rect=(0, 0, 1, 0.85))
