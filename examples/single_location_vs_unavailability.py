# -*- coding: utf-8 -*-
""" This example shows how to compute the attenuation vs. percentage of time
of the average year that values are exceeded for a single location.

The link is a space-to-Earth link between a ground station in Boston and a
satellite in GEO orbit (slot 77W). The link operates at 22.5 GHz and the
receiver antenna has a 1.2 m diameter.
"""
import itur
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Ground station coordinates (Boston)
lat_GS = 42.3601
lon_GS = -71.0942

# Satellite coordinates (GEO, 77 W)
lat_sat = 0
lon_sat = -77
h_sat = 35786 * itur.u.km

# Compute the elevation angle between satellite and ground station
el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat, lat_GS, lon_GS)

f = 22.5 * itur.u.GHz    # Link frequency
D = 1.2 * itur.u.m       # Antenna diameters

# Define unavailabilities vector in logarithmic scale
p = np.logspace(-1.5, 1.5, 100)

# Compute the attenuation values for different unavailabilities.
# The unavailability is the only parameter that is not vectorized in ITU-Rpy
A_g, A_c, A_r, A_s, A_t = \
    itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D,
                                            return_contributions=True)

# Plot the results using matplotlib
f, ax = plt.subplots(1, 1)
ax.semilogx(p, A_g.value, label='Gaseous attenuation')
ax.semilogx(p, A_c.value, label='Cloud attenuation')
ax.semilogx(p, A_r.value, label='Rain attenuation')
ax.semilogx(p, A_s.value, label='Scintillation attenuation')
ax.semilogx(p, A_t.value, label='Total atmospheric attenuation')

ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xlabel('Percentage of time attenuation value is exceeded [%]')
ax.set_ylabel('Attenuation [dB]')
ax.grid(which='both', linestyle=':')
plt.legend()
