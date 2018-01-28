# -*- coding: utf-8 -*-
""" This example shows how to compute the 'attenuation exceeded for 0.1 % of
time of the average year' vs. 'frequency' and 'elevation angle'
for a single location.

For the 'attenuation exceeded for 0.1 % of time of the average year' vs.
'frequency' case the link is assume to be a space-to-Earth link between
a ground station in Boston and a satellite in GEO orbit (slot 77W).

For the 'attenuation exceeded for 0.1 % of time of the average year' vs.
'elevation' case, the link operates at 22.5 GHz.

The receiver antenna has a 1.2 m diameter in both cases.
"""
import itur
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Ground station coordinates (Boston)
lat_GS = 42.3601
lon_GS = -71.0942

################################################
# First case: Attenuation vs. frequency        #
################################################

# Satellite coordinates (GEO, 77 W)
lat_sat = 0
lon_sat = -77
h_sat = 35786 * itur.u.km

# Compute the elevation angle between satellite and ground station
el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat, lat_GS, lon_GS)

f = 22.5 * itur.u.GHz    # Link frequency
D = 1.2 * itur.u.m       # Antenna diameters
p = 1

f = np.logspace(-0.2, 2, 100) * itur.u.GHz

Ag, Ac, Ar, As, A =\
    itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D,
                                            return_contributions=True)

# Plot the results
fig, ax = plt.subplots(1, 1)
ax.loglog(f, Ag, label='Gaseous attenuation')
ax.loglog(f, Ac, label='Cloud attenuation')
ax.loglog(f, Ar, label='Rain attenuation')
ax.loglog(f, As, label='Scintillation attenuation')
ax.loglog(f, A, label='Total atmospheric attenuation')

ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel('Atmospheric attenuation [dB]')
ax.grid(which='both', linestyle=':')
plt.legend()


################################################
# Second case: Attenuation vs. elevation angle #
################################################

f = 22.5 * itur.u.GHz
el = np.linspace(5, 90, 100)

Ag, Ac, Ar, As, A =\
    itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D,
                                            return_contributions=True)

# Plot the results
fig, ax = plt.subplots(1, 1)
ax.plot(el, Ag, label='Gaseous attenuation')
ax.plot(el, Ac, label='Cloud attenuation')
ax.plot(el, Ar, label='Rain attenuation')
ax.plot(el, As, label='Scintillation attenuation')
ax.plot(el, A, label='Total atmospheric attenuation')

ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.set_xlabel('Elevation angle [deg]')
ax.set_ylabel('Atmospheric attenuation [dB]')
ax.grid(which='both', linestyle=':')
plt.legend()
