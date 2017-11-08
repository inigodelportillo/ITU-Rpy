# -*- coding: utf-8 -*-
"""
"""
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import itur

# Ground station location
lat = 15
lon = 25
hs = itur.topographic_altitude(lat, lon)

# Elevation angle and frequency of the link
el = 35
p = 1

# Antenna parameters
D = 5 * u.m
eta = 0.65

rho = itur.models.itu836.surface_water_vapour_density(lat, lon, p, hs)
T = itur.models.itu1510.surface_mean_temperature(lat, lon)

# Unavailability probability
fs = np.linspace(5, 200, 200)

# Initalize result vectors
Ar = np.zeros_like(fs)
Ag = np.zeros_like(fs)
Ac = np.zeros_like(fs)
As = np.zeros_like(fs)
A = np.zeros_like(fs)

# Compute each of the attenuation components
for i, f in enumerate(fs):
    a_g, a_c, a_r, a_s, a =\
        itur.total_atmospheric_attenuation_slant_path(lat, lon, el, f, p, D,
                                                      hs=hs, eta=eta, rho=rho,
                                                      T=T,
                                                      return_contributions=True)
    Ar[i] = a_r.value
    Ag[i] = a_g.value
    Ac[i] = a_c.value
    As[i] = a_s.value
    A[i] = a.value

# Print the results
plt.figure()
plt.plot(fs, A, label='Total attenuation')
plt.plot(fs, Ar, label='Rain attenuation')
plt.plot(fs, Ag, label='Gaseous attenuation')
plt.plot(fs, Ac, label='Clouds attenuation')
plt.plot(fs, As, label='Scintillation attenuation')

# Format plot
plt.legend()
plt.grid(True)
plt.xlabel('Frequency [f]')
plt.ylabel('Attenuation [dB]')
plt.yscale('log')
plt.xlim([5, 200])
