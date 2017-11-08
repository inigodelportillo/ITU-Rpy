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

# Obtain the topographic altitude of the ground station
hs = itur.topographic_altitude(lat, lon)

# Elevation angle and frequency of the link
el = 35
f = 86 * u.GHz

# Antenna parameters
D = 5 * u.m
eta = 0.65

# Unavailability probability
ps = np.logspace(-2, 0, 50)

# Compute the atmospheric parameters
T = itur.surface_mean_temperature(lat, lon)
P = 1013 * u.hPa
rho = 7.5 * u.g / u.m**3

# Initalize result vectors
Ar = np.zeros_like(ps)
Ag = np.zeros_like(ps)
Ac = np.zeros_like(ps)
As = np.zeros_like(ps)
A = np.zeros_like(ps)

# Compute each of the attenuation components
for i, p in enumerate(ps):
    Ar[i] = itur.rain_attenuation(lat, lon, f, el, hs, p).value
    Ag[i] = itur.gaseous_attenuation_slant_path(f, el, rho,
                                                  P, T, mode='approx').value
    Ac[i] = itur.cloud_attenuation(lat, lon, el, f, p).value
    As[i] = itur.scintillation_attenuation(lat, lon, f, el, p, D, eta).value

    # Compute the total attenuation
    A[i] = Ag[i] + np.sqrt((Ar[i] + Ac[i])**2 + As[i]**2)

# Print the results
plt.figure()
plt.plot(ps, A , label='Total attenuation')
plt.plot(ps, Ar, label='Rain attenuation')
plt.plot(ps, Ag, label='Gaseous attenuation')
plt.plot(ps, Ac, label='Clouds attenuation')
plt.plot(ps, As, label='Scintillation attenuation')

# Format plot
plt.legend()
plt.grid(True)
plt.xlabel('Unavailability [%]')
plt.ylabel('Attenuation [dB]')
plt.xlim([0.01, 1])
