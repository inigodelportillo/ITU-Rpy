# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import itur

# Ground station location
lat = 15
lon = 25
hs = itur.topographic_altitude(lat, lon)

# Elevation angle and unavailability of the link
el = 35   # Elevation angle = 35 deg
p = 1      # Unavailability = 1 %

# Antenna diameter and efficiency for scintillation computation
D = 5 * u.m
eta = 0.65  

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

	# Using return_contributions we can get the atmospheric 
	# attenuation due to gases, clouds, rain, and scintillation
    a_g, a_c, a_r, a_s, a =\
        itur.total_atmospheric_attenuation_slant_path(lat, lon, el, f, p, D,
                                                      hs=hs, eta=eta, 
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
