# -*- coding: utf-8 -*-
""" An example of how to use the itur module to compute the attenuation at a
single ground station."""
import astropy.units as u
import numpy as np
import itur

# Ground station location
lat = 47.9358    # Latitude of the Ground Station [deg N]
lon = -124.5619  # Longitude of the Ground Station [deg E]
# hs = 0           # Altitude of the Ground Station [m]

# Obtain the topographic altitude of the ground station
hs = itur.topographic_altitude(lat, lon)

# Elevation angle and frequency of the link
el = 90
f = 86 * u.GHz

# Unavailability probability
p = 0.1

# Antenna parameters
D = 5 * u.m
eta = 0.65

# Compute the atmospheric parameters
T = itur.surface_mean_temperature(lat, lon)
P = 1013 * u.hPa
rho = 7.5 * u.g / u.m**3

# Compute each of the attenuation components
Ar = itur.rain_attenuation(lat, lon, f, el, hs, p)
Ag = itur.gaseous_attenuation_slant_path(f, el, rho, P, T, mode='approx')
Ac = itur.cloud_attenuation(lat, lon, el, f, p)
As = itur.scintillation_attenuation(lat, lon, f, el, p, D, eta, T=T, P=P)

#As = itur.scintillation_attenuation(lat, lon, f, el, p, D, eta, T, H, P, hL)

# Compute the total attenuation
A = Ag + np.sqrt((Ar + Ac)**2 + As**2)

# Print the results
print(('Attenuation for {0:0.2f} unavailability for a \n' +
       'ground station in lat  {1:.2f} and lon {2:.2f} \n'
       '==========================================\n' +
       'Rain attenuation {3:0.2f} \n' +
      'Gaseous attenuation {4:0.2f} \n' +
      'Scintillation attenuation {5:0.2f} \n' +
      'Clouds attenuation {6:0.2f} \n\n' +
      'Total attenuation {7:0.2f}').format(p, lat, lon, Ar, Ag, As, Ac, A))
