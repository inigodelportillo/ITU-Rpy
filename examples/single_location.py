# -*- coding: utf-8 -*-
"""
In this example we compute different atmospheric parameters for a single ground
station.

The ground station is located in Boston (GPS coordinates 41.36, -71.05) and
the link operates at a frequency of 22.5 GHz, elevation angle of 60 degrees,
and with a receiver antenna of 1 m diameter.

Values exceeded during 0.1 % of the average year are reported.
"""

import itur

# Location of the receiver ground stations
lat = 41.39
lon = -71.05

# Link parameters
el = 60                # Elevation angle equal to 60 degrees
f = 22.5 * itur.u.GHz  # Frequency equal to 22.5 GHz
D = 1 * itur.u.m       # Receiver antenna diameter of 1 m
p = 0.1                # We compute values exceeded during 0.1 % of the average
                       # year

# Compute atmospheric parameters
hs = itur.topographic_altitude(lat, lon)
T = itur.surface_mean_temperature(lat, lon)
P = itur.models.itu835.pressure(lat, hs)
rho_p = itur.surface_water_vapour_density(lat, lon, p, hs)
rho_sa = itur.models.itu835.water_vapour_density(lat, hs)
T_sa = itur.models.itu835.temperature(lat, hs)
V = itur.models.itu836.total_water_vapour_content(lat, lon, p, hs)

print(('The ITU recommendations predict the following values for the point '
       'located at coordinates ({0}, {1})')
      .format(lat, lon))

print('  - Height above the sea level                  [ITU-R P.1511]  {0:.1f}'
      .format(hs.to(itur.u.m)))
print('  - Surface mean temperature                    [ITU-R P.1510]  {0:.1f}'
      .format(T.to(itur.u.Celsius, equivalencies=itur.u.temperature())))
print('  - Surface pressure                            [ITU-R P.835]   {0:.1f}'
      .format(P))
print('  - Standard surface temperature                [ITU-R P.835]   {0:.1f}'
      .format(T_sa.to(itur.u.Celsius, equivalencies=itur.u.temperature())))
print('  - Standard water vapour density               [ITU-R P.835]   {0:.1f}'
      .format(rho_sa))
print('  - Water vapor density (p={0}%)                [ITU-R P.836]   {1:.1f}'
      .format(p, rho_p))
print('  - Total water vapour content (p={0}%)         [ITU-R P.836]   {1:.1f}'
      .format(p, V))

# Compute rain and cloud-related parameters
R_prob = itur.models.itu618.rain_attenuation_probability(lat, lon, el, hs)
R_pct_prob = itur.models.itu837.rainfall_probability(lat, lon)
R001 = itur.models.itu837.rainfall_rate(lat, lon, p)
h_0 = itur.models.itu839.isoterm_0(lat, lon)
h_rain = itur.models.itu839.rain_height(lat, lon)
L_red = itur.models.itu840.columnar_content_reduced_liquid(lat, lon, p)
A_w = itur.models.itu676.zenit_water_vapour_attenuation(lat, lon, p, f, h=hs)

print('  - Rain attenuation probability                [ITU-R P.618]   {0:.1f}'
      .format(R_prob))
print('  - Rain percentage probability                 [ITU-R P.837]   {0:.1f}'
      .format(R_pct_prob))
print('  - Rainfall rate exceeded for p={0}%           [ITU-R P.837]   {1:.1f}'
      .format(p, R001))
print('  - 0 degree C isotherm height                  [ITU-R P.839]   {0:.1f}'
      .format(h_0))
print('  - Rain height                                 [ITU-R P.839]   {0:.1f}'
      .format(h_rain))
print('  - Columnar content of reduced liquid (p={0}%) [ITU-R P.840]   {1:.1f}'
      .format(p, L_red))
print('  - Zenit water vapour attenuation (p={0}%)     [ITU-R P.676]   {1:.1f}'
      .format(p, A_w))

# Compute attenuation values
A_g = itur.gaseous_attenuation_slant_path(f, el, rho_p, P, T)
A_r = itur.rain_attenuation(lat, lon, f, el, hs=hs, p=p)
A_c = itur.cloud_attenuation(lat, lon, el, f, p)
A_s = itur.scintillation_attenuation(lat, lon, f, el, p, D)
A_t = itur.atmospheric_attenuation_slant_path(lat, lon, f, el, p, D)

print(('\n\nAttenuation values exceeded for p={0}% of the average year '
       'for a link with el={1} deg, f={2}, \nD={3} and '
       'receiver ground station located at coordinates ({4}, {5})')
      .format(p, el, f, D, lat, lon))

print('  - Rain attenuation                            [ITU-R P.618]   {0:.1f}'
      .format(A_r))
print('  - Gaseous attenuation                         [ITU-R P.676]   {0:.1f}'
      .format(A_g))
print('  - Clouds attenuation                          [ITU-R P.840]   {0:.1f}'
      .format(A_c))
print('  - Scintillation attenuation                   [ITU-R P.618]   {0:.1f}'
      .format(A_s))
print('  - Total atmospheric attenuation               [ITU-R P.618]   {0:.1f}'
      .format(A_t))
