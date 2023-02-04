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
el = 60  # Elevation angle equal to 60 degrees
f = 22.5 * itur.u.GHz  # Frequency equal to 22.5 GHz
D = 1 * itur.u.m  # Receiver antenna diameter of 1 m
p = 0.1  # We compute values exceeded during 0.1 % of the average
# year

# Compute atmospheric parameters
hs = itur.topographic_altitude(lat, lon)
T = itur.surface_mean_temperature(lat, lon)
P = itur.models.itu835.pressure(lat, hs)
rho_p = itur.surface_water_vapour_density(lat, lon, p, hs)
rho_sa = itur.models.itu835.water_vapour_density(lat, hs)
T_sa = itur.models.itu835.temperature(lat, hs)
V = itur.models.itu836.total_water_vapour_content(lat, lon, p, hs)

print(
    f"The ITU recommendations predict the following values for the point located at coordinates ({lat}, {lon})"
)

print(
    f"  - Height above the sea level                  [ITU-R P.1511]  {hs.to(itur.u.m):.1f}"
)
T_C = T.to(itur.u.Celsius, equivalencies=itur.u.temperature())
print(f"  - Surface mean temperature                    [ITU-R P.1510]  {T_C:.1f}")
print(f"  - Surface pressure                            [ITU-R P.835]   {P:.1f}")
T_sa_C = T_sa.to(itur.u.Celsius, equivalencies=itur.u.temperature())
print(f"  - Standard surface temperature                [ITU-R P.835]   {T_sa_C:.1f}")
print(f"  - Standard water vapour density               [ITU-R P.835]   {rho_sa:.1f}")
print(f"  - Water vapor density (p={p}%)                [ITU-R P.836]   {rho_p:.1f}")
print(f"  - Total water vapour content (p={p}%)         [ITU-R P.836]   {V:.1f}")

# Compute rain and cloud-related parameters
R_prob = itur.models.itu618.rain_attenuation_probability(lat, lon, el, hs)
R_pct = itur.models.itu837.rainfall_probability(lat, lon)
R001 = itur.models.itu837.rainfall_rate(lat, lon, p)
h_0 = itur.models.itu839.isoterm_0(lat, lon)
h_rain = itur.models.itu839.rain_height(lat, lon)
L_red = itur.models.itu840.columnar_content_reduced_liquid(lat, lon, p)
A_w = itur.models.itu676.zenit_water_vapour_attenuation(lat, lon, p, f, h=hs)

print(f"  - Rain attenuation probability                [ITU-R P.618]   {R_prob:.1f}")
print(f"  - Rain percentage probability                 [ITU-R P.837]   {R_pct:.1f}")
print(f"  - Rainfall rate exceeded for p={p}%           [ITU-R P.837]   {R001:.1f}")
print(f"  - 0 degree C isotherm height                  [ITU-R P.839]   {h_0:.1f}")
print(f"  - Rain height                                 [ITU-R P.839]   {h_rain:.1f}")
print(f"  - Columnar content of reduced liquid (p={p}%) [ITU-R P.840]   {L_red:.1f}")
print(f"  - Zenit water vapour attenuation (p={p}%)     [ITU-R P.676]   {A_w:.1f}")

# Compute attenuation values
A_g = itur.gaseous_attenuation_slant_path(f, el, rho_p, P, T)
A_r = itur.rain_attenuation(lat, lon, f, el, hs=hs, p=p)
A_c = itur.cloud_attenuation(lat, lon, el, f, p)
A_s = itur.scintillation_attenuation(lat, lon, f, el, p, D)
A_t = itur.atmospheric_attenuation_slant_path(lat, lon, f, el, p, D)

print(
    f"\n\nAttenuation values exceeded for p={p}% of the average year "
    f"for a link with el={el} deg, f={f}, \nD={D} and "
    f"receiver ground station located at coordinates ({lat}, {lon})"
)

print(f"  - Rain attenuation                            [ITU-R P.618]   {A_r:.1f}")
print(f"  - Gaseous attenuation                         [ITU-R P.676]   {A_g:.1f}")
print(f"  - Clouds attenuation                          [ITU-R P.840]   {A_c:.1f}")
print(f"  - Scintillation attenuation                   [ITU-R P.618]   {A_s:.1f}")
print(f"  - Total atmospheric attenuation               [ITU-R P.618]   {A_t:.1f}")
