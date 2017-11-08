# -*- coding: utf-8 -*-
""" An example of how to use the itur module to compute the attenuation at a
lat-lon surface."""
import astropy.units as u
import numpy as np
import itur
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from iturutils import regular_lat_lon_grid, elevation_angle

# Elevation angle and frequency of the link
lat, lon =regular_lat_lon_grid(0.1, 0.1, lat_min=-50, lat_max=50,
                               lon_min=-30, lon_max=90)

el = elevation_angle(36000, 0, 34.5, lat, lon)
f = 86 * u.GHz

# Obtain the topographic altitude of the ground station
hs = itur.topographic_altitude(lat, lon)

# Elevation angle and frequency of the link
f = 86 * u.GHz

# Unavailability probability
p = 8

# Antenna parameters
D = 5 * u.m
eta = 0.65

# Compute the atmospheric parameters
T = itur.surface_mean_temperature(lat, lon)
P = 1013 * u.hPa
rho = itur.models.itu836.surface_water_vapour_density(lat, lon, p)

# Compute each of the attenuation components
Ar = itur.rain_attenuation(lat, lon, f, el, hs, p)
Ag = itur.gaseous_attenuation_slant_path(f, el, rho, P, T, mode='approx')
Ac = itur.cloud_attenuation(lat, lon, el, f, p)
As = itur.scintillation_attenuation(lat, lon, f, el, p, D, eta, T=T, P=P)

# Compute the total attenuation
A = Ag + np.sqrt((Ar + Ac)**2 + As**2)

# Do some plotting
fig = plt.figure(figsize=(17, 8))
ax1 = fig.add_subplot(121)
m = Basemap(ax=ax1, projection='cyl',llcrnrlat=-50,urcrnrlat=50,\
            llcrnrlon=-30,urcrnrlon=90,resolution='l')
m.drawcountries(linewidth=0.2,color='w')
m.drawcoastlines(linewidth=0.2,color='w')
#m.drawmapscale()
m.drawmeridians([-20, 0, 20, 40, 60, 80], labels=[1,0,0,1], color='w', linewidth=0.1, fontdict={'size':11})
m.drawparallels([-40, -20, 0, 20, 40], labels=[1,0,0,1], color='w', linewidth=0.1, fontdict={'size':11})
nx = int((m.xmax-m.xmin)/0.025)+1; ny = int((m.ymax-m.ymin)/0.025)+1
topodat = m.transform_scalar(np.flipud(A.value), np.arange(-30,90,0.1),
                             np.arange(-50, 50, 0.1), nx, ny)
im = m.imshow(topodat, cmap='inferno')
cbar = m.colorbar(im, location='bottom', pad="8%")
cbar.set_label('Atmospheric Attenuation [dB]')
ax1.set_title('Attenuation for f={0} GHz, unavailability p={1} %'.format(f.value,
                                                                      p))

ax2 = fig.add_subplot(122)
m = Basemap(ax=ax2, projection='cyl',llcrnrlat=-50,urcrnrlat=50,\
            llcrnrlon=-30,urcrnrlon=90,resolution='l')
m.drawcountries(linewidth=0.4,color='w')
m.drawcoastlines(linewidth=0.4,color='w')
#m.drawmapscale()
m.drawmeridians([-20, 0, 20, 40, 60, 80], labels=[1,0,0,1], color='w', linewidth=0.1, fontdict={'size':11})
m.drawparallels([-40, -20, 0, 20, 40], labels=[1,0,0,1], color='w', linewidth=0.1, fontdict={'size':11})
nx = int((m.xmax-m.xmin)/0.025)+1; ny = int((m.ymax-m.ymin)/0.025)+1
topodat = m.transform_scalar(np.flipud(A.value), np.arange(-30,90,0.1),
                             np.arange(-50, 50, 0.1), nx, ny)
im = m.contourf(lon, lat, A.value, [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55], cmap='Vega20')
cbar = m.colorbar(im, location='bottom', pad="8%")
cbar.set_label('Atmospheric Attenuation [dB]')
ax2.set_title('Regions attenuation for f={0} GHz, unavailability p={1} %'.format(f.value,

plt.savefig('att_f_{0}_unav_{1}.png'.format(f.value, p))