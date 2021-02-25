# -*- coding: utf-8 -*-
"""
This script produces a plot that displays how changes in risk affects the 
final atmospheric attenuation value that is prodced by ITU-Rpy.
The script examines the risk values of 10%, 50%, and 90%.

Created on Thu Dec 31 11:44:27 2020

@author: MAW32652
"""
import itur
from itur.models.itu678 import p_lt_from_risk
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

#set up variables that will be used in calcs
lat = -22.8778
lon = -47.2121
p_arr = np.linspace(10, 0.1, 100)
av_arr = 100 - p_arr
#av_arr = np.flip(av_arr)


#calculate the attenuations at risk = 10% and plot the reuslts
p_r10_arr = p_lt_from_risk(lat, lon, p_arr, 10, N = 10000)
p_r10_atts = itur.atmospheric_attenuation_slant_path(lat, lon, 20.0, 45.0,
                                               p_r10_arr, 3.0)
ax.plot(av_arr, p_r10_atts, color = 'red')

#calculate the attenuations at risk = 50% and plot the reuslts
p_atts = itur.atmospheric_attenuation_slant_path(lat, lon, 20.0, 45.0,
                                               p_arr, 3.0)
ax.plot(av_arr, p_atts)

#calculate the attenuations at risk = 90% and plot the reuslts
p_r90_arr = p_lt_from_risk(lat, lon, p_arr, 90, N = 10000)
p_r90_atts = itur.atmospheric_attenuation_slant_path(lat, lon, 20.0, 45.0,
                                               p_r90_arr, 3.0)
ax.plot(av_arr, p_r90_atts, color = 'green')


#Add standard plot features
textstr = "Location = O3b ground station in Hortolandia, Brazil\n\
Frequency = 20 Ghz\nElevation Angle = 45 degrees\nAntenna Diameter = 3 m\n\
Polarization Angle  = 45 Degrees"
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform = ax.transAxes, fontsize = 8, verticalalignment = 'top', bbox = props)
ax.legend(['Risk = 10%', 'Risk = 50%', 'Risk = 90%'], loc = 'center left')
plt.ylabel('Total Atmospheric Attenuation [dB]')
plt.xlabel('Availability Over Any One Year, [0-100%]')
plt.title('Risk Analysis: Atmospheric Loss at Different Risk Values')