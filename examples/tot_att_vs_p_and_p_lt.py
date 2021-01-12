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

#set up variables that will be used in calcs
lat = 33.915335
lon = -118.38257
p_arr = np.linspace(50, 0.1, 100)

#calculate the attenuations at risk = 50% and plot the reuslts
p_atts = itur.atmospheric_attenuation_slant_path(lat, lon, 20.0, 45.0,
                                               p_arr, 3.0)
plt.plot(p_arr, p_atts)

#calculate the attenuations at risk = 10% and plot the reuslts
p_r10_arr = p_lt_from_risk(lat, lon, p_arr, 10, N = 10000)
p_r10_atts = itur.atmospheric_attenuation_slant_path(lat, lon, 20.0, 45.0,
                                               p_r10_arr, 3.0)
plt.plot(p_arr, p_r10_atts, color = 'red')

#calculate the attenuations at risk = 90% and plot the reuslts
p_r90_arr = p_lt_from_risk(lat, lon, p_arr, 90, N = 10000)
p_r90_atts = itur.atmospheric_attenuation_slant_path(lat, lon, 20.0, 45.0,
                                               p_r90_arr, 3.0)
plt.plot(p_arr, p_r90_atts, color = 'green')


#Add standard plot features
plt.legend(labels = ['Risk = 50%', 'Risk = 10%', 'Risk = 90%'])
plt.ylabel('Total Atmospheric Attenuation [dB]')
plt.xlabel('Unavailablity over a period of time, p')
plt.title('Risk Analysis: Atmospheric Loss at Different Risk Values')