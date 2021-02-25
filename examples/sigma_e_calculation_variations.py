# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 11:08:56 2020

@author: MAW32652
"""


from itur.models.itu678 import variance_of_estimation_calc
import numpy as np
import matplotlib.pyplot as plt

#create a numpy list of the correct probability inputs
pListL = np.arange(1, 51, 1)
pListL = np.flip(pListL)
toAdd = np.arange(.1, 1, .1)
toAdd = np.flip(toAdd)
pListL = np.append(pListL, toAdd)
pListL = np.append(pListL, 0.01)
pListL = np.append(pListL, 0.001)
pListL = np.append(pListL, .0005)

NList = np.array([10000, 1000, 100, 525960])

outputL = []
for i in NList:
    sigmas = variance_of_estimation_calc(pListL / 100, N = i)
    outputL.append(sigmas)

plt.plot(outputL[3], pListL, color = 'black')    
plt.plot(outputL[0], pListL, 'r--')
plt.plot(outputL[1], pListL, 'b--')
plt.plot(outputL[2], pListL, color = 'green')



plt.legend(labels = ['525960', '10000', '1000', '100'])
plt.ylabel('Probability of Exceedance [%]')
plt.xlabel('Calculated Variance of Estimation')
plt.title('Variance of Esitmation: Calculated Values Results Differing N')