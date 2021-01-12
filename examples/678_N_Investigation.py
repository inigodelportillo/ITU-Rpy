# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:59:57 2020

@author: MAW32652
"""
from itur.models.itu678 import variance_of_estimation_calc, p_lt_from_risk
import numpy as np
import matplotlib.pyplot as plt
import time

actual_sigma = variance_of_estimation_calc(0.01)

start = time.time()
actual_p = p_lt_from_risk(33.915335, -118.38257,1, 10, N = 525960)
end = time.time()
actual_time = end - start
nList = np.linspace(1000, 530000, 101)
nList = nList.astype(int)


def func(L):
    sigmas, ps, times = [], [], []
    for i in range(len(L)):
        sigma = variance_of_estimation_calc(0.01, L[i])
        
        start = time.time()
        p = p_lt_from_risk(33.915335, -118.38257,1, 10, L[i])
        end = time.time()
        runtime = end - start
        print(runtime)
        
        sigmas.append(sigma)
        ps.append(p)
        times.append(runtime)
    return sigmas, ps, times

sigmas, ps, times = func(nList)
sigmaPE = (abs(sigmas - actual_sigma)/actual_sigma) * 100
pPE = (abs(ps - actual_p) / actual_p) * 100


plt.plot(nList, sigmaPE)
plt.xscale('log')
plt.xlabel('N')
plt.ylabel('Percent Error')
plt.title('N vs Percent Error in Sigma')

plt.plot(nList, pPE)
plt.xscale('log')
plt.xlabel('N')
plt.ylabel('Percent Error')
plt.title('N vs Percent Error in p_LT')