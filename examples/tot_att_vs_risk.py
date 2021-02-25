# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 11:14:53 2020

@author: MAW32652
"""


import itur
from itur.models.itu678 import p_lt_from_risk
import numpy as np

risk_arr = np.linspace(10, 99, 90)
p_lt_arr = p_lt_from_risk(33.915335, -118.38257, 1, risk_arr)
print(p_lt_arr)