# -*- coding: utf-8 -*-
"""
Recreation of min_squares.py
Actually creates a rectangular data table from 
arrays of probability and elevation angle as inputs.

Created on Tue Oct 27 13:40:36 2020

@author: MAW32652
"""

import itur
import numpy as np
###input variables 
lat = 39 #latitude
lon = 283 #longitude
hs = 1 #altitude [km]
f = 20 #frequency [GHz]
d = 3 #antenna diameter [m] 
tau = 45 #polarization tilt angle [degrees]


def tableCreator(pList, eleList, output = [], pOff = 0, eleOff = 0):
    #check to see if this is the first iteration. 
    #if it is, create an array of correct dimmensions for output
    if str(type(output)) == "<class 'list'>":
        output = np.ones((len(pList), len(eleList)))

    #Base case where the input lists are of the same length
    if len(pList) == len(eleList):
        #Calculate and append the results 
        result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList, d, hs, return_contributions =False, include_gas = True)
        
        #replace the correct data cells in the output array.
        output[pOff : pOff + len(eleList), eleOff : eleOff + len(pList)] = result
        
        #return the rectangle of resutls
        return output
    
    #case where the pList dimension is larger than the eleList dimension.
    elif len(pList) > len(eleList):
        #Calculate the results for pList - eleList and all of eleList
        result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList[:len(eleList)], d, hs, return_contributions =False, include_gas = True)        

        #replace the correct data cells in the output array. 
        #the p rows that havent beel filled plus the number of rows equal to len(eleList)
        #all the ele columns that havent already been filled. 
        output[pOff : pOff + len(eleList), eleOff:] = result
        
        
        #update the offset variables
        pOff += len(eleList)
        
        #recursive call
        return tableCreator(pList[pOff:], eleList, output, pOff, eleOff)
        
    #case where the eleList dimension is larger than the pList dimension.
    elif len(eleList) > len(pList):
        #Calculate and append the results for all of pList and eleList - pList.
        result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList[:len(pList)], pList, d, hs, return_contributions =False, include_gas = True)

        #replace the correct data cells in the output array. 
        #all of the p rows that havent been filled
        #the ele columns that havent been filled plus the number of columns equal to len(pList)
        output[pOff:, eleOff: eleOff + len(pList)] = result
        
    
        #update the offset variables
        eleOff += len(pList)

        
        #recursive call
        return tableCreator(pList, eleList[eleOff:], output, pOff, eleOff)
    
    
testOut = tableCreator(np.linspace(1, 5, 62), np.linspace(30, 60, 181))