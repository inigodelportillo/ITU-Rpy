# -*- coding: utf-8 -*-
"""
Function that creates an 2D data table with probability and elevation angle as the
varrying inputs using the square tiled rectangle approach. 

Created on Wed Oct 28 13:05:02 2020

Copywrite The Aerospace Corporation 2020

@author: MAW32652
"""

import itur
import numpy as np
from tiler import tiler

###Default input variables 
lat = 39 #latitude
lon = 283 #longitude
hs = 1 #altitude [km]
f = 20 #frequency [GHz]
d = 3 #antenna diameter [m] 
tau = 45 #polarization tilt angle [degrees]

def tiledTable(pList, eleList, lat = lat, lon = lon, hs = hs, f = f, d = d, tau = tau):
    #check to see if this is the first iteration. 
    #if it is, create an array of correct dimmensions for output
    output = np.zeros((len(pList), len(eleList)))
    
    #Variable initialization
    eleCount = 0
    pCount = 0
    
    dimList = tiler(eleList, pList, dimList = [])
    
    for dim in dimList:
        
        #Case where the input lists are of the same length
        if len(pList) == len(eleList):
            #Calculate the results 
            result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList, d, hs, return_contributions =False, include_gas = True)
            
            #replace the correct data cells in the output array.
            output[pCount : pCount + len(eleList), eleCount : eleCount + len(pList)] = result
            
        #case where the pList dimension is larger than the eleList dimension.
        elif len(pList) > len(eleList):
            #Calculate the results for pList - eleList and all of eleList
            result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList[:dim], d, hs, return_contributions =False, include_gas = True)        
    
            #replace the correct data cells in the output array. 
            #the p rows that havent beel filled plus the number of rows equal to dim
            #all the ele columns that havent already been filled. 
            output[pCount : pCount + dim, eleCount:] = result
            
            
            #update the offset variables
            pCount += dim
            pList = pList[dim:]
            
            
            #case where the eleList dimension is larger than the pList dimension.
        elif len(eleList) > len(pList):
            #Calculate and append the results for all of pList and eleList - pList.
            result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList[:dim], pList, d, hs, return_contributions =False, include_gas = True)
    
            #replace the correct data cells in the output array. 
            #all of the p rows that havent been filled
            #the ele columns that havent been filled plus the number of columns equal to dim
            output[pCount:, eleCount: eleCount + dim] = result
            
        
            #update the offset variables
            eleCount += dim
            eleList = eleList[dim:]      
    return output
