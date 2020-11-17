# -*- coding: utf-8 -*-
"""
Function that creates an 2D data table with probability and elevation angle as the
varrying inputs using either the square tiled rectangle approach, or the column by column
approach, depending on the input mode. 

Created on Tue Nov  3 12:45:37 2020

@author: MAW32652
"""

import itur
import numpy as np
import warnings
from itur.models.itu678 import p_LT_vs_Risk
from tiledTable import tiledTable

def tableCreator(pList, eleList, lat, lon, hs, f, d, tau, mode, risk = None):
    """creates an 2D data table with probability and elevation angle as the
    varrying inputs using either the square tiled rectangle approach, or the column by column
    approach, depending on the input mode. 
    
    Inputs:        
        pList: List of probabilities. 
            type == numpy.array
        eleList: List of elevation angles. 
            type == numpy.array
        lat: Geographic latitude of ground station. 
            type == float
        lon: Geographic longitude of ground station. 
            type == float
        hs: Altitude of ground station [km]. 
            type == float
        f: Signal frequency [GHz]. 
            type == float
        d: Antenna diameter [m]. 
            type == float
        tau: Polarization tilt angle [degree]. 
            type == float
        mode: specified morde to run in. 
            type == string            
            values == 'tiled' or 'cbc'
        risk : The probability that the yearly probability is exceeded. 
            type == float
            values ==  0 < risk < 1 
            
    Output:
        2D data table with probability and elevation angle as variable inputs. 
        """
    if risk == None:
        if mode == 'tiled' or mode == None:
            print()
            result = tiledTable(pList, eleList, lat = lat, lon = lon, hs = hs, f = f, d = d, tau = tau)
            return result
        
        elif mode == 'cbc':
            #Filling in the for loops "column-by-column" (Passing in elevation angles as an array of values for each p)
            result = []
            for i in pList:
                toAdd = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, i, d, hs, return_contributions =False)
                result.append(toAdd)
            return result
        
        else: 
            warnings.warn("Please specify either 'tiled' or 'cbc' mdoe'")
    
    #if risk is to be implemented
    else:
        print("risk is being implemented")
        #pList / 100 beacue AL calc uses %'s where risk calc use 0-1
        newPList = p_LT_vs_Risk(lat, lon, (pList / 100), risk, plot = False)
        #transform newPList to %'s again so AL can be determined
        newPList = newPList * 100
        
        #comput AL
        result = tiledTable(newPList, eleList, lat = lat, lon = lon, hs = hs, f = f, d = d, tau = tau)
        return result    
        
        
        