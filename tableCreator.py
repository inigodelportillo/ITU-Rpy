# -*- coding: utf-8 -*-
"""
Recreation of min_squares.py
Actually creates a rectangular data table from 
arrays of probability and elevation angle as inputs. 

Created on Mon Oct 19 19:39:26 2020

@author: MAW32652
"""
import itur
###input variables 
lat = 39 #latitude
lon = 283 #longitude
hs = 1 #altitude [km]
f = 20 #frequency [GHz]
d = 3 #antenna diameter [m] 
tau = 45 #polarization tilt angle [degrees]



def tableCreator(pList, eleList, sqList = []):
       
    #if the input lists are of the same length
    if len(pList) == len(eleList):
        #append the add the results to the list of results.
        result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList, d, hs, return_contributions =False, include_gas = True)
        sqList.append(result)
        print("Done!")
        print("Your ist of square results is: " + str(sqList))
        #return the lsit of results. 
        return sqList
    
    #case where the m dimension is larger than the n dimension
    elif len(pList) > len(eleList):
        #append the results for m - n and all of n
        result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList[:len(eleList)], d, hs, return_contributions =False, include_gas = True)
        sqList.append(result)
        
        #do minSq on the rest of the rectangle
        #start at len(n) to the end and all of n. 
        return tableCreator(pList[len(eleList):], eleList, sqList)
    
    #case where the n dimension is larger than the m dimension.
    elif len(eleList) > len(pList):
        #append the results for all of m and n - m.
        result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList[:len(pList)], pList, d, hs, return_contributions =False, include_gas = True)        
        sqList.append(result)
        
        #do minSq on the rest of the rectangle
        #start at all of m and at len(m)
        return tableCreator(pList, eleList[len(pList):], sqList)

