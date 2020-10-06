# -*- coding: utf-8 -*-
"""
Python Script that computes all the values in various data tables and produces
the total runtim of the script at the end of teh program

Created on Tue Jul 28 13:26:23 2020

@author: MAW32652
"""
import datetime
import itur
import numpy as np

###input variables lists
lat = 39 #latitude
lon = 283 #longitude
hs = 1 #altitude [km]
f = 20 #frequency [GHz]
d = 3 #antenna diameter [m] 
tau = 45 #polarization tilt angle [degrees]


#p values for full large 2D Table
#create a numpy list of the correct probability inputs
pListL = np.arange(1, 51, 1)
pListL = np.flip(pListL)
toAdd = np.arange(.1, 1, .1)
toAdd = np.flip(toAdd)
pListL = np.append(pListL, toAdd)
pListL = np.append(pListL, 0.01)
pListL = np.append(pListL, 0.001)
pListL = np.append(pListL, .0005)

#p values for small 2D data table 
pS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
pListS = np.array(pS)

#elevation angle values for full large 2D table
#create a numpy list of the correct elevation angle inputs
eleListL =  np.arange(start = 0.5, stop = 90.5, step = 0.5) 
eleListL = np.flip(eleListL)
eleListL = np.append(eleListL, 0.1)
eleListL = np.flip(eleListL)

#elevation angle values for an extra large 2D table
eleListXl = np.arange(start = 0.5, stop = 90.25, step = 0.25)

#Elevation angle  values for creating a very large 2D Table
bigList = np.arange(start = 0.1, stop = 90.05, step = 0.05)

#Elevation angle values for creating a 2D data table with 1M entries
vBigList = np.arange(start = 0.1, stop = 90.0005, step = 0.005)

#elevation angle values for small 2D data table
eleS = [1, 2, 3, 4, 5, 6, 7, 8]
eleListS = np.array(eleS)

#probability lists
pList8 = np.linspace(1, 5, 8)
pList62 = np.linspace(1, 5, 62)
pList100 = np.linspace(1, 5, 100)
pList10000 = np.linspace(1, 5, 10000)

#elevation angle lists
eleList8 = np.linspace(30, 60, 8)
eleList100 =np.linspace(30, 60, 100)
eleList181 = np.linspace(30, 60, 181)
eleList359 = np.linspace(30, 60, 359)
eleList1779 = np.linspace(30, 60, 1779)
eleList10000 = np.linspace(30, 60, 10000)
eleList17981 = np.linspace(30, 60, 17981)



def runtime(pList, eleList, mode):
    #compute the outputs
    begin_time = datetime.datetime.now()
    outputs = [] 
    count = 0
    
    if mode == "cbc" or mode == None:       
        #Filling in the for loops "column-by-column" (Passing in elevation angles as an array of values for each p)
        for i in pList:
            toAdd = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, i, d, hs, return_contributions =False)
            outputs.append(toAdd)
            count += 1
        
    elif mode == "arrs":
        outputs = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList, d, hs, return_contributions =False, include_gas = True)
        count += 1 
    
    elif mode == "nested":
        #Using nested for loops used to fill in the 2D Data table       
        for i in pList:
            for j in eleList:
                toAdd = itur.atmospheric_attenuation_slant_path(lat, lon, f, j, i, d, hs, return_contributions =False)
                outputs.append(toAdd)
                count += 1
    
    #print(outputs)    
    print("Total runtime: "  + str(datetime.datetime.now() - begin_time))

    print("Number of Columns: " + str(len(outputs)))
    if mode == "cbc" or mode == "arrs":   
        print("Number of Rows: " + str(len(outputs[0])))        
        #print("Number of Rows: " + str(count))
    
    



