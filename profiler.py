# -*- coding: utf-8 -*-
"""
Python Script that computes all the values in various data tables and produces
the total runtim of the script at the end of teh program

Created on Tue Jul 28 13:26:23 2020

@author: MAW32652
"""
import datetime
import cProfile
import itur
import numpy as np
from tiledTable import tiledTable


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

#elevation angle values for full large 2D table
#create a numpy list of the correct elevation angle inputs
eleListL =  np.arange(start = 0.5, stop = 90.5, step = 0.5) 
eleListL = np.flip(eleListL)
eleListL = np.append(eleListL, 0.1)
eleListL = np.flip(eleListL)

#Extra probability lists
pList8 = np.linspace(1, 5, 8)
pList62 = np.linspace(1, 5, 62)
pList100 = np.linspace(1, 5, 100)
pList1000 =np.linspace(1, 5, 1000)
pList10000 = np.linspace(1, 5, 10000)

#Extra elevation angle lists
eleList8 = np.linspace(30, 60, 8)
eleList100 = np.linspace(30, 60, 100)
eleList181 = np.linspace(30, 60, 181)
eleList359 = np.linspace(30, 60, 359)
eleList1000 = np.linspace(30, 60, 1000)
eleList1779 = np.linspace(30, 60, 1779)
eleList10000 = np.linspace(30, 60, 10000)
eleList17981 = np.linspace(30, 60, 17981)



def runtime(pList, eleList, mode, lat = lat, lon = lon, hs = hs, f = f, d = d, tau = tau):
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
        
    elif mode == "minsq":
        #sqSizes = [1,2,5,1,3,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8] #M Case
        #sqSizes = [2,2,6,26,26,4,10,10,31,31,62,62] #L Case
        #sqSizes = [31,31,62,62,62,62,10,20,16,16,10,32,30,32,30] #XL Case
        #sqSizes = [62, 62, 62, 62, 62, 26, 18, 18, 12, 24, 26, 12, 36, 26] #3XL Case
        sqSizes = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]
        
        for thisSq in sqSizes:
            pList = np.linspace(1,5,thisSq)
            eleList = np.linspace(30,60,thisSq)
            outputs = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList, d, hs, 
                                                  return_contributions =False, include_gas = True)
        count +=1
    
    elif mode == "tiled":
        result = tiledTable(pList, eleList, lat = lat, lon = lon, hs = hs, f = f, d = d, tau = tau)
        count += 1
        print("Total runtime: "  + str(datetime.datetime.now() - begin_time))

        print("Number of Columns: " + str(len(result)))
       
        print("Number of Rows: " + str(len(result[0])))        
            
        return result
    
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
        #print("Number of Rows: " + str(len(outputs[0])))        
        print("Number of Rows: " + str(count))
    