# -*- coding: utf-8 -*-
"""
Recreation of min_squares.py
Actually creates a rectangular data table from 
arrays of probability and elevation angle as inputs. 

Created on Mon Oct 19 19:39:26 2020

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

def appendLeft():
    print("I appended to what was to the left.")
    
def appendAbove():
    print("I appended to what was above.")


def tableCreator(pList, eleList, output = [], toAdd = [], prev = ''):
    print("output is:")
    print(output)
    print(np.shape(output))
    print()
    
    print("toAdd is:")
    print(toAdd)
    print(np.shape(toAdd))
    print()
    #if the input lists are of the same length
    if len(pList) == len(eleList):
        #Calculate and append the results to the list of results.
        result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList, d, hs, return_contributions =False, include_gas = True)
        
        #if this is the first iteration
        if prev == '':
            output = result
            
        #if the previous dimensions had ele > p
        elif prev == 'e>p':
            
            #check if toAdd is emmpty
            if str(type(toAdd)) == "<class 'list'>":
                toAdd = result
                
            #if its not, append result to it
            else:
                np.append(toAdd, result, axis = 0)
                
            #check if toAdd and output are the correct shape to be appended and append if they are
            if np.shape(toAdd)[0] == np.shape(output)[0]:               
                output = np.append(output, toAdd, axis = 1)
                toAdd = [] #clearing toAdd so that it can be used again

        
        #if the previous dimensions had p > ele        
        elif prev == 'p>e':
            
             #check if toAdd is emmpty
            if str(type(toAdd)) == "<class 'list'>":
                toAdd = result
            #if its not, append result to it
            else:
                np.append(toAdd, result, axis = 1)
                
            #check if toAdd and output are the correct shape to be appended
            if np.shape(toAdd)[1] == np.shape(output)[1]:               
                output = np.append(output, toAdd, axis = 0)
                toAdd = [] #clearing toAdd so that it can be used again

        #return the lsit of results. 
        return output
    
    #case where the pList dimension is larger than the eleList dimension
    elif len(pList) > len(eleList):
        #Calculate the results for pList - eleList and all of eleList
        result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList, pList[:len(eleList)], d, hs, return_contributions =False, include_gas = True)
        
        ### Append then result dependent on values from previous calls. 
        
        #if this is the first iteration
        if prev == '':           
            output = result
            print("result is:")
            print(result)
            print(np.shape(result))
            print()
            return tableCreator(pList[len(eleList):], eleList, output, toAdd, prev = 'p>e')
             
        #if the previous iterations dimensions had ele > p
        elif prev == 'e>p':
            
            #check if toAdd is emmpty
            if str(type(toAdd)) == "<class 'list'>":
                toAdd = result
            #if its not, append result to it
            else:
                np.append(toAdd, result, axis = 0)
                
            #check if toAdd and output are the correct shape to be appended and append if they are
            if np.shape(toAdd)[0] == np.shape(output)[0]:               
                output = np.append(output, toAdd, axis = 1)
                toAdd = [] #clearing toAdd so that it can be used again
            
        
        #if the previous iterations dimensions had p > ele        
        elif prev == 'p>e':
            
            #check if toAdd is emmpty
            if str(type(toAdd)) == "<class 'list'>":         
                toAdd = result
            #if its not, append result to it
            else:
                np.append(toAdd, result, axis = 1)
                
            #check if toAdd and output are the correct shape to be appended
            if np.shape(toAdd)[1] == np.shape(output)[1]:               
                output = np.append(output, toAdd, axis = 0)
                toAdd = [] #clearing toAdd so that it can be used again

        #do tableCreator on the rest of the rectangle
        #start at len(eleList) to the end and all of eleList. 
        print("result is:")
        print(result)
        print(np.shape(result))
        print()
        return tableCreator(pList[len(eleList):], eleList, output, toAdd, prev = 'p>e')
    
    #case where the eleList dimension is larger than the pList dimension.
    elif len(eleList) > len(pList):
        #Calculate and append the results for all of pList and eleList - pList.
        result = itur.atmospheric_attenuation_slant_path(lat, lon, f, eleList[:len(pList)], pList, d, hs, return_contributions =False, include_gas = True)        
       
        ### Append then result dependent on values from previous calls. 
        
        #if this is the first iteration
        if prev == '':           
            output = result
            print("result is:")
            print(result)
            print(np.shape(result))
            print()
            return tableCreator(pList, eleList[len(pList):], output, toAdd, prev = 'e>p')
             
        #if the previous iterations dimensions had ele > p
        elif prev == 'e>p':
            
             #check if toAdd is emmpty
            if str(type(toAdd)) == "<class 'list'>":
                toAdd = result
                
            #if its not, append result to it
            else:
                np.append(toAdd, result, axis = 0)
                
            #check if toAdd and output are the correct shape to be appended
            if np.shape(toAdd)[0] == np.shape(output)[0]:               
                output = np.append(output, toAdd, axis = 1)
                toAdd = [] #clearing toAdd so that it can be used again

                
        #if the previous iterations dimensions had p > ele        
        elif prev == 'p>e':
            
             #check if toAdd is emmpty
            if str(type(toAdd)) == "<class 'list'>":
                toAdd = result
                
            #if its not, append result to it
            else:
                np.append(toAdd, result, axis = 1)
                
            #check if toAdd and output are the correct shape to be appended
            if np.shape(toAdd)[1] == np.shape(output)[1]:               
                output = np.append(output, toAdd, axis = 0)
                toAdd = [] #clearing toAdd so that it can be used again

        
        #do tableCreator on the rest of the rectangle
        #start at all of pList and at len(pList)
        print("result is:")
        print(result)
        print(np.shape(result))
        print()
        return tableCreator(pList, eleList[len(pList):], output, toAdd, prev = 'e>p')

