# -*- coding: utf-8 -*-
"""
Gamma Validation Scirpt
Created on THU Jun 18 3:13:08 2020

This script pulls data from ITU-R P.676 validation data and calculates the gamma_o and gamma_w using the itur package.
These are the interemediate values that are used to calculated the gasoues attenuation.

The script then comapres these claculated values to their expected values and determines the error. 
@author: MAW32652
"""

import itur
from itur.utils import prepare_quantity, prepare_output_array,\
    prepare_input_array, load_data, dataset_dir, memory
from itur.models.itu676 import gamma0_exact, gammaw_exact
import os
import numpy as np
import xlrd

def gamma_intermediate_value_validation():
    
    ########## VARIABLE INITIALIZATION ###########
    
    ###excel variables
    
    workbook = xlrd.open_workbook("CG-3M3J-13-ValEx-Rev5_0.xlsx")
    sheet = workbook.sheet_by_name("P676-12 SpAtt")
    
    ### ITU_rpy variables
    
    # Validation data V=variables
    gOData = [] #gamma_o expected values
    gWData = [] #gamma_w expected values
    gTData = [] #Gamma_tot expected values
    
    # Input variables
    p = 1013.25 #Dry Air pressure 
    e = 9.97288879 #Water vapor partial pressure [hPA]
    T = 288.15000000 #Temperature [K]
    rho = 7.50000000 #density [g/m^3]
    fList = [] #Frequency, iterable, populated later
    
    # Output variables
    
    #Oxygen
    siOList = [] #Fucntion parameter
    fiOList = [] #Function parameter
    dOList = [] #Debye spectrum width parameter
    NppO = [] #Oxygen frequency dependent complex refractivity
    gOList = [] #Calculated gamma_o values
    gOErrorList = [] #Error for calculated gamma_o values
    gOPEList = [] #gamma_O percent Error
    
    #Water Vapor
    siWList = [] #Fucntion parameter
    fiWList = [] #Function parameter
    dwList = [] #Debye spectrum width parameter
    NppW = [] #Water vapor frequency dependent complex refractivity
    gWList = [] #Calculated gamme_w values
    gWErrorList = [] #Error for calculated gamma_o values
    gWPEList = [] #Gamma_W percent Error. 
    
    ########## COMPUTATION ###########
    
    #Compute the gamma_o and gamma_w for frequencies 1-350 GHz
    for i in range(350):
        #data entries start on row 34
        fList.append(sheet.cell_value(i + 33, 2)) #Frequency Data in column C
        gOData.append(sheet.cell_value(i + 33, 3)) #Gamma_o data in column D
        gWData.append(sheet.cell_value(i + 33, 4))#Gamma_w data in column E
        gTData.append(sheet.cell_value(i + 33, 5)) #Gamma_tot data in column F
        
        #Compoute the gamma values
        gamma_o = gamma0_exact(fList[-1], p , rho , T )
        gamma_w = gammaw_exact(fList[-1], p, rho, T)
        
        #Add the gamme values to the end of the output list
        
        #Calculate the error for the determined errors and add them to the list of errors. 
        gOError = gOData[-1] - gamma_o.value
        gOPE = (gamma_o.value - gOData[-1]) / gOData[-1]
        gWError = gWData[-1] - gamma_w.value
        gWPE = (gamma_w.value - gWData[-1]) / gWData[-1]
        
        gOErrorList.append(gOError)
        gOPEList.append(gOPE)
        gWErrorList.append(gWError)
        gWPEList.append(gWPE)
           
        
    ########### ERROR ANALYSIS ##########
    
    gOAvg = sum(gOErrorList) / len(gOErrorList)
    gWAvg = sum(gWErrorList) / len(gWErrorList)
    
    gOMax = max(list(map(abs, gOErrorList)))
    gWMax = max(list(map(abs, gWErrorList)))
    
    gOAvgPE = sum(gOPEList)/len(gOPEList)
    gWAvgPE = sum(gWPEList)/len(gWPEList)
    
    
    print("Error Analysis: Expected value - Calculated value")
    
    print()
    print ("Gamma_O:")
    print("Average: " + str(gOAvg) + " dB/km")
    print("Max: " + str(gOMax) + " dB/km")    
    print("Average Percent Error: " + str(gOAvgPE))
    
    print()
    print ("Gamma_W:")
    print("Average: " + str(gWAvg) + " dB/km")
    print("Max: " + str(gWMax) + " dB/km")
    print("Average Percent Error: " + str(gWAvgPE))        
    
    
    

    




