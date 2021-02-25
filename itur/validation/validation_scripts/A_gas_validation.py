# -*- coding: utf-8 -*-
"""
A_gas_valiudation

Created on Tue Jul 14 11:30:38 2020

Validated the outputs computed with itur to the values expected in the A_gas
sheet in the validation data. 
@author: MAW32652
"""
import itur
from itur.models.itu676 import gamma0_approx
from itur.models.itu676 import _ITU676_12_
from itur.models.itu676 import gamma0_exact, gammaw_exact
import os
import numpy as np
import xlrd


def gaseous_attenuation_validation():
    ########## VARIABLE INITIALIZATION ###########
    
    ###excel variables
    
    workbook = xlrd.open_workbook("CG-3M3J-13-ValEx-Rev5_0.xlsx")
    sheet = workbook.sheet_by_name("P676-12 A_Gas ")
    
    ### ITU-Rpy variables
    
    #input paramters pulled from the validation data
    latList = [] #Latitude
    lonList = [] #Longitude
    hsList = [] #ground station altitude
    satLonList = [] #satelite longitude
    pList = [] #probability
    rhoList = [] #rho as a function of p (P.836)
    eleList = [] #elevation angle
    freqList = [] #frequency
    tList = [] #Temperature (Kelvin)
    eList = [] # e as a function of p 
    dryList = [] #dry air pressure. 
    vPList = [] 
    
    #Validation Data Expected Outputs (Includes intermediate values)
    
    #Oxygen
    gOData = [] #gamma_O intermediate value
    hOData = [] #Oxygen equivalent height intermediate value
    aOData = [] #Attenuation due to Oxygen
    aOsinData = [] #Attenuation due to Oxygen multipled by 1/sin(ele)
    
    #Water Vapor
    rhoData = [] #rho ref intermediate value
    tRefData = [] #Temperateure red intermediate value
    eData = [] #e ref intermediate value
    dryData = [] #Dry air pressure intermediate value
    gWData = [] #gamma_w intermediate value
    gWRefData = [] #gamma_w ref intermediate value
    aData = [] #a intermediate value
    bData = [] #b intermediate vlaue
    abData = [] #(a*h**b + 1) intermediate value
    funcData = [] #0.0176*(yw/yw_ref)*(a*h**b + 1) intermediate value
    aWData = [] #Attenuation due to water vapor
    aWsinData = [] #Attenuation due to water vapor multipled by 1/sin(ele)
    
    aGData = [] #gaseous attenuation
    
    # itur output lists
    
    #Oxygen
    gOList = [] #gamma_O intermediate value
    hOList = [] #Oxygen equivalent height intermediate value
    aOList = [] #Attenuation due to Oxygen
    aOsinList = [] #Attenuation due to Oxygen multipled by 1/sin(ele)
    
    #Water Vapor
    rhoList = [] #rho ref intermediate value
    tRefList = [] #Temperateure ref intermediate value
    eList = [] #e ref intermediate value
    dryList = [] #Dry air pressure intermediate value
    gWList = [] #gamma_w intermediate value
    gWRefList = [] #gamma_w ref intermediate value
    aList = [] #a intermediate value
    bList = [] #b intermediate vlaue
    abList = [] #(a*h**b + 1) intermediate value
    funcList = [] #0.0176*(yw/yw_ref)*(a*h**b + 1) intermediate value
    aWList = [] #Attenuation due to water vapor
    aWsinList = [] #Attenuation due to water vapor multipled by 1/sin(ele)
    
    aGList = [] #Total gaseous attenuation
    
    #Error Lists
    gOErrorList = [] 
    gOPEList = [] 
    hOErrorList = [] 
    hOPEList =[] 
    aOErrorList =[]
    aOPEList =[] 
    aOsinErrorList = [] 
    aOsinPEList = [] 
    rhoErrorList = []
    rhoPEList = [] 
    tRefErrorList = [] 
    tRefPEList = [] 
    aErrorList = [] 
    aPEList =[] 
    bErrorList = [] 
    bPEList = [] 
    abErrorList = [] 
    abPEList = []
    aWErrorList = [] 
    aWPEList =[] 
    aWsinErrorList = [] 
    aWsinPEList = []
    aGErrorList = [] 
    aGPEList = [] 
    
    ########## CALCULATION ##########
    
    #For lopp iterates through all 64 of the validation test case
    for i in range(64):
        
        #populate the variable lists
        #note that data entry starts on row 22
        
        #inputs
        latList.append(sheet.cell_value(i + 21, 2)) #latitude data  in column C
        lonList.append(sheet.cell_value(i + 21, 3)) #longitude data  in column D
        hsList.append(sheet.cell_value(i + 21, 4)) #ground station altitude data  in column E
        satLonList.append(sheet.cell_value(i + 21, 5)) #satelite longitude data  in column F
        pList.append(sheet.cell_value(i + 21, 6)) #probability data  in column G
        rhoList.append(sheet.cell_value(i + 21, 7)) #rho as a funciton of p data  in column H
        eleList.append(sheet.cell_value(i + 21, 8)) #elevation andle data  in column I
        freqList.append(sheet.cell_value(i + 21, 9)) #requency data  in column J
        tList.append(sheet.cell_value(i + 21, 10)) #Temperature data  in column K
        eList.append(sheet.cell_value(i + 21, 11)) #e as a function of p data  in column L
        dryList.append(sheet.cell_value(i + 21, 12)) #dry air pressure data  in column M
        vPList.append(sheet.cell_value(i + 21, 20)) 
        
        #Oxygen outputs including intermediate values
        gOData.append(sheet.cell_value(i + 21, 13)) #gamma_o data  in column N
        hOData.append(sheet.cell_value(i + 21, 14)) #hO data  in column O
        aOData.append(sheet.cell_value(i + 21, 15)) #Attenuation due to oxygen data  in column P
        aOsinData.append(sheet.cell_value(i + 21, 16)) #Attenuation due to oxygen (sin factor) data  in column Q
        
        #Water vapor outputs including intermediate values
        rhoData.append(sheet.cell_value(i + 21, 21)) #rho ref data  in column V
        tRefData.append(sheet.cell_value(i + 21, 22)) #T ref data  in column W
        eData.append(sheet.cell_value(i + 21, 23)) #e ref data  in column X
        dryData.append(sheet.cell_value(i + 21, 24)) #dry aior pressure data  in column Y
        gWData.append(sheet.cell_value(i + 21, 25)) #gamma_w data  in column Z
        gWRefData.append(sheet.cell_value(i + 21, 26)) #gamma_w ref data  in column AA
        aData.append(sheet.cell_value(i + 21, 27)) #a data  in column AB
        bData.append(sheet.cell_value(i + 21, 28)) #b data  in column AC
        abData.append(sheet.cell_value(i + 21, 29)) #ab func  data  in column AD
        funcData.append(sheet.cell_value(i + 21, 30)) #0.0176*(yw/yw_ref)*(a*h**b + 1) data  in column AE
        aWData.append(sheet.cell_value(i + 21, 31)) #Attemuation due to water vapor data  in column AF
        aWsinData.append(sheet.cell_value(i + 21, 32)) #attenuation due to water vapor (sin factor) data  in column AG
        
        #expected total gaseoues attenuation
        aGData.append(sheet.cell_value(i + 21, 34)) #Gaseous attenuation data  in column AI
        
        #Compute the Oxygen intermediate values and Attenuation due to Oxygen
        #Error Computations for Oxygen Values are also done in this section
        gamma_o = gamma0_exact(freqList[-1], dryList[-1], rhoList[-1], tList[-1])
        gOList.append(gamma_o.value)
        gOError = gOData[-1] - gOList[-1]
        gOPE = (gOList[-1] - gOData[-1]) / gOData[-1]
        gOErrorList.append(gOError)
        gOPEList.append(gOPE)
        
        hO = _ITU676_12_.slant_inclined_path_equivalent_height(freqList[-1], dryList[-1],  rhoList[-1], tList[-1])
        hOList.append(hO[0])
        hOError = hOData[-1] - hOList[-1]
        hOPE = (hOList[-1] - hOData[-1]) / hOData[-1]
        hOErrorList.append(hOError)
        hOPEList.append(hOPE)
        
        aO = gamma_o * hO[0]
        aOList.append(aO.value)
        aOError = aOData[-1] - aOList[-1]
        aOPE = (aOList[-1] - aOData[-1]) / aOData[-1]
        aOErrorList.append(aOError)
        aOPEList.append(aOPE)
        
        aOsin = aO/np.sin(np.deg2rad(eleList[-1]))
        aOsinList.append(aOsin.value)
        aOsinError = aOsinData[-1] - aOsinList[-1]
        aOsinPE = (aOsinList[-1] - aOsinData[-1]) / aOsinData[-1]
        aOsinErrorList.append(aOsinError)
        aOsinPEList.append(aOsinPE)
        
        #Compute the watervapor intermediate values and Attenuation due to water vapor
        #Error Computations for water vapor Values are also done in this section
        rho_ref, t_ref, a, b, ab, aW = _ITU676_12_.zenit_water_vapour_attenuation(latList[-1], lonList[-1], pList[-1], freqList[-1], vPList[-1], hsList[-1])
        
        rhoList.append(rho_ref)
        rhoError = rhoData[-1] - rhoList[-1]
        rhoPE = (rhoList[-1] - rhoData[-1]) / rhoData[-1]
        rhoErrorList.append(rhoError)
        rhoPEList.append(rhoPE)
        
        tRefList.append(t_ref)
        tRefError = tRefData[-1] -tRefList[-1]
        tRefPE = (tRefList[-1] - tRefData[-1]) / tRefData[-1]
        tRefErrorList.append(tRefError)
        tRefPEList.append(tRefPE)
        
        aWList.append(aW)
        aWError = aWData[-1] - aWList[-1]
        aWPE = (aWList[-1] - aWData[-1]) / aWData[-1]
        aWErrorList.append(aWError)
        aWPEList.append(aWPE)
        
        aWsin = aW / np.sin(np.deg2rad(eleList[-1]))
        aWsinList.append(aWsin)
        aWsinError = aWsinData[-1] - aWsinList[-1]
        aWsinPE = (aWsinList[-1] - aWsinData[-1]) / aWsinData[-1]
        aWsinErrorList.append(aWsinError)
        aWsinPEList.append(aWsinPE)
        
        #Compute and calculate the error for the total gaseous attenuation
        aG = aOsinList[-1] + aWsinList[-1]
        aGList.append(aG)
        aGError = aGData[-1] - aGList[-1]
        aGPE = (aGList[-1] - aGData[-1]) / aGData[-1]
        aGErrorList.append(aGError)
        aGPEList.append(aGPE)
    
    
    gOAvg = sum(gOErrorList) / len(gOErrorList)
    hOAvg = sum(hOErrorList) / len(hOErrorList)
    aOAvg = sum(aOErrorList) / len(aOErrorList)
    aOsinAvg = sum(aOsinErrorList) / len(aOsinErrorList)
    aWAvg = sum(aWErrorList) / len(aWErrorList)
    aWsinAvg = sum(aWsinErrorList) / len(aWsinErrorList)
    rhoRefAvg = sum(rhoErrorList) / len(rhoErrorList)
    tRefAvg = sum(tRefErrorList) / len(tRefErrorList)
    aGAvg = sum(aGErrorList) / len(aGErrorList)
    
    
    gOMax = max(list(map(abs, gOErrorList)))
    hOMax = max(list(map(abs, hOErrorList)))
    aOMax = max(list(map(abs, aOErrorList)))
    aOsinMax = max(list(map(abs, aOsinErrorList)))
    aWMax = max(list(map(abs, aWErrorList)))
    aWsinMax = max(list(map(abs, aWsinErrorList)))
    rhoRefMax = max(list(map(abs, rhoErrorList)))
    tRefMax = max(list(map(abs, tRefErrorList)))
    aGMax = max(list(map(abs, aGErrorList)))
    
    
    gOAvgPE = sum(gOPEList) / len(gOPEList)
    hOAvgPE = sum(hOPEList) / len(hOPEList)
    aOAvgPE = sum(aOPEList) / len(aOPEList)
    aOsinAvgPE = sum(aOsinPEList) / len(aOsinPEList)
    aWAvgPE = sum(aWPEList) / len(aWPEList)
    aWsinAvgPE = sum(aWsinPEList) / len(aWsinPEList)
    rhoRefAvgPE = sum(rhoPEList) / len(rhoPEList)
    tRefAvgPE = sum(tRefPEList) / len(tRefPEList)
    aGAvgPE = sum(aGPEList) / len(aGPEList)
    
    
    ########## ERROR ANALYSIS ##########
    print("Oxygen")
    
    print()
    print ("Gamma O: ")
    print("Average Error: " + '{:0.2e}'.format(gOAvg) + " dB/km")
    print("Max Error: " + '{:0.2e}'.format(gOMax) + " dB/km")
    print("Average Percent Error: " + '{:0.2e}'.format(gOAvgPE))    
    
    print()
    print ("Oxygen Equivalent Height (hO): ")
    print("Average Error: " + '{:0.2e}'.format(hOAvg) + " km")
    print("Max Error: " + '{:0.2e}'.format(hOMax) + " km")
    print("Average Percent Error: " + '{:0.2e}'.format(hOAvgPE))    
    
    print()
    print ("Attenuation due to Oxygen (A_O): ")
    print("Average Error: " + '{:0.2e}'.format(aOAvg) + " dB")
    print("Max Error: " + '{:0.2e}'.format(aOMax) + " dB")
    print("Average Percent Error: " + '{:0.2e}'.format(aOAvgPE)) 
    
    print()
    print ("A_O / sin(Elevation Angle): ")
    print("Average Error: " + '{:0.2e}'.format(aOsinAvg) + " dB")
    print("Max Error: " + '{:0.2e}'.format(aOsinMax) + " dB")
    print("Average Percent Error: " + '{:0.2e}'.format(aOsinAvgPE))     
    
    print()
    print("-----------------------------------------------------------------------")
    print()
    print("Water Vapour")
    
    print()
    print ("Attenuation due to Water Vapour (A_W):")
    print("Average Error: " + '{:0.2e}'.format(aWAvg) + " dB")
    print("Max Error: " + '{:0.2e}'.format(aWMax) + " dB")
    print("Average Percent Error: " + '{:0.2e}'.format(aWAvgPE))   
    
    print()
    print ("A_W / sin(Elevation Angle): ")
    print("Average Error: " + '{:0.2e}'.format(aWsinAvg) + " dB")
    print("Max Error: " + '{:0.2e}'.format(aWsinMax) + " dB")
    print("Average Percent Error: " + '{:0.2e}'.format(aWsinAvgPE))     
    
        
    print()
    print ("Rho_ref: ")
    print("Average Error: " + '{:0.2e}'.format(rhoRefAvg) + " g/m^3")
    print("Max Error: " + '{:0.2e}'.format(rhoRefMax) + " g/m^3")
    print("Average Percent Error: " + '{:0.2e}'.format(rhoRefAvgPE))  
      
    print()
    print ("t_ref: ")
    print("Note: The error in t_ref comes from the +3 constant in the algorithm.")
    print("It does not seem like the validation data includes this constant.")
    print("Average Error: " + '{:0.2e}'.format(tRefAvg) + " degrees C")
    print("Max Error: " + '{:0.2e}'.format(tRefMax) + " degrees C")
    print("Average Percent Error: " + '{:0.2e}'.format(tRefAvgPE))  
    
    print()
    print("-----------------------------------------------------------------------")
    print()
    print("Gaseous Attenuation: ")
    
    print()
    print ("Total Gaseous Attenuation ")
    print("Average Error: " + '{:0.2e}'.format(aGAvg) + " dB")
    print("Max Error: " + '{:0.2e}'.format(aGMax) + " dB")
    print("Average Percent Error: " + '{:0.2e}'.format(aGAvgPE)) 
        
        
        