# -*- coding: utf-8 -*-
"""
Validation Script
Created on Wed Jun 10 15:46:09 2020

This script pulls input data for the 64 cases in validation_results_comparison.xlsx

The script will then determine the total attenuation for each case.

Finally the script will determine the error between the publised value and the calculated value using ITU-Rpy

@author: MAW32652
"""
import datetime
import itur
import xlrd
import xlwt #not yet used but may be helpful for data presentation later

def total_attenuation_validation():
    begin_time = datetime.datetime.now()
    ########## VARIABLE INITIALIZATION ###########
    
    ###excel variables
    
    workbook = xlrd.open_workbook("CG-3M3J-13-ValEx-Rev5_0.xlsx")
    sheet = workbook.sheet_by_name("P618-13 Att_Tot")
    
    ### ITU-Rpy variables
    
    # Validation data variables
    pubData = [] #The published total attenuation
    agData = [] #The published gaseous attenuation
    acData = [] #The published cloud attenuation
    arData= [] #The published rain attenuation
    asData = [] #The published scintillation attenuation
    
    #input lists
    latList = [] #Latitude
    lonList = [] #Longitude
    altList = [] #Ground Station Altitude
    satLonList = [] #Satalite longitude
    freqList = [] #Signal Frequency
    eleList = [] #Elevation Angle
    diaList = [] #Antenna Diameter
    tauList = [] #Polarization tilt angle
    probList = [] #probablity
    
    #output lists
    attenList = [] #Total attenuation and componented attenuations. 2D List
    agList = [] #Gaseous attenuation
    acList = [] #Cload attenuation
    arList = [] #Rain attenuation
    asList = [] #Scintilation attenuation
    atList= [] #Total Attenuation
    errorList = [] #Error (Published - Calculated)
    atPEList = [] #total attenuation percent error
    agErrorList = [] #Erros of gaseous attenuation
    agPEList = [] #gaseous attenuation percent error
    acErrorList = [] #Errors of cloud attenuation
    acPEList = [] #cloud attenuation percent error
    arErrorList= [] #Errors of rain attenuation
    arPEList = [] #rain attenuation percent error
    asErrorList = [] #Errors  of scintillation attenuation
    asPEList = [] #scintilattion attenuation percent error
    largeErrors = [] #Errors greater than 0.1
    caseList = [] #Case's with large errors
    
    ########## CALCULATION ##########
    
    #For lopp iterates through all 64 of the validation test cases
    for i in range(64):
        
        #populate the variable lists
        #note that data entry starts on row 22
        pubData.append(sheet.cell_value(i + 21, 18)) #Published data in column S
        agData.append(sheet.cell_value(i + 21, 13)) #Gas attebnuation data in column M
        acData.append(sheet.cell_value(i + 21, 15)) #Cloud attenuation data in column O
        arData.append(sheet.cell_value(i + 21, 16)) #Rain attenuation data in column Q
        asData.append(sheet.cell_value(i + 21, 17)) #Scintalation attenuation datain column R
        latList.append(sheet.cell_value(i + 21, 2)) #Latitude data in column C
        lonList.append(sheet.cell_value(i + 21, 3)) #Longitude data in column D
        altList.append(sheet.cell_value(i + 21, 4)) #Ground station altitude data in column E
        satLonList.append(sheet.cell_value(i + 21, 5)) #Satalite longitude data in column F
        freqList.append(sheet.cell_value(i + 21, 6)) #Frequency data in column G
        eleList.append(sheet.cell_value(i + 21, 7)) #Elevation angle data in column H
        diaList.append(sheet.cell_value(i + 21, 8)) #Antenna diameter data in column I
        tauList.append(sheet.cell_value(i + 21, 10)) #Polarization tilt angle data in column K
        probList.append(sheet.cell_value(i + 21, 11)) #Probability data in column L
        
        #Calculate the attenuation and put the values into the output list
        a_g, a_c, a_r, a_s, a_t = itur.atmospheric_attenuation_slant_path(latList[-1], lonList[-1],
                                                                                  freqList[-1], 
                                                                                  eleList[-1], 
                                                                                  probList[-1], 
                                                                                  diaList[-1],
                                                                                  tau = tauList[-1],
                                                                                  hs = altList[-1],
                                                                                  return_contributions=True)
        attenuations = [a_g.value, a_c.value, a_r.value, a_s.value, a_t.value]
        attenList.append(attenuations)
        agList.append(a_g.value)
        acList.append(a_c.value)
        arList.append(a_r.value)
        asList.append(a_s.value)
        atList.append(a_t.value)
            
        #determine the error for each calculation including the attenuation components
        error = pubData[-1] - attenList[-1][-1]
        atPE = (attenList[-1][-1] - pubData[-1]) / pubData[-1]
        
        agError = agData[-1] - agList[-1]
        agPE = (a_g.value - agData[-1]) / agData[-1]
        
        acError = acData[-1] - acList[-1]
        acPE = (a_c.value - acData[-1]) / acData[-1]
        
        arError = arData[-1] - arList[-1]
        arPE = (a_r.value - arData[-1]) / arData[-1]
        
        asError = asData[-1] - asList[-1]
        asPE = (a_s.value - asData[-1]) / asData[-1]
        
        #place athe erros for each attenuation into a list
        agErrorList.append(agError)
        agPEList.append(agPE)
        
        acErrorList.append(acError)
        acPEList.append(acPE)
        
        arErrorList.append(arError)
        arPEList.append(arPE)
        
        asErrorList.append(asError)
        asPEList.append(asPE)
        
        errorList.append(error)
        atPEList.append(atPE)
        
        #if the error is >= 0.1 cl
        if abs(error) >= 0.1:
            largeErrors.append(error)
            caseList.append("Case " + str(i))
            
        print("Case " + str(i) + " had an error of " + str(round(error, 5)) + " dB")
    
    
    
    agAvg = sum(agErrorList) / len(agErrorList)
    acAvg = sum(acErrorList) / len(acErrorList)
    arAvg = sum(arErrorList) / len(arErrorList)
    asAvg = sum(asErrorList) / len(asErrorList)
    atAvg = sum(errorList) / len(errorList)
    
    agMax = max(list(map(abs, agErrorList)))
    acMax = max(list(map(abs, acErrorList)))
    arMax = max(list(map(abs, arErrorList)))
    asMax = max(list(map(abs, asErrorList)))
    atMax = max(list(map(abs, errorList)))
    
    agAvgPE = sum(agPEList) / len(agPEList)
    acAvgPE = sum(acPEList) / len(acPEList)
    arAvgPE = sum(arPEList) / len(arPEList)
    asAvgPE = sum(asPEList) / len(asPEList)
    atAvgPE = sum(atPEList) / len(atPEList)
    
        
    ########## ERROR ANALYSIS ##########
    #displays the number of large errors, and the cases that had large errors
    print()
    print("There was a total of "  + str(len(largeErrors)) + " large errors")
    print("These are the magnitudes of each of the errors [dB]")
    
    for i in range(len(largeErrors)):
        print(caseList[i] + " had an error of " + str(largeErrors[i]))
    
    #displays the average error and the max arror for each attenuation component
    print()
    print ("Gaseous Attenuation:")
    print("Average: " + str(agAvg) + " dB")
    print("Max: " + str(agMax) + " dB")
    print("Average Percent Error: " + str(agAvgPE))    
    
    print()
    print ("Cloud Attenuation:")
    print("Average: " + str(acAvg) + " dB")
    print("Max: " + str(acMax) + " dB")
    print("Average Percent Error: " + str(acAvgPE))
    
    print()
    print ("Rain Attenuation:")
    print("Average: " + str(arAvg) + " dB")
    print("Max: " + str(arMax) + " dB")
    print("Average Percent Error: " + str(arAvgPE))
    
    print()
    print ("Scictillation Attenuation:")
    print("Average: " + str(asAvg) + " dB")
    print("Max: " + str(asMax) + " dB")    
    print("Average Percent Error: " + str(asAvgPE))
    
    print()
    print ("Total Attenuation:")
    print("Average: " + str(atAvg) + " dB")
    print("Max: " + str(atMax) + " dB")   
    print("Average Percent Error: " + str(atAvgPE))

    print()
    print()
    print("Total runtime: "  + str(datetime.datetime.now() - begin_time))
    


    
    
    

    

        

    

    


