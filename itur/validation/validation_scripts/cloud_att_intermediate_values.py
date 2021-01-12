# -*- coding: utf-8 -*-
"""
cloud_att_intermediate_values.py

Created on Tue Jun 30 8:53:09 2020

Determined the error between the published L_red value in the ITU-R validation data,
sheet P840-8 Lred, and the calculated L_red value using iturpy

@author: MAW32652
"""
import itur
from itur.models.itu840 import columnar_content_reduced_liquid, specific_attenuation_coefficients
from itur.utils import prepare_quantity, prepare_output_array,\
    prepare_input_array, load_data, dataset_dir, memory
import os
import numpy as np
import xlrd

def cloud_attenuation_validation():
    
    ########## VARIABLE INITIALIZATION ###########
    
    ###excel variables
    
    workbook = xlrd.open_workbook("CG-3M3J-13-ValEx-Rev5_0.xlsx")
    sheet = workbook.sheet_by_name("P840-8 A_Clouds")
    
    ### ITU_rpy variables
    
    # Validation data V=variables
    epData = [] #epislon prime expected values
    eppData= [] #epsilon prime prime expected values
    etaData = [] #eta expected values
    klData = [] #kl expected values
    LredData = [] #L_red expected values
    
    # Input variables
    latList = [] #Latitude
    lonList = [] #Longitude
    fList = [] #frequency
    eleList = [] #Elevation angle
    pList = [] #Probability
    
    #Output Variables
    
    epList = [] #epsilon prime calculated values
    epErrorList = [] #error for epsilon prime (Validation - calculated)
    epPEList = [] #Epsilon prime percent error
    eppList = [] #epsilon prime prime calculated values
    eppErrorList = [] #error for epsilon prime prime (Validation - calculated)
    eppPEList = [] #epsilon prime prime percent error
    etaList= [] #eta calculated values
    etaErrorList = []  #error for eta (Validation - calculated)
    etaPEList = [] #eta percent error
    
    klList = [] #kl calculated values
    klErrorList = []  #error for kl (Validation - calculated)
    klPEList = [] # kl percent error
    LredList = [] #L_red calculated values
    LredErrorList= [] #error for L_red (Validation - calculated)
    LredPEList = [] # L_red percent error
    
    ########## Computation ##########
    for i in range(63):
        
        #populate the input variable lists
        #data entries start on row 21
           
        #expected outputs
        epData.append(sheet.cell_value(i + 20, 8)) #Published epsilon prime data in column I
        eppData.append(sheet.cell_value(i + 20, 9)) #Published epsilon prime prime data in column J
        etaData.append(sheet.cell_value(i + 20, 10)) #Published eta data in column K
        klData.append(sheet.cell_value(i + 20, 11)) #Published kl data in column L
        LredData.append(sheet.cell_value(i + 20, 12)) #Published L_red data in column M
           
        #input data
        latList.append(sheet.cell_value(i + 20, 3)) #Latitdue inputs in column D
        lonList.append(sheet.cell_value(i + 20, 4)) #Longitude inputs in column E
        fList.append(sheet.cell_value(i + 20, 5)) #Frequency inputs in column F
        eleList.append(sheet.cell_value(i + 20, 6)) #Elevation angle inputs in column G
        pList.append(sheet.cell_value(i + 20, 7)) #probability inputs in column H
    
        
        ### kl ###
       
        #detrmine kl and the associated intermediate values using P.840 functions
        ep, epp, eta, kl = specific_attenuation_coefficients(fList[-1], T = 0)
        
        epList.append(ep)
        eppList.append(epp)
        etaList.append(eta)
        klList.append(kl)
        
        #determine the error of the new calculated values
        epError = epData [-1] - ep
        epPE = (ep - epData[-1]) / epData[-1]
        epErrorList.append(epError)
        epPEList.append(epPE)
        
        eppError = eppData [-1] - epp
        eppPE = (epp - eppData[-1]) / eppData[-1] 
        eppErrorList.append(eppError)
        eppPEList.append(eppPE)
        
        etaError = etaData [-1] - eta
        etaPE = (eta - etaData[-1]) / etaData[-1] 
        etaErrorList.append(etaError)
        etaPEList.append(etaPE)
        
        klError = klData[-1] - kl
        klPE = (kl - klData[-1]) / klData[-1] 
        klErrorList.append(klError)
        klPEList.append(klPE)
                
        ### L_red ###
        
        #Calculate L red and put the value into the output list
        L_red = columnar_content_reduced_liquid(latList[-1], lonList[-1], pList[-1])
        LredList.append(L_red.value)
        
        #determine the error for the new calculated value
        LredError = LredData[-1] - L_red.value
        LredPE = (L_red.value - LredData[-1])/ LredData[-1]
        LredErrorList.append(LredError)
        LredPEList.append(LredPE)
        
        
    ########### ERROR ANALYSIS ##########
    
    epAvg = sum(epErrorList)/len(epErrorList)
    eppAvg = sum(eppErrorList)/len(eppErrorList)
    etaAvg = sum(etaErrorList)/len(etaErrorList)
    klAvg = sum(klErrorList)/len(klErrorList)
    LredAvg = sum(LredErrorList)/len(LredErrorList)
    
    epMax = max(list(map(abs, epErrorList)))
    eppMax = max(list(map(abs, eppErrorList)))
    etaMax = max(list(map(abs, etaErrorList)))
    klMax = max(list(map(abs, klErrorList)))
    LredMax = max(list(map(abs, LredErrorList)))
    
    epAvgPE = sum(epPEList)/len(epPEList)
    eppAvgPE = sum(eppPEList)/len(eppPEList)
    etaAvgPE = sum(etaPEList)/len(etaPEList)
    klAvgPE = sum(klPEList)/len(klPEList)
    LredAvgPE = sum(LredPEList)/len(LredPEList)
    
    print()  
    print("Epsilon prime: ")
    print("Average: " + str(epAvg))
    print("Max: " + str(epMax))
    print("Average Percent Error: " + str(epAvgPE))
    
    print()  
    print("Epsilon prime prime: ")
    print("Average: " + str(eppAvg))
    print("Max: " + str(eppMax))
    print("Average Percent Error: " + str(eppAvgPE))
    
    print()  
    print("Eta: ")
    print("Average: " + str(etaAvg))
    print("Max: " + str(etaMax))
    print("Average Percent Error: " + str(etaAvgPE))
    
    print()
    print("Note: Kl is a function of epsilon prime, epsilon prime prime, and eta")
    print()  
    print("Kl: ")
    print("Average: " + str(klAvg))
    print("Max: " + str(klMax))
    print("Average Percent Error: " + str(klAvgPE))
    
    print ()
    print("Note: L_red is not related to Kl")
    print()  
    print("L_red: ")
    print("Average: " + str(LredAvg))
    print("Max: " + str(LredMax))
    print("Average Percent Error: " + str(LredAvgPE))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    