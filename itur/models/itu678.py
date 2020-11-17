# -*- coding: utf-8 -*-
"""
This file contains functions that follow the methods used to determine the Risk 
metric as described in the reccomendation ITU-R P.678

Created on Mon Nov  5 12:36:06 2020

@author: MAW32652
"""
import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from itur.utils import load_data, dataset_dir
from itur.models.itu1144 import bilinear_2D_interpolator, nearest_2D_interpolator
from itur.models.itu840 import _ITU840_7



def variance_of_estimation_calc(p):
    print(p[0])
    """
    Parameters
    ----------
    p : float
        Probability. The percent of time where the atmospheric loss will exceed
        the atmospheric loss predicted by the model. 

    Returns
    -------
    sigma_e: float
    
        The Variance of Estimation. An intermediate value used to determine 
        Risk. 
    """ 
    N = 525960
    arr = np.linspace((-N + 1), (N - 1), 2*N)
    delta_t  = 60
    a = 0.0265
    b1 = -.0396
    b2 = 0.286   
    b = b1 * (np.log(p)) + b2
    
    if type(p) == float:
        C = np.sum(cu_func(a, b, arr, delta_t))
        
        sigma_e = ((p * (1 - p)) / N )* C 
        return sigma_e
    else:
        sigma_e = np.array([])
        
        #I tried to use vectorized operations but could not get len(b) to broadcast
        #to the 1M len result from cu_func
        #ASK ANDY ABOUT THIS MAKES IT APPEAR SOMEWHAT SLOW
        for i in range(len(b)):
            C = np.sum(cu_func(a, b[i], arr, delta_t))
        
            toAdd = ((p[i] * (1 - p[i])) / N )* C 
            sigma_e= np.append(sigma_e, toAdd)
        
        
        #print("sigma_e is", sigma_e)
        return sigma_e
            
                
def cu_func(a, b, i, delta_t):
    """
    helper function used in variance_of_estimation_calc
    """
    return np.exp(-a * (np.absolute(i * delta_t)**b))


def climatic_ratio_calc(lat, lon):
    """
    

    Parameters
    ----------
    lat : float
        Geopgraphic Latitude
    lon : float
        Geographic Longitude

    Returns
    -------
    rc:
        Climatic Ratio. An intermediate value used to determine risk.

    """
    #load the data that you will be interpolating over
    lats = load_data(os.path.join(dataset_dir, '678/LAT.txt'))
    lons = load_data(os.path.join(dataset_dir, '678/LON.txt'))
    rcs = load_data(os.path.join(dataset_dir, '678/CLIMATIC_RATIO.txt'))
    
    func = bilinear_2D_interpolator(lats, lons, rcs)
    rc = func([-(lat), lon])[0] #lat negative otherwise you will get wrong lat
    return rc



def inter_annual_climatic_variance_calc(lat, lon, p):    
    """
    Parameters
    ----------
    lat : float
        Geographic latitude
    lon : float
        Geographic longitude
    p : float
        Probability. The percent of time where the atmospheric loss will exceed
        the atmospheric loss predicted by the model.     

    Returns
    -------
    sigma_c : float
        inter-annual climatic variance. An intermediate value used to determine Risk.
    """
    rc = climatic_ratio_calc(lat, lon)
    sigma_c = (rc * p)**2
    
    #print("sigma_c is", sigma_c)
    return sigma_c

def sigma_calc(lat, lon, p):
    """
    Parameters
    ----------
    lat : float
        Geographic latitude
    lon : float
        Geographic longitude
    p : float
        Probability. The percent of time where the atmospheric loss will exceed
        the atmospheric loss predicted by the model.     

    Returns
    -------
    sigma : float
        An intermediate value used to determine Risk.
    """
    sigma_e = variance_of_estimation_calc(p)
    sigma_c = inter_annual_climatic_variance_calc(lat, lon, p)
    
    sigma_squared = sigma_e + sigma_c
    sigma = np.sqrt(sigma_squared)
    
    #print("sigma is", sigma)
    return sigma
    

def risk_from_yearly_probability(lat, lon, p, p_yearly):
    """
    Parameters
    ----------
    lat : float
        Geographic latitude.
    lon : float
        Geographic longitude.
    p : float, np.array
        Probability. The percent of time where the atmospheric loss will exceed
        the atmospheric loss predicted by the model over a long term period.
    p_yearly : float
        Probability. The percent of time where the atmospheric loss will exceed
        the atmospheric loss predicted by the model in any given year.    

    Returns
    -------
    risk : float, np.array
        The probability that the yearly probability is exceeded. 
    """
    sigma = sigma_calc(lat, lon, p)
    
    risk = (scipy.special.erfc((p_yearly - p)/ (np.sqrt(2) * sigma))) / 2
    return risk


def p_LT_vs_Risk(lat, lon, p_yearly, target_risk, plot = False):
    """
    Parameters
    ----------
    lat : float
        Geographic latitude.
    lon : float
        Geographic longitude.
    p_yearly : float
        Probability. The percent of time where the atmospheric loss will exceed
        the atmospheric loss predicted by the model in any given year. 
    plot : Boolean
        Whether or not you want to produce a plot of p_LT vs Risk at p_yearly
        Does not work if p_yearly is an array of values
    
    Returns
    -------
    Plot of p_LT vs Risk
    p_at_target_risk : float
        Probability that is required to achieve your target risk
    """
    #create the pArr array based off of p_yearly
    
    #determine risk values
    if type(p_yearly) == float:
        print("p_yearly is a float") 
        pArr = np.linspace((p_yearly / 10), p_yearly, 50)
        riskArr = risk_from_yearly_probability(lat, lon, pArr, p_yearly)
        
        #determine the p that gives you your target risk
        func = scipy.interpolate.interp1d(riskArr, pArr)
        p_at_target_risk = func([target_risk])[0]
        
        if plot:        
        #plot p_LT vs risk
        plt.plot(pArr * 100, riskArr*100)
        plt.title("p_LT vs Risk to Achieve p_yearly = "+ str(p_yearly*100)+ "%")
        plt.xlabel("p_LT [%]")
        plt.ylabel("Risk [%]")
        plt.plot(p_at_target_risk*100, target_risk*100, marker = 'o', markersize = 3, color = 'red')
    
    else:
        print("p_yearly is an array")
        p_at_target_risk = np.array([])
        for i in range(len(p_yearly)):
            pArr = np.linspace((p_yearly[i] / 10), p_yearly[i], 50) 
            riskArr = risk_from_yearly_probability(lat, lon, pArr, p_yearly[i])
            
            #determine the p that gives you your target risk
            func = scipy.interpolate.interp1d(riskArr, pArr)
            toAdd = func([target_risk])[0]
            p_at_target_risk = np.append(p_at_target_risk, toAdd)           

    return p_at_target_risk
    
   
    