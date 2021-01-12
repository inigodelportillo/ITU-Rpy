# -*- coding: utf-8 -*-
"""
This file contains functions that follow the methods used to determine the Risk 
metric and other associated intermediate values as described in the 
reccomendation ITU-R P.678

Created on Mon Nov  5 12:36:06 2020

@author: MAW32652
"""
import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from itur.utils import load_data, dataset_dir
from itur.models.itu1144 import bilinear_2D_interpolator
from joblib import Memory
location = './cachedir'
memory = Memory(location, verbose=0)
import time

def variance_of_estimation_calc(p, N = 525960):
    """
    Parameters
    ----------
    p : float
        Probability. The percent of time where the atmospheric loss will exceed
        the atmospheric loss predicted by the model. 
        Values : 0 - 1
    Returns
    -------
    sigma_e: float
    
        The Variance of Estimation. An intermediate value used to determine 
        Risk. 
    """ 
    siny = 31557600
    #N = 525960
    arr = np.linspace((-N + 1), (N - 1), (2*N) - 1)
    delta_t  = siny / N
    a = 0.0265
    b1 = -.0396
    b2 = 0.286   
    b = b1 * (np.log(p)) + b2
    
    if type(p) == float:      
        C = np.sum(cu_func(a, b, arr, delta_t))
        sigma_e = ((p * (1 - p)) / N )* C 
        return sigma_e
   
    else:        
        b = b.reshape(len(b), 1) #reshape b so that it can broadcast correctly        
        result = cu_func(a, b, arr, delta_t)
        C = np.sum(result, axis = 1)
        sigma_e = ((p * (1 - p)) / N )* C 
        
    return sigma_e

def variance_of_estimation_calc2(p, N = 525960):
    N = 525960
    # the function is symmetric, so only need to do one side of the sum (half of eqn 2)
    arr = np.linspace(0, (N - 1), N).astype(np.float32)  # include zero
    delta_t  = 60
    a = 0.0265
    b1 = -.0396
    b2 = 0.286   
    b = ( b1 * (np.log(p)) + b2 ).astype(np.float32)
    
    if type(p) == float:
        C = cu_func(a, b, arr, delta_t)
        
        sigma_e = ((p * (1 - p)) / N )* C 
        return sigma_e
    else:
        sigma_e = np.array([])
        idt = np.absolute(arr * delta_t)
        
        
        # this is the slowest line driving the whole slowdown
        start = time.time()
        # for i = 0 to N-1
        expThis = -a * ( idt.reshape([len(idt),1]) ** b.reshape([1,len(b)]) )
        end = time.time()
        print('The power function took {:.2f} s to compute.'.format(end - start))
        
        # (os meaning one sided)
        osExp = np.exp( expThis ).astype(np.float64)
        
        # exp(i=0)=1 always.  So take twice the i>0 sum plus 1 to get eqn2 summation result C
        start = time.time()
        C = 1 + 2 * np.sum( osExp[1:,:], axis=0 )
        end = time.time()
        print('The sum function took {:.2f} s to compute.'.format(end - start))
        
        for i in range(np.shape(C)[0]):
            toAdd = ((p[i] * (1 - p[i])) / N )* C[i]
            sigma_e= np.append(sigma_e, toAdd)
            
        return sigma_e    
                
def cu_func(a, b, i, delta_t):
    """
    Helper function used in variance_of_estimation_calc
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
    lats = load_data(os.path.join(dataset_dir, '678/2015_LAT.txt'))
    lons = load_data(os.path.join(dataset_dir, '678/2015_LON.txt'))
    rcs = load_data(os.path.join(dataset_dir, '678/2015_CLIMATIC_RATIO.txt'))
    
    func = bilinear_2D_interpolator(lats, lons, rcs)
    rc = func([-(lat), lon])[0] #lat negative otherwise you will get wrong lat
    #rc = func((lat, lon))
    #print(rc)
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
        Values : 0 - 1
    Returns
    -------
    sigma_c : float
        inter-annual climatic variance. An intermediate value used to determine Risk.
    """
    climatic_ratio_calc_cached = memory.cache(climatic_ratio_calc)
    rc = climatic_ratio_calc_cached(lat, lon)
    sigma_c = (rc * p)**2
    
    return sigma_c

def sigma_calc(lat, lon, p, N = 525960):
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
        Values : 0 - 1
    Returns
    -------
    sigma : float
        An intermediate value used to determine Risk.
    """
    sigma_e = variance_of_estimation_calc(p, N)
    sigma_c = inter_annual_climatic_variance_calc(lat, lon, p)
    
    sigma_squared = sigma_e + sigma_c
    sigma = np.sqrt(sigma_squared)
    
    return sigma
    

def risk_from_yearly_probability(lat, lon, p, p_yearly, N = 525960):
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
        Values : 0 - 1
    p_yearly : float
        Probability. The percent of time where the atmospheric loss will exceed
        the atmospheric loss predicted by the model in any given year.    
        Values : 0 - 1
    Returns
    -------
    risk : float, np.array
        The probability that the yearly probability is exceeded. 
    """
    sigma = sigma_calc(lat, lon, p, N)
    
    risk = (scipy.special.erfc((p_yearly - p)/ (np.sqrt(2) * sigma))) / 2
    return risk


def p_lt_from_risk(lat, lon, p_yearly, target_risk, N = 525960, plot = False):
    """
    Parameters
    ----------
    lat : float
        Geographic latitude.
    lon : float
        Geographic longitude.
    p_yearly : float, np.array
        Probability. The percent of time where the atmospheric loss will exceed
        the atmospheric loss predicted by the model in any given year.
        Values : 0 - 100 %
    target_risk : float
        risk you want to input to determine p_longterm at
        Values : 0 - 100 %
    plot : Boolean
        Whether or not you want to produce a plot of p_LT vs Risk at p_yearly
        Does not work if p_yearly is an array of values
    
    Returns
    -------
    Plot of p_LT vs Risk
    p_lt : float
        Probability that is required to achieve your target risk
    """
    p_yearly = p_yearly / 100 #calc requires input to be in decimal
    target_risk = target_risk / 100
    #determine risk values
    if type(p_yearly) == float:
        pArr = np.linspace((p_yearly / 10), (p_yearly * 10), 50)
        riskArr = risk_from_yearly_probability(lat, lon, pArr, p_yearly, N)
        
        #determine the p that gives you your target risk
        func = scipy.interpolate.interp1d(riskArr, pArr)
        p_lt = func([target_risk])[0]
        
        if plot:        
            #plot p_LT vs risk
            plt.plot(pArr * 100, riskArr*100)
            plt.title("p_LT vs Risk to Achieve p_yearly = "+ str(p_yearly*100)+ "%")
            plt.xlabel("p_LT [%]")
            plt.ylabel("Risk [%]")
            plt.plot(p_lt*100, target_risk*100, marker = 'o', markersize = 3, color = 'red')
    
    else:
        p_lt = np.array([])
        for i in range(len(p_yearly)):
            pArr = np.linspace((p_yearly[i] / 10), (p_yearly[i] * 10), 50) 
            riskArr = risk_from_yearly_probability(lat, lon, pArr, p_yearly[i])
            
            #determine the p that gives you your target risk
            func = scipy.interpolate.interp1d(riskArr, pArr)
            toAdd = func([target_risk])[0]
            p_lt = np.append(p_lt, toAdd)           

    return p_lt * 100 #outputing as a percentage because rest of iturpy uses % not decimal
    