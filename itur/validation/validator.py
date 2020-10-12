# -*- coding: utf-8 -*-
"""
validator.py

This validator file is the master validation file that is capable of calling all
the validation scripts in /validation/validation_scripts

Each of the validation scripts determines the error for certain values.
This error is calculated by comparing the expected values puiblished by the ITU
to the values computed by the ITU-Rpy module for 64 validation cases. 
An error result will be computed for each of these cases and statistics such as 
averge error, max error, and average percent error will be printed to the console.

The validation scripts are called using functions.
The validation scripts that are currently supported are as follows:

-------------------------------------------------------------------------------
    total_attenuation_validation()
        
        This function determines the error for the total attenuation as well 
        as the error for each of the attenuation components: gaseous attenuation,
        cloud attenuation, rain attenuation, and scintillation attenuation.
        
        This function validates off of the P618-13 Att_Tot Sheet. 
        
        Outputs:
            
            Individually prints out the error for each of the 64 validation cases
            
            Prints the number of cases with large errors (error > 0.1) and which
            cases produced these errors
            
            Error statistics for gaseous attenuation, cloud attenuation, rain
            attenuation, scintillation attenuation, and total attenuation
-------------------------------------------------------------------------------
    gaseous_attenuation_validation()
    
        This function determines the error for the gaseous attenuation as well 
        as the error for intermediate values that are used to determine the
        gaseous attenuation.
        
        This function validates off of the P676-12 A_gas sheet.
        
        Outputs:
            
            Error statistics for attenuation due to oxygen and associated intermediate
            values: gamma_o, Oxygen equivalent height (h_o), attenuation due to
            oxygen, attenuation due to oxygen with sin(elevation angle) factor
            
            Error statistics for attenuation due to water vapor and associated intermediate
            values: Attenuation due to water vapor, attenuation due to watervapor with
            sin(elevation angle) factor, rho_ref, t_ref.
            
            Error statistics for the total gasoues attenuation
-------------------------------------------------------------------------------
    gamma_intermediate_value_validation()
    
        This function determines the error for the gamma_o and gamma_w intermediate
        values that are used to determine the total gaseous attenuation
        
        This function does not validate of the regular 64 validation cases.
        Instead, it validates off of the P676-12 SpAtt sheet. 
        
        Outputs:
            
            Error statistics for the gamma intermediate values used to determine
            the total gaseous attenuation: gamma_o, and gamma_w.
-------------------------------------------------------------------------------
    cloud_attenuation_validation()
        
        This function determines the error for the cloud attenuation as well as
        the error for the intermediate values that are used to determine the 
        cloud attenuation.
        
        This function validates off of the P840-8 A_Clouds sheet. 
        
        Outputs:
            
            Error statistics for the Kl intermediate value as well as the intermediate
            values that are used to determine Kl: epsilon prime, epsilon prime prime, 
            eta. 
            
            Error Statistics for the L_red intermediate value.
-------------------------------------------------------------------------------

Created on Fri Jul 16 15:15:29 2020

@author: MAW32652
"""

### importing the various validation scripts 
import itur
from itur.validation.validation_scripts.total_attenuation_validation import total_attenuation_validation
from itur.validation.validation_scripts.A_gas_validation import gaseous_attenuation_validation
from itur.validation.validation_scripts.gamma_validation_script import gamma_intermediate_value_validation
from itur.validation.validation_scripts.cloud_att_intermediate_values import cloud_attenuation_validation




