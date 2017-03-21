# -*- coding: utf-8 -*-
import numpy as np
from utils import prepare_input_array, prepare_quantity
from astropy import units as u
from itur837 import rainfall_rate
from itur838 import specific_attenuation, specific_attenuation_coefficients

class __ITU530():
    """Reference Standard Atmospheres
    
    Available versions:
    * P.530-16 (07/15) (Current version)
    
    Not available versions:
    * P.530-1 (08/94) (Superseded)
    * P.530-2 (08/97) (Superseded)
    * P.530-3 (10/99) (Superseded)
    * P.530-4 (03/05) (Superseded)
    
    The procedures to compute the reference standard atmosphere parameters
    pressented in these versions are identical to those included in version
    ITU_T P.530-5. Version 3 includes a dataset with vertical profiles for
    353 locations over the world using 10 years of radiosonde observations 
    . Version 4 includes another dataset with monthly vertical profiles worldwide
    (grid of 1.5 by 1.5 deg). None of these are currently implemente but TODO
    work.
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.530 recommendation.    
    def __init__(self, version = 16):
        if version == 16:
            self.instance = _ITU530_16()
        else:
            raise ValueError('Version ' + str(version) + ' is not implemented'+
            ' for the ITU-R P.530 model.')

    @property
    def __version__(self):
        return self.instance.__version__
   
        

class _ITU530_16():
    
    def __init__(self):
        self.__version__ = 16
        self.year = 2015
        self.month = 7
        self.link = 'https://www.itu.int/rec/R-REC-P.530-16-201507-S/en'   

    def multipath_loss_for_A(self, h_e, h_r, d, f, A):
        # h_e antena emitter
        # h_r antenna receiver
        # d pathlength (km)
        # f : frequency (GHz
        # A : Attenuation (dB)
        pass
        # Step 1: Estimate the geoclimatic factor K
        # DN1 point refractivity gradient in the lowest 65 m of the atmospher 
        # not exceeded for 1% of an average year
        # s_a is the area terrain roughness
        K = 10**(-4.4 - 0.0027 * dN1)*(10 + s_a)**(-0.46)
        
        # Step 2: Claculate the magnitude of the path inclination
        e_p = np.abs(h_r - h_e)/d #change to mrad
        
        #Step 3: For detailed link design applications calculate the percentage
        # of time (p_W) that fade depth A (dB) is exceeded in the average worst
        #month
        h_L = np.minimum(h_e, h_r)
        p_W = K *d**3.4 (1 + e_p)**-1.03 * f**0.8 * 10**(-0.00076*h_L- A/10)
        
        return p_W
        
    def multipath_loss(self, h_e, h_r, d, f, A):
        
        # Step 1: Using the method multipath_loss_for_A calculate the 
        # multipath occurrence factor, p0
        p0 = self.multipath_loss_for_A(h_e, h_r, d, f, 0)
        
        # Step 2: Calculate the value of fade depth, At, at which the transition
        # occurs between the deep-fading distribution and the shallow-fading
        # distribution 
        At = 25 + 1.2 * np.log10(p0)
        
        # Step 3: Calculate the percentage of time that A is exceeded in the 
        # average worst month:
        def step_3b(p_0, At, A):
            p_t = p_0* 10 **(-At / 10)
            qa_p = -20 * np.log10( -np.log((100 - p_t)/100))/At
            q_t = (qa_p - 2) / (1 + 0.3 * 10 ** (-At/20)*10**(-0.016*At)) +\
                         - 4.3 * (10**(-At/20) + At/800)
            q_a = 2 + (1 + 0.3 * 10**(-A / 20)) * (10**(-0.016 * A)) *\
                        (q_t + 4.3*( 10**(-A/20 + A/800)))
            p_W = 100 * ( 1 - np.exp(-10 ** (-q_a * A / 20)))
            return p_W
            
        p_W = np.where(A >= At, p0* 10 **(-A /10), step_3b(p0, At, A))
        return p_W
        
    def rain_attenuation(lat, lon, d, f, el, p, tau=45, R001 = None):
        # Step 1: Obtain the rain rate R0.01 exceeded for 0.01% of the time
        # (with an integration time of 1 min). 
        if R001 is None:
            R001 = rainfall_rate(lat, lon, 0.01)
        
        # Step 2: Compute the specific attenuation, gammar (dB/km) for the
        # frequency, polarization and rain rate of interest using 
        # Recommendation ITU-R P.838            
        gammar = specific_attenuation(R001, f, el, tau)
        _, alpha = specific_attenuation_coefficients(f, el, tau)
        
        # Step 3: Compute the effective path length, deff, of the link by 
        # multiplying the actual path length d by a distance factor r
        r = 1 / ( 0.477 * d **0.633 * R001 ** (0.073* alpha) * f**(0.123) - 10.579 * (1 - np.exp(-0.024*d)))
        deff = np.minimum(r, 2.5)
        
        # Step 4: An estimate of the path attenuation exceeded for 0.01% of 
        # the time is given by:
        A001 = gammar * deff
        
        # Step 5: The attenuation exceeded for other percentages of time p in
        # the range 0.001% to 1% may be deduced from the following power law
        C0 = np.where( f >= 10, 0.12 * 0.4 * ( np.log10(f / 10)**0.8), 0.12)
        C1 = (0.07**C0)*(0.12**(1-C0))
        C2 = 0.855*C0 + 0.546*(1-C0)
        C3 = 0.139 * C0 + 0.043 * (1 - C0)
        Ap = A001 * C1 * p ** (- (C2 + C3 * np.log10(p) ))
        return Ap
        
    def inverse_rain_attenuation(lat, lon, d, f, el, Ap, tau=45, R001 = None):
        # Step 1: Obtain the rain rate R0.01 exceeded for 0.01% of the time
        # (with an integration time of 1 min). 
        if R001 is None:
            R001 = rainfall_rate(lat, lon, 0.01)
        
        # Step 2: Compute the specific attenuation, gammar (dB/km) for the
        # frequency, polarization and rain rate of interest using 
        # Recommendation ITU-R P.838            
        gammar = specific_attenuation(R001, f, el, tau)
        _, alpha = specific_attenuation_coefficients(f, el, tau)
        
        # Step 3: Compute the effective path length, deff, of the link by 
        # multiplying the actual path length d by a distance factor r
        r = 1 / ( 0.477 * d **0.633 * R001 ** (0.073* alpha) * f**(0.123) - 10.579 * (1 - np.exp(-0.024*d)))
        deff = np.minimum(r, 2.5)
        
        # Step 4: An estimate of the path attenuation exceeded for 0.01% of 
        # the time is given by:
        A001 = gammar * deff
        
        # Step 5: The attenuation exceeded for other percentages of time p in
        # the range 0.001% to 1% may be deduced from the following power law
        C0 = np.where( f >= 10, 0.12 * 0.4 * ( np.log10(f / 10)**0.8), 0.12)
        C1 = (0.07**C0)*(0.12**(1-C0))
        C2 = 0.855*C0 + 0.546*(1-C0)
        C3 = 0.139 * C0 + 0.043 * (1 - C0)
        Ap = A001 * C1 * p ** (- (C2 + C3 * np.log10(p) ))
        return Ap        
        
        
    def rain_event_count(self, lat, lon, d, f, el, A, tau=45, R001 = None):
        # Compute the the percentage of time that the rain attenuation A(dB) 
        # exceeded in the average year. 
        p_A = self.inverse_rain_attenuation(lat, lon, d, f, el, A)
        
        # The number of fade events exceeding attenuation A for 10 s or longer
        N10s = 1 + 1313*p_A**0.945
        
        return N10s
        
        
        