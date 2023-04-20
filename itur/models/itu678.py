# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from astropy import units as u
from scipy.special import erfc

from itur.models.itu1144 import bilinear_2D_interpolator
from itur.utils import (prepare_input_array, prepare_output_array,
                        load_data_interpolator, get_input_type)

class __ITU678__():

    """Characterization of the variability of propagation phenomena and
    estimation of the risk associated with propagation margin.

    Not available versions:
    * P.678-1 (03/92) (Superseded)
    * P.678-2 (09/2013) (Superseded)
    Available versions include:
    * P.678-3 (07/2015) (Current version)
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.678 recommendation.

    def __init__(self, version=3):
        if version == 3:
            self.instance = _ITU678_3_()
#        elif version == 1:
#            self.instance = _ITU678_2()
#        elif version == 0:
#            self.instance = _ITU678_1()
        else:
            raise ValueError(
                f"Version {version} is not implemented for the ITU-R P.678 model."
            )

        self._climatic_ratio_data = {}

    @property
    def __version__(self):
        return self.instance.__version__

    def inter_annual_variability(self, p, lat, lon):
        # Abstract method to compute the inter-annual variability
        return self.instance.inter_annual_variability(p, lat, lon)

class _ITU678_3_():

    def __init__(self):
        self.__version__ = 3
        self.year = 2015
        self.month = 7
        self.link = 'https://www.itu.int/rec/R-REC-P.678-3-201507-I/en'

        self._climatic_ratio_data = {}
        
    def climatic_ratio(self, lat, lon):
        if not self._climatic_ratio_data:
            self._climatic_ratio_data = load_data_interpolator(
                '678/v3_lat.npz', '678/v3_lon.npz',
                '678/v3_climatic_ratio.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return self._climatic_ratio_data(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)
        
    def variance_of_estimation(self, p):
        # Eq. (4)
        b1 = -0.0396
        b2 = 0.286
        b = b1*np.log(p)+b2
        a = 0.0265
        
        # Eq. (3)
        N = 525960
        dt = 60
        i = np.arange(-N+1, N)[:, np.newaxis]
        cu = np.exp(-a*np.abs(i*dt)**b)
        
        # Eq. (2)
        C = np.sum(cu, axis=0)

        # Eq. (5)
        return p*(1-p)*C/N
        
    def inter_annual_climatic_variance(self, p, lat, lon):
        rc = self.climatic_ratio(lat, lon)
        
        # Eq. (6)
        return (rc*p)**2
        
    def inter_annual_variability(self, p, lat, lon):
        """
        The inter-annual variance.
        """
        # Eq. (1)
        return (self.inter_annual_climatic_variance(p, lat, lon)
            + self.variance_of_estimation(p)).squeeze()
    
    def risk_of_exceedance(self, p, pr, lat, lon):
        """
        The risk to exceed a probability of exceedance.
        """
        # Eq. (8)
        return (0.5*erfc((pr-p)/self.inter_annual_variability(p, lat, lon)
            / np.sqrt(2))).squeeze()

__model = __ITU678__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.878 recommendation currently being used.

    This function changes the model used for the ITU-R P.878 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          * 3: Activates recommendation ITU-R P.678-3 (07/2015) (Current version)

    """
    global __model
    __model = __ITU678__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.839 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__


def inter_annual_variability(p, lat, lon):
    """
    Estimate the variance of the inter-annual variability for a given time
    percentage, at a given location.


    Parameters
    ----------
    p : number
        Probability level (%).
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points (degrees north).
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points (degrees east).


    Returns
    -------
    sigma: numpy.ndarray
        Variance of the inter-annual variability (%).


    References
    ----------
    [1] Characterization of the variability of propagation phenomena and
    estimation of the risk associated with propagation margin:
    https://www.itu.int/rec/R-REC-P.678/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.inter_annual_variability(p*100, lat, lon) / 100
    return prepare_output_array(val, type_output) * u.pct

def risk_of_exceedance(p, pr, lat, lon):
    """
    Estimate the risk (or probability) that an attenuation threshold exceeded
    for p % of the time on a long-term basis is exceeded for more than pr % of
    the time on an annual basis.


    Parameters
    ----------
    p : number
        Long term probability of exceedance (%).
    pr : number
        Annual proportion of exceedance (%).
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points (degrees north).
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points (degrees east).


    Returns
    -------
    r: numpy.ndarray
        Risk (or probability) to exceed the proportion pr on a annual basis (%).


    References
    ----------
    [1] Characterization of the variability of propagation phenomena and
    estimation of the risk associated with propagation margin:
    https://www.itu.int/rec/R-REC-P.678/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.risk_of_exceedance(p*100, pr*100, lat, lon) / 100
    return prepare_output_array(val, type_output) * u.pct