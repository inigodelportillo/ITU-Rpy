# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import os
from astropy import units as u

from itur import utils
from itur.models.itu1144 import bilinear_2D_interpolator
from itur.utils import load_data, dataset_dir, prepare_input_array,\
    prepare_output_array, memory


class __ITU839():
    """Rain height model for prediction methods.

    Not available versions:
    * P.839-0 (03/92) (Superseded)
    * P.839-1 (83/97) (Superseded)
    Available versions include:
    * P.839-2 (10/99) (Superseded)
    * P.839-3 (02/01) (Superseded)
    * P.839-4 (09/2013) (Current version)
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.839 recommendation.

    def __init__(self, version=4):
        if version == 4:
            self.instance = _ITU839_4()
        elif version == 3:
            self.instance = _ITU839_3()
        elif version == 2:
            self.instance = _ITU839_2()
#        elif version == 1:
#            self.instance = _ITU839_1()
#        elif version == 0:
#            self.instance = _ITU839_0()
        else:
            raise ValueError(
                'Version {0} is not implemented for the ITU-R P.839 model.'
                .format(version))

        self._zero_isoterm_data = {}

    @property
    def __version__(self):
        return self.instance.__version__

    def rain_height(self, lat, lon):
        # Abstract method to compute the rain height
        return self.instance.rain_height(lat, lon)

    def isoterm_0(self, lat, lon):
        # Abstract method to compute the zero isoterm height
        return self.instance.isoterm_0(lat, lon)


class _ITU839_4():

    def __init__(self):
        self.__version__ = 4
        self.year = 2013
        self.month = 9
        self.link = 'https://www.itu.int/rec/R-REC-P.839/' +\
                    'recommendation.asp?lang=en&parent=R-REC-P.839-4-201309-I'

        self._zero_isoterm_data = {}

    def isoterm_0(self, lat, lon):
        if not self._zero_isoterm_data:
            vals = load_data(os.path.join(dataset_dir,
                                          '839/v4_ESA0HEIGHT.txt'))
            lats = load_data(os.path.join(dataset_dir, '839/v4_ESALAT.txt'))
            lons = load_data(os.path.join(dataset_dir, '839/v4_ESALON.txt'))
            self._zero_isoterm_data = bilinear_2D_interpolator(
                lats, lons, vals)

        return self._zero_isoterm_data(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def rain_height(self, lat_d, lon_d):
        """
        The rain height is computed as

        ..math:
            h_r = h_0 + 0.36 (km)
        """
        return self.isoterm_0(lat_d, lon_d) + 0.36


class _ITU839_3():

    def __init__(self):
        self.__version__ = 3
        self.year = 2001
        self.month = 2
        self.link = 'https://www.itu.int/rec/R-REC-P.839-3-200102-S/en'

        self._zero_isoterm_data = {}

    def isoterm_0(self, lat, lon):
        if not self._zero_isoterm_data:
            vals = load_data(os.path.join(dataset_dir,
                                          '839/v3_ESA0HEIGHT.txt'))
            lats = load_data(os.path.join(dataset_dir, '839/v3_ESALAT.txt'))
            lons = load_data(os.path.join(dataset_dir, '839/v3_ESALON.txt'))
            self._zero_isoterm_data = bilinear_2D_interpolator(
                lats, lons, vals)

        return self._zero_isoterm_data(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def rain_height(self, lat_d, lon_d):
        """
        The rain height is computed as

        ..math:
            h_r = h_0 + 0.36 (km)
        """
        return self.isoterm_0(lat_d, lon_d) + 0.36


class _ITU839_2():

    def __init__(self):
        self.__version__ = 2
        self.year = 1999
        self.month = 10
        self.link = 'https://www.itu.int/rec/R-REC-P.839-2-199910-S/en'

    def isoterm_0(self, lat_d, lon_d):
        """
        The 0C mean isotherm height can be approximated as

        """
        # TODO: Complete this with the equation in ITU-R P.839-2
        h0 = np.where(
            lat_d > 23, 5 - 0.075 * (lat_d - 23),
            np.where(
                np.logical_and(0 < lat_d, lat_d < 23),
                5, np.where(
                    np.logical_and(-21 < lat_d, lat_d < 0),
                    5, np.where(
                        np.logical_and(-71 < lat_d, lat_d < -21),
                        5 + 0.1 * (lat_d + 21),
                        0))))
        return h0

    def rain_height(self, lat_d, lon_d):
        """
        For areas of the world where no specific information is available,
        the mean rain height, may be approximated by the mean 0C isotherm
        height, and for for North America and for Europe west of 60° E
        longitude the mean rain height is approximated by

        ..math:
            h_r = 3.2 - 0.075 (\\lambda - 35) \\qquad for \\qquad
            35 \\le \\lambda \\le 70  (km)
        """
        h0 = self.isoterm_0(lat_d, lon_d)
        return np.where(np.logical_and(np.logical_and(35 < lat_d, lat_d < 70),
                                       lon_d < 60),
                        3.2 - 0.075 * (lat_d - 35), h0)


__model = __ITU839()


def change_version(new_version):
    """
    Change the version of the ITU-R P.839 recommendation currently being used.


    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
        * version 0: P.839-0 (03/92) (Superseded)
        * version 1: P.839-1 (83/97) (Superseded)
        * version 2: P.839-2 (10/99) (Superseded)
        * version 3: P.839-3 (02/01) (Superseded)
        * version 4: P.839-4 (09/2013) (Current version)
    """
    global __model
    __model = __ITU839(new_version)
    utils.memory.clear()


def get_version():
    """
    Obtain the version of the ITU-R P.839 recommendation currently being used.
    """
    global __model
    return __model.__version__


@memory.cache
def isoterm_0(lat, lon):
    """
    A method to estimate the zero isoterm height for propagation prediction.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points


    Returns
    -------
    zero_isoterm: numpy.ndarray
        Zero isoterm height (km)


    References
    ----------
    [1] Rain height model for prediction methods:
    https://www.itu.int/rec/R-REC-P.839/en

    """
    global __model
    type_output = type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.isoterm_0(lat, lon)
    return prepare_output_array(val, type_output) * u.km


@memory.cache
def rain_height(lat, lon):
    """
    A method to estimate the rain height for propagation prediction.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points


    Returns
    -------
    rain_height: numpy.ndarray
        Rain height (km)


    References
    ----------
    [1] Rain height model for prediction methods:
    https://www.itu.int/rec/R-REC-P.839/en
    """
    global __model
    type_output = type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.rain_height(lat, lon)
    return prepare_output_array(val, type_output) * u.km
