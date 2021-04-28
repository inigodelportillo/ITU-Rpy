# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np
from astropy import units as u

from itur.models.itu1144 import bilinear_2D_interpolator
from itur.utils import (prepare_input_array, prepare_output_array,
                        load_data_interpolator, get_input_type)


class __ITU839__():

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
            self.instance = _ITU839_4_()
        elif version == 3:
            self.instance = _ITU839_3_()
        elif version == 2:
            self.instance = _ITU839_2_()
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


class _ITU839_4_():

    def __init__(self):
        self.__version__ = 4
        self.year = 2013
        self.month = 9
        self.link = 'https://www.itu.int/rec/R-REC-P.839/' +\
                    'recommendation.asp?lang=en&parent=R-REC-P.839-4-201309-I'

        self._zero_isoterm_data = {}

    def isoterm_0(self, lat, lon):
        if not self._zero_isoterm_data:
            self._zero_isoterm_data = load_data_interpolator(
                '839/v4_esalat.npz', '839/v4_esalon.npz',
                '839/v4_esa0height.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return self._zero_isoterm_data(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def rain_height(self, lat_d, lon_d):
        """The rain height is computed as

        ..math:
            h_r = h_0 + 0.36 (km)
        """
        return self.isoterm_0(lat_d, lon_d) + 0.36


class _ITU839_3_():

    def __init__(self):
        self.__version__ = 3
        self.year = 2001
        self.month = 2
        self.link = 'https://www.itu.int/rec/R-REC-P.839-3-200102-S/en'

        self._zero_isoterm_data = {}

    def isoterm_0(self, lat, lon):
        if not self._zero_isoterm_data:
            self._zero_isoterm_data = load_data_interpolator(
                '839/v3_esalat.npz', '839/v3_esalon.npz',
                '839/v3_esa0height.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return self._zero_isoterm_data(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def rain_height(self, lat_d, lon_d):
        """
        The rain height is computed as

        ..math:
            h_r = h_0 + 0.36 (km)
        """
        return self.isoterm_0(lat_d, lon_d) + 0.36


class _ITU839_2_():

    def __init__(self):
        self.__version__ = 2
        self.year = 1999
        self.month = 10
        self.link = 'https://www.itu.int/rec/R-REC-P.839-2-199910-S/en'

    @staticmethod
    def isoterm_0(lat_d, lon_d):
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


__model = __ITU839__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.839 recommendation currently being used.

    This function changes the model used for the ITU-R P.839 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          * 4: Activates recommendation ITU-R P.839-4 (09/2013) (Current version)
          * 3: Activates recommendation ITU-R P.839-3 (02/01) (Superseded)
          * 2: Activates recommendation ITU-R P.839-2 (10/99) (Superseded)


    """
    global __model
    __model = __ITU839__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.839 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__


def isoterm_0(lat, lon):
    """
    Estimate the zero degree Celsius isoterm height for propagation prediction.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points


    Returns
    -------
    h0: numpy.ndarray
        Zero degree Celsius isoterm height (km)


    References
    ----------
    [1] Rain height model for prediction methods:
    https://www.itu.int/rec/R-REC-P.839/en

    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.isoterm_0(lat, lon)
    return prepare_output_array(val, type_output) * u.km


def rain_height(lat, lon):
    """
    Estimate the annual mean rain height for propagation prediction.

    The mean annual rain height above mean sea level, :math:`h_R`,
    may be obtained from the 0° C isotherm as:

    .. math::

      h_R = h_0 + 0.36 \\qquad \\text{km}


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points


    Returns
    -------
    hR: numpy.ndarray
        Annual mean rain height (km)


    References
    ----------
    [1] Rain height model for prediction methods:
    https://www.itu.int/rec/R-REC-P.839/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.rain_height(lat, lon)
    return prepare_output_array(val, type_output) * u.km
