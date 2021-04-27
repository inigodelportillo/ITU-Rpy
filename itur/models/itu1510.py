# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np
from astropy import units as u

from itur.models.itu1144 import (bilinear_2D_interpolator,
                                 bicubic_2D_interpolator)
from itur.utils import (prepare_input_array, prepare_output_array,
                        load_data_interpolator, get_input_type)


class __ITU1510__():

    """Annual mean surface temperature

    Available versions include:
    * P.1510-0 (02/01) (Current version)
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.1510 recommendation.

    def __init__(self, version=1):
        if version == 1:
            self.instance = _ITU1510_1_()
        elif version == 0:
            self.instance = _ITU1510_0_()
        else:
            raise ValueError('Version ' + str(version) + ' is not implemented'
                             ' for the ITU-R P.1510 model.')

    @property
    def __version__(self):
        return self.instance.__version__

    def surface_mean_temperature(self, lat, lon):
        """
        Method to compute the annual mean surface temperature (K).

        The temperature is computed at 2 m above the surface of the Earth.
        """
        return self.instance.temperature(lat, lon)

    def surface_month_mean_temperature(self, lat, lon, m):
        # Abstract method to compute the monthly surface mean temperature
        fcn = np.vectorize(self.instance.surface_month_mean_temperature,
                           excluded=[0, 1], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, m).tolist())


class _ITU1510_1_():

    def __init__(self):
        self.__version__ = 1
        self.year = 2017
        self.month = 6
        self.link = 'https://www.itu.int/rec/R-REC-P.1510/' +\
                    'recommendation.asp?lang=en&parent=R-REC-P.1510-1-201706-I'

        self.__months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        self._temperature = {}
        self._month_temperature = {}

    def temperature(self, lat, lon):
        if not self._temperature:
            self._temperature = load_data_interpolator(
                '1510/v1_lat.npz', '1510/v1_lon.npz',
                '1510/v1_t_annual.npz', bilinear_2D_interpolator)

        lon[lon > 180] = lon[lon > 180] - 360
        return self._temperature(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def month_temperature(self, lat, lon, m):
        if not self._month_temperature:
            for _m in self.__months:
                self._month_temperature[_m] = load_data_interpolator(
                    '1510/v1_lat.npz', '1510/v1_lon.npz',
                    '1510/v1_t_month{0:02d}.npz'.format(_m),
                    bilinear_2D_interpolator)

        lon[lon > 180] = lon[lon > 180] - 360
        return self._month_temperature[m](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def surface_mean_temperature(self, lat, lon):
        """
        Method to compute the annual mean surface temperature (K) at 2 m
        above the surface of the Earth
        """
        return self.temperature(lat, lon)

    def surface_month_mean_temperature(self, lat, lon, m):
        return self.month_temperature(lat, lon, m)


class _ITU1510_0_():

    def __init__(self):
        self.__version__ = 0
        self.year = 2001
        self.month = 2
        self.link = 'https://www.itu.int/rec/R-REC-P.1510/' +\
                    'recommendation.asp?lang=en&parent=R-REC-P.1510-0-200102-I'

        self._temperature = {}

    def temperature(self, lat, lon):
        if not self._temperature:
            self._temperature = load_data_interpolator(
                '1510/v0_lat.npz', '1510/v0_lon.npz',
                '1510/v0_temp.npz', bicubic_2D_interpolator)

        return self._temperature(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def surface_month_mean_temperature(self, lat, lon, m):
        raise NotImplementedError('The monthly mean temperature is not'
                                  'implemented in recomendation ITU-R P.1510'
                                  '-{0}'.format(self.__version__))


__model = __ITU1510__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.1510 recommendation currently being used.


    This function changes the model used for the ITU-R P.1510 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:

        * 1: Activates recommendation ITU-R P.1510-1 (06/17) (Current version)
        * 0: Activates recommendation ITU-R P.1510-0 (02/01) (Current version)
    """
    global __model
    __model = __ITU1510__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.1510 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__


def surface_mean_temperature(lat, lon):
    """
    Annual mean surface temperature (K) at 2 m above the surface of the Earth.

    A method to estimate the annual mean surface temperature (K) at 2 m
    above the surface of the Earth


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points


    Returns
    -------
    annual_temperature: numpy.ndarray
        Annual mean surface temperature (K). Same dimensions as lat and lon.


    References
    ----------
    [1] Annual mean surface temperature:
    https://www.itu.int/rec/R-REC-P.1510/en

    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.surface_mean_temperature(lat, lon)
    return prepare_output_array(val, type_output) * u.Kelvin


def surface_month_mean_temperature(lat, lon, m):
    """
    Monthly mean surface temperature (K) at 2 m above the surface of the Earth.

    A method to estimate the monthly mean surface temperature (K) at 2 m
    above the surface of the Earth


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    m : integer
      Index of the month (1=Jan, 2=Feb, 3=Mar, 4=Apr, ...)


    Returns
    -------
    monthly_temperature: numpy.ndarray
        Monthly mean surface temperature (K). Same dimensions as lat and lon.


    References
    ----------
    [1] Annual mean surface temperature:
    https://www.itu.int/rec/R-REC-P.1510/en

    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.surface_month_mean_temperature(lat, lon, m)
    return prepare_output_array(val, type_output) * u.Kelvin
