# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np
from astropy import units as u

from itur.models.itu1144 import bicubic_2D_interpolator
from itur.utils import (prepare_input_array, prepare_output_array,
                        load_data_interpolator, get_input_type)


class __ITU1511():
    """Topography for Earth-to-space propagation modelling. This model shall be
    used to obtain the height above mean sea level when no local data are
    available or when no data with a better spatial resolution is available.

    Available versions include:
    * P.1511-0 (02/01) (Superseded)
    * P.1511-1 (07/15) (Superseded)
    * P.1511-2 (08/19) (Current version)
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.1511 recommendation.

    def __init__(self, version=2):
        if version == 2:
            self.instance = _ITU1511_2_()
        elif version == 1:
            self.instance = _ITU1511_1_()
        elif version == 0:
            self.instance = _ITU1511_0_()
        else:
            raise ValueError('Version ' + str(version) + ' is not implemented'
                             ' for the ITU-R P.1511 model.')

    @property
    def __version__(self):
        """
        Version of the model (similar to version of the ITU Recommendation)
        """
        return self.instance.__version__

    def topographic_altitude(self, lat, lon):
        # Abstract method to compute the topographic altitude
        return self.instance.topographic_altitude(lat, lon)


class _ITU1511_2_():
    """
    The values of topographical height (km) above mean sea level of the surface
    of the Earth are  provided on a 0.5° grid in both latitude and longitude.
    For a location different from the gridpoints, the height above mean sea
    level at the desired location can be obtained by performing a bi-cubic
    interpolation.
    """

    def __init__(self):
        self.__version__ = 2
        self.year = 2019
        self.month = 8
        self.link = 'https://www.itu.int/rec/R-REC-P.1511/' +\
                    'recommendation.asp?lang=en&parent=R-REC-P.1511-2-201908-I'

        self._altitude = None
        self._wgs4_altitude = None

    def altitude(self, lat, lon):
        if not self._altitude:
            self._altitude = load_data_interpolator(
                '1511/v2_lat.npz', '1511/v2_lon.npz',
                '1511/v2_topo.npz', bicubic_2D_interpolator)

        return self._altitude(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def wgs4_altitude(self, lat, lon):
        if not self._wgs4_altitude:
            self._wgs4_altitude = load_data_interpolator(
                '1511/v2_lat.npz', '1511/v2_lon.npz',
                '1511/v2_egm2008.npz', bicubic_2D_interpolator)

        return self._wgs4_altitude(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def topographic_altitude(self, lat_d, lon_d):
        """
        Method to compute the values of topographical height (km) above mean
        sea level of the surface of the Earth.
        """

        # The new recommendation provides the output in meters and uses
        # a -180 to 180 longitude refernce
        lon_d = np.where(lon_d > 180, lon_d - 360, lon_d)
        return self.altitude(lat_d, lon_d) / 1000


class _ITU1511_1_():
    """
    The values of topographical height (km) above mean sea level of the surface
    of the Earth are  provided on a 0.5° grid in both latitude and longitude.
    For a location different from the gridpoints, the height above mean sea
    level at the desired location can be obtained by performing a bi-cubic
    interpolation.
    """

    def __init__(self):
        self.__version__ = 1
        self.year = 2015
        self.month = 7
        self.link = 'https://www.itu.int/rec/R-REC-P.1511/' +\
                    'recommendation.asp?lang=en&parent=R-REC-P.1511-1-201507-I'

        self._altitude = None

    def altitude(self, lat, lon):
        if not self._altitude:
            self._altitude = load_data_interpolator(
                '1511/v1_lat.npz', '1511/v1_lon.npz',
                '1511/v1_topo_0dot5.npz', bicubic_2D_interpolator)

        return self._altitude(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def topographic_altitude(self, lat_d, lon_d):
        """
        Method to compute the values of topographical height (km) above mean
        sea level of the surface of the Earth.
        """
        return self.altitude(lat_d, lon_d)


class _ITU1511_0_():
    """
    The values of topographical height (km) above mean sea level of the surface
    of the Earth are  provided on a 0.5° grid in both latitude and longitude.
    For a location different from the gridpoints, the height above mean sea
    level at the desired location can be obtained by performing a bi-cubic
    interpolation.
    """

    def __init__(self):
        self.__version__ = 0
        self.year = 2001
        self.month = 2
        self.link = 'https://www.itu.int/rec/R-REC-P.1511/' +\
                    'recommendation.asp?lang=en&parent=R-REC-P.1511-0-200102-I'

        self._altitude = None

    def altitude(self, lat, lon):
        if not self._altitude:
            self._altitude = load_data_interpolator(
                '1511/v1_lat.npz', '1511/v1_lon.npz',
                '1511/v1_topo_0dot5.npz', bicubic_2D_interpolator)

        return self._altitude(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def topographic_altitude(self, lat_d, lon_d):
        """
        Method to compute the values of topographical height (km) above mean
        sea level of the surface of the Earth.
        """
        return self.altitude(lat_d, lon_d)


__model = __ITU1511()


def change_version(new_version):
    """
    Change the version of the ITU-R P.1511 recommendation currently being used.


    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          * 1:  Activates recommendation ITU-R P.1511-1 (07/15) (Current version)
          * 0:  Activates recommendation ITU-R P.1511-0 (02/01) (Superseded)

    """
    global __model
    __model = __ITU1511(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.1511 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__


def topographic_altitude(lat, lon):
    """
    Topographical height (km) above mean sea level of the surface of the Earth.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points


    Returns
    -------
    altitude: numpy.ndarray
        Topographic altitude (km)


    References
    ----------
    [1] Topography for Earth-to-space propagation modelling:
    https://www.itu.int/rec/R-REC-P.1511/en

    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.topographic_altitude(lat, lon)
    val = np.maximum(val, 1e-9)
    return prepare_output_array(val, type_output) * u.km
