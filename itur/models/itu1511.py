# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import os
from astropy import units as u

from itur import utils
from itur.models.itu1144 import bicubic_2D_interpolator
from itur.utils import load_data, dataset_dir, prepare_input_array,\
    prepare_output_array, memory


class __ITU1511():
    """Topography for Earth-to-space propagation modelling. This model shall be
    used to obtain the height above mean sea level when no local data are
    available or when no data with a better spatial resolution is available.

    Available versions include:
    * P.1511-0 (02/01) (Superseded)
    * P.1511-1 (07/15) (Current version)
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.1511 recommendation.

    def __init__(self, version=1):
        if version == 1:
            self.instance = _ITU1511_1()
        elif version == 0:
            self.instance = _ITU1511_0()
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


class _ITU1511_1():
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

        self._altitude = {}

    def altitude(self, lat, lon):

        if not self._altitude:
            vals = load_data(os.path.join(dataset_dir,
                                          '1511/v1_TOPO_0DOT5.txt'))
            lats = load_data(os.path.join(dataset_dir, '1511/v1_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '1511/v1_Lon.txt'))
            self._altitude = bicubic_2D_interpolator(np.flipud(lats), lons,
                                                     np.flipud(vals))

        return self._altitude(
                np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def topographic_altitude(self, lat_d, lon_d):
        """
        Method to compute the values of topographical height (km) above mean
        sea level of the surface of the Earth.
        """
        return self.altitude(lat_d, lon_d)


class _ITU1511_0():
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

        self._altitude = {}

    def altitude(self, lat, lon):
        if not self._altitude:
            vals = load_data(os.path.join(dataset_dir,
                                          '1511/v1_TOPO_0DOT5.txt'))
            lats = load_data(os.path.join(dataset_dir, '1511/v1_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '1511/v1_Lon.txt'))
            self._altitude = bicubic_2D_interpolator(np.flipud(lats), lons,
                                                     np.flipud(vals))

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
        *0: P.1511-0 (02/01) (Superseded)
        *1: P.1511-1 (07/15) (Current version)
    """
    global __model
    __model = __ITU1511(new_version)
    utils.memory.clear()


def get_version():
    """
    Obtain the version of the ITU-R P.1511 recommendation currently being used.
    """
    global __model
    return __model.__version__


@memory.cache
def topographic_altitude(lat, lon):
    """
    The values of topographical height (km) above mean sea level of the surface
    of the Earth a



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
    global __model
    type_output = type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.topographic_altitude(lat, lon)
    val = np.maximum(val, 1e-7)
    return prepare_output_array(val, type_output) * u.km
