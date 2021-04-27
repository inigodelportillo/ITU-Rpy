# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import numpy as np


from astropy import units as u

from itur.models.itu1511 import topographic_altitude
from itur.models.itu1144 import (bilinear_2D_interpolator,
                                 bicubic_2D_interpolator)
from itur.utils import (prepare_input_array, prepare_output_array,
                        dataset_dir, prepare_quantity, get_input_type,
                        load_data_interpolator)


def __interpolator_836__(self, data, lat, lon, p, alt=None,
                         alt_res_fcn=topographic_altitude):
    lat_f = lat.flatten()
    lon_f = lon.flatten()

    available_p = np.array([0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10,
                            20, 30, 50, 60, 70, 80, 90, 95, 99])

    if p in available_p:
        p_below = p_above = p
        pExact = True
    else:
        pExact = False
        idx = available_p.searchsorted(p, side='right') - 1
        idx = np.clip(idx, 0, len(available_p) - 1)

        p_below = available_p[idx]
        idx = np.clip(idx + 1, 0, len(available_p) - 1)
        p_above = available_p[idx]

    R = -(lat_f - 90) // 1.125
    C = lon_f // 1.125

    lats = np.array([90 - R * 1.125, 90 - (R + 1) * 1.125,
                     90 - R * 1.125, 90 - (R + 1) * 1.125])

    lons = np.mod(np.array([C * 1.125, C * 1.125,
                            (C + 1) * 1.125, (C + 1) * 1.125]), 360)

    r = - (lat_f - 90) / 1.125
    c = lon_f / 1.125

    data_a = data(lats, lons, p_above)
    VSCH_a = self.VSCH(lats, lons, p_above)

    # Compute the altitude of the data point
    if alt_res_fcn is topographic_altitude:
        altitude_res = alt_res_fcn(lats, lons).value.reshape(lats.shape)
    else:
        altitude_res = alt_res_fcn(lats, lons)

    if alt is None:
        alt = altitude_res
    else:
        alt = alt.flatten()

    data_a = data_a * np.exp(- (alt - altitude_res) * 1.0 / (VSCH_a))

    data_a = (data_a[0, :] * ((R + 1 - r) * (C + 1 - c)) +
              data_a[1, :] * ((r - R) * (C + 1 - c)) +
              data_a[2, :] * ((R + 1 - r) * (c - C)) +
              data_a[3, :] * ((r - R) * (c - C)))

    if not pExact:
        data_b = data(lats, lons, p_below)
        VSCH_b = self.VSCH(lats, lons, p_below)
        data_b = data_b * np.exp(- (alt - altitude_res) / (VSCH_b))

        data_b = (data_b[0, :] * ((R + 1 - r) * (C + 1 - c)) +
                  data_b[1, :] * ((r - R) * (C + 1 - c)) +
                  data_b[2, :] * ((R + 1 - r) * (c - C)) +
                  data_b[3, :] * ((r - R) * (c - C)))

    # Compute the values of Lred_a
    if not pExact:
        rho = data_b + (data_a - data_b) * (np.log(p) - np.log(p_below)) / \
            (np.log(p_above) - np.log(p_below))
        return rho.reshape(lat.shape)
    else:
        return data_a.reshape(lat.shape)


class __ITU836():
    """Private class to model the ITU-R P.836 recommendations.

    Water vapour: surface density and total columnar content

    Available versions:
       * P.836-6 (12/17) (Current version)
       * P.836-5 (09/13) (Superseded)
       * P.836-4 (10/09) (Superseded)

    Not available versions:
       * P.836-0 (03/92) (Superseded)
       * P.836-1 (08/97) (Superseded)
       * P.836-2 (02/01) (Superseded)
       * P.836-3 (11/01) (Superseded)
    """

    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.836 recommendation.

    def __init__(self, version=6):
        if version == 6:
            self.instance = _ITU836_6()
        elif version == 5:
            self.instance = _ITU836_5()
        elif version == 4:
            self.instance = _ITU836_4()
        else:
            raise ValueError(
                'Version {0} is not implemented for the ITU-R P.836 model.'
                .format(version))

        self._V = {}
        self._VSCH = {}
        self._rho = {}
        self._topo_alt = None

    @property
    def __version__(self):
        return self.instance.__version__

    def surface_water_vapour_density(self, lat, lon, p, alt):
        fcn = np.vectorize(self.instance.surface_water_vapour_density,
                           excluded=[0, 1, 3], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, p, alt).tolist())

    def total_water_vapour_content(self, lat, lon, p, alt):
        fcn = np.vectorize(self.instance.total_water_vapour_content,
                           excluded=[0, 1, 3], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, p, alt).tolist())


class _ITU836_6():

    def __init__(self):
        self.__version__ = 6
        self.year = 2017
        self.month = 12
        self.link = 'https://www.itu.int/rec/R-REC-P.836-6-201712-I/en'

        self._V = {}
        self._VSCH = {}
        self._rho = {}
        self._topo_alt = None

    def V(self, lat, lon, p):
        if not self._V:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '836/v6_v_%s.npz')
            for p_loads in ps:
                self._V[float(p_loads)] = load_data_interpolator(
                    '836/v6_lat.npz', '836/v6_lon.npz',
                    d_dir % (str(p_loads).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._V[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def VSCH(self, lat, lon, p):
        if not self._VSCH:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '836/v6_vsch_%s.npz')
            for p_loads in ps:
                self._VSCH[float(p_loads)] = load_data_interpolator(
                    '836/v6_lat.npz', '836/v6_lon.npz',
                    d_dir % (str(p_loads).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._VSCH[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def rho(self, lat, lon, p):
        if not self._rho:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '836/v6_rho_%s.npz')
            for p_loads in ps:
                self._rho[float(p_loads)] = load_data_interpolator(
                    '836/v6_lat.npz', '836/v6_lon.npz',
                    d_dir % (str(p_loads).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._rho[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def topo_alt(self, lat, lon):
        if self._topo_alt is None:
            self._topo_alt = load_data_interpolator(
                '836/v6_topolat.npz', '836/v6_topolon.npz',
                '836/v6_topo_0dot5.npz', bicubic_2D_interpolator)

        return self._topo_alt(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def surface_water_vapour_density(self, lat, lon, p, alt=None):
        return __interpolator_836__(
            self, data=self.rho, lat=lat, lon=lon, p=p, alt=alt,
            alt_res_fcn=self.topo_alt)

    def total_water_vapour_content(self, lat, lon, p, alt=None):
        return __interpolator_836__(
            self, data=self.V, lat=lat, lon=lon, p=p, alt=alt,
            alt_res_fcn=self.topo_alt)


class _ITU836_5():

    def __init__(self):
        self.__version__ = 5
        self.year = 2013
        self.month = 9
        self.link = 'https://www.itu.int/rec/R-REC-P.836-5-201309-I/en'

        self._V = {}
        self._VSCH = {}
        self._rho = {}

    def V(self, lat, lon, p):
        if not self._V:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '836/v5_v_%s.npz')
            for p_loads in ps:
                self._V[float(p_loads)] = load_data_interpolator(
                    '836/v5_lat.npz', '836/v5_lon.npz',
                    d_dir % (str(p_loads).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._V[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def VSCH(self, lat, lon, p):
        if not self._VSCH:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '836/v5_vsch_%s.npz')
            for p_loads in ps:
                self._VSCH[float(p_loads)] = load_data_interpolator(
                    '836/v5_lat.npz', '836/v5_lon.npz',
                    d_dir % (str(p_loads).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._VSCH[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def rho(self, lat, lon, p):
        if not self._rho:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '836/v5_rho_%s.npz')
            for p_loads in ps:
                self._rho[float(p_loads)] = load_data_interpolator(
                    '836/v5_lat.npz', '836/v5_lon.npz',
                    d_dir % (str(p_loads).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._rho[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def surface_water_vapour_density(self, lat, lon, p, alt=None):
        return __interpolator_836__(
            self, data=self.rho, lat=lat, lon=lon, p=p, alt=alt,
            alt_res_fcn=topographic_altitude)

    def total_water_vapour_content(self, lat, lon, p, alt=None):
        return __interpolator_836__(
            self, data=self.V, lat=lat, lon=lon, p=p, alt=alt,
            alt_res_fcn=topographic_altitude)


class _ITU836_4():

    def __init__(self):
        self.__version__ = 4
        self.year = 2009
        self.month = 10
        self.link = 'https://www.itu.int/rec/R-REC-P.836-4-200910-S/en'

        self._V = {}
        self._VSCH = {}
        self._rho = {}

    def V(self, lat, lon, p):
        if not self._V:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '836/v4_v_%s.npz')
            for p_loads in ps:
                self._V[float(p_loads)] = load_data_interpolator(
                    '836/v4_lat.npz', '836/v4_lon.npz',
                    d_dir % (str(p_loads).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._V[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def VSCH(self, lat, lon, p):
        if not self._VSCH:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '836/v4_vsch_%s.npz')
            for p_loads in ps:
                self._VSCH[float(p_loads)] = load_data_interpolator(
                    '836/v4_lat.npz', '836/v4_lon.npz',
                    d_dir % (str(p_loads).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._VSCH[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def rho(self, lat, lon, p):
        if not self._rho:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '836/v4_rho_%s.npz')
            for p_loads in ps:
                self._rho[float(p_loads)] = load_data_interpolator(
                    '836/v4_lat.npz', '836/v4_lon.npz',
                    d_dir % (str(p_loads).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._rho[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    # The procedure to compute the surface water vapour density and the
    # total water vapour content is similar to the ones in recommendation
    # ITU-P R.836-5.
    def surface_water_vapour_density(self, lat, lon, p, alt=None):
        return __interpolator_836__(
            self, data=self.rho, lat=lat, lon=lon, p=p, alt=alt,
            alt_res_fcn=topographic_altitude)

    def total_water_vapour_content(self, lat, lon, p, alt=None):
        return __interpolator_836__(
            self, data=self.V, lat=lat, lon=lon, p=p, alt=alt,
            alt_res_fcn=topographic_altitude)


__model = __ITU836()


def change_version(new_version):
    """
    Change the version of the ITU-R P.836 recommendation currently being used.

    This function changes the model used for the ITU-R P.836 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          *  6: Activates recommendation ITU-R P.836-6 (12/17) (Current version)
          *  5: Activates recommendation ITU-R P.836-5 (09/13) (Superseded)
          *  4: Activates recommendation ITU-R P.836-4 (10/09) (Superseded)

    """
    global __model
    __model = __ITU836(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.836 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__


def surface_water_vapour_density(lat, lon, p, alt=None):
    """
    Compute the surface water vapour density along a path.

    This method computes the surface water vapour density along a path at a
    desired location on the surface of the Earth.

    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    p : number
        Percentage of time exceeded for p% of the average year
    alt : number, sequence, or numpy.ndarray
        Altitude of the receivers. If None, use the topographical altitude as
        described in recommendation ITU-R P.1511


    Returns
    -------
    rho: Quantity
       Surface water vapour density (g/m3)


    References
    ----------
    [1] Water vapour: surface density and total columnar content
    https://www.itu.int/rec/R-REC-P.836/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    alt = prepare_input_array(alt)
    alt = prepare_quantity(alt, u.km, 'Altitude of the receivers')
    val = __model.surface_water_vapour_density(lat, lon, p, alt)
    return prepare_output_array(val, type_output) * u.g / u.m**3


def total_water_vapour_content(lat, lon, p, alt=None):
    """
    Compute the total water vapour content along a path.

    This method computes the total water vapour content along a path at a
    desired location on the surface of the Earth.

    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    p : number
        Percentage of time exceeded for p% of the average year
    alt : number, sequence, or numpy.ndarray
        Altitude of the receivers. If None, use the topographical altitude as
        described in recommendation ITU-R P.1511


    Returns
    -------
    V: Quantity
       Total water vapour content (kg/m2)


    References
    ----------
    [1] Water vapour: surface density and total columnar content
    https://www.itu.int/rec/R-REC-P.836/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    alt = prepare_input_array(alt)
    alt = prepare_quantity(alt, u.km, 'Altitude of the receivers')
    val = __model.total_water_vapour_content(lat, lon, p, alt)
    return prepare_output_array(val, type_output) * u.kg / u.m**2
