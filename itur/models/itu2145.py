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

def __interpolator_2145__(
    self, data, scale_height_fcn, lat, lon, p, alt=None, alt_res_fcn=topographic_altitude
):
    lat_f = lat.flatten()
    lon_f = lon.flatten()

    available_p = np.array([
        0.01, 0.02, 0.03, 0.05,
        0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10,
        20, 30, 50, 60, 70, 80, 90, 95, 99
    ])

    if p in available_p or p == 'mean':
        p_below = p_above = p
        pExact = True
    else:
        pExact = False
        idx = available_p.searchsorted(p, side='right') - 1
        idx = np.clip(idx, 0, len(available_p) - 1)

        p_below = available_p[idx]
        idx = np.clip(idx + 1, 0, len(available_p) - 1)
        p_above = available_p[idx]

    latlon_step = 0.25

    R = (lat_f+90) // latlon_step
    C = (lon_f+180) // latlon_step

    lats = np.array([
        R * latlon_step - 90,
        (R + 1) * latlon_step - 90,
        R * latlon_step - 90,
        (R + 1) * latlon_step - 90
    ])

    lons = np.mod(
        np.array([
            (C * latlon_step),
            (C * latlon_step),
            ((C + 1) * latlon_step),
            ((C + 1) * latlon_step)
        ]),
        360
    ) - 180

    r = (lat_f+90) / latlon_step
    c = (lon_f+180) / latlon_step

    data_a = data(lats, lons, p_above)

    # Compute the altitude of the data point
    if alt_res_fcn is topographic_altitude:
        altitude_res = alt_res_fcn(lats, lons).value.reshape(lats.shape)
    else:
        altitude_res = alt_res_fcn(lats, lons)

    if alt is None:
        alt = alt_res_fcn(lat,lon)
    else:
        alt = alt.flatten()

    data_a = scale_height_fcn(data_a, lats, lons, alt, altitude_res)

    data_a = (
        data_a[0, :] * ((R + 1 - r) * (C + 1 - c)) +
        data_a[1, :] * ((r - R) * (C + 1 - c)) +
        data_a[2, :] * ((R + 1 - r) * (c - C)) +
        data_a[3, :] * ((r - R) * (c - C))
    )

    if not pExact:
        data_b = data(lats, lons, p_below)
        data_b = scale_height_fcn(data_b, lats, lons, alt, altitude_res)

        data_b = (
            data_b[0, :] * ((R + 1 - r) * (C + 1 - c)) +
            data_b[1, :] * ((r - R) * (C + 1 - c)) +
            data_b[2, :] * ((R + 1 - r) * (c - C)) +
            data_b[3, :] * ((r - R) * (c - C))
        )

    # Compute the values of Lred_a
    if not pExact:
        result = data_b + (data_a - data_b) * (
            np.log(p) - np.log(p_below)
        ) / (
            np.log(p_above) - np.log(p_below)
        )
        return result.reshape(lat.shape)
    else:
        return data_a.reshape(lat.shape)
    

class __ITU2145():
    """Private class to model the ITU-R P.2145 recommendations

    Calculates Pressure, Temperature, Water Vapour density, and Water Vapour Content

    Available versons:
        * P.2145-0 (08/22) (Current version) 
    """

    def __init__(self, version=0):
        if version == 0:
            self.instance = _ITU2145_0()
        else:
            raise ValueError(
                f"Version {version} is not implemented for the ITU-R P.2145 model."
            )
        
        self._P = {}
        self._T = {}
        self._V = {}
        self._rho = {}

    @property
    def __version__(self):
        return self.instance.__version__
    
    def barometric_surface_pressure(self, lat, lon, p, alt):
        fcn = np.vectorize(
            self.instance.barometric_surface_pressure,
            excluded=[0, 1, 3],
            otypes=[np.ndarray]
        )
        return np.array(fcn(lat, lon, p, alt).tolist())
    
    def surface_temperature(self, lat, lon, p, alt):
        fcn = np.vectorize(
            self.instance.surface_temperature,
            excluded=[0,1,3],
            otypes=[np.ndarray]
        )
        return np.array(fcn(lat, lon, p, alt).tolist())
    
    def surface_water_vapour_density(self, lat, lon, p, alt):
        fcn = np.vectorize(
            self.instance.surface_water_vapour_density,
            excluded=[0, 1, 3],
            otypes=[np.ndarray]
        )
        return np.array(fcn(lat, lon, p, alt).tolist())

    def total_water_vapour_content(self, lat, lon, p, alt):
        fcn = np.vectorize(
            self.instance.total_water_vapour_content,
            excluded=[0, 1, 3],
            otypes=[np.ndarray]
        )
        return np.array(fcn(lat, lon, p, alt).tolist())


class _ITU2145_0():

    def __init__(self):
        self.__version__ = 0
        self.year = 2022
        self.month = 8
        self.link = 'https://www.itu.int/rec/R-REC-P.2145/en'

        self._P = {}
        self._T = {}
        self._V = {}
        self._rho = {}

        self._VSCH = None
        self._PSCH = None
        self._TSCH = None

        self._topo_alt = None

    def P(self, lat, lon, p='mean'):
        if not self._P:
            ps = [
                'mean', 0.01, 0.02, 0.03, 0.05,
                0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                50, 60, 70, 80, 90, 95, 99
            ]
            d_dir = os.path.join(
                dataset_dir,
                '2145/v0_P_%s.npz'
            )
            for p_loads in ps:
                self._P[p_loads] = load_data_interpolator(
                    '2145/v0_lat.npz',
                    '2145/v0_lon.npz',
                    d_dir % (str(p_loads).replace('.','')),
                    bilinear_2D_interpolator,
                    flip_ud=False
                )
            
        return self._P[p](
            np.array([lat.ravel(), lon.ravel()]).T
        ).reshape(lat.shape)

    def T(self, lat, lon, p='mean'):
        if not self._T:
            ps = [
                'mean', 0.01, 0.02, 0.03, 0.05,
                0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                50, 60, 70, 80, 90, 95, 99
            ]
            d_dir = os.path.join(
                dataset_dir,
                '2145/v0_T_%s.npz'
            )
            for p_loads in ps:
                self._T[p_loads] = load_data_interpolator(
                    '2145/v0_lat.npz',
                    '2145/v0_lon.npz',
                    d_dir % (str(p_loads).replace('.','')),
                    bilinear_2D_interpolator,
                    flip_ud=False
                )
            
        return self._T[p](
            np.array([lat.ravel(), lon.ravel()]).T
        ).reshape(lat.shape)

    def rho(self, lat, lon, p='mean'):
        if not self._rho:
            ps = [
                'mean', 0.01, 0.02, 0.03, 0.05,
                0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                50, 60, 70, 80, 90, 95, 99
            ]
            d_dir = os.path.join(
                dataset_dir,
                '2145/v0_RHO_%s.npz'
            )
            for p_loads in ps:
                self._rho[p_loads] = load_data_interpolator(
                    '2145/v0_lat.npz',
                    '2145/v0_lon.npz',
                    d_dir % (str(p_loads).replace('.','')),
                    bilinear_2D_interpolator,
                    flip_ud=False
                )
            
        return self._rho[p](
            np.array([lat.ravel(), lon.ravel()]).T
        ).reshape(lat.shape)

    def V(self, lat, lon, p='mean'):
        if not self._V:
            ps = [
                'mean', 0.01, 0.02, 0.03, 0.05,
                0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                50, 60, 70, 80, 90, 95, 99
            ]
            d_dir = os.path.join(
                dataset_dir,
                '2145/v0_V_%s.npz'
            )
            for p_loads in ps:
                self._V[p_loads] = load_data_interpolator(
                    '2145/v0_lat.npz',
                    '2145/v0_lon.npz',
                    d_dir % (str(p_loads).replace('.','')),
                    bilinear_2D_interpolator,
                    flip_ud=False
                )
            
        return self._V[p](
            np.array([lat.ravel(), lon.ravel()]).T
        ).reshape(lat.shape)
        
    def VSCH_scale(self, data, lats, lons, alt, alt_res):
        if self._VSCH is None:
            d_dir = os.path.join(dataset_dir, '2145/v0_VSCH.npz')
            self._VSCH = load_data_interpolator(
                '2145/v0_lat.npz', '2145/v0_lon.npz',
                d_dir, bilinear_2D_interpolator, flip_ud = False
            )

        VSCHs = self._VSCH(
            np.array([lats.ravel(), lons.ravel()]).T
            ).reshape(lats.shape)
        
        return data * np.exp(
            -(alt - alt_res) * 1.0 / (VSCHs)
        )

    def PSCH_scale(self, data, lats, lons, alt, alt_res):
        if self._PSCH is None:
            d_dir = os.path.join(dataset_dir, '2145/v0_PSCH.npz')
            self._PSCH = load_data_interpolator(
                '2145/v0_lat.npz', '2145/v0_lon.npz',
                d_dir, bilinear_2D_interpolator, flip_ud = False
            )

        PSCHs = self._PSCH(
            np.array([lats.ravel(), lons.ravel()]).T
            ).reshape(lats.shape)
        
        return data * np.exp(
            -(alt - alt_res) * 1.0 / (PSCHs)
        )

    def TSCH_scale(self, data, lats, lons, alt, alt_res):
        if self._TSCH is None:
            d_dir = os.path.join(dataset_dir, '2145/v0_TSCH.npz')
            self._TSCH = load_data_interpolator(
                '2145/v0_lat.npz', '2145/v0_lon.npz',
                d_dir, bilinear_2D_interpolator, flip_ud = False
            )

        TSCHs = self._TSCH(
            np.array([lats.ravel(), lons.ravel()]).T
            ).reshape(lats.shape)
        
        return data + TSCHs * (alt-alt_res)

    def topo_alt(self, lat, lon):
        if self._topo_alt is None:
            self._topo_alt = load_data_interpolator(
                '2145/v0_Z_ground_lat.npz', '2145/v0_Z_ground_lon.npz',
                '2145/v0_Z_ground.npz', bicubic_2D_interpolator, flip_ud = False
            )
        
        return self._topo_alt(
            np.array([lat.ravel(),lon.ravel()]).T
        ).reshape(lat.shape)
    
    def barometric_surface_pressure(self, lat, lon, p, alt):
        if p is None:
            p = 'mean'
        return __interpolator_2145__(
            self, data=self.P, scale_height_fcn=self.PSCH_scale, lat=lat, lon=lon, 
            p=p, alt=alt, alt_res_fcn=self.topo_alt
        )

    def surface_temperature(self, lat, lon, p, alt):
        if p is None:
            p = 'mean'
        return __interpolator_2145__(
            self, data=self.T, scale_height_fcn=self.TSCH_scale, lat=lat, lon=lon, 
            p=p, alt=alt, alt_res_fcn=self.topo_alt
        )

    def surface_water_vapour_density(self, lat, lon, p, alt):
        if p is None:
            p = 'mean'
        return __interpolator_2145__(
            self, data=self.rho, scale_height_fcn=self.VSCH_scale, lat=lat, lon=lon, 
            p=p, alt=alt, alt_res_fcn=self.topo_alt
        )

    def total_water_vapour_content(self, lat, lon, p, alt):
        if p is None:
            p = 'mean'
        return __interpolator_2145__(
            self, data=self.V, scale_height_fcn=self.VSCH_scale, lat=lat, lon=lon,
            p=p, alt=alt, alt_res_fcn=self.topo_alt
        )
    

global __model
__model = __ITU2145()


def change_version(new_version):
    """
    Change the version of the ITU-R P.2145 recommendation currently being used.

    This function changes the model used for the ITU-R P.2145 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          *  0: Activates recommendation ITU-R P.2145-0 (08/22) (Current version)

    """
    global __model
    __model = __ITU2145(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.2145 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__

def barometric_surface_pressure(lat, lon, p, alt=None):
    """
    Compute the barometric surface pressure along a path.

    This method computes the barometric surface pressure density along a path at a
    desired location on the surface of the Earth.

    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    p : number
        Percentage of time exceeded for p% of the average year. if None or 'mean',
        find the mean density based on supplied coordinates
    alt : number, sequence, or numpy.ndarray
        Altitude of the receivers. If None, use the topographical altitude as
        described in recommendation ITU-R P.1511


    Returns
    -------
    P: Quantity
       Barometric Surface Pressure (hPa)


    References
    ----------
    [1] Digital maps related to the calculation of gaseous attenuation and related effects  
    https://www.itu.int/rec/R-REC-P.2145/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon+180, 360) - 180
    alt = prepare_input_array(alt)
    alt = prepare_quantity(alt, u.km, 'Altitude of the receivers')
    val = __model.barometric_surface_pressure(lat, lon, p, alt)
    return prepare_output_array(val, type_output) * u.hPa


def surface_temperature(lat, lon, p, alt=None):
    """
    Computes the surface temperature along a path.

    This method computes the surface temperature along a path at a
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
    T: Quantity
       Annual Surface Temperature at a designated location (K)


    References
    ----------
    [1] Digital maps related to the calculation of gaseous attenuation and related effects  
    https://www.itu.int/rec/R-REC-P.2145/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon+180, 360) - 180
    alt = prepare_input_array(alt)
    alt = prepare_quantity(alt, u.km, 'Altitude of the receivers')
    val = __model.surface_temperature(lat, lon, p, alt)
    return prepare_output_array(val, type_output) * u.Kelvin

    


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
        Percentage of time exceeded for p% of the average year. if None or 'mean',
        find the mean density based on supplied coordinates
    alt : number, sequence, or numpy.ndarray
        Altitude of the receivers. If None, use the topographical altitude as
        described in recommendation ITU-R P.1511


    Returns
    -------
    rho: Quantity
       Surface water vapour density (g/m3)


    References
    ----------
    [1] Digital maps related to the calculation of gaseous attenuation and related effects  
    https://www.itu.int/rec/R-REC-P.2145/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon+180, 360) - 180
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
    [1] Digital maps related to the calculation of gaseous attenuation and related effects  
    https://www.itu.int/rec/R-REC-P.2145/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon+180, 360) - 180
    alt = prepare_input_array(alt)
    alt = prepare_quantity(alt, u.km, 'Altitude of the receivers')
    val = __model.total_water_vapour_content(lat, lon, p, alt)
    return prepare_output_array(val, type_output) * u.kg / u.m**2

    