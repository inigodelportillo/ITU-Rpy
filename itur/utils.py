# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import os
import numbers

from tempfile import mkdtemp
from joblib import Memory
from astropy import units as u

from pyproj import Geod

dir_path = os.path.dirname(os.path.realpath(__file__))
dataset_dir = os.path.join(dir_path, 'data/')

# Create a memory cache to memoize results of some functions
cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)

__NUMERIC_TYPES__ = [numbers.Number, int, float, complex,
                     np.float, np.float16, np.float32, np.float64,
                     np.int, np.int8, np.int16, np.int32, np.int64]
__wgs84_geod__ = Geod(ellps='WGS84')


def load_data(path, is_text=False, **kwargs):
    """ Loads data files from /itur/data/


    Parameters
    ----------
    path : string
        Path of the data to load
    is_text : bool
        Indicates whether the data is numerical or text


    Returns
    -------
    data: numpy.ndarray
        Numpy-array with the data. Numerical data is returned as a float
    """
    if is_text:
        data = np.loadtxt(path, dtype=np.string_, delimiter=',', **kwargs)
    else:
        data = np.genfromtxt(path, dtype=float, delimiter=',', **kwargs)
    return data


def prepare_input_array(array):
    """ Formats an array to be a 2-D numpy-array
    """
    if array is None:
        return None

    return np.atleast_2d(array)


def prepare_output_array(array, type_input=None):
    """ Formats the output to have the same shape and type as the input
    """
    global output_quantity

    if isinstance(array, u.Quantity):
        value = array.value
        unit = array.unit
    else:
        value = array
        unit = None

    # Squeeze output array to remove singleton dimensions
    if type(value) in [np.ndarray, list]:
        value = np.array(value).squeeze()

    if (type_input in __NUMERIC_TYPES__ and
        (type(array) in __NUMERIC_TYPES__) or
        ((isinstance(array, np.ndarray) and array.size == 1) or
         (not type(array) not in __NUMERIC_TYPES__ and len(array) == 1))):
        value = float(value)
    elif type_input is list:
        if isinstance(value, np.ndarray):
            value = value.tolist()
        else:
            value = list(value)
    else:
        value = value

    if unit is not None:
        return value * unit
    else:
        return value


def prepare_quantity(value, units=None, name_val=None):
    """ The function verifys that a
    """
    if value is None:
        return None

    if isinstance(value, u.Quantity):
        if units in [u.K, u.deg_C, u.Kelvin, u.Celsius, u.imperial.deg_F]:
            return value.to(units, equivalencies=u.temperature()).value
        else:
            return value.to(units).value

    elif isinstance(value, numbers.Number) and units is not None:
        return value
    elif isinstance(value, np.ndarray) and units is not None:
        return value
    elif isinstance(value, list) and units is not None:
        return np.array([prepare_quantity(v, units, name_val) for v in value])
    elif isinstance(value, tuple) and units is not None:
        return np.array([prepare_quantity(v, units, name_val) for v in value])
    else:
        raise ValueError('%s has not the correct format. It must be a value,'
                         'sequence, array, or a Quantity with %s units' %
                         (name_val, str(units)))


def compute_distance_earth_to_earth(lat_p, lon_p, lat_grid, lon_grid,
                                    method=None):
    '''
    Compute the distance between a point (P) in (lat_s, lon_s) and a matrix of
    latitude and longitudes (lat_grid, lon_grid).

    If the number of elements in lat_grid is smaller than 100,000, uses the
    WGS84 method, otherwise, uses the harvesine formula.


    Parameters
    ----------
    lat_p : number
        latitude projection of the point P (degrees)
    lon_p : number
        longitude projection of the point P (degrees)
    lat_grid : number, sequence of np.ndarray
        Grid of latitude points to which compute the distance (degrees)
    lon_grid : number, sequence of np.ndarray
        Grid of longitude points to which compute the distance (degrees)


    Returns
    -------
    d : numpy.ndarray
        Distance between the point P and each point in (lat_grid, lon_grid)
        (km)

    '''
    if ((method == 'WGS84' and not(method is not None)) or
        (type(lat_p) in __NUMERIC_TYPES__) or
        (type(lat_grid) in __NUMERIC_TYPES__) or
        (len(lat_grid) < 10000) or
        (isinstance(lat_grid, np.ndarray) and lat_grid.size < 1e5)):
            return compute_distance_earth_to_earth_wgs84(
                    lat_p, lon_p, lat_grid, lon_grid)
    else:
            return compute_distance_earth_to_earth_haversine(
                    lat_p, lon_p, lat_grid, lon_grid)


def compute_distance_earth_to_earth_wgs84(lat_p, lon_p, lat_grid, lon_grid):
    '''
    Compute the distance between a point (P) in (lat_s, lon_s) and a matrix of
    latitude and longitudes (lat_grid, lon_grid) using the WGS84 inverse method


    Parameters
    ----------
    lat_p : number
        latitude projection of the point P (degrees)
    lon_p : number
        longitude projection of the point P (degrees)
    lat_grid : number, sequence of np.ndarray
        Grid of latitude points to which compute the distance (degrees)
    lon_grid : number, sequence of np.ndarray
        Grid of longitude points to which compute the distance (degrees)


    Returns
    -------
    d : numpy.ndarray
        Distance between the point P and each point in (lat_grid, lon_grid)
        (km)

    '''
    lat_p = lat_p * np.ones_like(lat_grid)
    lon_p = lon_p * np.ones_like(lon_grid)
    _a, _b, d = __wgs84_geod__.inv(lon_p, lat_p, lon_grid, lat_grid)
    return d/1e3


def compute_distance_earth_to_earth_haversine(lat_p, lon_p,
                                              lat_grid, lon_grid):
    '''
    Compute the distance between a point (P) in (lat_s, lon_s) and a matrix of
    latitude and longitudes (lat_grid, lon_grid) using the Haversine formula


    Parameters
    ----------
    lat_p : number
        latitude projection of the point P (degrees)
    lon_p : number
        longitude projection of the point P (degrees)
    lat_grid : number, sequence of np.ndarray
        Grid of latitude points to which compute the distance (degrees)
    lon_grid : number, sequence of np.ndarray
        Grid of longitude points to which compute the distance (degrees)


    Returns
    -------
    d : numpy.ndarray
        Distance between the point P and each point in (lat_grid, lon_grid)
        (km)


    References
    This is based on the Haversine formula
    '''
    RE = 6371.0  # Radius of the Earth, km

    lat1 = np.deg2rad(lat_grid)
    lat2 = np.deg2rad(lat_p)
    lon1 = np.deg2rad(lon_grid)
    lon2 = np.deg2rad(lon_p)

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Compute the distance
    a = np.clip((np.sin(dlat / 2.0))**2 + np.cos(lat1) * np.cos(lat2) *
                (np.sin(dlon / 2))**2, -1, 1)
    c = 2 * np.arcsin(np.sqrt(a))
    d = RE * c
    return d


def regular_lat_lon_grid(resolution_lat=1, resolution_lon=1, lon_start_0=False,
                         lat_min=-90, lat_max=90, lon_min=-180, lon_max=180):
    '''
    Build latitude and longitude coordinate matrix with resolution
    resolution_lat, resolution_lon


    Parameters
    ----------
    resolution_lat: number
        Resolution for the latitude axis (deg)
    resolution_lon: number
        Resolution for the longitude axis (deg)
    lon_start_0: boolean
        Indicates whether the longitude is indexed using a 0 - 360 scale (True)
        or using -180 - 180 scale (False). Default value is False


    Returns
    -------
    lat: numpy.ndarray
        Grid of coordinates of the latitude point
    lon: numpy.ndarray
        Grid of coordinates of the latitude point
    '''
    if lon_start_0:
        lon, lat = np.meshgrid(np.arange(lon_min + 180.0, lon_max + 180.0,
                                         resolution_lon),
                               np.arange(lat_max, lat_min, - resolution_lat))
    else:
        lon, lat = np.meshgrid(np.arange(lon_min, lon_max, resolution_lon),
                               np.arange(lat_max, lat_min, - resolution_lat))

    return lat, lon


def elevation_angle(h, lat_s, lon_s, lat_grid, lon_grid):
    '''
    Compute the elevation angle between a satellite located in an orbit
    at height h and located above coordinates (lat_s, lon_s) and a matrix of
    latitude and longitudes (lat_grid, lon_grid)


    Parameters
    ----------
    h : float
        Orbital altitude of the satellite (km)
    lat_s : float
        latitude of the projection of the satellite (degrees)
    lon_s : float
        longitude of the projection of the satellite (degrees)
    lat_grid :  number, sequence of np.ndarray
        Grid of latitude points to which compute the elevation angle (degrees)
    lon_grid :  number, sequence of np.ndarray
        Grid of longitude points to which compute the elevation angle (degrees)


    Returns
    -------
    elevation : numpy.ndarray
        Elevation angle between the satellite and each point in
        (lat_grid, lon_grid) (degrees)


    References
    [1] http://www.propagation.gatech.edu/ECE6390/notes/ASD5.pdf - Slides 3, 4
    '''
    h = prepare_quantity(h, u.km, name_val='Orbital altitude of the satellite')

    RE = 6371.0     # Radius of the Earth (km)
    rs = RE + h

    # Transform latitude_longitude values to radians
    lat1 = np.deg2rad(lat_grid)
    lat2 = np.deg2rad(lat_s)
    lon1 = np.deg2rad(lon_grid)
    lon2 = np.deg2rad(lon_s)

    # Compute the elevation angle as described in
    gamma = np.arccos(
        np.clip(np.sin(lat2) * np.sin(lat1) +
                np.cos(lat1) * np.cos(lat2) * np.cos(lon2 - lon1), -1, 1))
    elevation = np.arccos(np.sin(gamma) /
                          np.sqrt(1 + (RE / rs)**2 -
                                  2 * (RE / rs) * np.cos(gamma)))  # In radians

    return np.rad2deg(elevation)


def plot_in_map(data, lat=None, lon=None, lat_min=None, lat_max=None,
                lon_min=None, lon_max=None, cbar_text='', ax=None,
                **kwargs):
    '''
    Displays the value sin data in a map.

    Either {lat, lon} or {lat_min, lat_max, lon_min, lon_max} need to be
    provided as inputs.

    Parameters
    ----------
    data : np.ndarray
        Data values to be plotted.
    lat : np.ndarray
        Matrix with the latitudes for each point in data (deg N)
    lon : np.ndarray
        Matrix with the longitudes for each point in data (deg E)
    lat_min :  float
        Minimum latitude of the data (deg N)
    lat_max :  float
        Maximum latitude of the data (deg N)
    lon_min :  float
        Minimum longitude of the data (deg E)
    lat_max :  float
        Maximum longitude of the data (deg E)
    cbar_text : string
        Colorbar text caption.
    ax : Axes
        matplotlib axes where the data will be plotted.
    **kwargs: dict
        Key-value arguments that will be passed to the imshow function.


    Returns
    -------
    m : Basemap
        The map object generated by Basemap
    '''
    import matplotlib.pyplot as plt

    try:
        from mpl_toolkits.basemap import Basemap
    except BaseException:
        raise RuntimeError('Basemap is not installed and therefore plot_in_map'
                           ' cannot be used')

    if all([el is None for el in [lat, lon, lat_min, lon_min,
                                  lat_max, lon_max]]):
        raise ValueError('Either \{lat, lon\} or \{lat_min, lon_min, lat_max,'
                         'lon_max\} need to be provided')

    elif lat is not None and lon is not None:
        assert(np.shape(lat) == np.shape(lon) and
               np.shape(lat) == np.shape(data))
        lat_max = np.max(lat)
        lat_min = np.min(lat)
        lon_max = np.max(lon)
        lon_min = np.min(lon)

    if ax is None:
        ax = plt.subplot(111)

    m = Basemap(ax=ax, projection='cyl', llcrnrlat=lat_min,
                urcrnrlat=lat_max, llcrnrlon=lon_min, urcrnrlon=lon_max,
                resolution='l')

    m.drawcoastlines(color='grey', linewidth=0.8)
    m.drawcountries(color='grey', linewidth=0.8)
    parallels = np.arange(-80, 81, 20)
    m.drawparallels(parallels, labels=[1, 0, 0, 1], dashes=[2, 1],
                    linewidth=0.2, color='white')
    meridians = np.arange(0., 360., 30.)
    m.drawmeridians(meridians, labels=[1, 0, 0, 1], dashes=[2, 1],
                    linewidth=0.2, color='white')

    im = m.imshow(np.flipud(data), **kwargs)
    cbar = m.colorbar(im, location='bottom', pad="8%")
    cbar.set_label(cbar_text)
    return m
