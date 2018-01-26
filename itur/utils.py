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

dir_path = os.path.dirname(os.path.realpath(__file__))
dataset_dir = os.path.join(dir_path, './data/')

# Create a memory cache to memoize results of some functions
cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)

__NUMERIC_TYPES__ = [numbers.Number, int, float, complex,
                     np.float, np.float16, np.float32, np.float64,
                     np.int, np.int8, np.int16, np.int32, np.int64]


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

    if type_input in __NUMERIC_TYPES__ and \
       ((isinstance(array, np.ndarray) and array.size == 1) or
        (not isinstance(array, np.ndarray) and len(array) == 1)):
            value = float(value)
    elif type_input is list:
        if isinstance(value, np.ndarray):
            value = value.to_list()
        else:
            value = list(value)
    else:
        value = value

    # Squeeze output array to remove singleton dimensions
    if isinstance(value, np.ndarray):
        value = value.squeeze()

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
    else:
        raise ValueError('%s has not the correct format. It must be a value,'
                         'sequence, array, or a Quantity with %s units' %
                         (name_val, str(units)))


def compute_distance_earth_to_earth(lat_p, lon_p, lat_grid, lon_grid):
    '''
    Compute the distance between a point (P) in (lat_s, lon_s) and a matrix of
    latitude and longitudes (lat_grid, lon_grid)


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
                lon_min=None, lon_max=None, cbar_text='',
                **kwargs):

    import matplotlib.pyplot as plt

    try:
        from mpl_toolkits.basemap import Basemap
    except:
        raise RuntimeError('Basemap is not installed and therefore plot_in_map'
                           ' cannot be used')

    if all([el is None for el in [lat, lon, lat_min, lon_min, lat_max, lon_max]]):
        raise ValueError('Either \{lat, lon\} or \{lat_min, lon_min, lat_max,'
                         'lon_max\} need to be provided')

    elif lat is not None and lon is not None:
        assert(np.shape(lat) == np.shape(lon) and
               np.shape(lat) == np.shape(data))
        lat_max = np.max(lat)
        lat_min = np.min(lat)
        lon_max = np.max(lon)
        lon_min = np.min(lon)

    ax = plt.subplot(111)
    m = Basemap(ax=ax, projection='cyl', llcrnrlat=lat_min,
                urcrnrlat=lat_max, llcrnrlon=lon_min, urcrnrlon=lon_max,
                resolution='l')

    m.drawcoastlines(color='white', linewidth=0.5)
    m.drawcountries(color='grey', linewidth=0.3)
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
