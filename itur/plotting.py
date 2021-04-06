# -*- coding: utf-8 -*-
"""``itur.utils`` is a utilities library for ITU-Rpy.

This utility library for ITU-Rpy contains methods to:
    * Load data and build an interpolator object.
    * Prepare the input and output arrays, and handle unit transformations.
    * Compute distances and elevation angles between two points on Earth and
      or space.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cpf
    import cartopy.mpl.gridliner as cmpl
    plotting_installed = True
except BaseException:
    plotting_installed = False


def plot_in_map(data, lat=None, lon=None, lat_min=None, lat_max=None,
                lon_min=None, lon_max=None, cbar_text='', ax=None,
                figsize=(6, 4), **kwargs):
    """Plot the values in `data` in a map.

    The map uses an equidistant cylindrical projection. Either
    {`lat`, `lon`} or {`lat_min`, `lat_max`, `lon_min`, `lon_max`} need
    to be provided as inputs. This function requires that `basemap` and
    `matplotlib` are installed.

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
    figsize : tuple
        Dimensions of the Figure
    **kwargs: dict
        Key-value arguments that will be passed to the contourf function.

    Returns
    -------
    ax : Axes
        The matploltib axes object
    """
    import matplotlib.pyplot as plt

    if not plotting_installed:
        raise RuntimeError('Neither cartopy nor matplotlib are installed. '
                           'Therefore plot_in_map cannot be used. '
                           'To use this function you need to install '
                           'the cartopy and matplotlib libraries')

    if all([el is None for el in [lat, lon, lat_min, lon_min,
                                  lat_max, lon_max]]):
        raise ValueError('Either {{lat, lon}} or {{lat_min, lon_min, lat_max,'
                         'lon_max}} need to be provided')

    elif lat is not None and lon is not None:
        if not(np.shape(lat) == np.shape(lon) and
               np.shape(lat) == np.shape(data)):
            raise RuntimeError('Shape of latitude grid is not equal to shape'
                               'of longitude grid')
        lat_max = np.max(lat)
        lat_min = np.min(lat)
        lon_max = np.max(lon)
        lon_min = np.min(lon)

    if ax is None:
        fig = plt.figure(figsize=figsize)
        proj = ccrs.PlateCarree(central_longitude=0.0)
        ax = fig.add_subplot(111, projection=proj)

    ax.set_extent([lon_min, lon_max, lat_min, lat_max],
                  crs=ccrs.PlateCarree())
    ax.coastlines(color='grey', linewidth=0.8)
    ax.add_feature(cpf.BORDERS, edgecolor='grey', linestyle=':')

    parallels = np.arange(-80, 81, 20)
    meridians = np.arange(-180., 181., 30.)

    ax.gridlines(xlocs=meridians, ylocs=parallels, draw_labels=True,
                 xformatter=cmpl.LONGITUDE_FORMATTER,
                 # yformatter=cmpl.LATITUDE_FORMATTER,
                 color='white', linestyle=':', linewidth=0.2)

    im = ax.contourf(lon, lat, data, 100, transform=ccrs.PlateCarree(),
                     **kwargs)

    cbar = fig.colorbar(im, orientation='horizontal', fraction=0.046, pad=0.04)
    cbar.set_label(cbar_text)
    fig.show()
    return ax


def plot_in_map_basemap(data, lat=None, lon=None, lat_min=None,
                        lat_max=None, lon_min=None, lon_max=None,
                        cbar_text='', ax=None, figsize=(6, 4), **kwargs):
    """Plot the values in `data` in a map.

    The map uses an equidistant cylindrical projection. Either
    {`lat`, `lon`} or {`lat_min`, `lat_max`, `lon_min`, `lon_max`} need
    to be provided as inputs. This function requires that `basemap` and
    `matplotlib` are installed.

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
    figsize : tuple
        Dimensions of the Figure
    **kwargs: dict
        Key-value arguments that will be passed to the imshow function.

    Returns
    -------
    m : Basemap
        The map object generated by Basemap
    """
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
    except BaseException:
        raise RuntimeError('Basemap is not installed and therefore '
                           'plot_in_map_basemap cannot be used. To use this '
                           'function you need to install the basemap library')

    if all([el is None for el in [lat, lon, lat_min, lon_min,
                                  lat_max, lon_max]]):
        raise ValueError('Either {{lat, lon}} or {{lat_min, lon_min, lat_max,'
                         'lon_max}} need to be provided')

    elif lat is not None and lon is not None:
        if not(np.shape(lat) == np.shape(lon) and
               np.shape(lat) == np.shape(data)):
            raise RuntimeError('Shape of latitude grid is not equal to shape'
                               'of longitude grid')
        lat_max = np.max(lat)
        lat_min = np.min(lat)
        lon_max = np.max(lon)
        lon_min = np.min(lon)

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

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
