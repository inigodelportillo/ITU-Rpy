# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import numpy as np
from astropy import units as u

from itur.models.itu1144 import bilinear_2D_interpolator
from itur.utils import (prepare_input_array, prepare_quantity, load_data,
                        prepare_output_array, dataset_dir, get_input_type)


class __ITU453__():
    """ Private class to model the ITU-R P.453 recommendations.

    Implementation of the methods in Recommendation ITU-R P.453
    "The radio refractive index: its formula and refractivity data"

    Available versions:
       * P.453-13 (12/17)
       * P.453-12 (07/15)

    TODO: Implement version P.453-13

    Recommendation ITU-R P.453 provides methods to estimate the radio
    refractive index and its behaviour for locations worldwide; describes both
    surface and vertical profile characteristics; and provides global maps for
    the distribution of refractivity parameters and their statistical
    variation.
    """

    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.453 recommendation.

    def __init__(self, version=13):
        if version == 13:
            self.instance = _ITU453_13_()
        elif version == 12:
            self.instance = _ITU453_12_()
        else:
            raise ValueError(
                'Version {0} is not implemented for the ITU-R P.453 model.'
                .format(version))

    @property
    def __version__(self):
        return self.instance.__version__

    def wet_term_radio_refractivity(self, e, T):
        return self.instance.wet_term_radio_refractivity(e, T)

    def dry_term_radio_refractivity(self, Pd, T):
        return self.instance.dry_term_radio_refractivity(Pd, T)

    def radio_refractive_index(self, Pd, e, T):
        return self.instance.radio_refractive_index(Pd, e, T)

    def water_vapour_pressure(self, T, P, H, type_hydrometeor='water'):
        return self.instance.water_vapour_pressure(
            T, P, H, type_hydrometeor=type_hydrometeor)

    def saturation_vapour_pressure(self, T, P, type_hydrometeor='water'):
        return self.instance.saturation_vapour_pressure(
            T, P, type_hydrometeor=type_hydrometeor)

    def map_wet_term_radio_refractivity(self, lat, lon, p=50):
        fcn = np.vectorize(self.instance.map_wet_term_radio_refractivity,
                           excluded=[0, 1], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, p).tolist())

    def DN65(self, lat, lon, p):
        fcn = np.vectorize(self.instance.DN65, excluded=[0, 1],
                           otypes=[np.ndarray])
        return np.array(fcn(lat, lon, p).tolist())

    def DN1(self, lat, lon, p):
        fcn = np.vectorize(self.instance.DN1, excluded=[0, 1],
                           otypes=[np.ndarray])
        return np.array(fcn(lat, lon, p).tolist())


class _ITU453_13_():

    def __init__(self):
        self.__version__ = 13
        self.year = 2017
        self.month = 12
        self.link = 'https://www.itu.int/rec/R-REC-P.453-13-201712-I/en'

        self._N_wet = {}
        self._DN65 = {}
        self._DN1 = {}

    def DN65(self, lat, lon, p):
        if not self._DN65:
            ps = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80,
                  90, 95, 98, 99, 99.5, 99.8, 99.9]
            d_dir = os.path.join(dataset_dir, '453/v12_dn65m_%02dd%02d_v1.npz')
            lats = load_data(os.path.join(dataset_dir, '453/v12_lat0d75.npz'))
            lons = load_data(os.path.join(dataset_dir, '453/v12_lon0d75.npz'))
            for p_loads in ps:
                int_p = p_loads // 1
                frac_p = round((p_loads % 1.0) * 100)
                vals = load_data(d_dir % (int_p, frac_p))
                self._DN65[float(p_loads)] = bilinear_2D_interpolator(
                    lats, lons, vals)

        return self._DN65[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def DN1(self, lat, lon, p):
        if not self._DN1:
            ps = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80,
                  90, 95, 98, 99, 99.5, 99.8, 99.9]
            d_dir = os.path.join(dataset_dir, '453/v12_dn_%02dd%02d_v1.npz')
            lats = load_data(os.path.join(dataset_dir, '453/v12_lat0d75.npz'))
            lons = load_data(os.path.join(dataset_dir, '453/v12_lon0d75.npz'))
            for p_loads in ps:
                int_p = p_loads // 1
                frac_p = round((p_loads % 1.0) * 100)
                vals = load_data(d_dir % (int_p, frac_p))
                self._DN1[float(p_loads)] = bilinear_2D_interpolator(
                    lats, lons, vals)

        return self._DN1[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def N_wet(self, lat, lon, p):
        if not self._N_wet:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 60, 70, 80,
                  90, 95, 99]
            d_dir = os.path.join(dataset_dir, '453/v13_nwet_annual_%s.npz')
            lats = load_data(os.path.join(dataset_dir, '453/v13_lat_n.npz'))
            lons = load_data(os.path.join(dataset_dir, '453/v13_lon_n.npz'))
            for p_loads in ps:
                vals = load_data(d_dir % (str(p_loads).replace('.', '')))
                self._N_wet[float(p_loads)] = bilinear_2D_interpolator(
                    np.flipud(lats), lons, np.flipud(vals))

        lon[lon > 180] = lon[lon > 180] - 360
        return self._N_wet[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    @classmethod
    def wet_term_radio_refractivity(self, e, T):
        N_wet = (72 * e / (T + 273.15) + 3.75e5 * e / (T + 273.15)**2) * 1e-6
        return N_wet

    @classmethod
    def dry_term_radio_refractivity(self, Pd, T):
        N_dry = 77.6 * Pd / T  # Eq. 3
        return N_dry

    @classmethod
    def radio_refractive_index(self, Pd, e, T):
        N = 77.6 * Pd / T + 72 * e / T + 3.75e5 * e / T**2   # Eq. 2 [N-units]
        n = 1 + N * 1e-6   # Eq. 1
        return n

    @classmethod
    def water_vapour_pressure(self, T, P, H, type_hydrometeor='water'):
        e_s = self.saturation_vapour_pressure(T, P, type_hydrometeor)
        return H * e_s / 100   # Eq. 8

    @classmethod
    def saturation_vapour_pressure(self, T, P, type_hydrometeor='water'):

        if type_hydrometeor == 'water':
            EF = 1 + 1e-4 * (7.2 + P * (0.0320 + 5.9e-6 * T**2))
            a = 6.1121
            b = 18.678
            c = 257.14
            d = 234.5

        elif type_hydrometeor == 'ice':
            EF = 1 + 1e-4 * (2.2 + P * (0.0383 + 6.4e-6 * T**2))
            a = 6.1115
            b = 23.036
            c = 279.82
            d = 333.7

        e_s = EF * a * np.exp((b - T / d) * T / (T + c))
        return e_s

    def map_wet_term_radio_refractivity(self, lat, lon, p):
        # Fix lon because the data-set is now indexed -180 to 180 instead
        # of 0 to 360
        lon[lon > 180] = lon[lon > 180] - 360

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

        R = -(lat_f - 90) // 0.75
        C = (lon_f + 180) // 0.75

        lats = np.array([90 - R * 0.75, 90 - (R + 1) * 0.75,
                         90 - R * 0.75, 90 - (R + 1) * 0.75])

        lons = np.array([C * 0.75, C * 0.75,
                         (C + 1) * 0.75, (C + 1) * 0.75]) - 180

        r = - (lat_f - 90) / 0.75
        c = (lon_f + 180) / 0.75

        N_wet_a = self.N_wet(lats, lons, p_above)
        N_wet_a = (N_wet_a[0, :] * ((R + 1 - r) * (C + 1 - c)) +
                   N_wet_a[1, :] * ((r - R) * (C + 1 - c)) +
                   N_wet_a[2, :] * ((R + 1 - r) * (c - C)) +
                   N_wet_a[3, :] * ((r - R) * (c - C)))

        if not pExact:
            N_wet_b = self.N_wet(lats, lons, p_below)
            N_wet_b = (N_wet_b[0, :] * ((R + 1 - r) * (C + 1 - c)) +
                       N_wet_b[1, :] * ((r - R) * (C + 1 - c)) +
                       N_wet_b[2, :] * ((R + 1 - r) * (c - C)) +
                       N_wet_b[3, :] * ((r - R) * (c - C)))

        # Compute the values of Lred_a
        if not pExact:
            rho = N_wet_b + (N_wet_a - N_wet_b) * \
                (np.log(p) - np.log(p_below)) / \
                (np.log(p_above) - np.log(p_below))
            return rho.reshape(lat.shape)
        else:
            return N_wet_a.reshape(lat.shape)


class _ITU453_12_():

    def __init__(self):
        self.__version__ = 12
        self.year = 2016
        self.month = 9
        self.link = 'https://www.itu.int/rec/R-REC-P.453-12-201609-I/en'

        self._N_wet = {}
        self._DN65 = {}
        self._DN1 = {}

    def DN65(self, lat, lon, p):
        if not self._DN65:
            ps = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80,
                  90, 95, 98, 99, 99.5, 99.8, 99.9]
            d_dir = os.path.join(dataset_dir, '453/v12_dn65m_%02dd%02d_v1.npz')
            lats = load_data(os.path.join(dataset_dir, '453/v12_lat0d75.npz'))
            lons = load_data(os.path.join(dataset_dir, '453/v12_lon0d75.npz'))
            for p_loads in ps:
                int_p = p_loads // 1
                frac_p = round((p_loads % 1.0) * 100)
                vals = load_data(d_dir % (int_p, frac_p))
                self._DN65[float(p_loads)] = bilinear_2D_interpolator(
                    lats, lons, vals)

        return self._DN65[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def DN1(self, lat, lon, p):
        if not self._DN1:
            ps = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80,
                  90, 95, 98, 99, 99.5, 99.8, 99.9]
            d_dir = os.path.join(dataset_dir, '453/v12_dn_%02dd%02d_v1.npz')
            lats = load_data(os.path.join(dataset_dir, '453/v12_lat0d75.npz'))
            lons = load_data(os.path.join(dataset_dir, '453/v12_lon0d75.npz'))
            for p_loads in ps:
                int_p = p_loads // 1
                frac_p = round((p_loads % 1.0) * 100)
                vals = load_data(d_dir % (int_p, frac_p))
                self._DN1[float(p_loads)] = bilinear_2D_interpolator(
                    lats, lons, vals)

        return self._DN1[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def N_wet(self, lat, lon):
        if not self._N_wet:
            vals = load_data(os.path.join(dataset_dir, '453/v12_esanwet.npz'))
            lats = load_data(os.path.join(dataset_dir, '453/v12_esalat.npz'))
            lons = load_data(os.path.join(dataset_dir, '453/v12_esalon.npz'))
            self._N_wet = bilinear_2D_interpolator(lats, lons, vals)

        return self._N_wet(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    @staticmethod
    def wet_term_radio_refractivity(e, T):
        return _ITU453_13_.wet_term_radio_refractivity(e, T)

    @staticmethod
    def dry_term_radio_refractivity(Pd, T):
        return _ITU453_13_.dry_term_radio_refractivity(Pd, T)

    @staticmethod
    def radio_refractive_index(Pd, e, T):
        return _ITU453_13_.radio_refractive_index(Pd, e, T)

    @staticmethod
    def water_vapour_pressure(T, P, H, type_hydrometeor='water'):
        return _ITU453_13_.water_vapour_pressure(T, P, H, type_hydrometeor)

    @staticmethod
    def saturation_vapour_pressure(T, P, type_hydrometeor='water'):
        return _ITU453_13_.saturation_vapour_pressure(T, P, type_hydrometeor)

    def map_wet_term_radio_refractivity(self, lat, lon, p):
        return self.N_wet(lat, lon)


__model = __ITU453__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.453 recommendation currently being used.

    This function changes the model used for the ITU-R P.453 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          * 13: Activates recommendation ITU-R P.453-13 (12/17)
          * 12: Activates recommendation ITU-R P.453-12 (07/15)

    """
    global __model
    __model = __ITU453__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.453 recommendation currently being used.

    Returns
    -------
    version: int
       The version of the ITU-R P.453 recommendation being used.
    """
    return __model.__version__


def wet_term_radio_refractivity(e, T):
    """Determine the wet term of the radio refractivity.

    Parameters
    ----------
    e : number or Quantity
        Water vapour pressure  (hPa)
    T : number or Quantity
        Absolute temperature (K)


    Returns
    -------
    N_wet: Quantity
        Wet term of the radio refractivity (-)



    References
    ----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    e = prepare_quantity(e, u.hPa, 'Water vapour pressure ')
    T = prepare_quantity(T, u.K, 'Absolute temperature')
    val = __model.wet_term_radio_refractivity(e, T)
    return val * u.dimensionless_unscaled


def dry_term_radio_refractivity(Pd, T):
    """Determine the dry term of the radio refractivity.

    Parameters
    ----------
    Pd : number or Quantity
        Dry atmospheric pressure (hPa)
    T : number or Quantity
        Absolute temperature (K)


    Returns
    -------
    N_dry: Quantity
        Dry term of the radio refractivity (-)



    References
    ----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    Pd = prepare_quantity(Pd, u.hPa, 'Dry atmospheric pressure')
    T = prepare_quantity(T, u.K, 'Absolute temperature')
    val = __model.dry_term_radio_refractivity(Pd, T)
    return val * u.dimensionless_unscaled


def radio_refractive_index(Pd, e, T):
    """Compute the radio refractive index.

    Parameters
    ----------
    Pd : number or Quantity
        Dry atmospheric pressure (hPa)
    e : number or Quantity
        Water vapour pressure  (hPa)
    T : number or Quantity
        Absolute temperature (K)


    Returns
    -------
    n: Quantity
        Radio refractive index (-)



    References
    ----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    Pd = prepare_quantity(Pd, u.hPa, 'Dry atmospheric pressure')
    e = prepare_quantity(e, u.hPa, 'Water vapour pressure ')
    T = prepare_quantity(T, u.K, 'Absolute temperature')
    val = __model.radio_refractive_index(Pd, e, T)
    return val * u.dimensionless_unscaled


def water_vapour_pressure(T, P, H, type_hydrometeor='water'):
    """Determine the water vapour pressure.

    Parameters
    ----------
    T : number or Quantity
        Absolute temperature (C)
    P : number or Quantity
        Total atmospheric pressure (hPa)
    H : number or Quantity
        Relative humidity (%)
    type_hydrometeor : string
        Type of hydrometeor. Valid strings are 'water' and 'ice'


    Returns
    -------
    e: Quantity
        Water vapour pressure (hPa)



    References
    ----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    T = prepare_quantity(T, u.deg_C, 'Absolute temperature')
    P = prepare_quantity(P, u.hPa, 'Total atmospheric pressure')
    H = prepare_quantity(H, u.percent, 'Total atmospheric pressure')
    val = __model.water_vapour_pressure(T, P, H, type_hydrometeor)
    return val * u.hPa


def saturation_vapour_pressure(T, P, type_hydrometeor='water'):
    """Determine the saturation water vapour pressure.

    Parameters
    ----------
    T : number or Quantity
        Absolute temperature (C)
    P : number or Quantity
        Total atmospheric pressure (hPa)
    type_hydrometeor : string
        Type of hydrometeor. Valid strings are 'water' and 'ice'


    Returns
    -------
    e_s: Quantity
        Saturation water vapour pressure (hPa)



    References
    ----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    T = prepare_quantity(T, u.deg_C, 'Absolute temperature')
    P = prepare_quantity(P, u.hPa, 'Total atmospheric pressure')
    val = __model.saturation_vapour_pressure(T, P, type_hydrometeor)
    return val * u.hPa


def map_wet_term_radio_refractivity(lat, lon, p=50):
    """Determine the wet term of the radio refractivity.

    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points


    Returns
    -------
    N_wet: Quantity
        Wet term of the radio refractivity (-)



    References
    ----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.map_wet_term_radio_refractivity(lat, lon, p)
    return prepare_output_array(val, type_output) * u.g / u.m**3


def DN65(lat, lon, p):
    """Determine the statistics of the vertical gradient of radio
       refractivity in the lower 65 m from the surface of the Earth.

    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    p : number
        Percentage of time exceeded for p% of the average year


    Returns
    -------
    DN65_p: Quantity
        Vertical gradient of radio refractivity in the lowest 65 m from the
        surface of the Earth exceeded for p% of the average year



    References
    ----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.DN65(lat, lon, p)
    return prepare_output_array(val, type_output) * u.one


def DN1(lat, lon, p):
    """Determine the statistics of the vertical gradient of radio
       refractivity over 1 km layer from the surface.

    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    p : number
        Percentage of time exceeded for p% of the average year


    Returns
    -------
    DN1_p: Quantity
        Vertical gradient of radio refractivity over a 1 km layer from the
        surface exceeded for p% of the average year

    References
    ----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.DN1(lat, lon, p)
    return prepare_output_array(val, type_output) * u.one
