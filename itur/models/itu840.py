# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

import numpy as np
from astropy import units as u

from itur.models.itu1144 import bilinear_2D_interpolator
from itur.utils import (dataset_dir, prepare_input_array, prepare_output_array,
                        prepare_quantity, load_data_interpolator,
                        get_input_type)


def __fcn_columnar_content_reduced_liquid__(Lred, lat, lon, p):
    available_p = np.array(
        [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0,
         60.0, 70.0, 80.0, 90.0, 95.0, 99.0])

    if p in available_p:
        p_below = p_above = p
        pExact = True
    else:
        pExact = False
        idx = available_p.searchsorted(p, side='right') - 1
        idx = np.clip(idx, 0, len(available_p))

        p_below = available_p[idx]
        p_above = available_p[idx + 1]

    # Compute the values of Lred_a
    Lred_a = Lred(lat, lon, p_above)
    if not pExact:
        Lred_b = Lred(lat, lon, p_below)
        Lred = Lred_b + (Lred_a - Lred_b) * (np.log(p) - np.log(p_below)) \
            / (np.log(p_above) - np.log(p_below))
        return Lred
    else:
        return Lred_a


class __ITU840__():
    """Attenuation due to clouds and fog: This Recommendation provides methods
    to predict the attenuation due to clouds and fog on Earth-space paths.

    Available versions include:
    * P.840-4 (10/09) (Superseded)
    * P.840-5 (02/12) (Superseded)
    * P.840-6 (09/13) (Superseded)
    * P.840-7 (12/17) (Superseded)
    * P.840-8 (08/19) (Current version)

    Non-available versions include:
    * P.840-1 (08/94) (Superseded) - Tentative similar to P.840-4
    * P.840-2 (08/97) (Superseded) - Tentative similar to P.840-4
    * P.840-3 (10/99) (Superseded) - Tentative similar to P.840-4

    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.840 recommendation.

    def __init__(self, version=7):
        if version == 8:
            self.instance = _ITU840_8_()
        elif version == 7:
            self.instance = _ITU840_7_()
        elif version == 6:
            self.instance = _ITU840_6_()
        elif version == 5:
            self.instance = _ITU840_5_()
        elif version == 4:
            self.instance = _ITU840_4_()
        else:
            raise ValueError(
                'Version {0}  is not implemented for the ITU-R P.840 model.'
                .format(version))

    @property
    def __version__(self):
        return self.instance.__version__

    def specific_attenuation_coefficients(self, f, T):
        # Abstract method to compute the specific attenuation coefficients
        fcn = np.vectorize(self.instance.specific_attenuation_coefficients)
        return fcn(f, T)

    def columnar_content_reduced_liquid(self, lat, lon, p):
        # Abstract method to compute the columnar content of reduced liquid
        fcn = np.vectorize(__fcn_columnar_content_reduced_liquid__,
                           excluded=[0, 1, 2], otypes=[np.ndarray])
        return np.array(fcn(self.instance.Lred, lat, lon, p).tolist())

    def cloud_attenuation(self, lat, lon, el, f, p, Lred=None):
        # Abstract method to compute the cloud attenuation
        Kl = self.specific_attenuation_coefficients(f, T=0)
        if Lred is None:
            Lred = self.columnar_content_reduced_liquid(lat, lon, p)
        A = Lred * Kl / np.sin(np.deg2rad(el))

        return A

    def lognormal_approximation_coefficient(self, lat, lon):
        # Abstract method to compute the lognormal approximation coefficients
        return self.instance.lognormal_approximation_coefficient(lat, lon)


class _ITU840_8_():

    def __init__(self):
        self.__version__ = 8
        self.year = 2019
        self.month = 8
        self.link = 'https://www.itu.int/rec/R-REC-P.840-8-201908-I/en'

        self._Lred = {}
        self._M = None
        self._sigma = None
        self._Pclw = None

    # Note: The dataset used in recommendation 840-8 is the same as the
    # dataset use in recommendation 840-7. (The zip files included in
    # both recommendations are identical)

    def Lred(self, lat, lon, p):
        if not self._Lred:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '840/v7_lred_%s.npz')
            for p_load in ps:
                self._Lred[float(p_load)] = load_data_interpolator(
                    '840/v7_lat.npz', '840/v7_lon.npz',
                    d_dir % (str(p_load).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._Lred[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def M(self, lat, lon):
        if not self._M:
            self._M = load_data_interpolator(
                '840/v7_lat.npz', '840/v7_lon.npz',
                '840/v7_m.npz', bilinear_2D_interpolator, flip_ud=False)

        return self._M(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def sigma(self, lat, lon):
        if not self._sigma:
            self._sigma = load_data_interpolator(
                '840/v7_lat.npz', '840/v7_lon.npz',
                '840/v7_sigma.npz', bilinear_2D_interpolator, flip_ud=False)

        return self._sigma(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def Pclw(self, lat, lon):
        if not self._Pclw:
            self._Pclw = load_data_interpolator(
                '840/v7_lat.npz', '840/v7_lon.npz',
                '840/v7_pclw.npz', bilinear_2D_interpolator, flip_ud=False)

        return self._Pclw(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    @staticmethod
    def specific_attenuation_coefficients(f, T):
        """
        """
        return _ITU840_6_.specific_attenuation_coefficients(f, T)

    def lognormal_approximation_coefficient(self, lat, lon):
        m = self.M(lat, lon)
        sigma = self.sigma(lat, lon)
        Pclw = self.Pclw(lat, lon)

        return m, sigma, Pclw


class _ITU840_7_():

    def __init__(self):
        self.__version__ = 7
        self.year = 2017
        self.month = 12
        self.link = 'https://www.itu.int/rec/R-REC-P.840-7-201712-I/en'

        self._Lred = {}
        self._M = {}
        self._sigma = {}
        self._Pclw = {}

    def Lred(self, lat, lon, p):
        if not self._Lred:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '840/v7_lred_%s.npz')
            for p_load in ps:
                self._Lred[float(p_load)] = load_data_interpolator(
                    '840/v7_lat.npz', '840/v7_lon.npz',
                    d_dir % (str(p_load).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._Lred[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def M(self, lat, lon):
        if not self._M:
            self._M = load_data_interpolator(
                '840/v7_lat.npz', '840/v7_lon.npz',
                '840/v7_m.npz', bilinear_2D_interpolator, flip_ud=False)

        return self._M(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def sigma(self, lat, lon):
        if not self._sigma:
            self._sigma = load_data_interpolator(
                '840/v7_lat.npz', '840/v7_lon.npz',
                '840/v7_sigma.npz', bilinear_2D_interpolator, flip_ud=False)

        return self._sigma(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def Pclw(self, lat, lon):
        if not self._Pclw:
            self._Pclw = load_data_interpolator(
                '840/v7_lat.npz', '840/v7_lon.npz',
                '840/v7_pclw.npz', bilinear_2D_interpolator, flip_ud=False)

        return self._Pclw(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    @staticmethod
    def specific_attenuation_coefficients(f, T):
        """
        """
        return _ITU840_6_.specific_attenuation_coefficients(f, T)

    def lognormal_approximation_coefficient(self, lat, lon):
        # TODO: This is the wrong method, Need to update
        m = self.M(lat, lon)
        sigma = self.sigma(lat, lon)
        Pclw = self.Pclw(lat, lon)

        return m, sigma, Pclw


class _ITU840_6_():

    def __init__(self):
        self.__version__ = 6
        self.year = 2013
        self.month = 9
        self.link = 'https://www.itu.int/rec/R-REC-P.840-6-201202-I/en'

        self._Lred = {}
        self._M = {}
        self._sigma = {}
        self._Pclw = {}

    def Lred(self, lat, lon, p):
        if not self._Lred:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '840/v6_lred_%s.npz')
            for p_load in ps:
                self._Lred[float(p_load)] = load_data_interpolator(
                    '840/v6_lat.npz', '840/v6_lon.npz',
                    d_dir % (str(p_load).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._Lred[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def M(self, lat, lon):
        if not self._M:
            self._M = load_data_interpolator(
                '840/v6_lat.npz', '840/v6_lon.npz',
                '840/v6_m.npz', bilinear_2D_interpolator, flip_ud=False)

        return self._M(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def sigma(self, lat, lon):
        if not self._sigma:
            self._sigma = load_data_interpolator(
                '840/v6_lat.npz', '840/v6_lon.npz',
                '840/v6_sigma.npz', bilinear_2D_interpolator, flip_ud=False)

        return self._sigma(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def Pclw(self, lat, lon):
        if not self._Pclw:
            self._Pclw = load_data_interpolator(
                '840/v6_lat.npz', '840/v6_lon.npz',
                '840/v6_pclw.npz', bilinear_2D_interpolator, flip_ud=False)

        return self._Pclw(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    @staticmethod
    def specific_attenuation_coefficients(f, T):
        """
        """
        if np.any(f > 1000):
            raise ValueError('Frequency must be introduced in GHz and the '
                             'maximum range is 1000 GHz')

        T_kelvin = T + 273.15
        theta = 300.0 / T_kelvin                # Eq. 9

        # Compute the values of the epsilons
        epsilon0 = 77.66 + 103.3 * (theta - 1)  # Eq. 6
        epsilon1 = 0.0671 * epsilon0            # Eq. 7
        epsilon2 = 3.52                         # Eq. 8

        # Compute the principal and secondary relacation frequencies
        fp = 20.20 - 146 * (theta - 1) + 316.0 * (theta - 1)**2     # Eq. 10
        fs = 39.8 * fp                                              # Eq. 11

        # Compute the dielectric permitivity of water
        epsilonp = (epsilon0 - epsilon1) / (1 + (f / fp) ** 2) + \
            (epsilon1 - epsilon2) / (1 + (f / fs) ** 2) + epsilon2  # Eq. 5

        epsilonpp = f * (epsilon0 - epsilon1) / (fp * (1 + (f / fp)**2)) + \
            f * (epsilon1 - epsilon2) / (fs * (1 + (f / fs)**2))       # Eq. 4

        eta = (2 + epsilonp) / epsilonpp                    # Eq. 3
        Kl = (0.819 * f) / (epsilonpp * (1 + eta**2))       # Eq. 2

        return Kl       # Specific attenuation coefficient  (dB/km)/(g/m3)

    def lognormal_approximation_coefficient(self, lat, lon):
        m = self.M(lat, lon)
        sigma = self.sigma(lat, lon)
        Pclw = self.Pclw(lat, lon)

        return m, sigma, Pclw


class _ITU840_5_():

    def __init__(self):
        self.__version__ = 5
        self.year = 2012
        self.month = 2
        self.link = 'https://www.itu.int/rec/R-REC-P.840-5-201202-S/en'

        self._Lred = {}
        self._M = {}
        self._sigma = {}
        self._Pclw = {}

    def Lred(self, lat, lon, p):
        if not self._Lred:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '840/v4_esawred_%s.npz')
            for p_load in ps:
                self._Lred[float(p_load)] = load_data_interpolator(
                    '840/v4_lat.npz', '840/v4_lon.npz',
                    d_dir % (str(p_load).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._Lred[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def M(self, lat, lon):
        if not self._M:
            self._M = load_data_interpolator(
                '840/v6_lat.npz', '840/v6_lon.npz',
                '840/v4_wred_lognormal_mean.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return self._M(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def sigma(self, lat, lon):
        if not self._sigma:
            self._sigma = load_data_interpolator(
                '840/v6_lat.npz', '840/v6_lon.npz',
                '840/v4_wred_lognormal_stdev.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return self._sigma(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def Pclw(self, lat, lon):
        if not self._Pclw:
            self._Pclw = load_data_interpolator(
                '840/v6_lat.npz', '840/v6_lon.npz',
                '840/v4_wred_lognormal_pclw.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return self._Pclw(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    @staticmethod
    def specific_attenuation_coefficients(f, T):
        """
        """
        return _ITU840_4_.specific_attenuation_coefficients(f, T)

    def lognormal_approximation_coefficient(self, lat, lon):
        m = self.M(lat, lon)
        sigma = self.sigma(lat, lon)
        Pclw = self.Pclw(lat, lon)

        return m, sigma, Pclw


class _ITU840_4_():

    def __init__(self):
        self.__version__ = 4
        self.year = 2013
        self.month = 9
        self.link = 'https://www.itu.int/rec/R-REC-P.840-6-201202-I/en'

        self._Lred = {}
        self._M = {}
        self._sigma = {}
        self._Pclw = {}

    def Lred(self, lat, lon, p):
        if not self._Lred:
            ps = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30,
                  50, 60, 70, 80, 90, 95, 99]
            d_dir = os.path.join(dataset_dir, '840/v4_esawred_%s.npz')
            for p_load in ps:
                self._Lred[float(p_load)] = load_data_interpolator(
                    '840/v4_lat.npz', '840/v4_lon.npz',
                    d_dir % (str(p_load).replace('.', '')),
                    bilinear_2D_interpolator, flip_ud=False)

        return self._Lred[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def M(self, lat, lon):
        if not self._M:
            self._M = load_data_interpolator(
                '840/v6_lat.npz', '840/v6_lon.npz',
                '840/v4_wred_lognormal_mean.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return self._M(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def sigma(self, lat, lon):
        if not self._sigma:
            self._sigma = load_data_interpolator(
                '840/v6_lat.npz', '840/v6_lon.npz',
                '840/v4_wred_lognormal_stdev.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return self._sigma(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def Pclw(self, lat, lon):
        if not self._Pclw:
            self._Pclw = load_data_interpolator(
                '840/v6_lat.npz', '840/v6_lon.npz',
                '840/v4_wred_lognormal_pclw.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return self._Pclw(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    @staticmethod
    def specific_attenuation_coefficients(f, T):
        """
        """
        if np.any(f > 1000):
            raise ValueError(
                'Frequency must be introduced in GHz and the maximum range'
                ' is 1000 GHz')

        T_kelvin = T + 273.15
        theta = 300.0 / T_kelvin                # Eq. 9

        # Compute the values of the epsilons
        epsilon0 = 77.66 + 103.3 * (theta - 1)  # Eq. 6
        epsilon1 = 5.48                         # Eq. 7
        epsilon2 = 3.51                         # Eq. 8

        # Compute the principal and secondary relacation frequencies
        fp = 20.09 - 142 * (theta - 1) + 294.0 * (theta - 1)**2     # Eq. 10
        fs = 590 - 1500 * (theta - 1)                               # Eq. 11

        # Compute the dielectric permitivity of water
        epsilonp = (epsilon0 - epsilon1) / (1 + (f / fp) ** 2) + \
            (epsilon1 - epsilon2) / (1 + (f / fs) ** 2) + epsilon2  # Eq. 5

        epsilonpp = f * (epsilon0 - epsilon1) / (fp * (1 + (f / fp)**2)) + \
            f * (epsilon1 - epsilon2) / (fs * (1 + (f / fs)**2))       # Eq. 4

        eta = (2 + epsilonp) / epsilonpp                    # Eq. 3
        Kl = (0.819 * f) / (epsilonpp * (1 + eta**2))       # Eq. 2

        return Kl       # Specific attenuation coefficient  (dB/km)/(g/m3)

    def lognormal_approximation_coefficient(self, lat, lon):
        m = self.M(lat, lon)
        sigma = self.sigma(lat, lon)
        Pclw = self.Pclw(lat, lon)

        return m, sigma, Pclw


__model = __ITU840__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.840 recommendation currently being used.


    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          * 8: Activates recommendation ITU-R P.840-8 (08/19) (Current version)
          * 7: Activates recommendation ITU-R P.840-7 (12/17) (Superseded)
          * 6: Activates recommendation ITU-R P.840-6 (09/13) (Superseded)
          * 5: Activates recommendation ITU-R P.840-5 (02/12) (Superseded)
          * 4: Activates recommendation ITU-R P.840-4 (10/09) (Superseded)
    """
    global __model
    __model = __ITU840__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.840 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__


def specific_attenuation_coefficients(f, T):
    """
    Compute the specific attenuation coefficient for cloud attenuation.

    A method to compute the specific attenuation coefficient. The method is
    based on Rayleigh scattering, which uses a double-Debye model for the
    dielectric permittivity of water.

    This model can be used to calculate the value of the specific attenuation
    coefficient for frequencies up to 1000 GHz:


    Parameters
    ----------
    f : number
        Frequency (GHz)
    T : number
        Temperature (degrees C)


    Returns
    -------
    Kl: numpy.ndarray
        Specific attenuation coefficient (dB/km)


    References
    ----------
    [1] Attenuation due to clouds and fog:
    https://www.itu.int/rec/R-REC-P.840/en
    """
    f = prepare_quantity(f, u.GHz, 'Frequency')
    T = prepare_quantity(T, u.deg_C, 'Temperature')
    return __model.specific_attenuation_coefficients(f, T)


def columnar_content_reduced_liquid(lat, lon, p):
    """
    Compute the total columnar contents of reduced cloud liquid water.

    A method to compute the total columnar content of reduced cloud liquid
    water, Lred (kg/m2), exceeded for p% of the average year


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
    Lred: numpy.ndarray
        Total columnar content of reduced cloud liquid water, Lred (kg/m2),
        exceeded for p% of the average year



    References
    ----------
    [1] Attenuation due to clouds and fog:
    https://www.itu.int/rec/R-REC-P.840/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.columnar_content_reduced_liquid(lat, lon, p)
    return prepare_output_array(val, type_output) * u.kg / u.m**2


def cloud_attenuation(lat, lon, el, f, p, Lred=None):
    """
    Compute the cloud attenuation in a slant path.

    A method to estimate the attenuation due to clouds along slant paths for
    a given probability. If local measured data of the total columnar content
    of cloud liquid water reduced to a temperature of 273.15 K, Lred, is
    available from other sources, (e.g., from ground radiometric measurements,
    Earth observation products, or meteorological numerical products), the
    value should be used directly.

    The value of the cloud attenuation is computed as:

    .. math::
      A=\\frac{L_{red}(\\text{lat}, \\text{lon}, p, T) \\cdot K_l(f, T)}{\\sin(\\text{el})}


    where:
        * :math:`L_{red}` : total columnar content of liquid water reduced to a
          temperature of 273.15 K (kg/m2);
        * :math:`K_l` : specific attenuation coefficient ((dB/km)/(g/m3));
        * :math:`el` : path elevation angle (deg).
        * :math:`f` : frequency (GHz).
        * :math:`p` : Percentage of time exceeded for p% of the average year (%).
        * :math:`T` : temperature (K). Equal to 273.15 K.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    el : number, sequence, or numpy.ndarray
        Elevation angle of the receiver points (deg)
    f : number
        Frequency (GHz)
    p : number
         Percentage of time exceeded for p% of the average year
    Lred: number
        Total columnar contents of reduced cloud liquid water. (kg/m2)



    Returns
    -------
    A: numpy.ndarray
        Cloud attenuation, A (dB), exceeded for p% of the average year



    References
    ----------
    [1] Attenuation due to clouds and fog:
    https://www.itu.int/rec/R-REC-P.840/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    Lred = prepare_quantity(
        Lred, u.kg / u.m**2,
        'Total columnar contents of reduced cloud liquid water.')
    val = __model.cloud_attenuation(lat, lon, el, f, p, Lred)
    
    # The values of attenuation cannot be negative. The ITU models end up
    # giving out negative values for certain inputs
    val[val < 0] = 0

    return prepare_output_array(val, type_output) * u.dB


def lognormal_approximation_coefficient(lat, lon):
    """
    Total columnar contents of cloud liquid water distribution coefficients.

    The annual statistics of the total columnar content of reduced cloud
    liquid water content can be approximated by a log-normal distribution.
    This function computes the coefficients for the mean, :math:`m`,
    standard deviation, :math:`\\sigma`, and probability of non-zero reduced
    total columnar content of cloud liquid water, :math:`Pclw`, for such the
    log-normal distribution.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points


    Returns
    -------
    m: numpy.ndarray
        Mean of the lognormal distribution
    σ: numpy.ndarray
        Standard deviation of the lognormal distribution
    Pclw: numpy.ndarray
        Probability of cloud liquid water of the lognormal distribution



    References
    ----------
    [1] Attenuation due to clouds and fog:
    https://www.itu.int/rec/R-REC-P.840/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.lognormal_approximation_coefficient(lat, lon)
    u_adim = u.dimensionless_unscaled
    return (prepare_output_array(val[0], type_output) * u_adim,
            prepare_output_array(val[1], type_output) * u_adim,
            prepare_output_array(val[2], type_output) * u_adim)
