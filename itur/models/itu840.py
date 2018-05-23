# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import os
from astropy import units as u

from itur.models.itu1144 import bilinear_2D_interpolator
from itur.utils import load_data, dataset_dir, prepare_input_array, \
    prepare_output_array, prepare_quantity, memory


class __ITU840():
    """Attenuation due to clouds and fog: This Recommendation provides methods
    to predict the attenuation due to clouds and fog on Earth-space paths.

    Available versions include:
    * P.840-4 (10/09) (Superseded)
    * P.840-5 (02/12) (Superseded)
    * P.840-6 (09/13) (Superseded)
    * P.840-7 (12/17) (Current version)

    Non-available versions include:
    * P.840-1 (08/94) (Superseded) - Tentative similar to P.840-4
    * P.840-2 (08/97) (Superseded) - Tentative similar to P.840-4
    * P.840-3 (10/99) (Superseded) - Tentative similar to P.840-4

    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.840 recommendation.

    def __init__(self, version=7):
        if version == 7:
            self.instance = _ITU840_7()
        elif version == 6:
            self.instance = _ITU840_6()
        elif version == 5:
            self.instance = _ITU840_5()
        elif version == 4:
            self.instance = _ITU840_4()
        else:
            raise ValueError(
                'Version {0}  is not implemented for the ITU-R P.840 model.'
                .format(version))

        self._Lred = {}
        self._M = {}
        self._sigma = {}
        self._Pclw = {}

    @property
    def __version__(self):
        return self.instance.__version__

    def specific_attenuation_coefficients(self, f, T):
        # Abstract method to compute the specific attenuation coefficients
        fcn = np.vectorize(self.instance.specific_attenuation_coefficients)
        return fcn(f, T)

    def columnar_content_reduced_liquid(self, lat, lon, p):
        # Abstract method to compute the columnar content of reduced liquid
        fcn = np.vectorize(self.instance.columnar_content_reduced_liquid,
                           excluded=[0, 1], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, p).tolist())

    def cloud_attenuation(self, lat, lon, el, f, p):
        # Abstract method to compute the cloud attenuation
        Kl = self.specific_attenuation_coefficients(f, T=0)
        Lred = self.columnar_content_reduced_liquid(lat, lon, p)
        A = Lred * Kl / np.sin(np.deg2rad(el))

        return A

    def lognormal_approximation_coefficient(self, lat, lon):
        # Abstract method to compute the lognormal approximation coefficients
        return self.instance.lognormal_approximation_coefficient(lat, lon)


class _ITU840_7():

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
                  50, 60, 70, 80, 90, 95]
            d_dir = os.path.join(dataset_dir, '840/v7_Lred_%s.txt')
            lats = load_data(os.path.join(dataset_dir, '840/v7_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v7_Lon.txt'))
            for p_load in ps:
                vals = load_data(d_dir % (str(p_load).replace('.', '')))
                self._Lred[float(p_load)] = bilinear_2D_interpolator(
                    lats, lons, vals)

        return self._Lred[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def M(self, lat, lon):
        if not self._M:
            vals = load_data(os.path.join(dataset_dir, '840/v7_M.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v7_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v7_Lon.txt'))
            self._M = bilinear_2D_interpolator(lats, lons, vals)

        return self._M(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def sigma(self, lat, lon):
        if not self._sigma:
            vals = load_data(os.path.join(dataset_dir, '840/v7_sigma.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v7_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v7_Lon.txt'))
            self._sigma = bilinear_2D_interpolator(lats, lons, vals)

        return self._sigma(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def Pclw(self, lat, lon):
        if not self._Pclw:
            vals = load_data(os.path.join(dataset_dir, '840/v6_Pclw.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._Pclw = bilinear_2D_interpolator(lats, lons, vals)

        return self._Pclw(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def specific_attenuation_coefficients(self, f, T):
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

    def columnar_content_reduced_liquid(self, lat, lon, p):
        """
        """
        available_p = np.array(
            [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0,
             60.0, 70.0, 80.0, 90.0, 95.0])

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
        Lred_a = self.Lred(lat, lon, p_above)
        if not pExact:
            Lred_b = self.Lred(lat, lon, p_below)
            Lred = Lred_b + (Lred_a - Lred_b) * (np.log(p) - np.log(p_below)) \
                / (np.log(p_above) - np.log(p_below))
            return Lred
        else:
            return Lred_a

    def lognormal_approximation_coefficient(self, lat, lon):
        m = self.M(lat, lon)
        sigma = self.sigma(lat, lon)
        Pclw = self.Pclw(lat, lon)

        return m, sigma, Pclw


class _ITU840_6():

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
                  50, 60, 70, 80, 90, 95]
            d_dir = os.path.join(dataset_dir, '840/v6_Lred_%s.txt')
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            for p_load in ps:
                vals = load_data(d_dir % (str(p_load).replace('.', '')))
                self._Lred[float(p_load)] = bilinear_2D_interpolator(
                    lats, lons, vals)

        return self._Lred[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def M(self, lat, lon):
        if not self._M:
            vals = load_data(os.path.join(dataset_dir, '840/v6_M.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._M = bilinear_2D_interpolator(lats, lons, vals)

        return self._M(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def sigma(self, lat, lon):
        if not self._sigma:
            vals = load_data(os.path.join(dataset_dir, '840/v6_sigma.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._sigma = bilinear_2D_interpolator(lats, lons, vals)

        return self._sigma(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def Pclw(self, lat, lon):
        if not self._Pclw:
            vals = load_data(os.path.join(dataset_dir, '840/v6_Pclw.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._Pclw = bilinear_2D_interpolator(lats, lons, vals)

        return self._Pclw(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def specific_attenuation_coefficients(self, f, T):
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

    def columnar_content_reduced_liquid(self, lat, lon, p):
        """
        """
        available_p = np.array(
            [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0,
             60.0, 70.0, 80.0, 90.0, 95.0])

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
        Lred_a = self.Lred(lat, lon, p_above)
        if not pExact:
            Lred_b = self.Lred(lat, lon, p_below)
            Lred = Lred_b + (Lred_a - Lred_b) * (np.log(p) - np.log(p_below)) \
                / (np.log(p_above) - np.log(p_below))
            return Lred
        else:
            return Lred_a

    def lognormal_approximation_coefficient(self, lat, lon):
        m = self.M(lat, lon)
        sigma = self.sigma(lat, lon)
        Pclw = self.Pclw(lat, lon)

        return m, sigma, Pclw


class _ITU840_5():

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
                  50, 60, 70, 80, 90, 95]
            d_dir = os.path.join(dataset_dir, '840/v4_ESAWRED_%s.txt')
            lats = load_data(os.path.join(dataset_dir, '840/v4_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v4_Lon.txt'))
            for p_load in ps:
                vals = load_data(d_dir % (str(p_load).replace('.', '')))
                self._Lred[float(p_load)] = bilinear_2D_interpolator(
                    lats, lons, vals)

        return self._Lred[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def M(self, lat, lon):
        if not self._M:
            vals = load_data(os.path.join(dataset_dir,
                                          '840/v4_WRED_LOGNORMAL_MEAN.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._M = bilinear_2D_interpolator(lats, lons, vals)

        return self._M(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def sigma(self, lat, lon):
        if not self._sigma:
            vals = load_data(os.path.join(dataset_dir,
                                          '840/v4_WRED_LOGNORMAL_STDEV.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._sigma = bilinear_2D_interpolator(lats, lons, vals)

        return self._sigma(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def Pclw(self, lat, lon):
        if not self._Pclw:
            vals = load_data(os.path.join(dataset_dir,
                                          '840/v4_WRED_LOGNORMAL_PCLW.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._Pclw = bilinear_2D_interpolator(lats, lons, vals)

        return self._Pclw(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def specific_attenuation_coefficients(self, f, T):
        """
        """
        if np.any(f > 1000):
            raise ValueError(
                'Frequency must be introduced in GHz and the maximum range '
                'is 1000 GHz')

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

    def columnar_content_reduced_liquid(self, lat, lon, p):
        """
        """
        available_p = np.array(
            [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0,
             60.0, 70.0, 80.0, 90.0, 95.0])

        if p in available_p:
            p_below = p_above = p
            pExact = True
        else:
            pExact = False
            idx = available_p.searchsorted(p, side='right') - 1
            idx = np.clip(idx, 0, len(available_p))

            p_below = available_p[idx]
            p_above = available_p[idx + 1]

        Lred_a = self.Lred(lat, lon, p_above)
        if not pExact:
            Lred_b = self.Lred(lat, lon, p_below)
            Lred = Lred_b + (Lred_a - Lred_b) * (np.log(p) - np.log(p_below)) \
                / (np.log(p_above) - np.log(p_below))
            return Lred
        else:
            return Lred_a

    def lognormal_approximation_coefficient(self, lat, lon):
        m = self.M(lat, lon)
        sigma = self.sigma(lat, lon)
        Pclw = self.Pclw(lat, lon)

        return m, sigma, Pclw


class _ITU840_4():

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
                  50, 60, 70, 80, 90, 95]
            d_dir = os.path.join(dataset_dir, '840/v4_ESAWRED_%s.txt')
            lats = load_data(os.path.join(dataset_dir, '840/v4_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v4_Lon.txt'))
            for p_load in ps:
                vals = load_data(d_dir % (str(p_load).replace('.', '')))
                self._Lred[float(p_load)] = bilinear_2D_interpolator(
                    lats, lons, vals)

        return self._Lred[float(p)](
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def M(self, lat, lon):
        if not self._M:
            vals = load_data(os.path.join(dataset_dir,
                                          '840/v4_WRED_LOGNORMAL_MEAN.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._M = bilinear_2D_interpolator(lats, lons, vals)

        return self._M(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def sigma(self, lat, lon):
        if not self._sigma:
            vals = load_data(os.path.join(dataset_dir,
                                          '840/v4_WRED_LOGNORMAL_STDEV.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._sigma = bilinear_2D_interpolator(lats, lons, vals)

        return self._sigma(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def Pclw(self, lat, lon):
        if not self._Pclw:
            vals = load_data(os.path.join(dataset_dir,
                                          '840/v4_WRED_LOGNORMAL_PCLW.txt'))
            lats = load_data(os.path.join(dataset_dir, '840/v6_Lat.txt'))
            lons = load_data(os.path.join(dataset_dir, '840/v6_Lon.txt'))
            self._Pclw = bilinear_2D_interpolator(lats, lons, vals)

        return self._Pclw(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def specific_attenuation_coefficients(self, f, T):
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

    def columnar_content_reduced_liquid(self, lat, lon, p):
        """
        """
        available_p = np.array(
            [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0,
             60.0, 70.0, 80.0, 90.0, 95.0])

        if p in available_p:
            p_below = p_above = p
            pExact = True
        else:
            pExact = False
            idx = available_p.searchsorted(p, side='right') - 1
            idx = np.clip(idx, 0, len(available_p))

            p_below = available_p[idx]
            p_above = available_p[idx + 1]

        Lred_a = self.Lred(lat, lon, p_above)
        if not pExact:
            Lred_b = Lred_a = self.Lred(lat, lon, p_below)
            Lred = Lred_b + (Lred_a - Lred_b) * (np.log(p) - np.log(p_below)) \
                / (np.log(p_above) - np.log(p_below))
            return Lred
        else:
            return Lred_a

    def lognormal_approximation_coefficient(self, lat, lon):
        m = self.M(lat, lon)
        sigma = self.sigma(lat, lon)
        Pclw = self.Pclw(lat, lon)

        return m, sigma, Pclw


__model = __ITU840()


def change_version(new_version):
    """
    Change the version of the ITU-R P.840 recommendation currently being used.


    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
        * 4: P.840-4 (10/09) (Superseded)
        * 5: P.840-5 (02/12) (Superseded)
        * 6: P.840-6 (09/13) (Current version)
    """
    global __model
    __model = __ITU840(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.840 recommendation currently being used.
    """
    global __model
    return __model.__version__


def specific_attenuation_coefficients(f, T):
    """
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
    global __model
    f = prepare_quantity(f, u.GHz, 'Frequency')
    T = prepare_quantity(T, u.deg_C, 'Temperature')
    return __model.specific_attenuation_coefficients(f, T)


@memory.cache
def columnar_content_reduced_liquid(lat, lon, p):
    """
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
    global __model
    type_output = type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.columnar_content_reduced_liquid(lat, lon, p)
    return prepare_output_array(val, type_output) * u.kg / u.m**2


def cloud_attenuation(lat, lon, el, f, p):
    """
    A method to estimate the attenuation due to clouds along slant paths for
    a given probability.


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


    Returns
    -------
    p: numpy.ndarray
        Rainfall rate exceeded for p% of the average year



    References
    ----------
    [1] Attenuation due to clouds and fog:
    https://www.itu.int/rec/R-REC-P.840/en
    """
    global __model
    type_output = type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    val = __model.cloud_attenuation(lat, lon, el, f, p)
    return prepare_output_array(val, type_output) * u.dB


@memory.cache
def lognormal_approximation_coefficient(lat, lon):
    """
    A method to estimate the paramerts of the lognormla distribution used to
    approximate the total columnar content of cloud liquid water


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
    sigma: numpy.ndarray
        Standard deviation of the lognormal distribution
    Pclw: numpy.ndarray
        Probability of liquid water of the lognormal distribution



    References
    ----------
    [1] Attenuation due to clouds and fog:
    https://www.itu.int/rec/R-REC-P.840/en
    """
    type_output = type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.lognormal_approximation_coefficient(lat, lon)
    u_adim = u.dimensionless_unscaled
    return (prepare_output_array(val[0], type_output) * u_adim,
            prepare_output_array(val[1], type_output) * u_adim,
            prepare_output_array(val[2], type_output) * u_adim)
