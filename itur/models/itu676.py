# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import warnings


import numpy as np
from astropy import units as u

from itur.models.itu453 import radio_refractive_index
from itur.models.itu835 import (standard_pressure, standard_temperature,
                                standard_water_vapour_density)
from itur.models.itu836 import total_water_vapour_content
from itur.models.itu1511 import topographic_altitude
from itur.utils import (prepare_quantity, prepare_output_array, get_input_type,
                        prepare_input_array, load_data, dataset_dir)


def __gamma0_exact__(self, f, p, rho, T):
    # T in Kelvin
    # e : water vapour partial pressure in hPa (total barometric pressure
    # ptot = p + e)
    theta = 300 / T
    e = rho * T / 216.7

    f_ox = self.f_ox

    D_f_ox = self.a3 * 1e-4 * (p * (theta ** (0.8 - self.a4)) +
                               1.1 * e * theta)

    D_f_ox = np.sqrt(D_f_ox**2 + 2.25e-6)

    delta_ox = (self.a5 + self.a6 * theta) * 1e-4 * (p + e) * theta**0.8

    F_i_ox = f / f_ox * ((D_f_ox - delta_ox * (f_ox - f)) /
                         ((f_ox - f) ** 2 + D_f_ox ** 2) +
                         (D_f_ox - delta_ox * (f_ox + f)) /
                         ((f_ox + f) ** 2 + D_f_ox ** 2))

    Si_ox = self.a1 * 1e-7 * p * theta**3 * np.exp(self.a2 * (1 - theta))

    N_pp_ox = Si_ox * F_i_ox

    d = 5.6e-4 * (p + e) * theta**0.8

    N_d_pp = f * p * theta**2 * \
        (6.14e-5 / (d * (1 + (f / d)**2)) +
         1.4e-12 * p * theta**1.5 / (1 + 1.9e-5 * f**1.5))

    N_pp = N_pp_ox.sum() + N_d_pp

    gamma = 0.1820 * f * N_pp           # Eq. 1 [dB/km]
    return gamma


def __gammaw_exact__(self, f, p, rho, T):
    # T in Kelvin
    # e : water vapour partial pressure in hPa (total barometric pressure
    # ptot = p + e)
    theta = 300 / T
    e = rho * T / 216.7

    f_wv = self.f_wv

    D_f_wv = self.b3 * 1e-4 * (p * theta ** self.b4 +
                               self.b5 * e * theta ** self.b6)

    D_f_wv = 0.535 * D_f_wv + \
        np.sqrt(0.217 * D_f_wv**2 + 2.1316e-12 * f_wv**2 / theta)

    F_i_wv = f / f_wv * ((D_f_wv) / ((f_wv - f)**2 + D_f_wv**2) +
                         (D_f_wv) / ((f_wv + f)**2 + D_f_wv**2))

    Si_wv = self.b1 * 1e-1 * e * theta**3.5 * np.exp(self.b2 * (1 - theta))

    N_pp_wv = Si_wv * F_i_wv

    N_pp = N_pp_wv.sum()

    gamma = 0.1820 * f * N_pp           # Eq. 1 [dB/km]
    return gamma


class __ITU676__():
    """Attenuation by atmospheric gases.

    Available versions include:
       * P.676-9 (02/12) (Superseded)
       * P.676-10 (09/13) (Superseded)
       * P.676-11 (09/16) (Superseded)
       * P.676-11 (08/19) (Current version)

    Not available versions:
       * P.676-1 (03/92) (Superseded)
       * P.676-2 (10/95) (Superseded)
       * P.676-3 (08/97) (Superseded)
       * P.676-4 (10/99) (Superseded)
       * P.676-5 (02/01) (Superseded)
       * P.676-6 (03/05) (Superseded)
       * P.676-7 (02/07) (Superseded)
       * P.676-8 (10/09) (Superseded)

    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.676 recommendation.

    def __init__(self, version=12):
        if version == 12:
            self.instance = _ITU676_12_()
        elif version == 11:
            self.instance = _ITU676_11_()
        elif version == 10:
            self.instance = _ITU676_10_()
        elif version == 9:
            self.instance = _ITU676_9_()
#        elif version == 8:
#            self.instance = _ITU676_8()
#        elif version == 7:
#            self.instance = _ITU676_7()
#        elif version == 6:
#            self.instance = _ITU676_6()
#        elif version == 5:
#            self.instance = _ITU676_5()
#        elif version == 4:
#            self.instance = _ITU676_4()
#        elif version == 3:
#            self.instance = _ITU676_3()
#        elif version == 2:
#            self.instance = _ITU676_2()
#        elif version == 1:
#            self.instance = _ITU676_1()
        else:
            raise ValueError(
                'Version {0} is not implemented for the ITU-R P.676 model.'
                .format(version))

    @property
    def __version__(self):
        return self.instance.__version__

    def gaseous_attenuation_terrestrial_path(self, r, f, el, rho, P, T, mode):
        # Abstract method to compute the gaseous attenuation over a slant path
        fcn = np.vectorize(self.instance.gaseous_attenuation_terrestrial_path)
        return fcn(r, f, el, rho, P, T, mode)

    def gaseous_attenuation_inclined_path(
            self, f, el, rho, P, T, h1, h2, mode):
        # Abstract method to compute the gaseous attenuation over an inclined
        # path
        fcn = np.vectorize(self.instance.gaseous_attenuation_inclined_path)
        return fcn(f, el, rho, P, T, h1, h2, mode)

    def gaseous_attenuation_slant_path(self, f, el, rho, P, T, V_t, h, mode):
        # Abstract method to compute the gaseous attenuation over a slant path
        fcn = np.vectorize(self.instance.gaseous_attenuation_slant_path)
        return fcn(f, el, rho, P, T, V_t, h, mode)

    def slant_inclined_path_equivalent_height(self, f, p, rho, T):
        fcn = np.vectorize(self.instance.slant_inclined_path_equivalent_height,
                           excluded=[0], otypes=[np.ndarray])
        return np.array(fcn(f, p, rho, T).tolist())

    def zenit_water_vapour_attenuation(
            self, lat, lon, p, f, V_t=None, h=None):
        # Abstract method to compute the water vapour attenuation over the
        # slant path
        fcn = np.vectorize(self.instance.zenit_water_vapour_attenuation,
                           excluded=[0, 1, 4, 5], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, p, f, V_t, h).tolist())

    def gamma_exact(self, f, p, rho, t):
        # Abstract method to compute the specific attenuation using the
        # line-by-line method
        fcn = np.vectorize(self.instance.gamma_exact)
        return fcn(f, p, rho, t)

    def gammaw_exact(self, f, p, rho, t):
        # Abstract method to compute the specific attenuation due to water
        # vapour
        fcn = np.vectorize(self.instance.gammaw_exact)
        return fcn(f, p, rho, t)

    def gamma0_exact(self, f, p, rho, t):
        # Abstract method to compute the specific attenuation due to dry
        # atmoshere
        fcn = np.vectorize(self.instance.gamma0_exact)
        return fcn(f, p, rho, t)

    def gammaw_approx(self, f, p, rho, t):
        # Abstract method to compute the specific attenuation due to water
        # vapour
        fcn = np.vectorize(self.instance.gammaw_approx)
        with np.errstate(invalid='ignore'):
            return fcn(f, p, rho, t)

    def gamma0_approx(self, f, p, rho, t):
        # Abstract method to compute the specific attenuation due to dry
        # atmoshere
        fcn = np.vectorize(self.instance.gamma0_approx)
        with np.errstate(invalid='ignore'):
            return fcn(f, p, rho, t)


class _ITU676_12_():

    tmp = load_data(os.path.join(dataset_dir, '676/v12_lines_oxygen.txt'),
                    skip_header=1)
    f_ox = tmp[:, 0]
    a1 = tmp[:, 1]
    a2 = tmp[:, 2]
    a3 = tmp[:, 3]
    a4 = tmp[:, 4]
    a5 = tmp[:, 5]
    a6 = tmp[:, 6]

    tmp = load_data(os.path.join(dataset_dir,
                                 '676//v12_lines_water_vapour.txt'),
                    skip_header=1)
    f_wv = tmp[:, 0]
    b1 = tmp[:, 1]
    b2 = tmp[:, 2]
    b3 = tmp[:, 3]
    b4 = tmp[:, 4]
    b5 = tmp[:, 5]
    b6 = tmp[:, 6]

    # Coefficients in table 3
    t2_coeffs = [(0.1597, 118.750334),
                 (0.1066, 368.498246),
                 (0.1325, 424.763020),
                 (0.1242, 487.249273),
                 (0.0938, 715.392902),
                 (0.1448, 773.839490),
                 (0.1374, 834.145546)]

    # Coefficients in table 4
    hw_coeffs = [(22.23508, 1.52, 2.56),
                 (183.310087, 7.62, 10.2),
                 (325.152888, 1.56, 2.7),
                 (380.197353, 4.15, 5.7),
                 (439.150807, 0.2, 0.91),
                 (448.001085, 1.63, 2.46),
                 (474.689092, 0.76, 2.22),
                 (488.490108, 0.26, 2.49),
                 (556.935985, 7.81, 10),
                 (620.70087, 1.25, 2.35),
                 (752.033113, 16.2, 20),
                 (916.171582, 1.47, 2.58),
                 (970.315022, 1.36, 2.44),
                 (987.926764, 1.6, 1.86)]

    def __init__(self):
        self.__version__ = 12
        self.year = 2019
        self.month = 8
        self.link = 'https://www.itu.int/rec/R-REC-P.676-11-201712-S/en'

    def gammaw_approx(self, f, p, rho, T):
        warnings.warn(
            RuntimeWarning(
                'Recommendation ITU-R P.676-12 does not have an explicit '
                'method to approximate gamma_w. The exact method shall be '
                'used instead.'))
        return self.gamma_exact(f, p, rho, T)

    def gamma0_approx(self, f, p, rho, T):
        warnings.warn(
            RuntimeWarning(
                'Recommendation ITU-R P.676-12 does not have an explicit '
                'method to approximate gamma_w. The exact method shall be '
                'used instead.'))
        return self.gamma_exact(f, p, rho, T)

    @classmethod
    def gamma0_exact(self, f, p, rho, T):
        return __gamma0_exact__(self, f, p, rho, T)

    @classmethod
    def gammaw_exact(self, f, p, rho, T):
        return __gammaw_exact__(self, f, p, rho, T)

    @classmethod
    def gamma_exact(self, f, p, rho, T):
        return (self.gamma0_exact(f, p, rho, T) +
                self.gammaw_exact(f, p, rho, T))

    @classmethod
    def gaseous_attenuation_approximation(self, f, el, rho, P, T):
        if np.any(f > 350):
            warnings.warn(
                RuntimeWarning(
                    'The approximated method to computes '
                    'the gaseous attenuation in recommendation ITU-P 676-11 '
                    'is only recommended for frequencies below 350GHz'))

        if np.any(5 > el) or np.any(np.mod(el, 90) < 5):
            warnings.warn(
                RuntimeWarning(
                    'The approximated method to compute '
                    'the gaseous attenuation in recommendation ITU-P 676-11 '
                    'is only recommended for elevation angles between '
                    '5 and 90 degrees'))

        # Water vapour attenuation (gammaw) computation as in Section 1 of
        # Annex 2 of [1]
        gamma0 = self.gamma0_exact(f, P, rho, T)
        gammaw = self.gammaw_exact(f, P, rho, T)

        return gamma0, gammaw

    @classmethod
    def slant_inclined_path_equivalent_height(self, f, P, rho, T):
        """
        """
        e = rho * T / 216.7
        rp = (P + e) / 1013.25

        # Eq. 31 - 34
        t1 = 5.1040 / (1 + 0.066 * rp**-2.3) * \
            np.exp(-((f - 59.7) / (2.87 + 12.4 * np.exp(-7.9 * rp)))**2)

        t2 = sum([(ci * np.exp(2.12 * rp)) /
                  ((f - fi)**2 + 0.025 * np.exp(2.2 * rp))
                  for ci, fi in self.t2_coeffs])

        t3 = 0.0114 * f / (1 + 0.14 * rp**-2.6) * \
            (15.02 * f**2 - 1353 * f + 5.333e4) / \
            (f**3 - 151.3 * f**2 + 9629 * f - 6803)

        A = 0.7832 + 0.00709 * (T - 273.15)

        # Eq. 30
        h0 = 6.1 * A / (1 + 0.17 * rp**-1.1) * (1 + t1 + t2 + t3)

        h0 = np.where(f < 70,
                      np.minimum(h0, 10.7 * rp**0.3),
                      h0)

        # Eq. 36 - 38
        A = 1.9298 - 0.04166 * (T - 273.15) + 0.0517 * rho
        B = 1.1674 - 0.00622 * (T - 273.15) + 0.0063 * rho
        sigmaw = 1.013 / (1 + np.exp(-8.6 * (rp - 0.57)))

        # Eq. 35 b
        hw = A + B * sum([(ai * sigmaw) / ((f - fi)**2 + bi * sigmaw)
                          for fi, ai, bi in self.hw_coeffs])
        return h0, hw

    @classmethod
    def gaseous_attenuation_terrestrial_path(
            self, r, f, el, rho, P, T, mode='approx'):
        """
        """
        if mode == 'approx':
            gamma0, gammaw = self.gaseous_attenuation_approximation(
                f, el, rho, P, T)
            return (gamma0 + gammaw) * r
        else:
            gamma = self.gamma_exact(f, P, rho, T)
            return gamma * r

    @classmethod
    def gaseous_attenuation_slant_path(self, f, el, rho, P, T, V_t=None,
                                       h=None, mode='approx'):
        """
        """
        if mode == 'approx':
            gamma0, gammaw = self.gaseous_attenuation_approximation(
                f, el, rho, P, T)

            h0, hw = self.slant_inclined_path_equivalent_height(f, P, rho, T)

            # Use the zenit water-vapour method if the values of V_t
            # and h are provided
            if V_t is not None and h is not None:
                Aw = self.zenit_water_vapour_attenuation(None, None, None,
                                                         f, V_t, h)
            else:
                Aw = gammaw * hw

            A0 = gamma0 * h0
            return (A0 + Aw) / np.sin(np.deg2rad(el))

        else:
            delta_h = 0.0001 * \
                np.exp((np.arange(0, 922)) / 100)             # Eq. 14
            h_n = 0.0001 * ((np.exp(np.arange(0, 922) / 100.0) -
                             1.0) / (np.exp(1.0 / 100.0) - 1.0))             # Eq. 15
            T_n = standard_temperature(h_n).to(u.K).value
            press_n = standard_pressure(h_n).value
            rho_n = standard_water_vapour_density(h_n, rho_0=rho).value

            e_n = rho_n * T_n / 216.7
            n_n = radio_refractive_index(press_n, e_n, T_n).value
            n_ratio = n_n / np.pad(n_n[1:], (0, 1), mode='edge')
            r_n = 6371 + h_n

            b = np.pi / 2 - np.deg2rad(el)
            Agas = 0
            for t, press, rho, r, delta, n_r in zip(
                    T_n, press_n, rho_n, r_n, delta_h, n_ratio):
                a = - r * np.cos(b) + 0.5 * np.sqrt(
                    4 * r**2 * np.cos(b)**2 + 8 * r * delta + 4 * delta**2)  # Eq. 17
                a_cos_arg = np.clip((-a**2 - 2 * r * delta - delta**2) /
                                    (2 * a * r + 2 * a * delta), -1, 1)
                # Eq. 18a
                alpha = np.pi - np.arccos(a_cos_arg)
                gamma = self.gamma_exact(f, press, rho, t)
                Agas += a * gamma                                            # Eq. 13
                b = np.arcsin(np.sin(alpha) *
                              n_r)                           # Eq. 19a

            return Agas

    @classmethod
    def gaseous_attenuation_inclined_path(
            self, f, el, rho, P, T, h1, h2, mode='approx'):
        """
        """
        if h1 > 10 or h2 > 10:
            raise ValueError(
                'Both the transmitter and the receiver must be at'
                'altitude of less than 10 km above the sea level.'
                'Current altitude Tx: %.2f km, Rx: %.2f km' % (h1, h2))

        if mode == 'approx':
            rho = rho * np.exp(h1 / 2)
            gamma0, gammaw = self.gaseous_attenuation_approximation(
                f, el, rho, P, T)
        else:
            gamma0 = self.gamma0_exact(f, P, rho, T)
            gammaw = 0

        e = rho * T / 216.7
        h0, hw = self.slant_inclined_path_equivalent_height(f, P + e, rho, T)

        if 5 < el and el < 90:
            h0_p = h0 * (np.exp(-h1 / h0) - np.exp(-h2 / h0))
            hw_p = hw * (np.exp(-h1 / hw) - np.exp(-h2 / hw))
            return (gamma0 * h0_p + gammaw * hw_p) / np.sin(np.deg2rad(el))
        else:
            def F(x):
                return 1 / (0.661 * x + 0.339 * np.sqrt(x**2 + 5.51))

            el1 = el
            Re = 8500  # TODO: change to ITU-R P 834
            el2 = np.rad2deg(
                np.arccos(((Re + h1) / (Re + h2)) * np.cos(np.deg2rad(el1))))

            def xi(eli, hi):
                return np.tan(np.deg2rad(eli)) * np.sqrt((Re + hi) / h0)

            def xi_p(eli, hi):
                return np.tan(np.deg2rad(eli)) * np.sqrt((Re + hi) / hw)

            def eq_33(h_num, h_den, el, x):
                return np.sqrt(Re + h_num) * F(x) * \
                    np.exp(-h_num / h_den) / np.cos(np.deg2rad(el))

            A = gamma0 * np.sqrt(h0) * (eq_33(h1, h0, el1, xi(el1, h1)) -
                                        eq_33(h2, h0, el2, xi(el2, h2))) +\
                gammaw * np.sqrt(hw) * (eq_33(h1, hw, el1, xi_p(el1, h1)) -
                                        eq_33(h2, hw, el2, xi_p(el2, h2)))
            return A

    @classmethod
    def zenit_water_vapour_attenuation(
            self, lat, lon, p, f, V_t=None, h=None):
        f_ref = 20.6        # [GHz]
        p_ref = 845         # [hPa]

        if h is None:
            h = topographic_altitude(lat, lon).value

        if V_t is None:
            V_t = total_water_vapour_content(lat, lon, p, h).value

        rho_ref = V_t / 2.38
        t_ref = 14 * np.log(0.22 * V_t / 2.38) + 3    # [Celsius]

        a = (0.2048 * np.exp(- ((f - 22.43) / 3.097)**2) +
             0.2326 * np.exp(- ((f - 183.5) / 4.096)**2) +
             0.2073 * np.exp(- ((f - 325) / 3.651)**2) - 0.1113)

        b = 8.741e4 * np.exp(-0.587 * f) + 312.2 * f**(-2.38) + 0.723
        h = np.clip(h, 0, 4)

        gammaw_approx_vect = np.vectorize(self.gammaw_exact)

        Aw_term1 = (0.0176 * V_t *
                    gammaw_approx_vect(f, p_ref, rho_ref, t_ref + 273.15) /
                    gammaw_approx_vect(f_ref, p_ref, rho_ref, t_ref + 273.15))

        return np.where(f < 20, Aw_term1, Aw_term1 * (a * h ** b + 1))


class _ITU676_11_():

    tmp = load_data(os.path.join(dataset_dir, '676/v11_lines_oxygen.txt'),
                    skip_header=1)
    f_ox = tmp[:, 0]
    a1 = tmp[:, 1]
    a2 = tmp[:, 2]
    a3 = tmp[:, 3]
    a4 = tmp[:, 4]
    a5 = tmp[:, 5]
    a6 = tmp[:, 6]

    tmp = load_data(os.path.join(dataset_dir,
                                 '676//v11_lines_water_vapour.txt'),
                    skip_header=1)
    f_wv = tmp[:, 0]
    b1 = tmp[:, 1]
    b2 = tmp[:, 2]
    b3 = tmp[:, 3]
    b4 = tmp[:, 4]
    b5 = tmp[:, 5]
    b6 = tmp[:, 6]

    idx_approx = np.zeros_like(b1, dtype=bool).squeeze()
    asterisk_rows = [0, 3, 4, 5, 7, 12, 20, 24, 34]
    idx_approx[np.array(asterisk_rows)] = True

    def __init__(self):
        self.__version__ = 11
        self.year = 2017
        self.month = 12
        self.link = 'https://www.itu.int/rec/R-REC-P.676-11-201712-S/en'

    @classmethod
    def gammaw_approx(self, f, p, rho, T):
        # T in Kelvin
        # e : water vapour partial pressure in hPa (total barometric pressure
        # ptot = p + e)
        theta = 300 / T
        e = rho * T / 216.7

        f_wv = self.f_wv[self.idx_approx]
        b1 = self.b1[self.idx_approx]
        b2 = self.b2[self.idx_approx]
        b3 = self.b3[self.idx_approx]
        b4 = self.b4[self.idx_approx]
        b5 = self.b5[self.idx_approx]
        b6 = self.b6[self.idx_approx]

        D_f_wv = b3 * 1e-4 * (p * theta ** b4 +
                              b5 * e * theta ** b6)

        F_i_wv = f / f_wv * ((D_f_wv) / ((f_wv - f)**2 + D_f_wv**2) +
                             (D_f_wv) / ((f_wv + f)**2 + D_f_wv**2))

        Si_wv = b1 * 1e-1 * e * theta**3.5 * np.exp(b2 * (1 - theta))

        N_pp_wv = Si_wv * F_i_wv

        N_pp = N_pp_wv.sum()

        gamma = 0.1820 * f * N_pp           # Eq. 1 [dB/km]
        return gamma

    @classmethod
    def gamma0_approx(self, f, p, rho, T):
        # T in Kelvin
        # e : water vapour partial pressure in hPa (total barometric pressure
        # ptot = p + e)
        theta = 300 / T
        e = rho * T / 216.7

        f_ox = self.f_ox

        D_f_ox = self.a3 * 1e-4 * (p * (theta ** (0.8 - self.a4)) +
                                   1.1 * e * theta)

        delta_ox = (self.a5 + self.a6 * theta) * 1e-4 * (p + e) * theta**0.8

        F_i_ox = f / f_ox * ((D_f_ox - delta_ox * (f_ox - f)) /
                             ((f_ox - f) ** 2 + D_f_ox ** 2) +
                             (D_f_ox - delta_ox * (f_ox + f)) /
                             ((f_ox + f) ** 2 + D_f_ox ** 2))

        Si_ox = self.a1 * 1e-7 * p * theta**3 * np.exp(self.a2 * (1 - theta))

        N_pp_ox = Si_ox * F_i_ox

        d = 5.6e-4 * (p + e) * theta**0.8

        N_d_pp = f * p * theta**2 * \
            (6.14e-5 / (d * (1 + (f / d)**2)) +
             1.4e-12 * p * theta**1.5 / (1 + 1.9e-5 * f**1.5))

        N_pp = N_pp_ox.sum() + N_d_pp

        gamma = 0.1820 * f * N_pp           # Eq. 1 [dB/km]
        return gamma

    @classmethod
    def gamma0_exact(self, f, p, rho, T):
        return __gamma0_exact__(self, f, p, rho, T)

    @classmethod
    def gammaw_exact(self, f, p, rho, T):
        return __gammaw_exact__(self, f, p, rho, T)

    @classmethod
    def gamma_exact(self, f, p, rho, T):
        return (self.gamma0_exact(f, p, rho, T) +
                self.gammaw_exact(f, p, rho, T))

    @classmethod
    def gaseous_attenuation_approximation(self, f, el, rho, P, T):
        """
        T goes in Kelvin
        """
        if np.any(f > 350):
            warnings.warn(
                RuntimeWarning(
                    'The approximated method to computes '
                    'the gaseous attenuation in recommendation ITU-P 676-11 '
                    'is only recommended for frequencies below 350GHz'))

        if np.any(5 > el) or np.any(np.mod(el, 90) < 5):
            warnings.warn(
                RuntimeWarning(
                    'The approximated method to compute '
                    'the gaseous attenuation in recommendation ITU-P 676-11 '
                    'is only recommended for elevation angles between '
                    '5 and 90 degrees'))

        # Water vapour attenuation (gammaw) computation as in Section 1 of
        # Annex 2 of [1]
        gamma0 = self.gamma0_approx(f, P, rho, T)
        gammaw = self.gammaw_approx(f, P, rho, T)

        return gamma0, gammaw

    @classmethod
    def slant_inclined_path_equivalent_height(self, f, p, rho=None, T=None):
        """
        """
        rp = p / 1013.25
        t1 = 4.64 / (1 + 0.066 * rp**-2.3) * \
            np.exp(- ((f - 59.7) / (2.87 + 12.4 * np.exp(-7.9 * rp)))**2)
        t2 = (0.14 * np.exp(2.12 * rp)) / \
            ((f - 118.75)**2 + 0.031 * np.exp(2.2 * rp))
        t3 = 0.0114 / (1 + 0.14 * rp**-2.6) * f * \
            (-0.0247 + 0.0001 * f + 1.61e-6 * f**2) / \
            (1 - 0.0169 * f + 4.1e-5 * f**2 + 3.2e-7 * f**3)

        h0 = 6.1 / (1 + 0.17 * rp**-1.1) * (1 + t1 + t2 + t3)

        h0 = np.where(f < 70,
                      np.minimum(h0, 10.7 * rp**0.3),
                      h0)

        sigmaw = 1.013 / (1 + np.exp(-8.6 * (rp - 0.57)))
        hw = 1.66 * (1 + (1.39 * sigmaw) / ((f - 22.235)**2 + 2.56 * sigmaw) +
                     (3.37 * sigmaw) / ((f - 183.31)**2 + 4.69 * sigmaw) +
                     (1.58 * sigmaw) / ((f - 325.1)**2 + 2.89 * sigmaw))

        return h0, hw

    @classmethod
    def gaseous_attenuation_terrestrial_path(
            self, r, f, el, rho, P, T, mode='approx'):
        """
        """
        if mode == 'approx':
            gamma0, gammaw = self.gaseous_attenuation_approximation(
                f, el, rho, P, T)
            return (gamma0 + gammaw) * r
        else:
            gamma = self.gamma_exact(f, P, rho, T)
            return gamma * r

    @classmethod
    def gaseous_attenuation_slant_path(self, f, el, rho, P, T, V_t=None,
                                       h=None, mode='approx'):
        """
        """
        if mode == 'approx':
            gamma0, gammaw = self.gaseous_attenuation_approximation(
                f, el, rho, P, T)

            e = rho * T / 216.7
            h0, hw = self.slant_inclined_path_equivalent_height(f, P + e)

            # Use the zenit water-vapour method if the values of V_t
            # and h are provided
            if V_t is not None and h is not None:
                Aw = self.zenit_water_vapour_attenuation(None, None, None,
                                                         f, V_t, h)
            else:
                Aw = gammaw * hw

            A0 = gamma0 * h0
            return (A0 + Aw) / np.sin(np.deg2rad(el))

        else:
            delta_h = 0.0001 * \
                np.exp((np.arange(0, 922)) / 100)             # Eq. 21
            h_n = 0.0001 * ((np.exp(np.arange(0, 922) / 100.0) -
                             1.0) / (np.exp(1.0 / 100.0) - 1.0))
            T_n = standard_temperature(h_n).to(u.K).value
            press_n = standard_pressure(h_n).value
            rho_n = standard_water_vapour_density(h_n, rho_0=rho).value

            e_n = rho_n * T_n / 216.7
            n_n = radio_refractive_index(press_n, e_n, T_n).value
            n_ratio = n_n / np.pad(n_n[1:], (0, 1), mode='edge')
            r_n = 6371 + h_n

            b = np.pi / 2 - np.deg2rad(el)
            Agas = 0
            for t, press, rho, r, delta, n_r in zip(
                    T_n, press_n, rho_n, r_n, delta_h, n_ratio):
                a = - r * np.cos(b) + 0.5 * np.sqrt(
                    4 * r**2 * np.cos(b)**2 + 8 * r * delta + 4 * delta**2)  # Eq. 17
                a_cos_arg = np.clip((-a**2 - 2 * r * delta - delta**2) /
                                    (2 * a * r + 2 * a * delta), -1, 1)
                # Eq. 18
                alpha = np.pi - np.arccos(a_cos_arg)
                gamma = self.gamma_exact(f, press, rho, t)
                Agas += a * gamma                                            # Eq. 20
                b = np.arcsin(np.sin(alpha) *
                              n_r)                           # Eq. 19

            return Agas

    @classmethod
    def gaseous_attenuation_inclined_path(
            self, f, el, rho, P, T, h1, h2, mode='approx'):
        """
        """
        if h1 > 10 or h2 > 10:
            raise ValueError(
                'Both the transmitter and the receiver must be at'
                'altitude of less than 10 km above the sea level.'
                'Current altitude Tx: %.2f km, Rx: %.2f km' % (h1, h2))

        if mode == 'approx':
            rho = rho * np.exp(h1 / 2)
            gamma0, gammaw = self.gaseous_attenuation_approximation(
                f, el, rho, P, T)
        else:
            gamma0 = self.gamma0_exact(f, P, rho, T)
            gammaw = 0

        e = rho * T / 216.7
        h0, hw = self.slant_inclined_path_equivalent_height(f, P + e)

        if 5 < el and el < 90:
            h0_p = h0 * (np.exp(-h1 / h0) - np.exp(-h2 / h0))
            hw_p = hw * (np.exp(-h1 / hw) - np.exp(-h2 / hw))
            return (gamma0 * h0_p + gammaw * hw_p) / np.sin(np.deg2rad(el))
        else:
            def F(x):
                return 1 / (0.661 * x + 0.339 * np.sqrt(x**2 + 5.51))

            el1 = el
            Re = 8500  # TODO: change to ITU-R P 834
            el2 = np.rad2deg(
                np.arccos(((Re + h1) / (Re + h2)) * np.cos(np.deg2rad(el1))))

            def xi(eli, hi):
                return np.tan(np.deg2rad(eli)) * np.sqrt((Re + hi) / h0)

            def xi_p(eli, hi):
                return np.tan(np.deg2rad(eli)) * np.sqrt((Re + hi) / hw)

            def eq_33(h_num, h_den, el, x):
                return np.sqrt(Re + h_num) * F(x) * \
                    np.exp(-h_num / h_den) / np.cos(np.deg2rad(el))

            A = gamma0 * np.sqrt(h0) * (eq_33(h1, h0, el1, xi(el1, h1)) -
                                        eq_33(h2, h0, el2, xi(el2, h2))) +\
                gammaw * np.sqrt(hw) * (eq_33(h1, hw, el1, xi_p(el1, h1)) -
                                        eq_33(h2, hw, el2, xi_p(el2, h2)))
            return A

    @classmethod
    def zenit_water_vapour_attenuation(
            self, lat, lon, p, f, V_t=None, h=None):
        f_ref = 20.6        # [GHz]
        p_ref = 815         # [hPa]

        if h is None:
            h = topographic_altitude(lat, lon).value

        if V_t is None:
            V_t = total_water_vapour_content(lat, lon, p, h).value

        rho_ref = V_t / 3.67
        t_ref = 14 * np.log(0.22 * V_t / 3.67) + 3    # [Celsius]

        a = (0.2048 * np.exp(- ((f - 22.43) / 3.097)**2) +
             0.2326 * np.exp(- ((f - 183.5) / 4.096)**2) +
             0.2073 * np.exp(- ((f - 325) / 3.651)**2) - 0.113)

        b = 8.741e4 * np.exp(-0.587 * f) + 312.2 * f**(-2.38) + 0.723
        h = np.minimum(h, 4)

        gammaw_approx_vect = np.vectorize(self.gammaw_approx)

        Aw_term1 = (0.0176 * V_t *
                    gammaw_approx_vect(f, p_ref, rho_ref, t_ref + 273.15) /
                    gammaw_approx_vect(f_ref, p_ref, rho_ref, t_ref + 273.15))

        return np.where(f < 20, Aw_term1, Aw_term1 * (a * h ** b + 1))


class _ITU676_10_():

    tmp = load_data(os.path.join(dataset_dir, '676/v10_lines_oxygen.txt'),
                    skip_header=1)
    f_ox = tmp[:, 0]
    a1 = tmp[:, 1]
    a2 = tmp[:, 2]
    a3 = tmp[:, 3]
    a4 = tmp[:, 4]
    a5 = tmp[:, 5]
    a6 = tmp[:, 6]

    tmp = load_data(os.path.join(dataset_dir,
                                 '676//v10_lines_water_vapour.txt'),
                    skip_header=1)
    f_wv = tmp[:, 0]
    b1 = tmp[:, 1]
    b2 = tmp[:, 2]
    b3 = tmp[:, 3]
    b4 = tmp[:, 4]
    b5 = tmp[:, 5]
    b6 = tmp[:, 6]

    def __init__(self):
        self.__version__ = 10
        self.year = 2013
        self.month = 9
        self.link = 'https://www.itu.int/rec/R-REC-P.676-10-201309-S/en'

    @classmethod
    def gammaw_approx(self, f, P, rho, T):
        rp = P / 1013
        rt = 288 / (T)
        eta1 = 0.955 * rp * rt**0.68 + 0.006 * rho
        eta2 = 0.735 * rp * rt**0.50 + 0.0353 * rt**4 * rho

        def g(f, fi): return 1 + ((f - fi) / (f + fi))**2
        gammaw = (
            (3.98 * eta1 * np.exp(2.23 * (1 - rt))) /
            ((f - 22.235) ** 2 + 9.42 * eta1 ** 2) * g(f, 22.0) +
            (11.96 * eta1 * np.exp(0.70 * (1 - rt))) /
            ((f - 183.310) ** 2 + 11.14 * eta1 ** 2) +
            (0.081 * eta1 * np.exp(6.44 * (1 - rt))) /
            ((f - 321.226) ** 2 + 6.29 * eta1 ** 2) +
            (3.660 * eta1 * np.exp(1.60 * (1 - rt))) /
            ((f - 325.153) ** 2 + 9.22 * eta1 ** 2) +
            (25.37 * eta1 * np.exp(1.09 * (1 - rt))) / ((f - 380.000) ** 2) +
            (17.40 * eta1 * np.exp(1.46 * (1 - rt))) / ((f - 448.000) ** 2) +
            (844.6 * eta1 * np.exp(0.17 * (1 - rt))) / ((f - 557.000) ** 2) *
            g(f, 557.0) + (290.0 * eta1 * np.exp(0.41 * (1 - rt))) /
            ((f - 752.000) ** 2) * g(f, 752.0) +
            (8.3328e4 * eta2 * np.exp(0.99 * (1 - rt))) /
            ((f - 1780.00) ** 2) *
            g(f, 1780.0)) * f ** 2 * rt ** 2.5 * rho * 1e-4
        return gammaw

    @classmethod
    def gamma0_approx(self, f, P, rho, T):
        rp = P / 1013.0
        rt = 288.0 / (T)

        def phi(rp, rt, a, b, c, d): return (
            rp**a * np.power(rt, b) * np.exp(c * (1 - rp) + d * (1 - rt)))

        # Dry air attenuation (gamma0) computation as in Section 1 of Annex 2
        # of [1]
        delta = -0.00306 * phi(rp, rt, 3.211, -14.94, 1.583, -16.37)

        xi1 = phi(rp, rt, 0.0717, -1.8132, 0.0156, -1.6515)
        xi2 = phi(rp, rt, 0.5146, -4.6368, -0.1921, -5.7416)
        xi3 = phi(rp, rt, 0.3414, -6.5851, 0.2130, -8.5854)
        xi4 = phi(rp, rt, -0.0112, 0.0092, -0.1033, -0.0009)
        xi5 = phi(rp, rt, 0.2705, -2.7192, -0.3016, -4.1033)
        xi6 = phi(rp, rt, 0.2445, -5.9191, 0.0422, -8.0719)
        xi7 = phi(rp, rt, -0.1833, 6.5589, -0.2402, 6.131)

        gamma54 = 2.192 * phi(rp, rt, 1.8286, -1.9487, 0.4051, -2.8509)
        gamma58 = 12.59 * phi(rp, rt, 1.0045, 3.5610, 0.1588, 1.2834)
        gamma60 = 15.00 * phi(rp, rt, 0.9003, 4.1335, 0.0427, 1.6088)
        gamma62 = 14.28 * phi(rp, rt, 0.9886, 3.4176, 0.1827, 1.3429)
        gamma64 = 6.819 * phi(rp, rt, 1.4320, 0.6258, 0.3177, -0.5914)
        gamma66 = 1.908 * phi(rp, rt, 2.0717, -4.1404, 0.4910, -4.8718)

        def fcn_le_54():
            return (((7.2 * rt**2.8) / (f**2 + 0.34 * rp**2 * rt**1.6) +
                     (0.62 * xi3) / ((54 - f)**(1.16 * xi1) + 0.83 * xi2)) *
                    f**2 * rp**2 * 1e-3)

        def fcn_le_60():
            return (np.exp(np.log(gamma54) / 24.0 * (f - 58) * (f - 60) -
                           np.log(gamma58) / 8.0 * (f - 54) * (f - 60) +
                           np.log(gamma60) / 12.0 * (f - 54) * (f - 58)))

        def fcn_le_62():
            return (gamma60 + (gamma62 - gamma60) * (f - 60) / 2.0)

        def fcn_le_66():
            return (np.exp(np.log(gamma62) / 8.0 * (f - 64) * (f - 66) -
                           np.log(gamma64) / 4.0 * (f - 62) * (f - 66) +
                           np.log(gamma66) / 8.0 * (f - 62) * (f - 64)))

        def fcn_le_120():
            return ((3.02e-4 * rt**3.5 + (0.283 * rt**3.8) /
                     ((f - 118.75)**2 + 2.91 * rp**2 * rt**1.6) +
                     (0.502 * xi6 * (1 - 0.0163 * xi7 * (f - 66))) /
                     ((f - 66)**(1.4346 * xi4) + 1.15 * xi5)) *
                    f**2 * rp**2 * 1e-3)

        def fcn_rest():
            return (((3.02e-4) / (1 + 1.9e-5 * f**1.5) +
                     (0.283 * rt**0.3) / ((f - 118.75)**2 +
                                          2.91 * rp**2 * rt**1.6)) *
                    f**2 * rp**2 * rt**3.5 * 1e-3 + delta)

        with np.errstate(invalid='ignore'):
            gamma0 = \
                np.where(
                    f <= 54, fcn_le_54(),
                    np.where(
                        np.logical_and(54 < f, f <= 60), fcn_le_60(),
                        np.where(
                            np.logical_and(60 < f, f <= 62), fcn_le_62(),
                            np.where(
                                np.logical_and(62 < f, f <= 66), fcn_le_66(),
                                np.where(
                                    np.logical_and(66 < f, f <= 120),
                                    fcn_le_120(),
                                    fcn_rest())))))

        return gamma0

    @classmethod
    def gamma0_exact(self, f, p, rho, T):
        return __gamma0_exact__(self, f, p, rho, T)

    @classmethod
    def gammaw_exact(self, f, p, rho, T):
        return __gammaw_exact__(self, f, p, rho, T)

    @classmethod
    def gamma_exact(self, f, p, rho, T):
        return (self.gamma0_exact(f, p, rho, T) +
                self.gammaw_exact(f, p, rho, T))

    @classmethod
    def gaseous_attenuation_approximation(self, f, el, rho, P, T):
        """
        T goes in Kelvin
        """
        if np.any(f > 350):
            warnings.warn(
                RuntimeWarning(
                    'The approximated method to computes '
                    'the gaseous attenuation in recommendation ITU-P 676-11 '
                    'is only recommended for frequencies below 350GHz'))

        if np.any(5 > el) or np.any(np.mod(el, 90) < 5):
            warnings.warn(
                RuntimeWarning(
                    'The approximated method to compute '
                    'the gaseous attenuation in recommendation ITU-P 676-11 '
                    'is only recommended for elevation angles between '
                    '5 and 90 degrees'))

        # Water vapour attenuation (gammaw) computation as in Section 1 of
        # Annex 2 of [1]
        gamma0 = self.gamma0_approx(f, P, rho, T)
        gammaw = self.gammaw_approx(f, P, rho, T)

        return gamma0, gammaw

    @classmethod
    def slant_inclined_path_equivalent_height(self, f, p, rho=None, T=None):
        """
        """
        rp = p / 1013.0
        t1 = (4.64) / (1 + 0.066 * rp**-2.3) * \
            np.exp(- ((f - 59.7) / (2.87 + 12.4 * np.exp(-7.9 * rp)))**2)
        t2 = (0.14 * np.exp(2.21 * rp)) / \
            ((f - 118.75)**2 + 0.031 * np.exp(2.2 * rp))
        t3 = (0.0114) / (1 + 0.14 * rp**-2.6) * f * \
             (-0.0247 + 0.0001 * f + 1.61e-6 * f**2) / \
             (1 - 0.0169 * f + 4.1e-5 * f**2 + 3.2e-7 * f**3)

        h0 = (6.1) / (1 + 0.17 * rp**-1.1) * (1 + t1 + t2 + t3)

        h0 = np.where(f < 70,
                      np.minimum(h0, 10.7 * rp**0.3),
                      h0)

        sigmaw = (1.013) / (1 + np.exp(-8.6 * (rp - 0.57)))
        hw = 1.66 * (1 + (1.39 * sigmaw) / ((f - 22.235)**2 + 2.56 * sigmaw) +
                     (3.37 * sigmaw) / ((f - 183.31)**2 + 4.69 * sigmaw) +
                     (1.58 * sigmaw) / ((f - 325.1)**2 + 2.89 * sigmaw))

        return h0, hw

    @classmethod
    def gaseous_attenuation_terrestrial_path(
            self, r, f, el, rho, P, T, mode='approx'):
        """
        """
        if mode == 'approx':
            gamma0, gammaw = self.gaseous_attenuation_approximation(
                f, el, rho, P, T)
            return (gamma0 + gammaw) * r
        else:
            gamma = self.gamma_exact(f, P, rho, T)
            return gamma * r

    @classmethod
    def gaseous_attenuation_slant_path(self, f, el, rho, P, T, V_t=None,
                                       h=None, mode='approx'):
        """
        """
        if mode == 'approx':
            gamma0, gammaw = self.gaseous_attenuation_approximation(
                f, el, rho, P, T)
            e = rho * T / 216.7
            h0, hw = self.slant_inclined_path_equivalent_height(f, P + e)

            # Use the zenit water-vapour method if the values of V_t
            # and h are provided
            if V_t is not None and h is not None:
                Aw = self.zenit_water_vapour_attenuation(None, None, None,
                                                         f, V_t, h)
            else:
                Aw = gammaw * hw

            A0 = gamma0 * h0
            return (A0 + Aw) / np.sin(np.deg2rad(el))

        else:
            delta_h = 0.0001 * \
                np.exp((np.arange(0, 922)) / 100)             # Eq. 21
            h_n = 0.0001 * ((np.exp(np.arange(0, 922) / 100.0) -
                             1.0) / (np.exp(1.0 / 100.0) - 1.0))
            T_n = standard_temperature(h_n).to(u.K).value
            press_n = standard_pressure(h_n).value
            rho_n = standard_water_vapour_density(h_n, rho_0=rho).value

            e_n = rho_n * T_n / 216.7
            n_n = radio_refractive_index(press_n, e_n, T_n).value
            n_ratio = n_n / np.pad(n_n[1:], (0, 1), mode='edge')
            r_n = 6371 + h_n

            b = np.pi / 2 - np.deg2rad(el)
            Agas = 0
            for t, press, rho, r, delta, n_r in zip(
                    T_n, press_n, rho_n, r_n, delta_h, n_ratio):
                a = - r * np.cos(b) + 0.5 * np.sqrt(
                    4 * r**2 * np.cos(b)**2 + 8 * r * delta + 4 * delta**2)  # Eq. 17
                a_cos_arg = np.clip((-a**2 - 2 * r * delta - delta**2) /
                                    (2 * a * r + 2 * a * delta), -1, 1)
                # Eq. 18
                alpha = np.pi - np.arccos(a_cos_arg)
                gamma = self.gamma_exact(f, press, rho, t)
                Agas += a * gamma                                            # Eq. 20
                b = np.arcsin(np.sin(alpha) *
                              n_r)                           # Eq. 19

            return Agas

    @classmethod
    def gaseous_attenuation_inclined_path(
            self, f, el, rho, P, T, h1, h2, mode='approx'):
        """
        """
        if h1 > 10 or h2 > 10:
            raise ValueError(
                'Both the transmitter and the receiver must be at'
                'altitude of less than 10 km above the sea level.'
                'Current altitude Tx: %.2f km, Rx: %.2f km' % (h1, h2))

        if mode == 'approx':
            rho = rho * np.exp(h1 / 2)
            gamma0, gammaw = self.gaseous_attenuation_approximation(
                f, el, rho, P, T)
        else:
            gamma0 = self.gamma_exact(f, P, rho, T)
            gammaw = 0

        e = rho * T / 216.7
        h0, hw = self.slant_inclined_path_equivalent_height(f, P + e)

        if 5 < el and el < 90:
            h0_p = h0 * (np.exp(-h1 / h0) - np.exp(-h2 / h0))
            hw_p = hw * (np.exp(-h1 / hw) - np.exp(-h2 / hw))
            return (gamma0 * h0_p + gammaw * hw_p) / np.sin(np.deg2rad(el))
        else:
            def F(x):
                return 1 / (0.661 * x + 0.339 * np.sqrt(x**2 + 5.51))

            el1 = el
            Re = 8500  # TODO: change to ITU-R P 834
            el2 = np.rad2deg(
                np.arccos(((Re + h1) / (Re + h2)) * np.cos(np.deg2rad(el1))))

            def xi(eli, hi):
                return np.tan(np.deg2rad(eli)) * np.sqrt((Re + hi) / h0)

            def xi_p(eli, hi):
                return np.tan(np.deg2rad(eli)) * np.sqrt((Re + hi) / hw)

            def eq_33(h_num, h_den, el, x):
                return np.sqrt(Re + h_num) * F(x) * \
                    np.exp(-h_num / h_den) / np.cos(np.deg2rad(el))

            A = gamma0 * np.sqrt(h0) * (eq_33(h1, h0, el1, xi(el1, h1)) -
                                        eq_33(h2, h0, el2, xi(el2, h2))) +\
                gammaw * np.sqrt(hw) * (eq_33(h1, hw, el1, xi_p(el1, h1)) -
                                        eq_33(h2, hw, el2, xi_p(el2, h2)))
            return A

    @classmethod
    def zenit_water_vapour_attenuation(
            self, lat, lon, p, f, V_t=None, h=None):
        f_ref = 20.6        # [GHz]
        p_ref = 780         # [hPa]
        if V_t is None:
            V_t = total_water_vapour_content(lat, lon, p, h).value
        rho_ref = V_t / 4     # [g/m3]
        t_ref = 14 * np.log(0.22 * V_t / 4) + 3    # [Celsius]

        gammaw_approx_vect = np.vectorize(self.gammaw_approx)
        return (0.0173 * V_t *
                gammaw_approx_vect(f, p_ref, rho_ref, t_ref + 273) /
                gammaw_approx_vect(f_ref, p_ref, rho_ref, t_ref + 273))


class _ITU676_9_():

    tmp = load_data(os.path.join(dataset_dir, '676//v9_lines_oxygen.txt'),
                    skip_header=1)
    f_ox = tmp[:, 0]
    a1 = tmp[:, 1]
    a2 = tmp[:, 2]
    a3 = tmp[:, 3]
    a4 = tmp[:, 4]
    a5 = tmp[:, 5]
    a6 = tmp[:, 6]

    tmp = load_data(os.path.join(dataset_dir,
                                 '676//v9_lines_water_vapour.txt'),
                    skip_header=1)
    f_wv = tmp[:, 0]
    b1 = tmp[:, 1]
    b2 = tmp[:, 2]
    b3 = tmp[:, 3]
    b4 = tmp[:, 4]
    b5 = tmp[:, 5]
    b6 = tmp[:, 6]

    def __init__(self):
        self.__version__ = 9
        self.year = 2012
        self.month = 2
        self.link = 'https://www.itu.int/rec/R-REC-P.676-9-201202-S/en'

    # Recommendation ITU-P R.676-9 has most of the methods similar to those
    # in Recommendation ITU-P R.676-10.
    @staticmethod
    def gammaw_approx(*args, **kwargs):
        return _ITU676_10_.gammaw_approx(*args, **kwargs)

    @staticmethod
    def gamma0_approx(*args, **kwargs):
        return _ITU676_10_.gamma0_approx(*args, **kwargs)

    @staticmethod
    def gaseous_attenuation_inclined_path(*args, **kwargs):
        return _ITU676_10_.gaseous_attenuation_inclined_path(*args, **kwargs)

    @staticmethod
    def zenit_water_vapour_attenuation(*args, **kwargs):
        return _ITU676_10_.zenit_water_vapour_attenuation(*args, **kwargs)

    @staticmethod
    def gaseous_attenuation_approximation(*args, **kwargs):
        return _ITU676_10_.gaseous_attenuation_approximation(*args, **kwargs)

    @staticmethod
    def slant_inclined_path_equivalent_height(*args, **kwargs):
        return _ITU676_10_.slant_inclined_path_equivalent_height(*args,
                                                                 **kwargs)

    @staticmethod
    def gaseous_attenuation_terrestrial_path(*args, **kwargs):
        return _ITU676_10_.gaseous_attenuation_terrestrial_path(*args,
                                                                **kwargs)

    @staticmethod
    def gaseous_attenuation_slant_path(*args, **kwargs):
        return _ITU676_10_.gaseous_attenuation_slant_path(*args, **kwargs)

    def gamma0_exact(self, f, p, rho, T):
        return __gamma0_exact__(self, f, p, rho, T)

    def gammaw_exact(self, f, p, rho, T):
        return __gammaw_exact__(self, f, p, rho, T)

    def gamma_exact(self, f, p, rho, T):
        return (self.gamma0_exact(f, p, rho, T) +
                self.gammaw_exact(f, p, rho, T))


__model = __ITU676__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.676 recommendation currently being used.

    This function changes the model used for the ITU-R P.676 recommendation
    to a different version.


    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          *  12: Activates recommendation ITU-R P.676-12 (08/19) (Current version)
          *  11: Activates recommendation ITU-R P.676-11 (09/16) (Superseded)
          *  10: Activates recommendation ITU-R P.676-10 (09/13) (Superseded)
          *  9: Activates recommendation ITU-R P.676-9 (02/12) (Superseded)

    """
    global __model
    __model = __ITU676__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.676 recommendation currently being used.

    Returns
    -------
    version: int
       The version of the ITU-R P.530 recommendation being used.
    """
    return __model.__version__


def gaseous_attenuation_terrestrial_path(r, f, el, rho, P, T, mode):
    """
    Estimate the attenuation of atmospheric gases on terrestrial paths.
    This function operates in two modes, 'approx', and 'exact':
      * 'approx': a simplified approximate method to estimate gaseous attenuation
        that is applicable in the frequency range 1-350 GHz.
      * 'exact': an estimate of gaseous attenuation computed by summation of
        individual absorption lines that is valid for the frequency
        range 1-1,000 GHz


    Parameters
    ----------
    r : number or Quantity
        Path length (km)
    f : number or Quantity
        Frequency (GHz)
    el : sequence, number or Quantity
        Elevation angle (degrees)
    rho : number or Quantity
        Water vapor density (g/m**3)
    P : number or Quantity
        Atmospheric pressure (hPa)
    T : number or Quantity
        Absolute temperature (K)
    mode : string, optional
        Mode for the calculation. Valid values are 'approx', 'exact'. If
        'approx' Uses the method in Annex 2 of the recommendation (if any),
        else uses the method described in Section 1. Default, 'approx'


    Returns
    -------
    attenuation: Quantity
        Terrestrial path attenuation (dB)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en
    """
    type_output = get_input_type(el)
    r = prepare_quantity(r, u.km, 'Path Length')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(prepare_input_array(el), u.deg, 'Elevation angle')
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapor density')
    P = prepare_quantity(P, u.hPa, 'Atospheric pressure')
    T = prepare_quantity(T, u.K, 'Temperature')
    val = __model.gaseous_attenuation_terrestrial_path(
        r, f, el, rho, P, T, mode)
    return prepare_output_array(val, type_output) * u.dB


def gaseous_attenuation_slant_path(f, el, rho, P, T, V_t=None, h=None,
                                   mode='approx'):
    """
    Estimate the attenuation of atmospheric gases on slant paths. This function
    operates in two modes, 'approx', and 'exact':
      * 'approx': a simplified approximate method to estimate gaseous attenuation
        that is applicable in the frequency range 1-350 GHz.
      * 'exact': an estimate of gaseous attenuation computed by summation of
        individual absorption lines that is valid for the frequency
        range 1-1,000 GHz


    Parameters
    ----------
    f : number or Quantity
        Frequency (GHz)
    el : sequence, number or Quantity
        Elevation angle (degrees)
    rho : number or Quantity
        Water vapor density (g/m3)
    P : number or Quantity
        Atmospheric pressure (hPa)
    T : number or Quantity
        Absolute temperature (K)
    V_t: number or Quantity (kg/m2)
        Integrated water vapour content from: a) local radiosonde or
        radiometric data or b) at the required percentage of time (kg/m2)
        obtained from the digital maps in Recommendation ITU-R P.836 (kg/m2).
        If None, use general method to compute the wet-component of the
        gaseous attenuation. If provided, 'h' must be also provided. Default
        is None.
    h : number, sequence, or numpy.ndarray
        Altitude of the receivers. If None, use the topographical altitude as
        described in recommendation ITU-R P.1511. If provided, 'V_t' needs to
        be also provided. Default is None.
    mode : string, optional
        Mode for the calculation. Valid values are 'approx', 'exact'. If
        'approx' Uses the method in Annex 2 of the recommendation (if any),
        else uses the method described in Section 1. Default, 'approx'


    Returns
    -------
    attenuation: Quantity
        Slant path attenuation (dB)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en
    """
    type_output = get_input_type(el)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(prepare_input_array(el), u.deg, 'Elevation angle')
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapor density')
    P = prepare_quantity(P, u.hPa, 'Atospheric pressure')
    T = prepare_quantity(T, u.K, 'Temperature')
    V_t = prepare_quantity(V_t, u.kg / u.m**2,
                           'Integrated water vapour content')
    h = prepare_quantity(h, u.km, 'Altitude')
    val = __model.gaseous_attenuation_slant_path(
        f, el, rho, P, T, V_t, h, mode)
    
    # The values of attenuation cannot be negative. The ITU models end up
    # giving out negative values for certain inputs
    val[val < 0] = 0

    return prepare_output_array(val, type_output) * u.dB


def gaseous_attenuation_inclined_path(f, el, rho, P, T, h1, h2, mode='approx'):
    """
    Estimate the attenuation of atmospheric gases on inclined paths between two
    ground stations at heights h1 and h2. This function operates in two modes,
    'approx', and 'exact':
      * 'approx': a simplified approximate method to estimate gaseous attenuation
        that is applicable in the frequency range 1-350 GHz.
      * 'exact': an estimate of gaseous attenuation computed by summation of
        individual absorption lines that is valid for the frequency
        range 1-1,000 GHz


    Parameters
    ----------
    f : number or Quantity
        Frequency (GHz)
    el : sequence, number or Quantity
        Elevation angle (degrees)
    rho : number or Quantity
        Water vapor density (g/m3)
    P : number or Quantity
        Atmospheric pressure (hPa)
    T : number or Quantity
        Absolute temperature (K)
    h1 : number or Quantity
        Height of ground station 1 (km)
    h2 : number or Quantity
        Height of ground station 2 (km)
    mode : string, optional
        Mode for the calculation. Valid values are 'approx', 'exact'. If
        'approx' Uses the method in Annex 2 of the recommendation (if any),
        else uses the method described in Section 1. Default, 'approx'


    Returns
    -------
    attenuation: Quantity
        Inclined path attenuation (dB)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en
    """
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    type_output = get_input_type(el)
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapor density')
    P = prepare_quantity(P, u.hPa, 'Atospheric pressure')
    T = prepare_quantity(T, u.K, 'Temperature')
    h1 = prepare_quantity(h1, u.km, 'Height of Ground Station 1')
    h2 = prepare_quantity(h2, u.km, 'Height of Ground Station 2')
    val = __model.gaseous_attenuation_inclined_path(
        f, el, rho, P, T, h1, h2, mode=mode)
    return prepare_output_array(val, type_output) * u.dB


def slant_inclined_path_equivalent_height(f, p, rho=7.5, T=298.15):
    """ Computes the equivalent height to be used for oxygen and water vapour
    gaseous attenuation computations.

    Parameters
    ----------
    f : number or Quantity
        Frequency (GHz)
    p : number
        Percentage of the time the gaseous attenuation value is exceeded.
    rho : number or Quantity
        Water vapor density (g/m3)
    P : number or Quantity
        Atmospheric pressure (hPa)
    T : number or Quantity
        Absolute temperature (K)

    Returns
    -------
    ho, hw : Quantity
        Equivalent height for oxygen and water vapour (m)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en

    """
    type_output = get_input_type(f)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapor density')
    T = prepare_quantity(T, u.K, 'Temperature')

    val = __model.slant_inclined_path_equivalent_height(f, p, rho, T)
    return prepare_output_array(val, type_output) * u.m


def zenit_water_vapour_attenuation(lat, lon, p, f, V_t=None, h=None):
    """
    An alternative method may be used to compute the slant path attenuation by
    water vapour, in cases where the integrated water vapour content along the
    path, ``V_t``, is known.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    p : number
        Percentage of the time the zenit water vapour attenuation value is
        exceeded.
    f : number or Quantity
        Frequency (GHz)
    V_t : number or Quantity, optional
        Integrated water vapour content along the path (kg/m2 or mm).
        If not provided this value is estimated using Recommendation
        ITU-R P.836. Default value None
    h : number, sequence, or numpy.ndarray
        Altitude of the receivers. If None, use the topographical altitude as
        described in recommendation ITU-R P.1511


    Returns
    -------
    A_w : Quantity
        Water vapour attenuation along the slant path (dB)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    V_t = prepare_quantity(
        V_t,
        u.kg / u.m**2,
        'Integrated water vapour content along the path')
    h = prepare_quantity(h, u.km, 'Altitude')
    val = __model.zenit_water_vapour_attenuation(
        lat, lon, p, f, V_t=V_t, h=h)
    return prepare_output_array(val, type_output) * u.dB


def gammaw_approx(f, P, rho, T):
    """
    Method to estimate the specific attenuation due to water vapour using the
    approximate method descibed in Annex 2.


    Parameters
    ----------
    f : number or Quantity
        Frequency (GHz)
    P : number or Quantity
        Atmospheric pressure (hPa)
    rho : number or Quantity
        Water vapor density (g/m3)
    T : number or Quantity
        Absolute temperature (K)


    Returns
    -------
    gamma_w : Quantity
        Water vapour specific attenuation (dB/km)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en
    """
    type_output = get_input_type(f)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    P = prepare_quantity(P, u.hPa, 'Atmospheric pressure ')
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapour density')
    T = prepare_quantity(T, u.K, 'Temperature')
    val = __model.gammaw_approx(f, P, rho, T)
    return prepare_output_array(val, type_output) * u.dB / u.km


def gamma0_approx(f, P, rho, T):
    """
    Method to estimate the specific attenuation due to dry atmosphere using the
    approximate method descibed in Annex 2.

    Parameters
    ----------
    f : number or Quantity
        Frequency (GHz)
    P : number or Quantity
        Atmospheric pressure (hPa)
    rho : number or Quantity
        Water vapor density (g/m3)
    T : number or Quantity
        Absolute temperature (K)


    Returns
    -------
    gamma_w : Quantity
        Dry atmosphere specific attenuation (dB/km)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en
    """
    type_output = get_input_type(f)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    P = prepare_quantity(P, u.hPa, 'Atmospheric pressure')
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapour density')
    T = prepare_quantity(T, u.K, 'Temperature')
    val = __model.gamma0_approx(f, P, rho, T)
    return prepare_output_array(val, type_output) * u.dB / u.km


def gammaw_exact(f, P, rho, T):
    """
    Method to estimate the specific attenuation due to water vapour using
    the line-by-line method described in Annex 1 of the recommendation.


    Parameters
    ----------
    f : number or Quantity
        Frequency (GHz)
    P : number or Quantity
        Atmospheric pressure (hPa)
    rho : number or Quantity
        Water vapor density (g/m3)
    T : number or Quantity
        Absolute temperature (K)


    Returns
    -------
    gamma_w : Quantity
        Water vapour specific attenuation (dB/km)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en
    """
    type_output = get_input_type(f)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    P = prepare_quantity(P, u.hPa, 'Atmospheric pressure ')
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapour density')
    T = prepare_quantity(T, u.K, 'Temperature')
    val = __model.gammaw_exact(f, P, rho, T)
    return prepare_output_array(val, type_output) * u.dB / u.km


def gamma0_exact(f, P, rho, T):
    """
    Method to estimate the specific attenuation due to dry atmosphere using
    the line-by-line method described in Annex 1 of the recommendation.

    Parameters
    ----------
    f : number or Quantity
        Frequency (GHz)
    P : number or Quantity
        Atmospheric pressure (hPa)
    rho : number or Quantity
        Water vapor density (g/m3)
    T : number or Quantity
        Absolute temperature (K)


    Returns
    -------
    gamma_w : Quantity
        Dry atmosphere specific attenuation (dB/km)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en
    """
    type_output = get_input_type(f)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    P = prepare_quantity(P, u.hPa, 'Atmospheric pressure')
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapour density')
    T = prepare_quantity(T, u.K, 'Temperature')
    val = __model.gamma0_exact(f, P, rho, T)
    return prepare_output_array(val, type_output) * u.dB / u.km


def gamma_exact(f, P, rho, T):
    """
    Method to estimate the specific attenuation using the line-by-line method
    described in Annex 1 of the recommendation.


    Parameters
    ----------
    f : number or Quantity
        Frequency (GHz)
    P : number or Quantity
        Atmospheric pressure (hPa)
    rho : number or Quantity
        Water vapor density (g/m3)
    T : number or Quantity
        Absolute temperature (K)

    Returns
    -------
    gamma : Quantity
        Specific attenuation (dB/km)

    References
    --------
    [1] Attenuation by atmospheric gases:
    https://www.itu.int/rec/R-REC-P.676/en
    """
    f = prepare_quantity(f, u.GHz, 'Frequency')
    type_output = get_input_type(f)

    P = prepare_quantity(P, u.hPa, 'Atmospheric pressure ')
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapour density')
    T = prepare_quantity(T, u.K, 'Temperature')
    val = __model.gamma_exact(f, P, rho, T)
    return prepare_output_array(val, type_output) * u.dB / u.km
