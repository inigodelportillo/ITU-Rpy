# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import scipy.stats as stats
import scipy.special
import scipy.integrate
from astropy import units as u

from itur.models.itu453 import water_vapour_pressure,\
    wet_term_radio_refractivity, map_wet_term_radio_refractivity
from itur.models.itu837 import rainfall_rate, rainfall_probability
from itur.models.itu838 import rain_specific_attenuation
from itur.models.itu839 import rain_height
from itur.models.itu1511 import topographic_altitude
from itur.utils import prepare_input_array, prepare_output_array,\
    prepare_quantity, compute_distance_earth_to_earth, get_input_type, EPSILON

import warnings


def __CDF_bivariate_normal__(alpha_x, alpha_y, rho):
    # This function calculates the complementary bivariate normal
    # distribution with limits alpha_x, alpha_y and correlation factor rho
    def CDF_bivariate_normal_fcn(x, y, rho):
        return np.exp(- (x**2 - 2 * rho * x * y + y**2) /
                      (2. * (1 - rho**2)))

    def CDF_bivariate_normal_int(alpha, y, rho):
        return scipy.integrate.quad(
            CDF_bivariate_normal_fcn, alpha, np.inf, args=(y, rho))[0]

    return 1 / (2 * np.pi * np.sqrt(1 - rho**2)) * scipy.integrate.quad(
        lambda y: CDF_bivariate_normal_int(alpha_x, y, rho),
        alpha_y,
        np.inf)[0]


class _ITU618():
    """
    Propagation data and prediction methods required for the design of
    Earth-space telecommunication systems.

    Available versions include:
       * P.618-13 (12/17) (Current version)
       * P.618-12 (07/15) (Superseded)

    Versions that need to be implemented
       * P.618-11
       * P.618-10
       * P.618-09
       * P.618-08
       * P.618-07
       * P.618-06
       * P.618-05
       * P.618-04
       * P.618-03
       * P.618-02
       * P.618-01

    Recommendation ITU-R P.618 provides methods to estimate the propagation
    loss on an Earth-space path, relative to the free-space loss. This value
    is the sum of different contributions as follows:
      * attenuation by atmospheric gases;
      * attenuation by rain, other precipitation and clouds;
      * focusing and defocusing;
      * decrease in antenna gain due to wave-front incoherence;
      * scintillation and multipath effects;
      * attenuation by sand and dust storms.

    Each of these contributions has its own characteristics as a function of
    frequency, geographic location and elevation angle. As a rule, at elevation
    angles above 10°, only gaseous attenuation, rain and cloud attenuation and
    possibly scintillation will be significant, depending on propagation
    conditions.
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.618 recommendation.

    def __init__(self, version=13):
        if version == 13:
            self.instance = _ITU618_13()
        elif version == 12:
            self.instance = _ITU618_12()
#        elif version == 11:
#            self.instance = _ITU618_11()
#        elif version == 10:
#            self.instance = _ITU618_10()
#        elif version == 9:
#            self.instance = _ITU618_9()
#        elif version == 8:
#            self.instance = _ITU618_8()
#        elif version == 7:
#            self.instance = _ITU618_7()
#        elif version == 6:
#            self.instance = _ITU618_6()
#        elif version == 5:
#            self.instance = _ITU618_5()
#        elif version == 4:
#            self.instance = _ITU618_4()
#        elif version == 3:
#            self.instance = _ITU618_3()
#        elif version == 2:
#            self.instance = _ITU618_2()
#        elif version == 1:
#            self.instance = _ITU618_1()
        else:
            raise ValueError(('Version {0} is not implemented'
                              ' for the ITU-R P.618 model.').format(version))

    @property
    def __version__(self):
        return self.instance.__version__

    def rain_attenuation(self, lat, lon, f, el, hs=None, p=0.01, R001=None,
                         tau=45, Ls=None):
        fcn = np.vectorize(self.instance.rain_attenuation,
                           excluded=[0, 1, 3, 4, 6], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, f, el, hs, p, R001, tau, Ls).tolist())

    def rain_attenuation_probability(self, lat, lon, el, hs, Ls, P0=None):
        fcn = np.vectorize(self.instance.rain_attenuation_probability,
                           excluded=[0, 1, 2], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, el, hs, Ls, P0).tolist())

    def rain_cross_polarization_discrimination(self, Ap, f, el, p, tau):
        fcn = np.vectorize(
            self.instance.rain_cross_polarization_discrimination)
        return fcn(Ap, f, el, p, tau)

    def scintillation_attenuation(self, lat, lon, f, el, p, D, eta,
                                  T, H, P, hL):
        fcn = np.vectorize(self.instance.scintillation_attenuation,
                           excluded=[0, 1, 3, 7, 8, 9], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, f, el, p, D, eta, T, H, P, hL).tolist())

    def scintillation_attenuation_sigma(self, lat, lon, f, el, p, D, eta,
                                        T, H, P, hL):
        fcn = np.vectorize(self.instance.scintillation_attenuation_sigma,
                           excluded=[0, 1, 3, 7, 8, 9], otypes=[np.ndarray])
        return np.array(fcn(lat, lon, f, el, p, D, eta, T, H, P, hL).tolist())

    def fit_rain_attenuation_to_lognormal(self, lat, lon, f, el, hs, P_k, tau):
        fcn = np.vectorize(self.instance.fit_rain_attenuation_to_lognormal)
        return fcn(lat, lon, f, el, hs, P_k, tau)

    def site_diversity_rain_outage_probability(self, lat1, lon1, a1, el1,
                                               lat2, lon2, a2, el2, f, tau=45,
                                               hs1=None, hs2=None):
        fcn = np.vectorize(
            self.instance.site_diversity_rain_outage_probability)
        return np.array(fcn(lat1, lon1, a1, el1,
                            lat2, lon2, a2, el2,
                            f, tau, hs1, hs2).tolist())


class _ITU618_13():

    def __init__(self):
        self.__version__ = 13

    @classmethod
    def rain_attenuation(self, lat, lon, f, el, hs=None, p=0.01, R001=None,
                         tau=45, Ls=None):
        if np.logical_or(p < 0.001, p > 5).any():
            warnings.warn(
                RuntimeWarning('The method to compute the rain attenuation in '
                               'recommendation ITU-P 618-12 is only valid for '
                               'unavailability values between 0.001 and 5'))

        Re = 8500   # Efective radius of the Earth (8500 km)

        if hs is None:
            hs = topographic_altitude(lat, lon).to(u.km).value

        # Step 1: Compute the rain height (hr) based on ITU - R P.839
        hr = rain_height(lat, lon).value

        # Step 2: Compute the slant path length
        if Ls is None:
            Ls = np.where(
                el >= 5, (hr - hs) / (np.sin(np.deg2rad(el))),         # Eq. 1
                2 * (hr - hs) / (((np.sin(np.deg2rad(el)))**2 +
                                  2 * (hr - hs) / Re)**0.5 + (np.sin(np.deg2rad(el)))))  # Eq. 2

        # Step 3: Calculate the horizontal projection, LG, of the
        # slant-path length
        Lg = np.abs(Ls * np.cos(np.deg2rad(el)))

        # Obtain the raingall rate, exceeded for 0.01% of an average year,
        # if not provided, as described in ITU-R P.837.
        if R001 is None:
            R001 = rainfall_rate(lat, lon, 0.01).to(u.mm / u.hr).value + EPSILON

        # Step 5: Obtain the specific attenuation gammar using the frequency
        # dependent coefficients as given in ITU-R P.838
        # https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.838-3-200503-I!!PDF-E.pdf
        gammar = rain_specific_attenuation(
            R001, f, el, tau).to(
            u.dB / u.km).value

        # Step 6: Calculate the horizontal reduction factor, r0.01,
        # for 0.01% of the time:
        r001 = 1. / (1 + 0.78 * np.sqrt(Lg * gammar / f) -
                     0.38 * (1 - np.exp(-2 * Lg)))

        # Step 7: Calculate the vertical adjustment factor, v0.01,
        # for 0.01% of the time:
        eta = np.rad2deg(np.arctan2(hr - hs, Lg * r001))

        Delta_h = np.where(hr - hs <= 0, EPSILON, (hr - hs))
        Lr = np.where(eta > el, Lg * r001 / np.cos(np.deg2rad(el)),
                      Delta_h / np.sin(np.deg2rad(el)))

        xi = np.where(np.abs(lat) < 36, 36 - np.abs(lat), 0)

        v001 = 1. / (1 + np.sqrt(np.sin(np.deg2rad(el))) *
                     (31 * (1 - np.exp(-(el / (1 + xi)))) *
                      np.sqrt(Lr * gammar) / f**2 - 0.45))

        # Step 8: calculate the effective path length:
        Le = Lr * v001   # (km)

        # Step 9: The predicted attenuation exceeded for 0.01% of an average
        # year
        A001 = gammar * Le   # (dB)

        # Step 10: The estimated attenuation to be exceeded for other
        # percentages of an average year
        if p >= 1:
            beta = np.zeros_like(A001)
        else:
            beta = np.where(np.abs(lat) >= 36,
                            np.zeros_like(A001),
                            np.where((np.abs(lat) < 36) & (el > 25),
                                     -0.005 * (np.abs(lat) - 36),
                                     -0.005 * (np.abs(lat) - 36) + 1.8 -
                                     4.25 * np.sin(np.deg2rad(el))))

        A = A001 * (p / 0.01)**(
            -(0.655 + 0.033 * np.log(p) - 0.045 * np.log(A001) -
              beta * (1 - p) * np.sin(np.deg2rad(el))))

        return A

    @classmethod
    def rain_attenuation_probability(self, lat, lon, el, hs=None,
                                     Ls=None, P0=None):

        Re = 8500
        if hs is None:
            hs = topographic_altitude(lat, lon).to(u.km).value

        # Step 1: Estimate the probability of rain, at the earth station either
        # from Recommendation ITU-R P.837 or from local measured rainfall
        # rate data
        if P0 is None:
            P0 = rainfall_probability(lat, lon).\
                to(u.dimensionless_unscaled).value

        # Step 2: Calculate the parameter alpha using the inverse of the
        # Q-function alpha = Q^{-1}(P0)  -> Q(alpha) = P0
        alpha = stats.norm.ppf(1 - P0)

        # Step 3: Calculate the spatial correlation function, rho:
        hr = rain_height(lat, lon).value

        if Ls is None:
            Ls = np.where(
                el >= 5, (hr - hs) / (np.sin(np.deg2rad(el))),         # Eq. 1
                2 * (hr - hs) / (((np.sin(np.deg2rad(el)))**2 +
                                  2 * (hr - hs) / Re)**0.5 + (np.sin(np.deg2rad(el)))))  # Eq. 2

        d = Ls * np.cos(np.deg2rad(el))
        rho = 0.59 * np.exp(-abs(d) / 31) + 0.41 * np.exp(-abs(d) / 800)

        # Step 4: Calculate the complementary bivariate normal distribution
        biva_fcn = np.vectorize(__CDF_bivariate_normal__)
        c_B = biva_fcn(alpha, alpha, rho)

        # Step 5: Calculate the probability of rain attenuation on the slant
        # path:
        P = 1 - (1 - P0) * ((c_B - P0**2) / (P0 * (1 - P0)))**P0
        return P

    @classmethod
    def fit_rain_attenuation_to_lognormal(self, lat, lon, f, el, hs, P_k, tau):
        # Performs the log-normal fit of rain attenuation vs. probability of
        # occurrence for a particular path

        # Step 1: Construct the set of pairs [Pi, Ai] where Pi (% of time) is
        # the probability the attenuation Ai (dB) is exceeded where Pi < P_K
        p_i = np.array([0.01, 0.02, 0.03, 0.05,
                        0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10])
        Pi = np.array([p for p in p_i if p < P_k], dtype=float)
        Ai = np.array([0 for p in p_i if p < P_k], dtype=float)

        for i, p in enumerate(Pi):
            Ai[i] = self.rain_attenuation(lat, lon, f, el, hs, p, tau=tau)

        # Step 2: Transform the set of pairs [Pi, Ai] to [Q^{-1}(Pi/P_k),
        # ln(Ai)]
        Q = stats.norm.ppf(1 - (Pi / P_k))
        lnA = np.log(Ai)

        # Step 3: Determine the variables sigma_lna, m_lna by performing a
        # least-squares fit to lnAi = sigma_lna Q^{-1}(Pi/P_k) + m_lna
        m_lna, sigma_lna = np.linalg.lstsq(np.vstack([np.ones(len(Q)), Q]).T,
                                           lnA, rcond=None)[0]

        return sigma_lna, m_lna

    @classmethod
    def site_diversity_rain_outage_probability(self, lat1, lon1, a1, lat2,
                                               lon2, a2, f, el1, el2, tau=45,
                                               hs1=None, hs2=None):
        # The diversity prediction method assumes a log-normal distribution of
        # rain intensity and rain attenuation. This method predicts
        # Pr(A1 > a1, A2 > a2), the joint probability (%) that the attenuation
        # on the path to the first site is greater than a1 and the attenuation
        # on the path to the second site is greater than a2.
        d = compute_distance_earth_to_earth(lat1, lon1, lat2, lon2)
        rho_r = 0.7 * np.exp(-d / 60) + 0.3 * np.exp(-(d / 700)**2)

        P_1 = rainfall_probability(lat1, lon1).\
            to(u.dimensionless_unscaled).value

        P_2 = rainfall_probability(lat2, lon2).\
            to(u.dimensionless_unscaled).value

        R_1 = stats.norm.ppf(1 - P_1)
        R_2 = stats.norm.ppf(1 - P_2)
        biva_fcn = np.vectorize(__CDF_bivariate_normal__)
        P_r = biva_fcn(R_1, R_2, rho_r)

        sigma_lna1, m_lna1 = self.fit_rain_attenuation_to_lognormal(
            lat1, lon1, f, el1, hs1, P_1 * 100, tau)

        sigma_lna2, m_lna2 = self.fit_rain_attenuation_to_lognormal(
            lat2, lon2, f, el2, hs2, P_2 * 100, tau)

        rho_a = 0.94 * np.exp(-d / 30) + 0.06 * np.exp(-(d / 500)**2)
        lim_1 = (np.log(a1) - m_lna1) / sigma_lna1
        lim_2 = (np.log(a2) - m_lna2) / sigma_lna2

        P_a = biva_fcn(lim_1, lim_2, rho_a)

        return 100 * P_r * P_a

    @classmethod
    def rain_cross_polarization_discrimination(self, Ap, f, el, p, tau):
        # Frequency reuse by means of orthogonal polarizations is often used to
        # increase the capacity of space telecommunication systems. This
        # technique is restricted, however, by depolarization on atmospheric
        # propagation paths. Various depolarization mechanisms, especially
        # hydrometeor effects, are important in the troposphere

        # The method described below to calculate cross-polarization
        # discrimination (XPD) statistics from rain attenuation statistics for
        # the same path is valid for 6 < f < 55 GHz and el < 60°.
        if f < 4 or f > 55:
            warnings.warn(
                RuntimeWarning(
                    'The method to compute the cross '
                    'polarization discrimination in recommendation '
                    'ITU-P 618-12 is only valid for frequency values between'
                    ' 4 and 55 GHz'))

        if el > 60:
            warnings.warn(
                RuntimeWarning(
                    'The method to compute thecross '
                    'polarization discrimination in recommendation ITU-P '
                    '618-12 is only valid for elevation angle values below '
                    '60 degrees'))

        # In case that the frequency is comprised between 4 and 6 GHz, scaling
        # is necessary
        scale_to_orig_f = False
        if 4 <= f < 6:
            f_orig = f
            f = 6
            scale_to_orig_f = True

        # Step 1: Calculate the frequency-dependent term:
        if 6 <= f < 9:
            C_f = 60 * np.log10(f) - 28.3
        elif 9 <= f < 36:
            C_f = 26 * np.log10(f) + 4.1
        elif 36 <= f <= 55:
            C_f = 35.9 * np.log10(f) - 11.3

        # Step 2: Calculate the rain attenuation dependent term:
        if 6 <= f < 9:
            V = 30.8 * f**-0.21
        elif 9 <= f < 20:
            V = 12.8 * f**0.19
        elif 20 <= f < 40:
            V = 22.6
        elif 40 <= f <= 55:
            V = 13.0 * f**0.15

        C_a = V * np.log10(Ap)

        # Step 3: Calculate the polarization improvement factor:
        C_tau = -10 * np.log10(1 - 0.484 * (1 + np.cos(np.deg2rad(4 * tau))))

        # Step 4: Calculate the elevation angle-dependent term:
        C_theta = -40 * np.log10(np.cos(np.deg2rad(el)))

        # Step 5: Calculate the canting angle dependent term:
        if p <= 0.001:
            C_sigma = 0.0053 * 15**2
        elif p <= 0.01:
            C_sigma = 0.0053 * 10**2
        elif p <= 0.1:
            C_sigma = 0.0053 * 5**2
        else:
            C_sigma = 0

        # Step 6: Calculate rain XPD not exceeded for p% of the time:
        XPD_rain = C_f - C_a + C_tau + C_theta + C_sigma

        # Step 7: Calculate the ice crystal dependent term:
        C_ice = XPD_rain * (0.3 + 0.1 * np.log10(p)) / 2

        # Step 8: Calculate the XPD not exceeded for p% of the time,
        # including the effects of ice:
        XPD_p = XPD_rain - C_ice

        if scale_to_orig_f:
            # Long-term XPD statistics obtained at one frequency and
            # polarization tilt angle can be scaled to another frequency and
            # polarization tilt angle using the semi-empirical formula:
            XPD_p = XPD_p - 20 * np.log10(
                f_orig * np.sqrt(1 - 0.484 * (1 + np.cos(np.deg2rad(4 * tau)))) /
                (f * np.sqrt(1 - 0.484 * (1 + np.cos(np.deg2rad(4 * tau))))))
        return XPD_p

    @classmethod
    def scintillation_attenuation_sigma(cls, lat, lon, f, el, p, D, eta=0.5,
                                        T=None, H=None, P=None, hL=1000):
        # Step 1: For the value of t, calculate the saturation water vapour
        # pressure, es, (hPa), as specified in Recommendation ITU-R P.453.
        if T is not None and H is not None and P is not None:
            e = water_vapour_pressure(T, P, H).value

            # Step 2: Compute the wet term of the radio refractivity, Nwet,
            # corresponding to es, t and H as given in Recommendation ITU-R
            # P.453.
            N_wet = wet_term_radio_refractivity(e, T).value
        else:
            N_wet = map_wet_term_radio_refractivity(lat, lon, 50).value

        # Step 3: Calculate the standard deviation of the reference signal
        # amplitude:
        sigma_ref = 3.6e-3 + 1e-4 * N_wet  # Eq. 43   [dB]

        # Step 4: Calculate the effective path length L:
        L = 2 * hL / (np.sqrt(np.sin(np.deg2rad(el))**2 + 2.35e-4) +
                      np.sin(np.deg2rad(el)))  # Eq. 44   [m]

        # Step 5: Estimate the effective antenna diameter, Deff
        D_eff = np.sqrt(eta) * D  # Eq. 45   [m]

        # Step 6: Step 6: Calculate the antenna averaging factor
        x = 1.22 * D_eff**2 * f / L
        g = np.where(x >= 7.0, 0,
                     np.sqrt(3.86 * (x**2 + 1)**(11. / 12) *
                             np.sin(11. / 6 * np.arctan2(1, x)) -
                             7.08 * x**(5. / 6)))  # Eq. 46    [-]

        # Step 7: Calculate the standard deviation of the signal for the
        # applicable period and propagation path:
        sigma = sigma_ref * f**(7. / 12) * g / np.sin(np.deg2rad(el))**1.2
        return sigma

    @classmethod
    def scintillation_attenuation(cls, lat, lon, f, el, p, D, eta=0.5, T=None,
                                  H=None, P=None, hL=1000):
        # Step 1 - 7: Calculate the standard deviation of the signal for the
        # applicable period and propagation path:
        sigma = cls.scintillation_attenuation_sigma(lat, lon, f, el, p,
                                                     D, eta, T, H, P, hL)
        # Step 8: Calculate the time percentage factor, a(p), for the time
        # percentage, p, in the range between 0.01% < p < 50%:
        a = -0.061 * np.log10(p)**3 + 0.072 * \
            np.log10(p)**2 - 1.71 * np.log10(p) + 3

        # Step 9: Calculate the fade depth, A(p), exceeded for p% of the time:
        A_s = a * sigma  # Eq. 49   [dB]

        return A_s


class _ITU618_12():

    def __init__(self):
        self.__version__ = 12

    @classmethod
    def rain_attenuation(self, lat, lon, f, el, hs=None, p=0.01, R001=None,
                         tau=45, Ls=None):

        if p < 0.001 or p > 5:
            warnings.warn(
                RuntimeWarning('The method to compute the rain attenuation in '
                               'recommendation ITU-P 618-12 is only valid for '
                               'unavailability values between 0.001% and 5%'))

        Re = 8500   # Efective radius of the Earth (8500 km)

        if hs is None:
            hs = topographic_altitude(lat, lon).to(u.km).value

        # Step 1: Compute the rain height (hr) based on ITU - R P.839
        hr = rain_height(lat, lon).value

        # Step 2: Compute the slant path length
        if Ls is None:
            Ls = np.where(
                el >= 5, (hr - hs) / (np.sin(np.deg2rad(el))),  # Eq. 1
                2 * (hr - hs) / (((np.sin(np.deg2rad(el)))**2 +
                                  2 * (hr - hs) / Re)**0.5 +
                                 (np.sin(np.deg2rad(el)))))  # Eq. 2

        # Step 3: Calculate the horizontal projection, LG, of the
        # slant-path length
        Lg = np.abs(Ls * np.cos(np.deg2rad(el)))

        # Obtain the raingall rate, exceeded for 0.01% of an average year,
        # if not provided, as described in ITU-R P.837.
        if R001 is None:
            R001 = rainfall_rate(lat, lon, 0.01).to(u.mm / u.hr).value + EPSILON

        # Step 5: Obtain the specific attenuation gammar using the frequency
        # dependent coefficients as given in ITU-R P.838
        gammar = rain_specific_attenuation(
            R001, f, el, tau).to(
            u.dB / u.km).value

        # Step 6: Calculate the horizontal reduction factor, r0.01,
        # for 0.01% of the time:
        r001 = 1. / (1 + 0.78 * np.sqrt(Lg * gammar / f) -
                     0.38 * (1 - np.exp(-2 * Lg)))

        # Step 7: Calculate the vertical adjustment factor, v0.01,
        # for 0.01% of the time:
        eta = np.rad2deg(np.arctan2(hr - hs, Lg * r001))
        
        Delta_h = np.where(hr - hs <= 0, EPSILON, (hr - hs))
        Lr = np.where(eta > el, Lg * r001 / np.cos(np.deg2rad(el)),
                      Delta_h / np.sin(np.deg2rad(el)))

        xi = np.where(np.abs(lat) < 36, 36 - np.abs(lat), 0)

        v001 = 1. / (1 + np.sqrt(np.sin(np.deg2rad(el))) *
                     (31 * (1 - np.exp(-(el / (1 + xi)))) *
                      np.sqrt(Lr * gammar) / f**2 - 0.45))

        # Step 8: calculate the effective path length:
        Le = Lr * v001   # (km)

        # Step 9: The predicted attenuation exceeded for 0.01% of an average
        # year
        A001 = gammar * Le   # (dB)

        # Step 10: The estimated attenuation to be exceeded for other
        # percentages of an average year
        if p >= 1:
            beta = np.zeros_like(A001)
        else:
            beta = np.where(np.abs(lat) >= 36,
                            np.zeros_like(A001),
                            np.where((np.abs(lat) < 36) & (el > 25),
                                     -0.005 * (np.abs(lat) - 36),
                                     -0.005 * (np.abs(lat) - 36) + 1.8 -
                                     4.25 * np.sin(np.deg2rad(el))))

        A = A001 * (p / 0.01)**(-(0.655 + 0.033 * np.log(p) -
                                  0.045 * np.log(A001) -
                                  beta * (1 - p) * np.sin(np.deg2rad(el))))

        return A

    @classmethod
    def rain_attenuation_probability(self, *args, **kwargs):
        return _ITU618_13.rain_attenuation_probability(*args, **kwargs)

    @classmethod
    def fit_rain_attenuation_to_lognormal(self, *args, **kwargs):
        return _ITU618_13.fit_rain_attenuation_to_lognormal(*args, **kwargs)

    @classmethod
    def site_diversity_rain_outage_probability(self, *args, **kwargs):
        return _ITU618_13.site_diversity_rain_outage_probability(*args,
                                                                 **kwargs)

    @classmethod
    def rain_cross_polarization_discrimination(self, *args, **kwargs):
        return _ITU618_13.rain_cross_polarization_discrimination(*args,
                                                                 **kwargs)

    @classmethod
    def scintillation_attenuation(self, *args, **kwargs):
        return _ITU618_13.scintillation_attenuation(*args, **kwargs)


__model = _ITU618()


def change_version(new_version):
    """
    Change the version of the ITU-R P.618 recommendation currently being used.

    This function changes the model used for the ITU-R P.618 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          *  13: Activates recommendation ITU-R P.618-13 (12/17) (Current version)
          *  12: Activates recommendation ITU-R P.618-12 (07/15) (Superseded)
    """
    global __model
    __model = _ITU618(new_version)


def get_version():
    """ The version of the current model for the ITU-R P.618 recommendation.

    Obtain the version of the ITU-R P.618 recommendation currently being used.

    Returns
    -------
    version: int
       The version of the ITU-R P.618 recommendation being used.
    """
    return __model.__version__


def rain_attenuation(lat, lon, f, el, hs=None, p=0.01, R001=None,
                     tau=45, Ls=None):
    """
    Calculation of long-term rain attenuation statistics from point rainfall
    rate.
    The following procedure provides estimates of the long-term statistics of
    the slant-path rain attenuation at a given location for frequencies up
    to 55 GHz.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    f : number
        Frequency (GHz)
    el : sequence, or number
        Elevation angle (degrees)
    hs : number, sequence, or numpy.ndarray, optional
        Heigh above mean sea level of the earth station (km). If local data for
        the earth station height above mean sea level is not available, an
        estimate is obtained from the maps of topographic altitude
        given in Recommendation ITU-R P.1511.
    p : number, optional
        Percentage of the time the rain attenuation value is exceeded.
    R001: number, optional
        Point rainfall rate for the location for 0.01% of an average year
        (mm/h).
        If not provided, an estimate is obtained from Recommendation
        Recommendation ITU-R P.837. Some useful values:
            * 0.25 mm/h : Drizzle
            *  2.5 mm/h : Light rain
            * 12.5 mm/h : Medium rain
            * 25.0 mm/h : Heavy rain
            * 50.0 mm/h : Downpour
            * 100  mm/h : Tropical
            * 150  mm/h : Monsoon
    tau : number, optional
        Polarization tilt angle relative to the horizontal (degrees)
        (tau = 45 deg for circular polarization). Default value is 45
    Ls :number, optional
        Slant path length (km). If not provided, it will be computed using the
        rain height and the elevation angle. The ITU model does not require
        this parameter as an input.


    Returns
    -------
    attenuation: Quantity
        Attenuation due to rain (dB)

    References
    --------
    [1] Propagation data and prediction methods required for the design of
    Earth-space telecommunication systems:
    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.618-12-201507-I!!PDF-E.pdf
    """
    type_output = get_input_type(lat)

    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)

    lon = np.mod(lon, 360)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(prepare_input_array(el), u.deg, 'Elevation angle')
    hs = prepare_quantity(
        hs, u.km, 'Heigh above mean sea level of the earth station')
    R001 = prepare_quantity(R001, u.mm / u.hr, 'Point rainfall rate')
    tau = prepare_quantity(tau, u.one, 'Polarization tilt angle')
    Ls = prepare_quantity(Ls, u.km, 'Slant path length')

    val = __model.rain_attenuation(lat, lon, f, el, hs=hs, p=p,
                                   R001=R001, tau=tau, Ls=Ls)
    
    # The values of attenuation cannot be negative. The ITU models end up
    # giving out negative values for certain inputs
    val[val < 0] = 0
    
    return prepare_output_array(val, type_output) * u.dB


def rain_attenuation_probability(lat, lon, el, hs=None, Ls=None, P0=None):
    """
    The following procedure computes the probability of non-zero rain
    attenuation on a given slant path Pr(Ar > 0).


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    el : sequence, or number
        Elevation angle (degrees)
    hs : number, sequence, or numpy.ndarray, optional
        Heigh above mean sea level of the earth station (km). If local data for
        the earth station height above mean sea level is not available, an
        estimate is obtained from the maps of topographic altitude
        given in Recommendation ITU-R P.1511.
    Ls : number, sequence, or numpy.ndarray, optional
        Slant path length from the earth station to the rain height (km). If
        data about the rain height is not available, this value is estimated
        automatically using Recommendation ITU-R P.838
    P0 : number, sequence, or numpy.ndarray, optional
        Probability of rain at the earth station, (0 ≤ P0 ≤ 1)



    Returns
    -------
    p: Quantity
         Probability of rain attenuation on the slant path (%)


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    Earth-space telecommunication systems:
    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.618-12-201507-I!!PDF-E.pdf
    """
    type_output = get_input_type(lat)

    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)

    lon = np.mod(lon, 360)
    el = prepare_quantity(prepare_input_array(el), u.deg, 'Elevation angle')
    hs = prepare_quantity(
        hs, u.km, 'Heigh above mean sea level of the earth station')
    Ls = prepare_quantity(
        Ls, u.km, 'Heigh above mean sea level of the earth station')
    P0 = prepare_quantity(P0, u.pct, 'Point rainfall rate')

    val = __model.rain_attenuation_probability(lat, lon, el, hs, Ls, P0)

    return prepare_output_array(val, type_output) * 100 * u.pct


def site_diversity_rain_outage_probability(lat1, lon1, a1, el1, lat2,
                                           lon2, a2, el2, f, tau=45, hs1=None,
                                           hs2=None):
    """
    Calculate the link outage probability in a diversity based scenario (two
    ground stations) due to rain attenuation. This method is valid for
    frequencies below 20 GHz, as at higher frequencies other impairments might
    affect affect site diversity performance.

    This method predicts Pr(A1 > a1, A2 > a2), the joint probability (%) that
    the attenuation on the path to the first site is greater than a1 and the
    attenuation on the path to the second site is greater than a2.


    Parameters
    ----------
    lat1 : number or Quantity
        Latitude of the first ground station (deg)
    lon1 : number or Quantity
        Longitude of the first ground station (deg)
    a1 : number or Quantity
        Maximum admissible attenuation of the first ground station (dB)
    el1 : number or Quantity
        Elevation angle to the first ground station (deg)
    lat2 : number or Quantity
        Latitude of the second ground station (deg)
    lon2 : number or Quantity
        Longitude of the second ground station (deg)
    a2 : number or Quantity
        Maximum admissible attenuation of the second ground station (dB)
    el2 : number or Quantity
        Elevation angle to the second ground station (deg)
    f : number or Quantity
        Frequency (GHz)
    tau : number, optional
        Polarization tilt angle relative to the horizontal (degrees)
        (tau = 45 deg for circular polarization). Default value is 45
    hs1 : number or Quantity, optional
        Altitude over the sea level of the first ground station (km). If not
        provided, uses Recommendation ITU-R P.1511 to compute the toporgraphic
        altitude
    hs2 : number or Quantity, optional
        Altitude over the sea level of the first ground station (km). If not
        provided, uses Recommendation ITU-R P.1511 to compute the toporgraphic
        altitude


    Returns
    -------
    probability: Quantity
        Joint probability (%) that the attenuation on the path to the first
        site is greater than a1 and the attenuation on the path to the second
        site is greater than a2


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    Earth-space telecommunication systems:
    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.618-12-201507-I!!PDF-E.pdf
    """
    type_output = get_input_type(lat1)
    lon1 = np.mod(lon1, 360)
    lat1 = prepare_quantity(lat1, u.deg, 'Latitude in ground station 1')
    lon1 = prepare_quantity(lon1, u.deg, 'Longitude in ground station 1')
    a1 = prepare_quantity(a1, u.dB, 'Attenuation margin in ground station 1')

    lon2 = np.mod(lon2, 360)
    lat2 = prepare_quantity(lat2, u.deg, 'Latitude in ground station 2')
    lon2 = prepare_quantity(lon2, u.deg, 'Longitude in ground station 2')
    a2 = prepare_quantity(a2, u.dB, 'Attenuation margin in ground station 2')

    f = prepare_quantity(f, u.GHz, 'Frequency')
    tau = prepare_quantity(tau, u.one, 'Polarization tilt angle')
    el1 = prepare_quantity(el1, u.deg, 'Elevation angle in ground station 1')
    el2 = prepare_quantity(el2, u.deg, 'Elevation angle in ground station 2')
    hs1 = prepare_quantity(
        hs1, u.km, 'Altitude over the sea level for ground station 1')
    hs2 = prepare_quantity(
        hs2, u.km, 'Altitude over the sea level for ground station 2')

    val = __model.site_diversity_rain_outage_probability(
        lat1, lon1, a1, lat2, lon2, a2, f, el1, el2, tau=tau, hs1=hs1, hs2=hs2)

    return prepare_output_array(val, type_output) * u.pct


def scintillation_attenuation(lat, lon, f, el, p, D, eta=0.5, T=None,
                              H=None, P=None, hL=1000):
    """
    Calculation of monthly and long-term statistics of amplitude scintillations
    at elevation angles greater than 5° and frequencies up to 20 GHz.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    f : number or Quantity
        Frequency (GHz)
    el : sequence, or number
        Elevation angle (degrees)
    p : number
        Percentage of the time the scintillation attenuation value is exceeded.
    D: number
        Physical diameter of the earth-station antenna (m)
    eta: number, optional
        Antenna efficiency. Default value 0.5 (conservative estimate)
    T: number, sequence, or numpy.ndarray, optional
        Average surface ambient temperature (°C) at the site. If None, uses the
        ITU-R P.453 to estimate the wet term of the radio refractivity.
    H: number, sequence, or numpy.ndarray, optional
        Average surface relative humidity (%) at the site. If None, uses the
        ITU-R P.453 to estimate the wet term of the radio refractivity.
    P: number, sequence, or numpy.ndarray, optional
        Average surface pressure (hPa) at the site. If None, uses the
        ITU-R P.453 to estimate the wet term of the radio refractivity.
    hL : number, optional
        Height of the turbulent layer (m). Default value 1000 m


    Returns
    -------
    attenuation: Quantity
        Attenuation due to scintillation (dB)


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    Earth-space telecommunication systems:
    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.618-12-201507-I!!PDF-E.pdf
    """
    type_output = get_input_type(lat)

    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)

    lon = np.mod(lon, 360)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(prepare_input_array(el), u.deg, 'Elevation angle')
    D = prepare_quantity(D, u.m, 'Antenna diameter')
    eta = prepare_quantity(eta, u.one, 'Antenna efficiency')
    T = prepare_quantity(T, u.deg_C, 'Average surface temperature')
    H = prepare_quantity(H, u.percent, 'Average surface relative humidity')
    P = prepare_quantity(P, u.hPa, 'Average surface pressure')
    hL = prepare_quantity(hL, u.m, 'Height of the turbulent layer')

    val = __model.scintillation_attenuation(
        lat, lon, f, el, p, D, eta=eta, T=T, H=H, P=P, hL=hL)
    
    # The values of attenuation cannot be negative. The ITU models end up
    # giving out negative values for certain inputs
    val[val < 0] = 0

    return prepare_output_array(val, type_output) * u.dB


def scintillation_attenuation_sigma(lat, lon, f, el, p, D, eta=0.5, T=None,
                                    H=None, P=None, hL=1000):
    """
    Calculation of the standard deviation of the amplitude of the
    scintillations attenuation at elevation angles greater than 5° and
    frequencies up to 20 GHz.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    f : number or Quantity
        Frequency (GHz)
    el : sequence, or number
        Elevation angle (degrees)
    p : number
        Percentage of the time the scintillation attenuation value is exceeded.
    D: number
        Physical diameter of the earth-station antenna (m)
    eta: number, optional
        Antenna efficiency. Default value 0.5 (conservative estimate)
    T: number, sequence, or numpy.ndarray, optional
        Average surface ambient temperature (°C) at the site. If None, uses the
        ITU-R P.453 to estimate the wet term of the radio refractivity.
    H: number, sequence, or numpy.ndarray, optional
        Average surface relative humidity (%) at the site. If None, uses the
        ITU-R P.453 to estimate the wet term of the radio refractivity.
    P: number, sequence, or numpy.ndarray, optional
        Average surface pressure (hPa) at the site. If None, uses the
        ITU-R P.453 to estimate the wet term of the radio refractivity.
    hL : number, optional
        Height of the turbulent layer (m). Default value 1000 m


    Returns
    -------
    attenuation: Quantity
        Attenuation due to scintillation (dB)


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    Earth-space telecommunication systems:
    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.618-12-201507-I!!PDF-E.pdf
    """
    type_output = get_input_type(lat)

    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)

    lon = np.mod(lon, 360)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(prepare_input_array(el), u.deg, 'Elevation angle')
    D = prepare_quantity(D, u.m, 'Antenna diameter')
    eta = prepare_quantity(eta, u.one, 'Antenna efficiency')
    T = prepare_quantity(T, u.deg_C, 'Average surface temperature')
    H = prepare_quantity(H, u.percent, 'Average surface relative humidity')
    P = prepare_quantity(P, u.hPa, 'Average surface pressure')
    hL = prepare_quantity(hL, u.m, 'Height of the turbulent layer')

    val = __model.scintillation_attenuation_sigma(
        lat, lon, f, el, p, D, eta=eta, T=T, H=H, P=P, hL=hL)

    return prepare_output_array(val, type_output) * u.dB


def rain_cross_polarization_discrimination(Ap, f, el, p, tau=45):
    """
    Calculation of the cross-polarization discrimination (XPD) statistics from
    rain attenuation statistics. The following procedure provides estimates of
    the long-term statistics of the cross-polarization discrimination (XPD)
    statistics for frequencies up to 55 GHz and elevation angles lower than 60
    deg.


    Parameters
    ----------
    Ap : number, sequence, or numpy.ndarray
        Rain attenuation (dB) exceeded for the required percentage of time, p,
        for the path in question, commonly called co-polar attenuation (CPA)
    f : number
        Frequency
    el : number, sequence, or numpy.ndarray
        Elevation angle (degrees)
    p : number
        Percentage of the time the XPD is exceeded.
    tau : number, optional
        Polarization tilt angle relative to the horizontal (degrees)
        (tau = 45 deg for circular polarization). Default value is 45


    Returns
    -------
    attenuation: Quantity
        Cross-polarization discrimination (dB)


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    Earth-space telecommunication systems:
    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.618-12-201507-I!!PDF-E.pdf
    """
    type_output = get_input_type(Ap)
    Ap = prepare_input_array(Ap)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    tau = prepare_quantity(tau, u.one, 'Polarization tilt angle')
    val = __model.rain_cross_polarization_discrimination(Ap, f, el, p, tau=tau)
    return prepare_output_array(val, type_output) * u.dB


def fit_rain_attenuation_to_lognormal(lat, lon, f, el, hs, P_k, tau):
    """
    Compute the log-normal fit of rain attenuation vs. probability of
    occurrence.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    f : number or Quantity
        Frequency (GHz)
    el : sequence, or number
        Elevation angle (degrees)
    hs : number, sequence, or numpy.ndarray, optional
        Heigh above mean sea level of the earth station (km). If local data for
        the earth station height above mean sea level is not available, an
        estimate is obtained from the maps of topographic altitude
        given in Recommendation ITU-R P.1511.
    P_k : number
        Rain probability
    tau : number, optional
        Polarization tilt angle relative to the horizontal (degrees)
        (tau = 45 deg for circular polarization). Default value is 45

    Returns
    -------
    sigma_lna:
        Standar deviation of the lognormal distribution
    m_lna:
        Mean of the lognormal distribution

    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    Earth-space telecommunication systems:
    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.618-12-201507-I!!PDF-E.pdf

    """
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(prepare_input_array(el), u.deg, 'Elevation angle')
    hs = prepare_quantity(
        hs, u.km, 'Heigh above mean sea level of the earth station')
    tau = prepare_quantity(tau, u.one, 'Polarization tilt angle')
    sigma_lna, m_lna = __model.fit_rain_attenuation_to_lognormal(
        lat, lon, f, el, hs, P_k, tau)
    return sigma_lna, m_lna
