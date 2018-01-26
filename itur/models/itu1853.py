# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import scipy.stats as stats
from scipy.signal import lfilter

from itu618 import fit_rain_attenuation_to_lognormal
from itu840 import lognormal_approximation_coefficient, specific_attenuation_coefficients
from itu836 import total_water_vapour_content
from itu837 import rain_percentage_probability
from itu676 import zenit_water_vapour_attenuation
from itu1510 import surface_mean_temperature
from itur.utils import load_data, dataset_dir, prepare_input_array,\
    prepare_output_array, prepare_quantity

from astropy import units as u


class __ITU1853():
    """Tropospheric attenuation time series synthesis

    Available versions include:
    * P.1853-1 (08/94) (Superseded)
    * P.1853-2 (10/99) (Superseded)
    * P.1853-3 (02/01) (Superseded)
    * P.1853-4 (04/03) (Superseded)
    * P.1853-5 (08/07) (Superseded)
    * P.1853-6 (02/12) (Current version)
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.1853 recommendation.

    def __init__(self, version=1):
        if version == 1:
            self.instance = _ITU1853_1()
        elif version == 0:
            self.instance = _ITU1853_0()
        else:
            raise ValueError(
                'Version ' +
                str(version) +
                ' is not implemented' +
                ' for the ITU-R P.1853 model.')

    @property
    def __version__(self):
        return self.instance.__version__

    def rain_attenuation_synthesis(self, lat, lon, f, el, hs, T, Ts=1, n=None):
        return self.instance.rain_attenuation_synthesis(lat, lon, f,
                                                        el, hs, T, Ts=Ts, n=n)


class _ITU1853_1():

    def __init__(self):
        self.__version__ = 1
        self.year = 2012
        self.month = 2
        self.link = 'https://www.itu.int/rec/R-REC-P.1853-1-201202-I/en'

    def rain_attenuation_synthesis(self, lat, lon, f, el, hs, T, Ts=1, n=None):
        """
        For Earth-space paths, the time series synthesis method is valid for
        frequencies between 4 GHz and 55 GHz and elevation angles between
        5 deg and 90 deg.
        """

        P_rain = rain_percentage_probability(lat, lon).\
            to(u.dimensionless_unscaled).value

        sigma, m = fit_rain_attenuation_to_lognormal(
            lat, lon, f, el, hs, P_rain * 100)

        # Step B: Set the low-pass filter parameter
        beta = 2e-5
        # Step C: compute the attenuation offset
        A_offset = np.exp(m + sigma * stats.norm.ppf(P_rain))
        # Step D: Time series synthesis
        # D1: Synthesize a white Gaussian noise time series
        if n is None:
            n = np.random.normal(0, 1, (T * Ts + 2e5))

        # D2, D3 : Filter the noise time series with a recursive low-pass
        # filter
        rho = np.exp(-beta * Ts)
        X = lfilter([np.sqrt(1 - rho**2)], [1, -rho], n, 0)
        # D4: Compute Y_rain
        Y_rain = np.exp(m + sigma * X)
        # D5: Compute Arain
        A_rain = np.maximum(Y_rain - A_offset, 0)

        # D6: Discard the first 200 000 samples from the synthesized
        A_rain = A_rain[200000::Ts]

        return A_rain

    def fftnoise(f):
        f = np.array(f, dtype='complex')
        Np = (len(f) - 1) // 2
        phases = np.random.rand(Np) * 2 * np.pi
        phases = np.cos(phases) + 1j * np.sin(phases)
        f[1:Np + 1] *= phases
        f[-1:-1 - Np:-1] = np.conj(f[1:Np + 1])
        return np.fft.ifft(f).real

    def scintillation_synthesis(self, T, f_c=0.1, Ts=1, f_c2=0.01):
        """
        For Earth-space paths, the time series synthesis method is valid for
        frequencies between 4 GHz and 55 GHz and elevation angles between
        5 deg and 90 deg.
        """
        freqs = np.abs(np.fft.fftfreq(2 * int(T + 2e5), 1 / Ts))
        H_f = np.where(freqs <= f_c, 1, 10 **
                       ((np.log10(freqs) - np.log10(f_c)) * (-8 / 3)))
        H_f = H_f[0:int(T + 2e5)]
        sci = self.fftnoise(np.fft.fftshift(H_f))
        return sci

    def cloud_liquid_water_synthesis(self, lat, lon, T, Ts=1, n=None):
        """
        """

        # Step A: Estimation of m, sigma and Pcwl
        m, sigma, Pcwl = lognormal_approximation_coefficient(lat, lon)
        m = m.value
        sigma = sigma.value
        Pcwl = Pcwl.value

        # Step B: Low pass filter parameters
        beta_1 = 7.17e-4
        beta_2 = 2.01e-5
        gamma_1 = 0.349
        gamma_2 = 0.830

        # Step C: Truncation threshold
        alpha = stats.norm.ppf(1 - Pcwl)

        # Step D: Time series synthesis
        # Step D1: Synthesize a white Gaussian noise time series
        if n is None:
            n = np.random.rand((T + 5e5) // Ts)
        # Step D3: Filter the noise time series, with two recursive low-pass
        # filters
        rho_1 = np.exp(-beta_1 * Ts)
        X_1 = lfilter([np.sqrt(1 - rho_1**2)], [1, -rho_1], n, 0)
        rho_2 = np.exp(-beta_2 * Ts)
        X_2 = lfilter([np.sqrt(1 - rho_2**2)], [1, -rho_2], n, 0)
        # Step D4: Compute Gc(kTs),
        G_c = gamma_1 * X_1 + gamma_2 * X_2
        # Step D5: Compute L(kTs) (dB)
        L = np.where(G_c > alpha, np.exp(m + sigma * stats.norm.ppf(
            1 - 1 / Pcwl * stats.norm.sf(G_c))), 0)

        # D6: Discard the first 500 000 samples from the synthesized
        L = L[5e5:]

        return L

    def integrated_water_vapour_synthesis(self, lat, lon, T, Ts=1, n=None):
        # A Estimation of κ and λ
        ps = np.array([0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50])
        Vi = np.array([total_water_vapour_content(lat, lon, p_i)
                       for p_i in ps])

        ln_lnPi = np.log(- np.log(ps))
        ln_Vi = np.log(Vi)

        a, b = np.linalg.lstsq(np.vstack([ln_Vi, np.ones(len(ln_Vi))]).T,
                               ln_lnPi)[0]
        kappa = a
        lambd = np.exp(-b / a)

        # B Low-pass filter parameter
        beta_V = 3.24e-6

        # Step C: Time series synthesis
        # Step C1: Synthesize a white Gaussian noise time series
        if n is None:
            n = np.random.rand((T + 5e6) // Ts)
        # Step C3: Filter the noise time series, with two recursive low-pass
        # filters
        rho = np.exp(-beta_V * Ts)
        G_v = lfilter([np.sqrt(1 - rho**2)], [1, -rho], n, 0)
        # Step C4: Compute Compute V(kTs),
        V = lambd * (- np.log10(stats.norm.sf(G_v)))**(1 / kappa)
        # Step C5: Discard the first 5 000 000 samples from the synthesized
        V = V[5e6:]

        return V

    def total_attenuation_synthesis(self, lat, lon, f, el, p, hs, T, Ts=1):
        # Step A Correlation coefficients:
        C_RC = 1
        C_CV = 0.8

        # Step B Scintillation polynomials
        def a_Fade(p):
            return -0.061 * np.log10(p)**3 + 0.072 * \
                np.log10(p)**2 - 1.71 * np.log10(p) + 3

        def a_Enhanc(p):
            return -0.0597 * np.log10(p)**3 - 0.0835 * \
                np.log10(p)**2 - 1.258 * np.log10(p) + 2.672

        # Step C1-C3:
        n_R = np.random.rand((T + 5e6) // Ts)
        n_L0 = np.random.rand((T + 5e6) // Ts)
        n_V0 = np.random.rand((T + 5e6) // Ts)
        # Step C4-C5:
        n_L = C_RC * n_R + np.sqrt(1 - C_RC**2) * n_L0
        n_V = C_CV * n_L + np.sqrt(1 - C_CV**2) * n_V0

        # Step C6: Compute the rain attenuation time series
        Ar = self.rain_attenuation_synthesis(
            lat, lon, f, el, hs, T + 5e6, n=n_R)
        Ar = Ar[5e6:]

        # Step C7: Compute the cloud integrated liquid water content time
        # series
        L = self.cloud_liquid_water_synthesis(lat, lon, T + 5e6, n=n_L)
        Ac = L * \
            specific_attenuation_coefficients(f, T=0) / np.sin(np.deg2rad(el))

        # Step C9: Identify time stamps where A_R > 0 L > 1
        idx = np.where(Ar > 0 & L > 1)[0]
        idx_no = np.where(Ar <= 0 | L <= 1)[0]
        # Step C10: Discard the previous values of Ac and re-compute them by
        # linear interpolation vs. time starting from the non-discarded cloud
        # attenuations values
        Ac[idx] = np.interp(idx, idx_no, Ac[idx_no])

        # Step C11: Compute the integrated water vapour content time series
        V = self.integrated_water_vapour_synthesis(lat, lon, T, n=n_V)

        # Step C12: Convert the integrated water vapour content time series
        # V into water vapour attenuation time series AV(kTs)
        Av = zenit_water_vapour_attenuation(lat, lon, p, f, V_t=V)

        # Step C13: Compute the mean annual temperature Tm for the location of
        # interest using experimental values if available.
        Tm = surface_mean_temperature(lat, lon)

        # Step C14: Convert the mean annual temperature Tm into mean annual
        # oxygen attenuation AO following the method recommended in
        # Recommendation ITU-R P.676.
        Ao = 1

        # Step C15: Synthesize unit variance scintillation time series
        sci_0 = self.scintillation_synthesis(T)
        # Step C16: Compute the correction coefficient time series Cx(kTs) in
        # order to distinguish between scintillation fades and enhancements:
        Q_sci = 100 * stats.norm.sf(sci_0)
        C_x = np.where(sci_0 > 0, a_Fade(Q_sci) / a_Enhanc(Q_sci), 1)

        # Step C17: Transform the integrated water vapour content time series
        # V(kTs) into the Gamma distributed time series Z(kTs) as follows:
        Z = scipy.stats.gamma.ppf(1 - np.exp(-(V / lambd)**kappa), 10, 0.1)

        # Step C18: Compute the scintillation standard deviation σ following
        # the method recommended in Recommendation ITU-R P.618.
        sigma = 2

        # Step C19: Compute the scintillation time series sci:
        sci = np.where(Ar > 1, sigma * sci_0 * C_x * Z * Ar ** (5 / 12),
                       sigma * sci_0 * C_x * Z)

        # Step C20: Compute total tropospheric attenuation time series A(kTs)
        # as follows:
        A = Ar + Ac + Av + Ao + sci

        return A


__model = __ITU1853()


def change_version(new_version):
    """
    Change the version of the ITU-R P.1853 recommendation currently being used.


    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
        * P.1853-1 (08/94) (Superseded)
        * P.1853-2 (10/99) (Superseded)
        * P.1853-3 (02/01) (Superseded)
        * P.1853-4 (04/03) (Superseded)
        * P.1853-5 (08/07) (Superseded)
        * P.1853-6 (02/12) (Current version)
    """
    global __model
    __model = __ITU1853(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.1853 recommendation currently being used.
    """
    global __model
    return __model.__version__


def rain_attenuation_synthesis(lat, lon, f, el, hs, T, Ts=1, n=None):
    """
    A method to compute the rainfall rate exceeded for p% of the average year


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
    T : int
        Number of samples
    Ts : int
        Time step between consecutive samples (seconds)
    n : list, np.array, optional
        Additive White Gaussian Noise used as input for the


    Returns
    -------
    rain_att: numpy.ndarray
        Synthesized rain attenuation time series (dB)


    References
    ----------
    [1] Characteristics of precipitation for propagation modelling
    https://www.itu.int/rec/R-REC-P.1853/en
    """
    global __model

    lon = np.mod(lon, 360)
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    hs = prepare_quantity(
        hs, u.km, 'Heigh above mean sea level of the earth station')
    Ts = prepare_quantity(f, u.second, 'Time step between samples')
    val = __model.rain_attenuation_synthesis(lat, lon, f, el, hs, T, Ts, n)
    return val * u.dB


def unavailability_from_rainfall_rate(lat, lon, R):
    """
    A method to estimate the percentage of time of the average year that a given
    rainfall rate (R) is exceeded. This method calls successively to the
    `rainfall_rate` method and interpolates its value.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    R : number, sequence, or numpy.ndarray
        Rainfall rate (mm/h)


    Returns
    -------
    p: numpy.ndarray
        Rainfall rate exceeded for p% of the average year


    References
    ----------
    [1] Characteristics of precipitation for propagation modelling
    https://www.itu.int/rec/R-REC-P.1853/en
    """
    global __model

    # TODO: write this function
