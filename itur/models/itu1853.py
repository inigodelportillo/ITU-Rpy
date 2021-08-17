# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import scipy.stats as stats
from scipy.signal import lfilter

from itur.models.itu618 import rain_attenuation, scintillation_attenuation_sigma
from itur.models.itu676 import gamma0_exact, \
    slant_inclined_path_equivalent_height
from itur.models.itu840 import lognormal_approximation_coefficient, \
    specific_attenuation_coefficients
from itur.models.itu835 import standard_pressure, standard_water_vapour_density
from itur.models.itu836 import total_water_vapour_content
from itur.models.itu837 import rainfall_probability
from itur.models.itu676 import zenit_water_vapour_attenuation
from itur.models.itu1510 import surface_mean_temperature
from itur.models.itu1511 import topographic_altitude
from itur.utils import prepare_quantity

from astropy import units as u


class __ITU1853():

    """Tropospheric attenuation time series synthesis

    Available versions include:
    * P.1853-0 (10/09) (Superseded)
    * P.1853-1 (02/12) (Current version)
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
                'Version {0} is not implemented for the ITU-R P.1853 model.'
                .format(version))

    @staticmethod
    def set_seed(seed):
        np.random.seed(seed)

    @property
    def __version__(self):
        return self.instance.__version__

    def rain_attenuation_synthesis(self, lat, lon, f, el, hs, Ns,
                                   Ts=1, tau=45, n=None):
        return self.instance.rain_attenuation_synthesis(
            lat, lon, f, el, hs, Ns, Ts=Ts, tau=tau, n=n)

    def total_attenuation_synthesis(self, lat, lon, f, el, p, D, Ns, Ts=1,
                                    tau=45, hs=None, eta=0.65, rho=None,
                                    H=None, P=None, hL=1000,
                                    return_contributions=False):
        return self.instance.total_attenuation_synthesis(
            lat, lon, f, el, p, D, Ns, Ts=Ts, tau=tau, hs=hs, eta=eta, rho=rho,
            H=H, P=P, hL=hL, return_contributions=return_contributions)

    def scintillation_attenuation_synthesis(self, Ns, f_c=0.1, Ts=1):
        return self.instance.scintillation_attenuation_synthesis(
            Ns, f_c=f_c, Ts=Ts)

    def cloud_liquid_water_synthesis(self, lat, lon, Ns, Ts=1, n=None):
        return self.instance.cloud_liquid_water_synthesis(
            lat, lon, Ns, Ts=Ts, n=n)

    def integrated_water_vapour_synthesis(self, lat, lon, Ns, Ts=1, n=None):
        return self.instance.integrated_water_vapour_synthesis(
            lat, lon, Ns, Ts=Ts, n=n)


class _ITU1853_1():

    def __init__(self):
        self.__version__ = 1
        self.year = 2012
        self.month = 2
        self.link = 'https://www.itu.int/rec/R-REC-P.1853-1-201202-I/en'

    @staticmethod
    def rain_attenuation_synthesis(
            lat, lon, f, el, hs, Ns, tau=45, Ts=1, n=None):
        """
        For Earth-space paths, the time series synthesis method is valid for
        frequencies between 4 GHz and 55 GHz and elevation angles between
        5 deg and 90 deg.
        """

        # Step A1: Determine Prain (% of time), the probability of rain on the
        # path. Prain can be well approximated as P0(lat, lon)
        P_rain = rainfall_probability(lat, lon).\
            to(u.dimensionless_unscaled).value

        # Step A2: Construct the set of pairs [Pi, Ai] where Pi (% of time) is
        # the probability the attenuation Ai (dB) is exceeded where Pi < P_K
        p_i = np.array([0.01, 0.02, 0.03, 0.05,
                        0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10])
        Pi = np.array([p for p in p_i if p < P_rain * 100], dtype=float)
        Ai = np.array([0 for p in p_i if p < P_rain * 100], dtype=float)

        for i, p in enumerate(Pi):
            Ai[i] = rain_attenuation(lat, lon, f, el, hs, p, tau=tau).value

        # Step A3: Transform the set of pairs [Pi, Ai] to [Q^{-1}(Pi/P_k),
        # ln(Ai)]
        Q = stats.norm.ppf((Pi / 100))
        lnA = np.log(Ai)

        # Step A4: Determine the variables sigma_lna, m_lna by performing a
        # least-squares fit to lnAi = sigma_lna Q^{-1}(Pi/P_k) + m_lna
        m, sigma = np.linalg.lstsq(np.vstack([np.ones(len(Q)), Q]).T,
                                   lnA, rcond=None)[0]

        # Step B: Set the low-pass filter parameter
        beta = 2e-4
        # Step C: compute the attenuation offset
        A_offset = np.exp(m + sigma * stats.norm.ppf(P_rain))
        # Step D: Time series synthesis
        # D1: Synthesize a white Gaussian noise time series
        if n is None:
            n = np.random.normal(0, 1, int(Ns * Ts + 2e5))[::Ts]
            discard_samples = True
        else:
            discard_samples = False

        # D2, D3 : Filter the noise time series with a recursive low-pass
        # filter
        rho = np.exp(-beta * Ts)
        X = lfilter([np.sqrt(1 - rho**2)], [1, -rho], n, 0)
        # D4: Compute Y_rain
        Y_rain = np.exp(m + sigma * X)
        # D5: Compute Arain
        A_rain = np.maximum(Y_rain - A_offset, 0)

        # D6: Discard the first 200 000 samples from the synthesized
        if discard_samples:
            A_rain = A_rain[np.ceil(200000 / Ts).astype(int):]

        return A_rain.flatten()

    @classmethod
    def fftnoise(cls, f):
        f = np.array(f, dtype='complex')
        Np = (len(f) - 1) // 2
        phases = np.random.rand(Np) * 2 * np.pi
        phases = np.cos(phases) + 1j * np.sin(phases)
        f[1:Np + 1] *= phases
        f[-1:-1 - Np:-1] = np.conj(f[1:Np + 1])
        return np.fft.ifft(f).real

    @staticmethod
    def scintillation_attenuation_synthesis(Ns, f_c=0.1, Ts=1):
        """
        For Earth-space paths, the time series synthesis method is valid for
        frequencies between 4 GHz and 55 GHz and elevation angles between
        5 deg and 90 deg.
        """
        freqs = np.abs(np.fft.fftfreq(2 * int(Ns + 2e5), 1 / Ts))
        H_f = np.where(freqs <= f_c, 1, 10 **
                       ((np.log10(freqs) - np.log10(f_c)) * (-8 / 3)))
        H_f = H_f[0:int(Ns + 2e5)]
        sci = _ITU1853_1.fftnoise(np.fft.fftshift(H_f))
        return sci[200000:].flatten()

    @staticmethod
    def cloud_liquid_water_synthesis(lat, lon, Ns, Ts=1, n=None):
        """
        """

        # Step A: Estimation of m, sigma and Pcwl
        m, sigma, Pcwl = lognormal_approximation_coefficient(lat, lon)
        m = m.value
        sigma = sigma.value
        Pcwl = Pcwl.value / 100

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
            n = np.random.normal(0, 1, int(Ns * Ts + 5e5))[::Ts]
            discard_samples = True
        else:
            discard_samples = False

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
        if discard_samples:
            L = L[np.ceil(500000 / Ts).astype(int):]

        return L.flatten()

    @staticmethod
    def integrated_water_vapour_coefficients(lat, lon):
        # A Estimation of κ and λ
        ps = np.array([0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50])
        Vi = np.array([total_water_vapour_content(lat, lon, p_i).value
                       for p_i in ps])

        ln_lnPi = np.log(- np.log(ps / 100))
        ln_Vi = np.log(Vi)

        a, b = np.linalg.lstsq(np.vstack([ln_Vi, np.ones(len(ln_Vi))]).T,
                               ln_lnPi, rcond=None)[0]
        kappa = a
        lambd = np.exp(-b / a)
        return kappa, lambd

    def integrated_water_vapour_synthesis(self, lat, lon, Ns, Ts=1, n=None):
        # A Estimation of κ and λ
        kappa, lambd = self.integrated_water_vapour_coefficients(lat, lon)

        # B Low-pass filter parameter
        beta_V = 3.24e-6

        # Step C: Time series synthesis
        # Step C1: Synthesize a white Gaussian noise time series
        if n is None:
            n = np.random.normal(0, 1, (Ns * Ts + 5000000))[::Ts]
            discard_samples = True
        else:
            discard_samples = False

        # Step C3: Filter the noise time series, with two recursive low-pass
        # filters
        rho = np.exp(-beta_V * Ts)
        G_v = lfilter([np.sqrt(1 - rho**2)], [1, -rho], n, 0)
        # Step C4: Compute Compute V(kTs),
        V = lambd * (- np.log10(stats.norm.sf(G_v)))**(1 / kappa)
        # Step C5: Discard the first 5 000 000 samples from the synthesized
        if discard_samples:
            V = V[np.ceil(5000000 / Ts).astype(int):]

        return V.flatten()

    def total_attenuation_synthesis(self, lat, lon, f, el, p, D, Ns, Ts=1,
                                    hs=None, tau=45, eta=0.65, rho=None,
                                    H=None, P=None, hL=1000,
                                    return_contributions=False):
        t_disc = int(5e6)
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
        n_R = np.random.normal(0, 1, int((Ns * Ts + t_disc)))
        n_L0 = np.random.normal(0, 1, int((Ns * Ts + t_disc)))
        n_V0 = np.random.normal(0, 1, int((Ns * Ts + t_disc)))

        # Step C4-C5:
        n_L = C_RC * n_R + np.sqrt(1 - C_RC**2) * n_L0
        n_V = C_CV * n_L + np.sqrt(1 - C_CV**2) * n_V0

        # Step C6: Compute the rain attenuation time series
        if hs is None:
            hs = topographic_altitude(lat, lon)
        Ar = self.rain_attenuation_synthesis(
            lat, lon, f, el, hs, Ns, Ts=1, tau=tau, n=n_R)
        Ar = Ar[t_disc:]

        # Step C7: Compute the cloud integrated liquid water content time
        # series
        L = self.cloud_liquid_water_synthesis(lat, lon, Ns, Ts=1, n=n_L)
        L = L[t_disc:]
        Ac = L * \
            specific_attenuation_coefficients(f, T=0) / np.sin(np.deg2rad(el))
        Ac = Ac.flatten()

        # Step C9: Identify time stamps where A_R > 0 L > 1
        idx = np.where(np.logical_and(Ar > 0, L > 1))[0]
        idx_no = np.where(np.logical_not(
            np.logical_and(Ar > 0, L > 1)))[0]

        # Step C10: Discard the previous values of Ac and re-compute them by
        # linear interpolation vs. time starting from the non-discarded cloud
        # attenuations values
        Ac[idx] = np.interp(idx, idx_no, Ac[idx_no])

        # Step C11: Compute the integrated water vapour content time series
        V = self.integrated_water_vapour_synthesis(lat, lon, Ns, Ts=1, n=n_V)
        V = V[t_disc:]

        # Step C12: Convert the integrated water vapour content time series
        # V into water vapour attenuation time series AV(kTs)
        Av = zenit_water_vapour_attenuation(lat, lon, p, f, V_t=V).value

        # Step C13: Compute the mean annual temperature Tm for the location of
        # interest using experimental values if available.
        Tm = surface_mean_temperature(lat, lon).value

        # Step C14: Convert the mean annual temperature Tm into mean annual
        # oxygen attenuation AO following the method recommended in
        # Recommendation ITU-R P.676.
        if P is None:
            P = standard_pressure(hs).value

        if rho is None:
            rho = standard_water_vapour_density(hs).value

        e = Tm * rho / 216.7
        go = gamma0_exact(f, P, rho, Tm).value
        ho, _ = slant_inclined_path_equivalent_height(f, P + e, rho).value
        Ao = ho * go * np.ones_like(Ar)

        # Step C15: Synthesize unit variance scintillation time series
        sci_0 = self.scintillation_attenuation_synthesis(Ns * Ts, Ts=1)

        # Step C16: Compute the correction coefficient time series Cx(kTs) in
        # order to distinguish between scintillation fades and enhancements:
        Q_sci = 100 * stats.norm.sf(sci_0)
        C_x = np.where(sci_0 > 0, a_Fade(Q_sci) / a_Enhanc(Q_sci), 1)

        # Step C17: Transform the integrated water vapour content time series
        # V(kTs) into the Gamma distributed time series Z(kTs) as follows:
        kappa, lambd = self.integrated_water_vapour_coefficients(lat, lon)
        Z = stats.gamma.ppf(1 - np.exp(-(V / lambd)**kappa), 10, 0.1)

        # Step C18: Compute the scintillation standard deviation σ following
        # the method recommended in Recommendation ITU-R P.618.
        sigma = scintillation_attenuation_sigma(lat, lon, f, el, p, D, eta, Tm,
                                                H, P, hL).value

        # Step C19: Compute the scintillation time series sci:
        As = np.where(Ar > 1, sigma * sci_0 * C_x * Z * Ar ** (5 / 12),
                      sigma * sci_0 * C_x * Z)

        # Step C20: Compute total tropospheric attenuation time series A(kTs)
        # as follows:
        A = Ar + Ac + Av + Ao + As

        if return_contributions:
            return (Ao + Av)[::Ts], Ac[::Ts], Ar[::Ts], As[::Ts], A[::Ts]
        else:
            return A[::Ts]


class _ITU1853_0():

    def __init__(self):
        self.__version__ = 0
        self.year = 2009
        self.month = 10
        self.link = 'https://www.itu.int/rec/R-REC-P.1853-0-200910-I/en'

    @staticmethod
    def rain_attenuation_synthesis(*args, **kwargs):
        return _ITU1853_1.rain_attenuation_synthesis(*args, **kwargs)

    @staticmethod
    def scintillation_attenuation_synthesis(*args, **kwargs):
        return _ITU1853_1.scintillation_attenuation_synthesis(*args, **kwargs)

    @staticmethod
    def cloud_liquid_water_synthesis(*args, **kwargs):
        raise NotImplementedError(
            "Recommendation ITU-R P.1853 does not specify a method to compute "
            "time series for the cloud liquid water content.")

    @staticmethod
    def integrated_water_vapour_synthesis(*args, **kwargs):
        raise NotImplementedError(
            "Recommendation ITU-R P.1853 does not specify a method to compute "
            "time series for the water vapour content.")

    @staticmethod
    def total_attenuation_synthesis(*args, **kwargs):
        raise NotImplementedError(
            "Recommendation ITU-R P.1853 does not specify a method to compute "
            "time series for the total atmospheric attenuation.")


__model = __ITU1853()


def change_version(new_version):
    """
    Change the version of the ITU-R P.1853 recommendation currently being used.


    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          * 1:  Activates recommendation ITU-R P.1853-1 (02/12) (Current version)
          * 0:  Activates recommendation ITU-R P.1853-0 (10/09) (Superseded)
    """
    global __model
    __model = __ITU1853(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.1853 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.

    """
    global __model
    return __model.__version__


def set_seed(seed):
    """
    Set the seed used to generate random numbers.

    Parameters
    ----------
    seed : int
        Seed used to generate random numbers
    """
    __model.set_seed(seed)


def rain_attenuation_synthesis(lat, lon, f, el, hs, Ns, Ts=1, tau=45, n=None):
    """
    A method to generate a synthetic time series of rain attenuation values.


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
    Ns : int
        Number of samples
    Ts : int
        Time step between consecutive samples (seconds)
    tau : number, optional
        Polarization tilt angle relative to the horizontal (degrees)
        (tau = 45 deg for circular polarization). Default value is 45
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
    Ts = prepare_quantity(Ts, u.second, 'Time step between samples')
    val = __model.rain_attenuation_synthesis(lat, lon, f, el, hs, Ns,
                                             Ts=Ts, tau=tau, n=n)
    return val * u.dB


def scintillation_attenuation_synthesis(Ns, f_c=0.1, Ts=1):
    """
    A method to generate a synthetic time series of scintillation attenuation
    values.


    Parameters
    ----------
    Ns : int
        Number of samples
    f_c : float
        Cut-off frequency for the low pass filter
    Ts : int
        Time step between consecutive samples (seconds)

    Returns
    -------
    sci_att: numpy.ndarray
        Synthesized scintilation attenuation time series (dB)


    References
    ----------
    [1] Characteristics of precipitation for propagation modelling
    https://www.itu.int/rec/R-REC-P.1853/en
    """
    global __model

    val = __model.scintillation_attenuation_synthesis(Ns, f_c, Ts)
    return val * u.dB


def integrated_water_vapour_synthesis(lat, lon, Ns, Ts=1, n=None):
    """ The time series synthesis method generates a time series that
    reproduces the spectral characteristics and the distribution of water
    vapour content.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    Ns : int
        Number of samples
    Ts : int
        Time step between consecutive samples (seconds)
    n : list, np.array, optional
        Additive White Gaussian Noise used as input for the

    Returns
    -------
    L: numpy.ndarray
        Synthesized water vapour content time series (kg/m2)


    References
    ----------
    [1] Characteristics of precipitation for propagation modelling
    https://www.itu.int/rec/R-REC-P.1853/en
    """
    global __model

    lon = np.mod(lon, 360)
    val = __model.integrated_water_vapour_synthesis(lat, lon, Ns, Ts, n)
    return val * u.kg / u.m**2


def cloud_liquid_water_synthesis(lat, lon, Ns, Ts=1, n=None):
    """ The time series synthesis method generates a time series that
    reproduces the spectral characteristics, rate of change and duration
    statistics of cloud liquid content events.

    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    Ns : int
        Number of samples
    Ts : int
        Time step between consecutive samples (seconds)
    n : list, np.array, optional
        Additive White Gaussian Noise used as input for the

    Returns
    -------
    V: numpy.ndarray
        Synthesized cloud liquid water time series (mm)


    References
    ----------
    [1] Characteristics of precipitation for propagation modelling
    https://www.itu.int/rec/R-REC-P.1853/en
    """
    global __model

    lon = np.mod(lon, 360)
    val = __model.cloud_liquid_water_synthesis(lat, lon, Ns, Ts, n)
    return val * u.mm


def total_attenuation_synthesis(lat, lon, f, el, p, D, Ns, Ts=1, hs=None,
                                tau=45, eta=0.65, rho=None, H=None, P=None,
                                hL=1000, return_contributions=False):
    """ The time series synthesis method generates a time series that
    reproduces the spectral characteristics, rate of change and duration
    statistics of the total atmospheric attenuation events.

    The time series is obtained considering the contributions of gaseous,
    cloud, rain, and scintillation attenuation.

    Parameters
    ----------
    lat : number
        Latitudes of the receiver points
    lon : number
        Longitudes of the receiver points
    f : number or Quantity
        Frequency (GHz)
    el : number
        Elevation angle (degrees)
    p : number
        Percentage of the time the rain attenuation value is exceeded.
    D: number or Quantity
        Physical diameter of the earth-station antenna (m)
    Ns : int
        Number of samples
    Ts : int
        Time step between consecutive samples (seconds)
    tau : number, optional
        Polarization tilt angle relative to the horizontal (degrees)
        (tau = 45 deg for circular polarization). Default value is 45
    hs : number, sequence, or numpy.ndarray, optional
        Heigh above mean sea level of the earth station (km). If local data for
        the earth station height above mean sea level is not available, an
        estimate is obtained from the maps of topographic altitude
        given in Recommendation ITU-R P.1511.
    eta: number, optional
        Antenna efficiency. Default value 0.5 (conservative estimate)
    rho : number or Quantity, optional
        Water vapor density (g/m3). If not provided, an estimate is obtained
        from Recommendation Recommendation ITU-R P.836.
    H: number, sequence, or numpy.ndarray, optional
        Average surface relative humidity (%) at the site. If None, uses the
        ITU-R P.453 to estimate the wet term of the radio refractivity.
    P: number, sequence, or numpy.ndarray, optional
        Average surface pressure (hPa) at the site. If None, uses the
        ITU-R P.453 to estimate the wet term of the radio refractivity.
    hL : number, optional
        Height of the turbulent layer (m). Default value 1000 m
    return_contributions: bool, optional
        Determines whether individual contributions from gases, rain, clouds
        and scintillation are returned in addition ot the total attenuation
        (True), or just the total atmospheric attenuation (False).
        Default is False

    Returns
    ---------
    A : Quantity
        Synthesized total atmospheric attenuation time series (dB)

    Ag, Ac, Ar, As, A : tuple
        Synthesized Gaseous, Cloud, Rain, Scintillation contributions to total
        attenuation time series, and synthesized total attenuation time seires
        (dB).

    References
    ----------
    [1] Characteristics of precipitation for propagation modelling
    https://www.itu.int/rec/R-REC-P.1853/en
    """
    global __model

    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    D = prepare_quantity(D, u.m, 'Antenna diameter')
    hs = prepare_quantity(
        hs, u.km, 'Heigh above mean sea level of the earth station')
    eta = prepare_quantity(eta, u.one, 'Antenna efficiency')
    rho = prepare_quantity(rho, u.g / u.m**3, 'Water vapor density')
    H = prepare_quantity(H, u.percent, 'Average surface relative humidity')
    P = prepare_quantity(P, u.hPa, 'Average surface pressure')
    hL = prepare_quantity(hL, u.m, 'Height of the turbulent layer')

    val = __model.total_attenuation_synthesis(
        lat, lon, f, el, p, D, Ns, Ts=Ts, tau=tau, hs=hs, eta=eta, rho=rho,
        H=H, P=P, hL=hL, return_contributions=return_contributions)
    if return_contributions:
        return tuple([v * u.dB for v in val])
    else:
        return val * u.dB
