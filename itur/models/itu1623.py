# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import warnings
import numpy as np
from astropy import units as u
from scipy.optimize import fsolve
from scipy.special import erf as erf

from itur.utils import get_input_type, prepare_quantity, prepare_output_array, prepare_input_array


def Qfunc(z):
    """Tail distribution function of the standard normal distribution.

    Q(z) is the probability that a normal (Gaussian) random variable will
    a value larger than z standard deviations

    The Q-function can be expressed in terms of the error function as

    .. math::

        Q(z) = \\frac{1}{2} \\left(1 - erf\\left(\\frac{z}{\\sqrt{2}}\\right)\\right)

    Parameters
    ----------
    z: float
        Value to evaluate Q at.

    Returns
    -------
    q : float
        Value of the Q function evaluated at z.
    """
    return 0.5 * (1 - erf(z / np.sqrt(2)))


class __ITU1623__:
    """Prediction method of fade dynamics on Earth-space paths.

    Available versions:
    * P.1623-0 (03/92) (Superseded)
    * P.1623-1 (83/97) (Current version)

    """

    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.1623 recommendation.

    def __init__(self, version=1):
        if version == 1:
            self.instance = _ITU1623_1_()
        elif version == 0:
            self.instance = _ITU1623_0_()
        else:
            raise ValueError(
                "Version {0} is not implemented for the ITU-R P.1623 model."
                .format(version))

        self._zero_isoterm_data = {}

    @property
    def __version__(self):
        return self.instance.__version__

    def fade_duration(self, D_arr, A, el, f, T_tot):
        # Abstract method to compute the fade duration
        return self.instance.fade_duration(D_arr, A, el, f, T_tot)

    def fade_slope(self, z, A, f_B, delta_t):
        # Abstract method to compute the fade slope
        return self.instance.fade_slope(z, A, f_B, delta_t)

    def fade_depth(self, N_target, D_target, A, PofA, el, f):
        # Abstract method to compute the fade depth
        return self.instance.fade_depth(N_target, D_target, A, PofA, el, f)


class _ITU1623_1_:
    def __init__(self):
        self.__version__ = 1
        self.year = 2005
        self.month = 3
        self.link = "https://www.itu.int/rec/R-REC-P.1623-1-200503-I/en"

    @classmethod
    def fade_duration(self, D_arr, A, el, f, T_tot):
        if np.any(f < 10) or np.any(f > 50):
            warnings.warn(
                RuntimeWarning(
                    'The method to compute fade duration parameters '
                    'in recommendation ITU-P 1623-11 is only '
                    'recommended for frequencies in the 10-50GHz range'))

        if np.any(el < 5) or np.any(el > 60):
            warnings.warn(
                RuntimeWarning(
                    'The method to compute fade duration parameters '
                    'in recommendation ITU-P 1623-11 is only '
                    'recommended for elevation angles in the 5-60deg range'))

        # Step 1: Calculate the mean duration D0 of the log-normal
        # distribution of the fraction of fading time due to fades of long
        # duration, given that the attenuation is greater than A:
        D_0 = 80 * (el ** (-0.4)) * (f ** (1.4)) * (A ** (-0.39))  # seconds

        # Step 2: Calculate the standard deviation σ of the lognormal
        # distribution of the fraction of fading time due to fades of long
        # duration:
        sigma = 1.85 * (f ** (-0.05)) * (A ** (-0.027))

        # Step 3: Calculate the exponent γ of the power-law distribution of
        # the fraction of fading time due to fades of short duration:
        gamma = 0.055 * (f ** 0.65) * (A ** (-0.003))

        # Step 4: Calculate the boundary between short and long fade
        # durations, Dt:
        p_1 = (0.885 * gamma) - 0.814
        p_2 = (-1.05 * (gamma ** 2)) + (2.23 * gamma) - 1.61
        D_t = D_0 * np.exp(p_1 * sigma ** 2 + p_2 * sigma - 0.39)

        # Step 5: Calculate the mean duration D2 of the log-normal distribution
        # of the probability of occurrence of fading events of long duration:
        D_2 = D_0 * np.exp(-(sigma ** 2))

        # Step 6: Calculate the fraction of time k due to fades of duration
        # less than Dt:
        Q_1 = Qfunc((np.log(D_t) - np.log(D_0)) / sigma)
        Q_2 = Qfunc((np.log(D_t) - np.log(D_2)) / sigma)
        k = 1. / (1 + ((np.sqrt(D_0 * D_2) *
                        (1 - gamma) * Q_1) / (D_t * gamma * Q_2)))

        # Step 7: Calculate the probability of occurrence of fade events
        # duration d longer than D given that attenuation a is greater than A:
        p = np.zeros_like(D_arr)  # initializes p for indexing ops.
        Q_ratio_p = (Qfunc(np.log(D_arr / D_2) / sigma) /
                     Qfunc(np.log(D_t / D_2) / sigma))

        p = np.where(
            np.logical_and(D_arr >= 1, D_arr <= D_t),
            D_arr ** -gamma,
            (D_t ** -gamma) * Q_ratio_p,
        )

        # Step 8: Calculate the cumulative probability of exceedance, i.e. the
        # total fraction of fade time due to fades of duration d longer than D:
        F = np.zeros_like(D_arr)
        Q_ratio_F = (Qfunc(np.log(D_arr / D_0) / sigma) /
                     Qfunc(np.log(D_t / D_0) / sigma))  # or divide by Q_2

        F = np.where(
            np.logical_and(D_arr >= 1, D_arr <= D_t),
            (1 - (k * (D_arr / D_t) ** (1 - gamma))),
            ((1 - k) * Q_ratio_F),
        )

        # Step 9: Compute N(D,A), The total number of fades of duration d
        # longer than D for a given threshold A.

        # Step 9a. Compute Ntot, using  the Ttot(A) parameter
        N_tot = T_tot * (k / gamma) * ((1 - gamma) / (D_t ** (1 - gamma)))

        # Compute number of fades N(D,A)
        N = N_tot * p

        # Compute T(d > D|a > A), total fading time due to fades of duration d
        # longer than D for the threshold A:
        T = T_tot * F

        return np.array([p, F, N, T])

    @classmethod
    def fade_slope(self, z, A, f_B, delta_t):
        # Step 1:   Calculate F
        b = 2.3
        F = np.sqrt(
            (2 * np.pi) ** 2 / ((1 / f_B ** b) + (2 * delta_t) ** b) ** (1 / b)
        )   # eq. 18

        # Step 2:   Calculate STD of the conditional fade slope
        # s is a parameter which depends on climate and elevation angle; an
        # overall average value in Europe and the United States of America,
        # at elevations between 10° and 50°, is s = 0.01
        s = 0.01
        sigma_z = s * F * A

        # Step 3a:  Calculate the conditional probability
        p = 2 / (np.pi * sigma_z * (1 + (z / (sigma_z) ** 2) ** 2))  # eq. 20

        # Step 3b : If required, calculate p(ζ |A), the conditional
        # probability (complementary cumulative distribution function) that
        # the fade slope ζ is exceeded for a given attenuation value, A:
        z_over_sigmaz = z / sigma_z
        abs_z_over_sigmaz = np.abs(z) / sigma_z
        P = (
            0.5
            - z_over_sigmaz / (np.pi * (1 + z_over_sigmaz ** 2))
            - (np.arctan(z_over_sigmaz) / np.pi)
        )

        # Step 3b2 : calculate p(ζ |A), the conditional probability that the
        # absolute value of the fade slope ζ is exceeded for a given
        # attenuation value, A:
        P2 = (
            1
            - 2 * abs_z_over_sigmaz / (np.pi * (1 + abs_z_over_sigmaz ** 2))
            - (2 * np.arctan(abs_z_over_sigmaz) / np.pi)
        )

        return p, P, P2, sigma_z

    @classmethod
    def fade_depth(self, N_target, D_target, A, PofA, el, f):

        d_target = np.atleast_1d(D_target)

        def delta_N_events(x):
            P_it = 10 ** (np.interp(x, A, np.log10(PofA)))
            T_tot_it = (P_it / 100) * 365.25 * 86400
            _, _, N_it, _ = self.fade_duration(d_target, x, el, f, T_tot_it)
            delta = N_target - N_it
            return delta

        a_min = fsolve(delta_N_events, 1)  # a_min should have shape (1,)
        return a_min.item()


class _ITU1623_0_:
    def __init__(self):
        self.__version__ = 0
        self.year = 2003
        self.month = 4
        self.link = "https://www.itu.int/rec/R-REC-P.1623-0-200304-S/en"

    @staticmethod
    def fade_duration(*args, **kwargs):
        return _ITU1623_1_.fade_duration(*args, **kwargs)

    @staticmethod
    def fade_slope(*args, **kwargs):
        return _ITU1623_1_.fade_slope(*args, **kwargs)

    @staticmethod
    def fade_depth(*args, **kwargs):
        return _ITU1623_1_.fade_depth(*args, **kwargs)


__model = __ITU1623__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.1623 recommendation currently being used.

    This function changes the model used for the ITU-R P.1623 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          * 1: Activates recommendation ITU-R P.1623-1 (03/2005) (Current version)
          * 0: Activates recommendation ITU-R P.1623-0 (04/2003) (Superseded)
    """
    global __model
    __model = __ITU1623__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.1623 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__


def fade_duration_probability(D, A, el, f):
    """Compute the probability of occurrence of fades of duration longer than D.

    Compute the probability of occurrence of fades of duration d longer than
    D (s), given that the attenuation a is greater than A (dB).

    This probability can be estimated from the ratio of the number of fades
    of duration longer than D to the total number of fades observed,
    given that the threshold A is exceeded.

    Parameters
    ----------
    D: number, sequence, or numpy.ndarray
        Event durations, array, (s)
    A: number
        Attenuation threshold, scalar, (dB)
    el: number
        Elevation angle towards the satellite, deg (5 - 60)
    f:  number
        Frequency, GHz (between 10 and 50 GHz)

    Returns
    -------
    p: number, sequence, or numpy.ndarray
        Probability of occurence of fade events of duration d longer than D
        given a>A, P(d > D|a > A)


    References
    ----------
    [1] Prediction method of fade dynamics on Earth-space paths:
    https://www.itu.int/rec/R-REC-P.1623/en

    """
    type_output = get_input_type(D)

    A = prepare_quantity(A, u.dB / u.s, 'Attenuation threshold')
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    f = prepare_quantity(f, u.GHz, 'Frequency')

    val = __model.fade_duration(D, A, el, f, 1)[0]
    return prepare_output_array(val, type_output) * u.dimensionless_unscaled


def fade_duration_cummulative_probability(D, A, el, f):
    """Compute the cumulative probability of exceedance of fades of duration
    longer than D.

    Compute the cummulative exceedance probability F(d > D|a > A),
    the total fraction (between 0 and 1) of fade time due to fades of duration
    d longer than D (s), given that the attenuation a is greater than A (dB).

    Parameters
    ----------
    D: number, sequence, or numpy.ndarray
        Event durations, array, (s)
    A: number
        Attenuation threshold, scalar, (dB)
    el: number
        Elevation angle towards the satellite, deg (5 - 60)
    f:  number
        Frequency, GHz (between 10 and 50 GHz)

    Returns
    -------

    F:  number, sequence, or numpy.ndarray
        Cumulative probability of exceedance, total fraction of fade time
        due to fades of d > D

    References
    ----------
    [1] Prediction method of fade dynamics on Earth-space paths:
    https://www.itu.int/rec/R-REC-P.1623/en

    """
    type_output = get_input_type(D)

    A = prepare_quantity(A, u.dB / u.s, 'Attenuation threshold')
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    f = prepare_quantity(f, u.GHz, 'Frequency')

    val = __model.fade_duration(D, A, el, f, 1)[1]
    return prepare_output_array(val, type_output) * u.dimensionless_unscaled


def fade_duration_number_fades(D, A, el, f, T_tot):
    """Compute the number of fades of duration longer than D.

    For a given reference period, the number of fades of duration longer
    D is estimated by multiplying the probability of occurrence P(d > D|a > A)
    by the total number of fades exceeding the threshold, Ntot(A).

    Parameters
    ----------
    D: number, sequence, or numpy.ndarray
        Event durations, array, (s)
    A: number
        Attenuation threshold, scalar, (dB)
    el: number
        Elevation angle towards the satellite, deg (5 - 60)
    f:  number
        Frequency, GHz (between 10 and 50 GHz)
    T_tot: number
        Total fade time from cumulative distribution (P(A)/100)*Reference time
        period. T_tot should be obtained from local data. If this long-term
        statistic is not available, an estimate can be calculated from
        Recommendation ITU-R P.618. In this case the procedure consists in
        calculating the CDF of total attenuation, deriving the percentage of
        time the considered attenuation threshold A is exceeded and then the
        associated total exceedance time T_tot for the reference period
        considered.

        For a reference period of a year,
        T_tot = ((100-availability_in_pctg)/100)*365.25*24*3600   [s]


    Returns
    -------
    N:  Total number of fades of duration d longer than D, for a given
        threshold A

    References
    ----------
    [1] Prediction method of fade dynamics on Earth-space paths:
    https://www.itu.int/rec/R-REC-P.1623/en

    """
    type_output = get_input_type(D)
    D = prepare_input_array(D)
    A = prepare_quantity(A, u.dB / u.s, 'Attenuation threshold')
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    f = prepare_quantity(f, u.GHz, 'Frequency')

    val = __model.fade_duration(D, A, el, f, T_tot)[2]
    return prepare_output_array(val, type_output) * u.dimensionless_unscaled


def fade_duration_total_exceedance_time(D, A, el, f, T_tot):
    """Compute the total exceedance time of fades of duration longer than D.

    The total exceedance time due to fade events of duration longer than D is
    obtained by multiplying the fraction of time F(d > D|a > A) by the total
    time that the threshold is exceeded, Ttot(A).

    Parameters
    ----------
    D: number, sequence, or numpy.ndarray
        Event durations, array, (s)
    A: number
        Attenuation threshold, scalar, (dB)
    el: number
        Elevation angle towards the satellite, deg (5 - 60)
    f:  number
        Frequency, GHz (between 10 and 50 GHz)
    T_tot: number
        Total fade time from cumulative distribution (P(A)/100)*Reference time
        period. T_tot should be obtained from local data. If this long-term
        statistic is not available, an estimate can be calculated from
        Recommendation ITU-R P.618. In this case the procedure consists in
        calculating the CDF of total attenuation, deriving the percentage of
        time the considered attenuation threshold A is exceeded and then the
        associated total exceedance time T_tot for the reference period
        considered.

        For a reference period of a year,
        T_tot = ((100-availability_in_pctg)/100)*365.25*24*3600   [s]

    Returns
    -------
    T:    Total fading time due to fades of d > D for A threshold.


    References
    ----------
    [1] Prediction method of fade dynamics on Earth-space paths:
    https://www.itu.int/rec/R-REC-P.1623/en

    """
    type_output = get_input_type(D)

    A = prepare_quantity(A, u.dB / u.s, 'Attenuation threshold')
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    f = prepare_quantity(f, u.GHz, 'Frequency')

    val = __model.fade_duration(D, A, el, f, T_tot)[3]
    return prepare_output_array(val, type_output) * u.s


def fade_duration(D, A, el, f, T_tot):
    """Compute the probability of occurrence of fades of duration longer than D.

    Compute the probability of occurrence of fades of duration d longer than
    D (s), given that the attenuation a is greater than A (dB) and
    F(d > D|a > A), the cumulative exceedance probability, or, equivalently,
    the total fraction (between 0 and 1) of fade time due to fades of duration
    d longer than D (s), given that the attenuation a is greater than A (dB).

    The function also returns other parameters associated to the fade duration
    prediction method. See ITU-R P.1623 Annex 1 Section 2.2

    Parameters
    ----------
    D: number, sequence, or numpy.ndarray
        Event durations, array, (s)
    A: number
        Attenuation threshold, scalar, (dB)
    el: number
        Elevation angle towards the satellite, deg (5 - 60)
    f:  number
        Frequency, GHz (between 10 and 50 GHz)
    T_tot: number
        Total fade time from cumulative distribution (P(A)/100)*Reference time
        period. T_tot should be obtained from local data. If this long-term
        statistic is not available, an estimate can be calculated from
        Recommendation ITU-R P.618. In this case the procedure consists in
        calculating the CDF of total attenuation, deriving the percentage of
        time the considered attenuation threshold A is exceeded and then the
        associated total exceedance time T_tot for the reference period
        considered.

        For a reference period of a year,
        T_tot = ((100-availability_in_pctg)/100)*365.25*24*3600   [s]

    Returns
    -------
    p:    probability of occurence of fade events of
             duration d longer than D given a>A, P(d > D|a > A)
    F:    cumulative probability of exceedance, total
          fraction of fade time due to fades of d > D
    N:    total number of fades of duration d longer than D, for a given
        threshold A
    T:    total fading time due to fades of d > D for A threshold


    References
    ----------
    [1] Prediction method of fade dynamics on Earth-space paths:
    https://www.itu.int/rec/R-REC-P.1623/en

    """
    get_input_type(D)

    A = prepare_quantity(A, u.dB / u.s, 'Attenuation threshold')
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    f = prepare_quantity(f, u.GHz, 'Frequency')

    val = __model.fade_duration(D, A, el, f, T_tot)
    return val


def fade_slope(z, A, f_B, delta_t):
    """Compute the probability of exceeding a valueo f fade slope.

    Fade slope is defined as the rate of change of attenuation with time
    information about the expected fade slope is essential to assess the
    required minimum tracking rate of a fade mitigation system.
    The model is valid for the following ranges of parameters:
        * frequencies from 10 to 30 GHz
        * elevation angles from 10° to 50°.

    See ITU-R P.1623 Annex 1 Section 3.2

    Parameters
    ----------
    z: number, sequence, or numpy.ndarray
        array of fade slope values (dB/s)
    A: number
        attenuation threshold, scalar, dB (range 0 - 20 dB)
    f_B: number
        3 dB cut-off frequency of the low pass filter (Hz, range 0.001 - 1)
        used to remove tropospheric scintillation and rapid variations of rain
        attenuation from the signal. Experimental results show that a 3 dB
        cut-off frequency of 0.02 Hz allows scintillation and rapid variations
        of rain attenuation to be filtered out adequately.
    delta_t: number
        Time interval length over which fade slope is calculated (s), 2-200 s

    Returns
    -------
    p: conditional  probability (probability  density  function)
       that the fade slope is equal to the fade slope for
       a given attenuation value, A
    P: conditional probability (complementary cumulative
       distribution function)that the fade slope is exceeded
       for a given attenuation value, A
    P2: conditional probability that the absolute value of
        the fade slope is exceeded for a given attenuation
        value, A
    sigma_z: standard deviation of the conditional fade slope

    Remark
    ------
    The output is an array of 4 elements.

    Example
    -------

    .. code-block:: python

        import itur.models.itu1623 as itu1623

        z = np.linspace(-2,2,100)
        A = 10
        f_B = 0.02
        delta_t = 1
        p, P, P2, sigma_z = itu1623.fade_slope(z, A, f_B, delta_t)

    References
    ----------
    [1] Prediction method of fade dynamics on Earth-space paths:
    https://www.itu.int/rec/R-REC-P.1623/en
    """
    get_input_type(z)
    z = prepare_quantity(z, u.dB / u.s, 'Fade slope values')
    A = prepare_quantity(A, u.dB / u.s, 'Attenuation threshold')
    delta_t = prepare_quantity(delta_t, u.s, 'Time interval')
    f_B = prepare_quantity(f_B, u.GHz, 'Cut-off Frequency')

    val = __model.fade_slope(z, A, f_B, delta_t)
    return val


def fade_depth(N_target, D_target, A, PofA, el, f):
    """Compute the maximum fade a link must tolerate given a target outage
    intensity value (number of events) and a target duration of event.

    The fade depth is computed by numerical solution of the fade_duration
    problem.

    See ITU-R P.1623 Annex 1 Section 3.2

    Parameters
    ----------
    N_target : int
        Target outage intensity (scalar)
    D_target : int
        Event duration (scalar)
    A : number, sequence, or numpy.ndarray
        Attenuation distribution (CDF, A) for the link under analysis
    PofA : number, sequence, or numpy.ndarray
        Probability that A is exceeded (CDF, probability)
    el : number
        Elevation angle (deg)
    f : number
        Frequency (GHz)

    Returns
    -------
    a_min: number
        Minimum attenuation the link must tolerate to meet the OI target

    Remark
    ------
    This function uses scipy's fsolve as optimizer.

    Example
    -------
    .. code-block:: python

        import itur.models.itu1623 as itu1623

        N_target = 25
        D_target = 60
        PofA = np.array([50, 30, 20, 10, 5, 3, 2, 1, .5, .3, .2, .1, .05, .03,
                         .02, .01, .005, .003, .002, .001])
        A = np.array([0.4, 0.6, 0.8, 1.8, 2.70, 3.5, 4.20, 5.7, 7.4, 9, 10.60,
                      14, 18.3, 22.3, 25.8, 32.6, 40.1, 46.1, 50.8, 58.8])
        el = 38.5
        f = 28
        itu1623.fade_depth(N_target, D_target, A, PofA, el, f) # 21.6922280


    References
    ----------
    [1] Prediction method of fade dynamics on Earth-space paths:
    https://www.itu.int/rec/R-REC-P.1623/en
    """
    get_input_type(A)
    el = prepare_quantity(el, u.deg, 'Elevation angle')
    f = prepare_quantity(f, u.GHz, 'Frequency')

    val = __model.fade_depth(N_target, D_target, A, PofA, el, f)
    return val
