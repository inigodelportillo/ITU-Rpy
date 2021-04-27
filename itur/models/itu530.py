# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from astropy import units as u
from scipy.optimize import bisect

from itur.models.itu453 import DN65
from itur.models.itu837 import rainfall_rate
from itur.models.itu1144 import bilinear_2D_interpolator
from itur.models.itu838 import (rain_specific_attenuation,
                                rain_specific_attenuation_coefficients)
from itur.utils import (prepare_input_array, prepare_quantity,
                        get_input_type,
                        prepare_output_array, load_data_interpolator)


class __ITU530__():
    """Private class to model the ITU-R P.530 recommendations.

    Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems

    Available versions:
       * P.530-16 (07/15) (Current version)

    Not available versions:

    This recommendation includes prediction methods for the propagation effects
    that should be taken into account in the design of digital fixed
    line-of-sight links, both in clear-air and rainfall conditions. It also
    provides link design guidance in clear step-by-step procedures including
    the use of mitigation techniques to minimize propagation impairments. The
    final outage predicted is the base for other Recommendations addressing
    error performance and availability.
    """

    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.530 recommendation.

    def __init__(self, version=17):
        if version == 17:
            self.instance = _ITU530_17_()
        elif version == 16:
            self.instance = _ITU530_16()
        else:
            raise ValueError(
                'Version ' +
                str(version) +
                ' is not implemented' +
                ' for the ITU-R P.530 model.')

    @property
    def __version__(self):
        return self.instance.__version__

    def fresnel_ellipse_radius(self, d1, d2, f):
        return self.instance.fresnel_ellipse_radius(d1, d2, f)

    def diffraction_loss(self, d1, d2, h, f):
        return self.instance.diffraction_loss(d1, d2, h, f)

    def multipath_loss_for_A(self, lat, lon, h_e, h_r, d, f, A):
        return self.instance.multipath_loss_for_A(lat, lon, h_e, h_r, d, f, A)

    def multipath_loss(self, lat, lon, h_e, h_r, d, f, A):
        return self.instance.multipath_loss(lat, lon, h_e, h_r, d, f, A)

    def rain_attenuation(self, lat, lon, d, f, el, p, tau=45, R001=None):
        return self.instance.rain_attenuation(lat, lon, d, f, el, p, tau, R001)

    def inverse_rain_attenuation(self, lat, lon, d, f, el,
                                 Ap, tau=45, R001=None):
        return self.instance.inverse_rain_attenuation(lat, lon, d, f, el,
                                                      Ap, tau, R001)

    def rain_event_count(self, lat, lon, d, f, el, A, tau=45, R001=None):
        return self.instance.rain_event_count(lat, lon, d, f, el, A, tau, R001)

    def XPD_outage_clear_air(self, lat, lon, h_e, h_r,
                             d, f, XPD_g, C0_I, XPIF=0):
        return self.instance.XPD_outage_clear_air(lat, lon, h_e, h_r, d, f,
                                                  XPD_g, C0_I, XPIF)

    def XPD_outage_precipitation(self, lat, lon, d, f, el, C0_I, tau=45,
                                 U0=15, XPIF=0):
        return self.instance.XPD_outage_precipitation(lat, lon, d, f, el, C0_I,
                                                      tau, U0, XPIF)


class _ITU530_17_():

    _s_a = {}

    def __init__(self):
        self.__version__ = 17
        self.year = 2017
        self.month = 12
        self.link = 'https://www.itu.int/rec/R-REC-P.530-17-201712-S/en'

    @classmethod
    def s_a(self, lat, lon):
        """
        Standard deviation of terrain heights.

        Computes the Standard deviation of terrain heights (m) within a
        110 km × 110 km area with a 30 s resolution (e.g. the Globe “gtopo30”
        data).

        The value for the mid-path may be obtained from an area roughness
        with 0.5 × 0.5 degree resolution of geographical coordinates
        using bi-linear interpolation.
        """
        if not _ITU530_17_._s_a:
            _ITU530_17_._s_a = load_data_interpolator(
                '530/v16_lat.npz', '530/v16_lon.npz',
                '530/v16_gtopo_30.npz', bilinear_2D_interpolator,
                flip_ud=False)

        return _ITU530_17_._s_a(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    ###########################################################################
    #                               Section 2.2                               #
    ###########################################################################
    @classmethod
    def fresnel_ellipse_radius(self, d1, d2, f):
        """
        Compute the Fresnel ellipse radius at a given frequency.

        Implementation of 'fresnel_ellipse_radius' method for recommendation
        ITU-P R.530-16. See documentation for function
        'ITUR530.fresnel_ellipse_radius'
        """
        return 17.3 * np.sqrt(d1 * d2 / (f * (d1 + d2)))

    @classmethod
    def diffraction_loss(self, d1, d2, h, f):
        """
        Compute the diffraction losses at a given frequency.

        Implementation of 'diffraction_loss' method for recommendation
        ITU-P R.530-16. See documentation for function
        'ITUR530.diffraction_loss'
        """
        F1 = self.fresnel_ellipse_radius(d1, d2, f)       # Eq. 2 [m]
        Ad = -20 * h / F1 + 10                  # Eq. 3 [dB]
        return Ad

    ###########################################################################
    #                               Section 2.3                               #
    ###########################################################################
    @classmethod
    def multipath_loss_for_A(self, lat, lon, h_e, h_r, d, f, A):
        """
        Computes the percentage of time that a fade depth A is exceeded
        due to multi-path losses.

        Implementation of 'multipath_loss_for_A' method for recommendation
        ITU-P R.530-16. See documentation for function
        'ITUR530.multipath_loss_for_A'
        """
        # Step 1: Estimate the geoclimatic factor K
        # DN1 point refractivity gradient in the lowest 65 m of the atmosphere
        # not exceeded for 1% of an average year
        # s_a is the area terrain roughness
        s_a = self.s_a(lat, lon)
        dN1 = DN65(lat, lon, 1).value
        K = 10**(-4.4 - 0.0027 * dN1) * (10 + s_a)**(-0.46)       # Eq. 4 [-]

        # Step 2: Calculate the magnitude of the path inclination
        # Eq. 5 [mrad]
        e_p = np.abs(h_r - h_e) / d

        # Step 3: For detailed link design applications calculate the
        # percentage of time (p_W) that fade depth A (dB) is exceeded in the
        # average worst month
        h_L = np.minimum(h_e, h_r)
        p_W = K * d**3.4 * (1 + e_p)**-1.03 * f**0.8 * \
            10**(-0.00076 * h_L - A / 10)
        # Eq. 7 [%]
        return p_W

    @classmethod
    def multipath_loss(self, lat, lon, h_e, h_r, d, f, A):
        """
        Estimate the number of fade events exceeding attenuation 'A'
        for 10 seconds or longer.

        Implementation of 'multipath_loss' method for recommendation
        ITU-P R.530-16. See documentation for function
        'ITUR530.multipath_loss'
        """
        # Step 1: Using the method multipath_loss_for_A calculate the
        # multi-path occurrence factor, p0
        p0 = self.multipath_loss_for_A(
            lat, lon, h_e, h_r, d, f, 0)   # Eq. 10 [%]

        # Step 2: Calculate the value of fade depth, At, at which the
        # transition occurs between the deep-fading distribution and the
        # shallow-fading distribution
        At = 25 + 1.2 * np.log10(p0)                        # Eq. 12 [dB]

        # Step 3: Calculate the percentage of time that A is exceeded in the
        # average worst month:
        def step_3b(p_0, At, A):
            p_t = p_0 * 10 ** (-At / 10)
            qa_p = -20 * np.log10(-np.log((100 - p_t) / 100)) / At
            q_t = ((qa_p - 2) /
                   (1 + 0.3 * 10 ** (-At / 20) * 10 ** (-0.016 * At)) -
                   4.3 * (10**(-At / 20) + At / 800))
            q_a = 2 + (1 + 0.3 * 10**(-A / 20)) * (10**(-0.016 * A)) *\
                (q_t + 4.3 * (10**(-A / 20 + A / 800)))
            p_W = 100 * (1 - np.exp(-10 ** (-q_a * A / 20)))
            return p_W

        p_W = np.where(A >= At, p0 * 10 ** (-A / 10), step_3b(p0, At, A))
        # Eq. 13 and Eq. 18 [%]
        return p_W

    ###########################################################################
    #                               Section 2.4                               #
    ###########################################################################
    @classmethod
    def rain_attenuation(self, lat, lon, d, f, el, p, tau=45, R001=None):
        """
        Estimate long-term statistics of rain attenuation.

        Implementation of 'rain_attenuation' method for recommendation
        ITU-P R.530-16. See documentation for function
        'ITUR530.rain_attenuation'
        """
        # Step 1: Obtain the rain rate R0.01 exceeded for 0.01% of the time
        # (with an integration time of 1 min).
        if R001 is None:
            R001 = rainfall_rate(lat, lon, 0.01).value

        # Step 2: Compute the specific attenuation, 'gammar' (dB/km) for the
        # frequency, polarization and rain rate of interest using
        # Recommendation ITU-R P.838
        gammar = rain_specific_attenuation(R001, f, el, tau).value
        _, alpha = rain_specific_attenuation_coefficients(f, el, tau)

        # Step 3: Compute the effective path length, 'deff', of the link by
        # multiplying the actual path length d by a distance factor r
        r = 1 / (0.477 * d ** 0.633 * R001 ** (0.073 * alpha) *
                 f**(0.123) - 10.579 * (1 - np.exp(-0.024 * d)))  # Eq. 32 [-]
        deff = np.minimum(r, 2.5) * d

        # Step 4: An estimate of the path attenuation exceeded for 0.01% of
        # the time is given by:
        A001 = gammar * deff                                    # Eq. 33 [dB]

        # Step 5: The attenuation exceeded for other percentages of time p in
        # the range 0.001% to 1% may be deduced from the following power law
        C0 = np.where(f >= 10, 0.12 + 0.4 * (np.log10(f / 10)**0.8), 0.12)
        # Eq. 35a [-]
        C1 = (0.07**C0) * (0.12**(1 - C0))
        # Eq. 35b [-]
        C2 = 0.855 * C0 + 0.546 * (1 - C0)
        C3 = 0.139 * C0 + 0.043 * (1 - C0)                      # Eq. 35c [-]
        Ap = A001 * C1 * p ** (- (C2 + C3 * np.log10(p)))       # Eq. 34 [dB]
        return Ap

    @staticmethod
    def inverse_rain_attenuation(
            lat, lon, d, f, el, Ap, tau=45, R001=None):
        """
        Estimate the percentage of time a given attenuation is exceeded due to
        rain events.

        Implementation of 'inverse_rain_attenuation' method for
        recommendation ITU-P R.530-16. See documentation for function
        'ITUR530.inverse_rain_attenuation'
        """

        # Step 1: Obtain the rain rate R0.01 exceeded for 0.01% of the time
        # (with an integration time of 1 min).
        if R001 is None:
            R001 = rainfall_rate(lat, lon, 0.01).value

        # Step 2: Compute the specific attenuation, 'gammar' (dB/km) for the
        # frequency, polarization and rain rate of interest using
        # Recommendation ITU-R P.838
        gammar = rain_specific_attenuation(R001, f, el, tau).value
        _, alpha = rain_specific_attenuation_coefficients(f, el, tau)

        # Step 3: Compute the effective path length, 'deff', of the link by
        # multiplying the actual path length d by a distance factor r
        r = 1 / (0.477 * d ** 0.633 * R001 ** (0.073 * alpha) *
                 f**(0.123) - 10.579 * (1 - np.exp(-0.024 * d)))  # Eq. 32 [-]
        deff = np.minimum(r, 2.5) * d

        # Step 4: An estimate of the path attenuation exceeded for 0.01% of
        # the time is given by:
        A001 = gammar * deff                                    # Eq. 33 [dB]

        # Step 5: The attenuation exceeded for other percentages of time p in
        # the range 0.001% to 1% may be deduced from the following power law
        C0 = np.where(f >= 10, 0.12 + 0.4 * (np.log10(f / 10)**0.8), 0.12)
        # Eq. 35a [-]
        C1 = (0.07**C0) * (0.12**(1 - C0))
        # Eq. 35b [-]
        C2 = 0.855 * C0 + 0.546 * (1 - C0)
        C3 = 0.139 * C0 + 0.043 * (1 - C0)                      # Eq. 35c [-]

        def func_bisect(p):
            return A001 * C1 * p ** (- (C2 + C3 * np.log10(p))) - Ap

        return bisect(func_bisect, 0.000001, 100)

    @classmethod
    def rain_event_count(self, lat, lon, d, f, el, A, tau=45, R001=None):
        """
        Estimate the number of fade events exceeding attenuation 'A'
        for 10 seconds or longer.

        Implementation of 'rain_event_count' method for recommendation
        ITU-P R.530-16. See documentation for function
        'ITUR530.rain_event_count'
        """
        # Compute the the percentage of time that the rain attenuation A(dB)
        # exceeded in the average year.
        p_A = self.inverse_rain_attenuation(lat, lon, d, f, el, Ap=A,
                                            tau=tau, R001=R001)

        # The number of fade events exceeding attenuation A for 10 s or longer
        N10s = 1 + 1313 * p_A**0.945                               # Eq. 78 [-]

        return N10s

    ###########################################################################
    #                                Section 4                                #
    ###########################################################################
    @classmethod
    def XPD_outage_clear_air(self, lat, lon, h_e, h_r,
                             d, f, XPD_g, C0_I, XPIF=0):
        """
        Estimate the probability of outage due to cross-polar discrimination
        reduction due to clear-air effects, assuming that a target C0_I is
        required.

        Implementation of 'XPD_outage_clear_air' method for recommendation
        ITU-P R.530-16. See documentation for function
        'ITUR530.XPD_outage_clear_air'
        """
        # Step 1
        XPD_0 = np.where(XPD_g <= 35, XPD_g + 5, 40)            # Eq. 101

        # Step 2: Evaluate the multi-path activity parameter:
        P0 = self.multipath_loss_for_A(lat, lon, h_e, h_r, d, f, 0)
        eta = 1 - np.exp(-0.2 * P0**0.75)                       # Eq. 102

        # Step 3:
        kXP = 0.7                                               # Eq. 104
        Q = - 10 * np.log10(kXP * eta / P0)                     # Eq. 103

        # Step 4: Derive the parameter C:
        C = XPD_0 + Q                                           # Eq. 105

        # Step 5: Calculate the probability of outage PXP due to clear-air
        # cross-polarization:
        M_XPD = C - C0_I + XPIF
        P_XP = P0 * 10 ** (- M_XPD / 10)                       # Eq. 106 [%]
        return P_XP

    @classmethod
    def XPD_outage_precipitation(self, lat, lon, d, f, el, C0_I, tau=45,
                                 U0=15, XPIF=0):
        """
        Estimate the probability of outage due to cross-polar discrimination
        reduction due to clear-air effects, assuming that a target C0_I is
        required.

        Implementation of 'XPD_outage_precipitation' method for recommendation
        ITU-P R.530-16. See documentation for function
        'ITUR530.XPD_outage_precipitation'
        """
        # Step 1: Determine the path attenuation, A0.01 (dB), exceeded
        # for 0.01% of the time
        A001 = self.rain_attenuation(lat, lon, d, f, el, 0.01)

        # Step 2: Determine the equivalent path attenuation, Ap
        U = U0 + 30 * np.log10(f)
        V = np.where(f < 20, 12.8 * f**0.19, 22.6)
        Ap = 10 ** ((U - C0_I + XPIF) / V)                     # Eq. 112

        # Step 3: Determine parameters m and n
        m = min(23.26 * np.log10(Ap / (0.12 * A001)), 40)      # Eq. 113
        n = (-12.7 + np.sqrt(161.23 - 4 * m)) / 2              # Eq. 114

        # Step 4 : Determine the outage probability
        P_XPR = 10**(n - 2)                                    # Eq. 115 [%]
        return P_XPR


class _ITU530_16():

    def __init__(self):
        self.__version__ = 16
        self.year = 2015
        self.month = 7
        self.link = 'https://www.itu.int/rec/R-REC-P.530-16-201507-S/en'

        self._s_a = {}

    @staticmethod
    def s_a(self, *args, **kwargs):
        return _ITU530_17_.s_a(*args, **kwargs)

    ###########################################################################
    #                               Section 2.2                               #
    ###########################################################################
    @staticmethod
    def fresnel_ellipse_radius(*args, **kwargs):
        return _ITU530_17_.fresnel_ellipse_radius(*args, **kwargs)

    @staticmethod
    def diffraction_loss(*args, **kwargs):
        return _ITU530_17_.diffraction_loss(*args, **kwargs)

    ###########################################################################
    #                               Section 2.3                               #
    ###########################################################################
    @staticmethod
    def multipath_loss_for_A(*args, **kwargs):
        return _ITU530_17_.multipath_loss_for_A(*args, **kwargs)

    @staticmethod
    def multipath_loss(*args, **kwargs):
        return _ITU530_17_.multipath_loss(*args, **kwargs)

    ###########################################################################
    #                               Section 2.4                               #
    ###########################################################################
    @staticmethod
    def rain_attenuation(*args, **kwargs):
        return _ITU530_17_.rain_attenuation(*args, **kwargs)

    @staticmethod
    def inverse_rain_attenuation(*args, **kwargs):
        return _ITU530_17_.inverse_rain_attenuation(*args, **kwargs)

    @staticmethod
    def rain_event_count(*args, **kwargs):
        return _ITU530_17_.rain_event_count(*args, **kwargs)

    ###########################################################################
    #                                Section 4                                #
    ###########################################################################
    @staticmethod
    def XPD_outage_clear_air(*args, **kwargs):
        return _ITU530_17_.XPD_outage_clear_air(*args, **kwargs)

    @staticmethod
    def XPD_outage_precipitation(*args, **kwargs):
        return _ITU530_17_.XPD_outage_precipitation(*args, **kwargs)


__model = __ITU530__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.530 recommendation currently being used.

    This function changes the model used for the ITU-R P.530 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          * 16: Activates recommendation ITU-R P.530-16 (07/15) (Current version)
    """
    global __model
    __model = __ITU530__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.530 recommendation currently being used.

    Returns
    -------
    version: int
       The version of the ITU-R P.530 recommendation being used.
    """
    return __model.__version__


def fresnel_ellipse_radius(d1, d2, f):
    """Compute the radius of the first Fresnel ellipsoid.

    Parameters
    ----------
    d1 : number, sequence, or numpy.ndarray
        Distances from the first terminal to the path obstruction. [km]
    d2 : number, sequence, or numpy.ndarray
        Distances from the second terminal to the path obstruction. [km]
    f : number
        Frequency of the link [GHz]


    Returns
    -------
    F1: Quantity
       Radius of the first Fresnel ellipsoid [m]


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems: https://www.itu.int/rec/R-REC-P.530/en
    """
    type_output = get_input_type(d1)
    d1 = prepare_quantity(d1, u.km, 'Distance to the first terminal')
    d2 = prepare_quantity(d2, u.km, 'Distance to the second terminal')
    f = prepare_quantity(f, u.GHz, 'Frequency')

    val = __model.fresnel_ellipse_radius(d1, d2, f)
    return prepare_output_array(val, type_output) * u.m


def diffraction_loss(d1, d2, h, f):
    """Estimate the diffraction loss over average terrain.

    Diffraction loss over average terrain. This value is valid for losses
    greater than 15 dB.


    Parameters
    ----------
    d1 : number, sequence, or numpy.ndarray
        Distances from the first terminal to the path obstruction. [km]
    d2 : number, sequence, or numpy.ndarray
        Distances from the second terminal to the path obstruction. [km]
    h : number, sequence, or numpy.ndarray
        Height difference between most significant path blockage
        and the path trajectory. h is negative if the top of the obstruction
        of interest is above the virtual line-of-sight. [m]
    f : number
        Frequency of the link [GHz]


    Returns
    -------
    A_d: Quantity
        Diffraction loss over average terrain  [dB]


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems: https://www.itu.int/rec/R-REC-P.530/en
    """
    type_output = get_input_type(d1)
    d1 = prepare_quantity(d1, u.km, 'Distance to the first terminal')
    d2 = prepare_quantity(d2, u.km, 'Distance to the second terminal')
    h = prepare_quantity(h, u.m, 'Height difference')
    f = prepare_quantity(f, u.GHz, 'Frequency')

    val = __model.diffraction_loss(d1, d2, h, f)
    return prepare_output_array(val, type_output) * u.m


def multipath_loss_for_A(lat, lon, h_e, h_r, d, f, A):
    """Estimate the single-frequency (or narrow-band) fading distribution.

    Method for predicting the single-frequency (or narrow-band) fading
    distribution at large fade depths in the average worst month in any part
    of the world. Given a fade depth value 'A', determines the amount of time
    it will be exceeded during a year

    This method does not make use of the path profile and can be used for
    initial planning, licensing, or design purposes.

    This method is only valid for small percentages of time.

    Multi-path fading and enhancement only need to be calculated for path
    lengths longer than 5 km, and can be set to zero for shorter paths.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    h_e : number
        Emitter antenna height (above the sea level) [m]
    h_r : number
        Receiver antenna height (above the sea level) [m]
    d : number, sequence, or numpy.ndarray
        Distances between antennas [km]
    f : number
        Frequency of the link [GHz]
    A : number
         Fade depth [dB]


    Returns
    -------
    p_w: Quantity
         percentage of time that fade depth A is exceeded in the average
         worst month  [%]


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems: https://www.itu.int/rec/R-REC-P.530/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    h_e = prepare_quantity(
        h_e, u.m, 'Emitter antenna height (above sea level)')
    h_r = prepare_quantity(
        h_r, u.m, 'Receiver antenna height (above sea level)')
    d = prepare_quantity(d, u.km, 'Distance between antennas')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    A = prepare_quantity(A, u.dB, 'Fade depth')

    val = __model.multipath_loss_for_A(lat, lon, h_e, h_r, d, f, A)
    return prepare_output_array(val, type_output) * u.percent


def multipath_loss(lat, lon, h_e, h_r, d, f, A):
    """Estimate the percentage of time that any fade depth is exceeded.

    Method for predicting the percentage of time that any fade depth is
    exceeded. This method combines the deep fading distribution given in the
    multipath_loss_for_A' and an empirical interpolation procedure for shallow
    fading down to 0 dB.

    This method does not make use of the path profile and can be used for
    initial planning, licensing, or design purposes.

    Multi-path fading and enhancement only need to be calculated for path
    lengths longer than 5 km, and can be set to zero for shorter paths.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    h_e : number
        Emitter antenna height (above the sea level) [m]
    h_r : number
        Receiver antenna height (above the sea level) [m]
    d : number, sequence, or numpy.ndarray
        Distances between antennas [km]
    f : number
        Frequency of the link [GHz]
    A : number
         Fade depth [dB]


    Returns
    -------
    p_w: Quantity
         percentage of time that fade depth A is exceeded in the average
         worst month  [%]


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems: https://www.itu.int/rec/R-REC-P.530/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    h_e = prepare_quantity(
        h_e, u.m, 'Emitter antenna height (above sea level)')
    h_r = prepare_quantity(
        h_r, u.m, 'Receiver antenna height (above sea level)')
    d = prepare_quantity(d, u.km, 'Distance between antennas')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    A = prepare_quantity(A, u.dB, 'Fade depth')

    val = __model.multipath_loss(lat, lon, h_e, h_r, d, f, A)
    return prepare_output_array(val, type_output) * u.percent


def rain_attenuation(lat, lon, d, f, el, p, tau=45, R001=None):
    """Estimate long-term statistics of rain attenuation.

    Attenuation can also
    occur as a result of absorption and scattering by such hydro-meteors as
    rain, snow, hail and fog. Although rain attenuation can be ignored at
    frequencies below about 5 GHz, it must be included in design calculations
    at higher frequencies, where its importance increases rapidly.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    d : number, sequence, or numpy.ndarray
        Path length [km]
    f : number
        Frequency of the link [GHz]
    el : sequence, or number
        Elevation angle (degrees)
    p : number
        Percentage of the time the rain attenuation value is exceeded.
    R001: number, optional
        Point rainfall rate for the location for 0.01% of an average year
        (mm/h). If not provided, an estimate is obtained from Recommendation
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


    Returns
    -------
    A_r: Quantity
         Attenuation exceeded during p percent of the time  [dB]


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems: https://www.itu.int/rec/R-REC-P.530/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    d = prepare_quantity(d, u.km, 'Distance between antennas')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(el, u.deg, 'Elevation Angle')
    R001 = prepare_quantity(R001, u.mm / u.hr, 'Rainfall Rate')

    val = __model.rain_attenuation(lat, lon, d, f, el, p, tau, R001)
    return prepare_output_array(val, type_output) * u.dB


def inverse_rain_attenuation(lat, lon, d, f, el, Ap, tau=45, R001=None):
    """Estimate the percentage of time a given attenuation is exceeded.

    Estimate the percentage of time a given attenuation is exceeded due to
    rain events.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    d : number, sequence, or numpy.ndarray
        Path length [km]
    f : number
        Frequency of the link [GHz]
    el : sequence, or number
        Elevation angle (degrees)
    Ap : number
        Fade depth
    R001: number, optional
        Point rainfall rate for the location for 0.01% of an average year
        (mm/h). If not provided, an estimate is obtained from Recommendation
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


    Returns
    -------
    p: Quantity
        Percentage of time that the attenuation A is exceeded.


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems: https://www.itu.int/rec/R-REC-P.530/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    d = prepare_quantity(d, u.km, 'Distance between antennas')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(el, u.deg, 'Elevation Angle')
    Ap = prepare_quantity(Ap, u.dB, 'Fade depth')
    R001 = prepare_quantity(R001, u.mm / u.hr, 'Rainfall Rate')

    val = __model.inverse_rain_attenuation(lat, lon, d, f, el, Ap, tau, R001)
    return prepare_output_array(val, type_output) * u.percent


def rain_event_count(lat, lon, d, f, el, A, tau=45, R001=None):
    """Estimate the number of fade events exceeding attenuation 'A'.

    Estimate the number of fade events exceeding attenuation 'A'
    for 10 seconds or longer.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    d : number, sequence, or numpy.ndarray
        Path length [km]
    f : number
        Frequency of the link [GHz]
    el : sequence, or number
        Elevation angle (degrees)
    A : number
        Fade depth
    R001: number, optional
        Point rainfall rate for the location for 0.01% of an average year
        (mm/h). If not provided, an estimate is obtained from Recommendation
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


    Returns
    -------
    p: Quantity
        Percentage of time that the attenuation A is exceeded.


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems: https://www.itu.int/rec/R-REC-P.530/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    d = prepare_quantity(d, u.km, 'Distance between antennas')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    el = prepare_quantity(el, u.deg, 'Elevation Angle')
    A = prepare_quantity(A, u.dB, 'Fade depth')
    R001 = prepare_quantity(R001, u.mm / u.hr, 'Rainfall Rate')

    val = __model.rain_event_count(lat, lon, d, f, el, A, tau, R001)
    return prepare_output_array(val, type_output) * u.percent


def XPD_outage_clear_air(lat, lon, h_e, h_r, d, f, XPD_g, C0_I, XPIF=0):
    """Estimate the probability of outage due to cross-polar discrimination.

    Estimate the probability of outage due to cross-polar discrimination
    reduction due to clear-air effects, assuming that a target C0_I is
    required.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    h_e : number
        Emitter antenna height (above the sea level) [m]
    h_r : number
        Receiver antenna height (above the sea level) [m]
    d : number, sequence, or numpy.ndarray
        Distances between antennas [km]
    f : number
        Frequency of the link [GHz]
    XPD_g : number
        Manufacturer's guaranteed minimum XPD at boresight for both the
        transmitting and receiving antennas [dB]
    C0_I : number
         Carrier-to-interference ratio for a reference BER [dB]
    XPIF : number, optional
        Laboratory-measured cross-polarization improvement factor that gives
        the difference in cross-polar isolation (XPI) at sufficiently large
        carrier-to-noise ratio (typically 35 dB) and at a specific BER for
        systems with and without cross polar interference canceler (XPIC).
        A typical value of XPIF is about 20 dB. Default value 0 dB (no XPIC)
        [dB]


    Returns
    -------
    p_XP: Quantity
        Probability of outage due to clear-air cross-polarization


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems: https://www.itu.int/rec/R-REC-P.530/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    h_e = prepare_quantity(
        h_e, u.m, 'Emitter antenna height (above sea level)')
    h_r = prepare_quantity(
        h_r, u.m, 'Receiver antenna height (above sea level)')
    d = prepare_quantity(d, u.km, 'Distance between antennas')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    XPD_g = prepare_quantity(XPD_g, u.dB, 'Manufacturers minimum XPD')
    C0_I = prepare_quantity(C0_I, u.dB, 'Carrier-to-interference ratio')
    XPIF = prepare_quantity(
        XPIF, u.dB, 'Cross-polarization improvement factor')

    val = __model.XPD_outage_clear_air(
        lat, lon, h_e, h_r, d, f, XPD_g, C0_I, XPIF)
    return prepare_output_array(val, type_output) * u.percent


def XPD_outage_precipitation(lat, lon, d, f, el, C0_I, tau=45,
                             U0=15, XPIF=0):
    """Estimate the probability of outage due to cross-polar discrimination.

    Estimate the probability of outage due to cross-polar discrimination
    reduction due to clear-air effects, assuming that a target C0_I is
    required.


    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    d : number, sequence, or numpy.ndarray
        Distances between antennas [km]
    f : number
        Frequency of the link [GHz]
    el : sequence, or number
        Elevation angle (degrees)
    C0_I : number
         Carrier-to-interference ratio for a reference BER [dB]
    tau : number, optional
        Polarization tilt angle relative to the horizontal (degrees)
        (tau = 45 deg for circular polarization). Default value is 45
    U0 : number, optional
        Coefficient for the cumulative distribution of the co-polar attenuation
        (CPA) for rain. Default 15 dB.
    XPIF : number, optional
        Laboratory-measured cross-polarization improvement factor that gives
        the difference in cross-polar isolation (XPI) at sufficiently large
        carrier-to-noise ratio (typically 35 dB) and at a specific BER for
        systems with and without cross polar interference canceler (XPIC).
        A typical value of XPIF is about 20 dB. Default value 0 dB (no XPIC)
        [dB]


    Returns
    -------
    p_XP: Quantity
        Probability of outage due to clear-air cross-polarization


    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    terrestrial line-of-sight systems: https://www.itu.int/rec/R-REC-P.530/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    d = prepare_quantity(d, u.km, 'Distance between antennas')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    C0_I = prepare_quantity(C0_I, u.dB, 'Carrier-to-interference ratio')
    U0 = prepare_quantity(U0, u.dB, 'Coefficient for the CPA')
    XPIF = prepare_quantity(
        XPIF, u.dB, 'Cross-polarization improvement factor')

    val = __model.XPD_outage_precipitation(lat, lon, d, f, C0_I, tau, U0, XPIF)
    return prepare_output_array(val, type_output) * u.percent
