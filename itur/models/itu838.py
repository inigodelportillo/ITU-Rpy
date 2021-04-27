# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from astropy import units as u
from itur.utils import prepare_quantity


class __ITU838__():
    """Specific attenuation model for rain for use in prediction methods

    Available versions include:
    * P.838-0 (03/92) (Superseded)
    * P.838-1 (10/99) (Superseded)
    * P.838-2 (04/03) (Superseded)
    * P.838-3 (03/05) (Current version)
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.838 recommendation.

    def __init__(self, version=3):
        if version == 3:
            self.instance = _ITU838_3_()
        elif version == 2:
            self.instance = _ITU838_2_()
        elif version == 1:
            self.instance = _ITU838_1_()
        elif version == 0:
            self.instance = _ITU838_0_()
        else:
            raise ValueError(
                'Version ' +
                str(version) +
                ' is not implemented' +
                ' for the ITU-R P.838 model.')

    @property
    def __version__(self):
        return self.instance.__version__

    def rain_specific_attenuation_coefficients(self, f, el, tau):
        # Abstract method to compute the rain height
        fcn = np.vectorize(self.instance.rain_specific_attenuation_coefficients,
                           excluded=[1], otypes=[np.ndarray])
        return np.array(fcn(f, el, tau).tolist())

    def rain_specific_attenuation(self, R, f, el, tau):
        # Abstract method to compute the zero isoterm height
        k, alpha = self.rain_specific_attenuation_coefficients(f, el, tau)
        return k * R**alpha


class _ITU838_3_():

    def __init__(self):
        self.__version__ = 3
        self.year = 2005
        self.month = 3
        self.link = 'https://www.itu.int/rec/R-REC-P.838-3-200503-I/en'

    @staticmethod
    def rain_specific_attenuation_coefficients(f, el, tau):
        """

        """
        kh = {'aj': [-5.33980, -0.35351, -0.23789, -0.94158],
              'bj': [-0.10008, 1.2697, 0.86036, 0.64552],
              'cj': [1.13098, 0.454, 0.15354, 0.16817],
              'mk': -0.18961,
              'ck': 0.71147}

        kv = {'aj': [-3.80595, -3.44965, -0.39902, 0.50167],
              'bj': [0.56934, -0.22911, 0.73042, 1.07319],
              'cj': [0.81061, 0.51059, 0.11899, 0.27195],
              'mk': -0.16398,
              'ck': 0.63297}

        alphah = {'aj': [-0.14318, 0.29591, 0.32177, -5.37610, 16.1721],
                  'bj': [1.82442, 0.77564, 0.63773, -0.96230, -3.29980],
                  'cj': [-0.55187, 0.19822, 0.13164, 1.47828, 3.4399],
                  'ma': 0.67849,
                  'ca': -1.95537}

        alphav = {'aj': [-0.07771, 0.56727, -0.20238, -48.2991, 48.5833],
                  'bj': [2.3384, 0.95545, 1.1452, 0.791669, 0.791459],
                  'cj': [-0.76284, 0.54039, 0.26809, 0.116226, 0.116479],
                  'ma': -0.053739,
                  'ca': 0.83433}

        def curve_fcn(f, a, b, c):
            return (a * np.exp(-((np.log10(f) - b) / c)**2))

        KH = np.power(10, sum([curve_fcn(f, kh['aj'][j], kh['bj'][j], kh['cj'][j])
                               for j in range(4)]) + kh['mk'] * np.log10(f) + kh['ck'])
        KV = np.power(10, sum([curve_fcn(f, kv['aj'][j], kv['bj'][j], kv['cj'][j])
                               for j in range(4)]) + kv['mk'] * np.log10(f) + kv['ck'])

        alphaH = sum([curve_fcn(f, alphah['aj'][j], alphah['bj'][j], alphah['cj'][j])
                      for j in range(5)]) + alphah['ma'] * np.log10(f) + alphah['ca']
        alphaV = sum([curve_fcn(f, alphav['aj'][j], alphav['bj'][j], alphav['cj'][j])
                      for j in range(5)]) + alphav['ma'] * np.log10(f) + alphav['ca']

        k = (KH + KV + (KH - KV) * np.cos(np.deg2rad(el))**2
             * np.cos(np.deg2rad(2 * tau))) / 2.0

        alpha = (KH * alphaH + KV * alphaV + (KH * alphaH - KV * alphaV) *
                 np.cos(np.deg2rad(el))**2 * np.cos(np.deg2rad(2 * tau))) / (2.0 * k)

        return k, alpha


class _ITU838_2_():

    def __init__(self):
        self.__version__ = 2
        self.year = 2003
        self.month = 4
        self.link = 'https://www.itu.int/rec/R-REC-P.838-2-200304-S/en'

    @staticmethod
    def rain_specific_attenuation_coefficients(f, el, tau):
        """

        """
        kh = {'aj': [0.3364, 0.7520, -0.9466],
              'bj': [1.1274, 1.6644, 2.8496],
              'cj': [0.2916, 0.5175, 0.4315],
              'mk': 1.9925,
              'ck': -4.4123}

        kv = {'aj': [0.3023, 0.7790, -1.0022],
              'bj': [1.1402, 1.6723, 2.9400],
              'cj': [0.2826, 0.5694, 0.4823],
              'mk': 1.9710,
              'ck': -4.4535}

        alphah = {'aj': [0.5564, 0.2237, -0.1961, -0.02219],
                  'bj': [0.7741, 1.4023, 0.5769, 2.2959],
                  'cj': [0.4011, 0.3475, 0.2372, 0.2801],
                  'ma': -0.08016,
                  'ca': 0.8993}

        alphav = {'aj': [0.5463, 0.2158, -0.1693, -0.01895],
                  'bj': [0.8017, 1.4080, 0.6353, 2.3105],
                  'cj': [0.3657, 0.3636, 0.2155, 0.2938],
                  'ma': -0.07059,
                  'ca': 0.8756}

        def curve_fcn(f, a, b, c):
            return a * np.exp(-((np.log10(f) - b) / c)**2)

        KH = np.power(10, sum([curve_fcn(f, kh['aj'][j], kh['bj'][j], kh['cj'][j])
                               for j in range(3)]) + kh['mk'] * np.log10(f) + kh['ck'])
        KV = np.power(10, sum([curve_fcn(f, kv['aj'][j], kv['bj'][j], kv['cj'][j])
                               for j in range(3)]) + kv['mk'] * np.log10(f) + kv['ck'])

        alphaH = sum([curve_fcn(f, alphah['aj'][j], alphah['bj'][j], alphah['cj'][j])
                      for j in range(4)]) + alphah['ma'] * np.log10(f) + alphah['ca']
        alphaV = sum([curve_fcn(f, alphav['aj'][j], alphav['bj'][j], alphav['cj'][j])
                      for j in range(4)]) + alphav['ma'] * np.log10(f) + alphav['ca']

        k = (KH + KV + (KH - KV) * np.cos(np.deg2rad(el))**2
             * np.cos(np.deg2rad(2 * tau))) / 2.0

        alpha = (KH * alphaH + KV * alphaV + (KH * alphaH - KV * alphaV) *
                 np.cos(np.deg2rad(el))**2 *
                 np.cos(np.deg2rad(2 * tau))) / (2.0 * k)

        return k, alpha


class _ITU838_1_():

    def __init__(self):
        self.__version__ = 1
        self.year = 1999
        self.month = 10
        self.link = 'https://www.itu.int/rec/R-REC-P.838-1-199910-S/en'

    @staticmethod
    def rain_specific_attenuation_coefficients(f, el, tau):
        """
        The frequency-dependent coefficients k and α are given in Table 1 for
        linear polarizations (horizontal: H, vertical: V) and horizontal paths.
        Values of k and α at frequencies other than those in Table 1 can be
        obtained by interpolator using a logarithmic scale for frequency,
        a logarithmic scale for k and a linear scale for α.

        The values in Table 1 have been tested and found to be sufficiently
        accurate for attenuation predictions up to frequencies of 55 GHz.
        """
        _f = [1, 2, 4, 6, 7, 8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70,
              80, 90, 100, 120, 150, 200, 300, 400]

        _kH = [0.0000387, 0.000154, 0.00065, 0.00175, 0.00301, 0.00454, 0.0101,
               0.0188, 0.0367, 0.0751, 0.124, 0.187, 0.263, 0.35, 0.442, 0.536,
               0.707, 0.851, 0.975, 1.06, 1.12, 1.18, 1.31, 1.45, 1.36, 1.32]

        _kV = [0.0000352, 0.000138, 0.000591, 0.00155, 0.00265, 0.00395,
               0.00887, 0.0168, 0.0335, 0.0691, 0.113, 0.167, 0.233, 0.31,
               0.393, 0.479, 0.642, 0.784, 0.906, 0.999, 1.06, 1.13, 1.27,
               1.42, 1.35, 1.31]

        _alphaH = [0.912, 0.963, 1.121, 1.308, 1.332, 1.327, 1.276, 1.217,
                   1.154, 1.099, 1.061, 1.021, 0.979, 0.939, 0.903, 0.873,
                   0.826, 0.793, 0.769, 0.753, 0.743, 0.731, 0.71, 0.689,
                   0.688, 0.683]

        _alphaV = [0.88, 0.923, 1.075, 1.265, 1.312, 1.31, 1.264, 1.2, 1.128,
                   1.065, 1.03, 1, 0.963, 0.929, 0.897, 0.868, 0.824, 0.793,
                   0.769, 0.754, 0.744, 0.732, 0.711, 0.69, 0.689, 0.684]

        KH = np.exp(np.interp(np.log(f), np.log(_f), np.log(_kH)))
        KV = np.exp(np.interp(np.log(f), np.log(_f), np.log(_kV)))
        alphaH = np.interp(np.log(f), np.log(_f), _alphaH)
        alphaV = np.interp(np.log(f), np.log(_f), _alphaV)

        k = (KH + KV + (KH - KV) * np.cos(np.deg2rad(el))**2 *
             np.cos(np.deg2rad(2 * tau))) / 2.0

        alpha = (KH * alphaH + KV * alphaV + (KH * alphaH - KV * alphaV) *
                 np.cos(np.deg2rad(el))**2 *
                 np.cos(np.deg2rad(2 * tau))) / (2.0 * k)

        return k, alpha


class _ITU838_0_():

    def __init__(self):
        self.__version__ = 0
        self.year = 1992
        self.month = 8
        self.link = 'https://www.itu.int/rec/R-REC-P.838-0-199203-S/en'

    @staticmethod
    def rain_specific_attenuation_coefficients(*args, **kwargs):
        return _ITU838_1_.rain_specific_attenuation_coefficients(*args,
                                                                 **kwargs)


__model = __ITU838__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.838 recommendation currently being used.


    This function changes the model used for the ITU-R P.838 recommendation
    to a different version.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          *  3: Activates recommendation ITU-R P.838-3 (03/05) (Current version)
          *  2: Activates recommendation ITU-R P.838-2 (04/03) (Superseded)
          *  1: Activates recommendation ITU-R P.838-1 (10/99) (Superseded)
          *  0: Activates recommendation ITU-R P.838-0 (03/92) (Superseded)

    """
    global __model
    __model = __ITU838__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.838 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__


def rain_specific_attenuation_coefficients(f, el, tau):
    """Compute the values for the coefficients k and α.

    A method to compute the values for the coefficients k and α to compute
    the rain specific attenuation :math:`\\gamma_R` (dB/km) (dB/km)


    Parameters
    ----------
    f : number or Quantity
        Frequency (GHz)
    el : number, sequence, or numpy.ndarray
        Elevation angle of the receiver points
    tau : number, sequence, or numpy.ndarray
        Polarization tilt angle relative to the horizontal (degrees). Tau = 45
        deg for circular polarization)


    Returns
    -------
    k: number
        Coefficient k (non-dimensional)
    α: number
        Coefficient α (non-dimensional)


    References
    ----------
    [1] Rain height model for prediction methods:
    https://www.itu.int/rec/R-REC-P.838/en
    """
    f = prepare_quantity(f, u.GHz, 'Frequency')
    return __model.rain_specific_attenuation_coefficients(f, el, tau)


def rain_specific_attenuation(R, f, el, tau):
    """Compute the specific attenuation γ_R (dB/km) given the rainfall rate.

    A method to compute the specific attenuation γ_R (dB/km) from rain. The
    value is obtained from the rainfall rate R (mm/h) using a power law
    relationship.

    .. math::
        \\gamma_R = k R^\\alpha


    Parameters
    ----------
    R : number, sequence, numpy.ndarray or Quantity
        Rain rate (mm/h)
    f : number or Quantity
        Frequency (GHz)
    el : number, sequence, or numpy.ndarray
        Elevation angle of the receiver points
    tau : number, sequence, or numpy.ndarray
        Polarization tilt angle relative to the horizontal (degrees). Tau = 45
        deg for circular polarization)


    Returns
    -------
    γ_R: numpy.ndarray
        Specific attenuation from rain (dB/km)


    References
    ----------
    [1] Rain height model for prediction methods:
    https://www.itu.int/rec/R-REC-P.838/en
    """
    R = prepare_quantity(R, u.mm / u.hr, 'Rain rate')
    f = prepare_quantity(f, u.GHz, 'Frequency')
    return __model.rain_specific_attenuation(R, f, el, tau) * u.dB / u.km
