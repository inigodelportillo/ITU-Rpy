# -*- coding: utf-8 -*-
import numpy as np
from models.itu1144 import bilinear_2D_interpolator
from iturutils import prepare_input_array, prepare_quantity, load_data,\
                prepare_output_array, dataset_dir
from astropy import units as u


class __ITU453():
    """The radio refractive index: its formula and refractivity data

    Available versions:
    * P.453-12 (07/15) (Current version)

    Recommendation ITU-R P.453 provides methods to estimate the radio refractive
    index and its behaviour for locations worldwide; describes both surface and
    vertical profile characteristics; and provides global maps for the
    distribution of refractivity parameters and their statistical variation.
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.453 recommendation.
    def __init__(self, version = 12):
        if version == 12:
            self.instance = _ITU453_12()
        else:
            raise ValueError('Version ' + str(version) + ' is not implemented'+
            ' for the ITU-R P.453 model.')

    @property
    def __version__(self):
        return self.instance.__version__

    def wet_term_radio_refractivity(self, e, T):
        return self.instance.wet_term_radio_refractivity(e, T)

    def dry_term_radio_refractivity(self, Pd, T):
        return self.instance.dry_term_radio_refractivity(Pd, T)

    def radio_refractive_index(self, P, e, T):
        return self.instance.radio_refractive_index(P, e, T)

    def water_vapour_pressure(self, T, P, H, type_hydrometeor='water'):
        return self.instance.water_vapour_pressure(T, P, H, type_hydrometeor='water')

    def saturation_vapour_pressure(self, T, P, type_hydrometeor='water'):
        return self.instance.saturation_vapour_pressure(T, P, type_hydrometeor='water')

    def map_wet_term_radio_refractivity(self, lat, lon):
        return self.instance.map_wet_term_radio_refractivity(lat, lon)

    def DN65(self, lat, lon, p):
        return self.instance.DN65(lat, lon, p)

    def DN1(self, lat, lon, p):
        return self.instance.DN1(lat, lon, p)

class _ITU453_12():

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
            ps = [ 0.1, 0.2, 0.5,  1, 2, 3, 5, 10, 20, 30, 50, 60, 70, 80,
                  90, 95, 98, 99, 99.5, 99.8, 99.9]
            d_dir = dataset_dir + '453/DN65m_%02dd%02d_v1.txt'
            lats = load_data(dataset_dir + '453/lat0d75.txt')
            lons = load_data(dataset_dir + '453/lon0d75.txt')
            for p_loads in ps:
                int_p = p_loads//1
                frac_p = p_loads % 1
                vals = load_data(d_dir%(int_p, frac_p))
                self._DN65[float(p_loads)] = bilinear_2D_interpolator(lats, lons, vals)

        return self._DN65[float(p)](np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def DN1(self, lat, lon, p):
        if not self._DN1:
            ps = [ 0.1, 0.2, 0.5,  1, 2, 3, 5, 10, 20, 30, 50, 60, 70, 80,
                  90, 95, 98, 99, 99.5, 99.8, 99.9]
            d_dir = dataset_dir + '453/DN_%02dd%02d_v1.txt'
            lats = load_data(dataset_dir + '453/lat0d75.txt')
            lons = load_data(dataset_dir + '453/lon0d75.txt')
            for p_loads in ps:
                int_p = p_loads//1
                frac_p = p_loads % 1
                vals = load_data(d_dir%(int_p, frac_p))
                self._DN1[float(p_loads)] = bilinear_2D_interpolator(lats, lons, vals)

        return self._DN1[float(p)](np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)


    def N_wet(self, lat, lon):
        if not self._N_wet:
            vals = load_data(dataset_dir + '453/v12_ESANWET.txt')
            lats = load_data(dataset_dir + '453/v12_ESALAT.txt')
            lons = load_data(dataset_dir + '453/v12_ESALON.txt')
            self._N_wet = bilinear_2D_interpolator(lats, lons, vals)

        return self._N_wet(np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def wet_term_radio_refractivity(self, e, T):
        N_wet = (72 * e / (T + 273.15) + 3.75e5 * e / (T + 273.15)**2) * 1e-6
        return N_wet

    def dry_term_radio_refractivity(self, Pd, T):
        N_dry = 77.6 * Pd/T
        return N_dry

    def radio_refractive_index(self, P, e, T):
        N = 77.6 * P/T  - 5.6 * e/T + 3.75e5 * e/T**2     # Eq. 2    [N-units]
        n = 1 + N * 1e-6
        return n

    def water_vapour_pressure(self, T, P, H, type_hydrometeor='water'):
        e_s = self.saturation_vapour_pressure(T, P, type_hydrometeor)
        return H * e_s / 100

    def saturation_vapour_pressure(self, T, P, type_hydrometeor='water'):

        if type_hydrometeor == 'water':
            EF = 1 + 1e-4*(7.2 + P*(0.00320 + 5.9e-6 * T**2))
            a = 6.1121
            b = 18.678
            c = 257.14
            d = 234.5

        elif type_hydrometeor == 'ice':
            EF = 1 + 1e-4*(2.2 + P*(0.00382 + 6.4e-7 * T**2))
            a = 6.11215
            b = 23.036
            c = 279.82
            d = 333.7

        e_s = EF * a * np.exp((b - T/d)*T/(T+c))
        return e_s

    def map_wet_term_radio_refractivity(self, lat, lon):
        return self.N_wet(lat, lon)


__model = __ITU453()

def change_version(new_version):
    """
    Change the version of the ITU-R P.453 recommendation currently being used.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
        * P.453-12 (02/12) (Current version)
    """
    global __model
    __model = __ITU453(new_version)

def get_version():
    """
    Obtain the version of the ITU-R P.453 recommendation currently being used.
    """
    global __model
    return __model.__version__

def wet_term_radio_refractivity(e, T):
    """
    Method to determine the wet term of the radio refractivity

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

    References:
    -----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    global __model
    e = prepare_quantity(e, u.hPa, 'Water vapour pressure ')
    T = prepare_quantity(T, u.K, 'Absolute temperature')
    val = __model.wet_term_radio_refractivity(e, T)
    return val*u.dimensionless_unscaled

def dry_term_radio_refractivity(Pd, T):
    """
    Method to determine the dry term of the radio refractivity

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

    References:
    -----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    global __model
    Pd = prepare_quantity(Pd, u.hPa, 'Dry atmospheric pressure')
    T = prepare_quantity(T, u.K, 'Absolute temperature')
    val = __model.dry_term_radio_refractivity(Pd, T)
    return val*u.dimensionless_unscaled

def radio_refractive_index(P, e, T):
    """
    Method to compute the radio refractive index

    Parameters
    ----------
    P : number or Quantity
        Total atmospheric pressure (hPa)
    e : number or Quantity
        Water vapour pressure  (hPa)
    T : number or Quantity
        Absolute temperature (K)

    Returns
    -------
    n: Quantity
        Radio refractive index (-)

    References:
    -----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    global __model
    P = prepare_quantity(P, u.hPa, 'Total atmospheric pressure')
    e = prepare_quantity(e, u.hPa, 'Water vapour pressure ')
    T = prepare_quantity(T, u.K, 'Absolute temperature')
    val = __model.radio_refractive_index(P, e, T)
    return val*u.dimensionless_unscaled


def water_vapour_pressure(T, P, H, type_hydrometeor='water'):
    """
    Method to determine the water vapour pressure

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

    References:
    -----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    global __model
    T = prepare_quantity(T, u.C, 'Absolute temperature')
    P = prepare_quantity(P, u.hPa, 'Total atmospheric pressure')
    H = prepare_quantity(H, u.percent, 'Total atmospheric pressure')
    val = __model.water_vapour_pressure(T, P, H,type_hydrometeor)
    return val*u.hPa



def saturation_vapour_pressure(T, P, type_hydrometeor='water'):
    """
    Method to determine the saturation water vapour pressure

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

    References:
    -----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    global __model
    T = prepare_quantity(T, u.C, 'Absolute temperature')
    P = prepare_quantity(P, u.hPa, 'Total atmospheric pressure')
    val = __model.saturation_vapour_pressure(T, P, type_hydrometeor)
    return val*u.hPa



def map_wet_term_radio_refractivity(lat, lon):
    """
    Method to determine the wet term of the radio refractivity

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

    References:
    -----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en
    """
    global __model
    type_output = type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.map_wet_term_radio_refractivity(lat, lon)
    return prepare_output_array(val, type_output) * u.g / u.m**3


def DN65(lat, lon, p):
    """
    Method to determine the statistics of the vertical gradient of radio
    refractivity in the lowest 65 m from the surface of the Earth.

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

    References:
    -----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en

    """
    global __model
    type_output = type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.DN65(lat, lon, p)
    return prepare_output_array(val, type_output) * u.one


def DN1(lat, lon, p):
    """
     Method to determine the statistics of the vertical gradient of radio
    refractivity over a 1 km layer from the surface.

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

    References:
    -----------
    [1] The radio refractive index: its formula and refractivity data
    https://www.itu.int/rec/R-REC-P.453/en
    """
    global __model
    type_output = type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.DN1(lat, lon, p)
    return prepare_output_array(val, type_output) * u.one
