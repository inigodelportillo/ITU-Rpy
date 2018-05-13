# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from satcomm_utils import RE, LLA_to_ECEF


"""
Method to compute the antenna gain in angles different from the boresight. 
The method is based on the 2nd one found in Recommendation ITU-R S.1528-0. 

References
--------
[1] Satellite antenna radiation patterns for non-geostationary orbit satellite
antennas operating in the fixed-satellite service below 30 GHz: 
https://www.itu.int/rec/R-REC-S.1528/en
"""


def convert_nparray(arr):
    if not isinstance(arr, np.ndarray):
        if isinstance(arr, list):
            arr = np.array(arr)
        else:
            arr = np.array([arr])
    return arr


def calculate_psi(lat_boresight, lon_boresight, lat_sat, lon_sat,
                  lat, lon, a):
    """
    Calculates the angle (psi, in deg) from the direction of the beam boresight
    to the points defined by (lat, lon) as seen from the satellite defined by
    (lat_sat, lon_sat, a)


    Parameters
    -----------
    lat_boresight: float
        Latitude coordinate of the beam boresight
    lon_boresight: float
        Longitude coordinate of the beam boresight
    lat_sat: float
        Latitude coordinate of the satellite
    lon_sat: float
        Longitude coordinate of the satellite
    lat: numpy.ndarray or list or float
        Latitude coordinates of the points for which psi will be calculated
    lon: numpy.ndarray or list or float
        Longitude coordinates of the points for which psi will be calculated
    a: float
        Semi-major axis of the satellite (km)


    Returns
    --------
    psi: numpy.ndarray
        Angles from the boresight to the points (deg)
    """
    lat = convert_nparray(lat)
    lon = convert_nparray(lon)

    boresight = LLA_to_ECEF([lon_boresight, lat_boresight, 0])
    sat = LLA_to_ECEF([lon_sat, lat_sat, (a - RE)*1000])
    points = LLA_to_ECEF([lon, lat, 0])  # [0]*len(lon))
    vec_b = boresight - sat
    vec_p = points - sat
    psi = np.arccos(np.dot(vec_b, vec_p) /
                    (np.linalg.norm(vec_b)*np.linalg.norm(vec_p)))
    return np.rad2deg(psi)


def calculate_gain(lat_boresight, lon_boresight, lat_sat, lon_sat, lat, lon, a,
                   G_m, D=None, l=None, L_n=-25, L_f=3, z=1, psi_b=None,
                   psi=None):
    """
    Calculates the gain pattern of a multiple-beam satellite antenna having
    circular or elliptical beams. Method #2 of the recommendation.


    Parameters
    -----------
    lat_boresight: float
        Latitude coordinate of the beam boresight
    lon_boresight: float
        Longitude coordinate of the beam boresight
    lat_sat: float
        Latitude coordinate of the satellite
    lon_sat: float
        Longitude coordinate of the satellite
    lat: numpy.ndarray or list or float
        Latitude coordinates of the points for which psi will be calculated
    lon: numpy.ndarray or list or float
        Longitude coordinates of the points for which psi will be calculated
    a: float
        Semi-major axis of the satellite (km)
    G_m: float
        Maximum gain in the boresight, the main lobe (dB)
    D: float
        Diameter of the antenna (m), used to calculate psi_b if not given
    l: float
        Wavelength of the lowest band edge of interest(m), used to calculate
        psi_b if not given
    L_n: float
        near-in-side-lobe level (dB) relative to the peak gain required
        by the system design. L_n should be in [-15, -20, -25, -30] if z != 1
    L_f: float
        0 dBi far-out side-lobe level (dBi)
    z: float
        (major axis/minor axis) for the radiated beam
    psi_b: float
        One-half the 3 dB beamwidth (3 dB below G_m) (deg)
    psi: np.array
        Angles  (deg) from the boresight to the points on which gain is to be
        calculated. If it is given, lats, lons and 'a' are ignored.


    Returns
    --------
    gain: numpy.ndarray
        Antenna gain in the direction specified by the input parameters (dB)
    """

    if psi is None:
        psi = calculate_psi(lat_boresight, lon_boresight, lat_sat, lon_sat,
                            lat, lon, a)
    else:
        psi = convert_nparray(psi)

    assert all(0 <= angle <= 90 for angle in psi)

    if psi_b is None:
        psi_b = np.sqrt(1200)/(D/l)

    if z == 1:
        a = 2.58
    elif L_n == -15:
        a = 2.58*np.sqrt(1-1.4*np.log10(z))
    elif L_n == -20:
        a = 2.58*np.sqrt(1-np.log10(z))
    elif L_n == -25:
        a = 2.58*np.sqrt(1-0.6*np.log10(z))
    elif L_n == -30:
        a = 2.58*np.sqrt(1-0.4*np.log10(z))
    else:
        raise ValueError('L_n should be in [-15, -20, -25, -30] if z != 1')
    b = 6.32
    alpha = 1.5

    X = G_m + L_n + 25*np.log10(b*psi_b)
    Y = b*psi_b*np.power(10, 0.04*(G_m + L_n - L_f))

    return np.where(psi <= a*psi_b, G_m - 3*(psi/psi_b)**alpha,
           np.where(psi <= 0.5*b*psi_b, G_m + L_n + 20*np.log10(z),
           np.where(psi <= b*psi_b, G_m + L_n,
           np.where(psi <= Y, X - 25*np.log10(psi), L_f))))


if __name__ == '__main__':
    lat_lon_boresight = [2, 0]
    lat_sat = 0
    lon_sat = 1
    a = 2*RE
    # lat_lon = [20.1, 9.9]
    lat = 0
    lon = 3
    epsilon = 0.4
    freq = 30e9
    G_m = 35
    print(calculate_gain(*lat_lon_boresight, lat_sat, lon_sat,
                         lat, lon, a, G_m, psi_b=epsilon))

    psi = np.linspace(0, 12)
    G = calculate_gain(0, 0, 0, 0, 0, 0, 0, G_m, psi_b=epsilon, psi=psi)
    plt.plot(psi, G)
    plt.title('Radiation pattern envelope function')
    plt.xlabel('Off-axis angle (deg)')
    plt.ylabel('Gain (dB)')
    plt.grid()
    plt.show()

# psi = np.rad2deg(psi)
# psi_b = 1.6  # deg
# G_m = 35
# L_s = -12
# Y = 2*psi_b  # Y = psi_b*(-L_s/3)^0.5
# # L_f = 3
# Z = 20  # deg  Y*np.power(10, 0.04*(G_m + L_s - L_f))
# return np.where(psi <= Y, G_m - 3*(psi/psi_b)**2,
#                 (np.where(psi <= Z, G_m + L_s - 25*np.log10(psi/Y), 3)))
