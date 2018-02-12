# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from satcomm_utils import RE

"""
Method to compute the antenna gain in angles different from the boresight. 
The method is based on the 3rd one found in Recommendation ITU-R S.1528-0. 

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


def calculate_psi(lat_boresight, lon_boresight, lat, lon, a):
    lat = convert_nparray(lat)
    lon = convert_nparray(lon)

    lat_boresight = np.deg2rad(lat_boresight)
    lon_boresight = np.deg2rad(lon_boresight)
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    radii = a / RE
    delta_lat = lat_boresight - lat
    delta_lon = lon_boresight - lon
    beta = np.arccos(np.cos(delta_lat) * np.cos(delta_lon))
    ret = np.arctan2(np.sin(beta), (radii - np.cos(beta)))
    return np.rad2deg(ret)


def calculate_gain(lat_boresight, lon_boresight, lat, lon, a, D, l, G_m,
                   L_n=-25, L_f=3, z=1, psi_b=None, psi=None):
    """
    Calculates the gain of a multiple-beam satellite antenna having either
    circular or elliptical beams. Method #2 of the recommendation.


    Parameters
    -----------
    lat: numpy.ndarray or list or float
        Latitude coordinates of
    lon: numpy.ndarray or list or float
        Longitude coordinates of


    Returns
    --------
    gain: numpy.ndarray
        Antenna gain...
    """

    if psi is None:
        psi = calculate_psi(lat_boresight, lon_boresight, lat, lon, a)
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
    lat_boresight = [20, 18]
    lon_boresight = [10, 12]
    lat = 20
    lon = 10
    a = 8000
    epsilon = 0.4
    freq = 30e9
    l = 299792458/freq
    D = 70*l/(2*epsilon)
    G_m = 35
    print(calculate_gain(lat_boresight, lon_boresight, lat, lon, a, D, l, G_m,
                         psi_b = epsilon))

    psi = np.linspace(0, 12)  # 28.8
    G = calculate_gain(0, 0, 0, 0, 0, D, l, G_m, psi_b=epsilon, psi=psi)
    plt.plot(psi, G)
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