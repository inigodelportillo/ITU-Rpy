# -*- coding: utf-8 -*-
"""
ITU-RPy is a python implementation of the ITU-P R Recommendations.

ITU-Rpy can be used to compute atmospheric attenuation for Earth-to-space
and horizontal paths, for frequencies in the GHz range.

The propagation loss on an Earth-space path and a horizontal-path, relative to
the free-space loss, is the sum of different contributions, namely:
 * attenuation by atmospheric gases;
 * attenuation by rain, other precipitation and clouds;
 * scintillation and multipath effects;
 * attenuation by sand and dust storms.

Each of these contributions has its own characteristics as a function of
frequency, geographic location and elevation angle. ITU-Rpy allows for fast,
vectorial computation of the different contributions to the atmospheric
attenuation.
"""
from __future__ import absolute_import, division, print_function

__all__ = ["utils", "plotting"]
import warnings

import astropy.units as u
import numpy as np

import itur.plotting
import itur.utils

from .__version__ import __version__
from .models.itu618 import rain_attenuation, scintillation_attenuation
from .models.itu676 import (
    gaseous_attenuation_inclined_path,
    gaseous_attenuation_slant_path,
    gaseous_attenuation_terrestrial_path,
)
from .models.itu835 import standard_pressure
from .models.itu836 import surface_water_vapour_density, total_water_vapour_content
from .models.itu840 import cloud_attenuation
from .models.itu1510 import surface_mean_temperature
from .models.itu1511 import topographic_altitude

# Ignore divide by zero errors
np.seterr(divide="ignore")

AUTHORS = "Inigo del Portillo"


def atmospheric_attenuation_slant_path(
    lat,
    lon,
    f,
    el,
    p,
    D,
    hs=None,
    rho=None,
    R001=None,
    eta=0.5,
    T=None,
    H=None,
    P=None,
    hL=1e3,
    Ls=None,
    tau=45,
    V_t=None,
    mode="approx",
    return_contributions=False,
    include_rain=True,
    include_gas=True,
    include_scintillation=True,
    include_clouds=True,
):
    """
    Calculate long-term atmospheric attenuation statistics for slant paths.

    This function provides estimates of the long-term statistics of
    the slant-path atmospheric attenuation at a given location, for
    frequencies up to 55 GHz and percentages of time 0.001% < `p` < 50%.

    The model used is based on the guidelines provided in Section 2 of
    ITU-R P.618. If optional values are not provided they will be
    automatically computed using the procedures described in other ITU-R P.
    recommendations.

    Parameters
    ----------
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points
    f : number or Quantity
        Frequency (GHz)
    el : sequence, number or Quantity
        Elevation angle (degrees)
    p : number
        Percentage of the time the rain attenuation value is exceeded.
    D: number or Quantity
        Physical diameter of the earth-station antenna (m)
    hs : number, sequence, or numpy.ndarray, optional
        Height above mean sea level of the earth station (km). If local data for
        the earth station height above mean sea level is not available, an
        estimate is obtained from the maps of topographic altitude
        given in Recommendation ITU-R P.1511.
    rho : number or Quantity, optional
        Water vapor density (g/m3). If not provided, an estimate is obtained
        from Recommendation Recommendation ITU-R P.836.
    R001: number  or Quantity, optional
        Point rainfall rate for the location for 0.01% of an average year \
        (mm/h). If not provided, an estimate is obtained from Recommendation
        ITU-R P.837. Some useful values:
            * 0.25 mm/h : Drizzle
            *  2.5 mm/h : Light rain
            * 12.5 mm/h : Medium rain
            * 25.0 mm/h : Heavy rain
            * 50.0 mm/h : Downpour
            * 100  mm/h : Tropical
            * 150  mm/h : Monsoon
    eta: number, optional
        Antenna efficiency. Default value 0.5 (conservative estimate)
    T: number, sequence, or numpy.ndarray, optional
        Average surface ambient temperature (°C) at the site. If None, uses the
        ITU-R P.1510 to estimate the surface ambient temperature.
    H: number, sequence, or numpy.ndarray, optional
        Average surface relative humidity (%) at the site. If None, uses the
        ITU-R P.836 to estimate the wet term of the surface relative humidity.
    P: number, sequence, or numpy.ndarray, optional
        Average surface pressure (hPa) at the site. If None, uses the
        ITU-R P.835 to estimate the average surface pressure.
    hL : number, optional
        Height of the turbulent layer (m). Default value 1000 m
    Ls :number, optional
        Slant path length (km). If not provided, it will be computed using the
        rain height and the elevation angle. The ITU model does not require
        this parameter as an input.
    tau : number, optional
        Polarization tilt angle relative to the horizontal (degrees)
        (tau = 45 deg for circular polarization). Default value is 45
    V_t : number or Quantity, optional
        Integrated water vapour content along the path (kg/m2 or mm).
        If not provided this value is estimated using Recommendation
        ITU-R P.836. Default value None
    mode : string, optional
        Mode for the calculation of gaseous attenuation. Valid values are
        'approx', 'exact'. If 'approx' Uses the method in Annex 2 of
        Recommendation ITU-R P.676, else uses the method described in
        Section 1. Default, 'approx'
    return_contributions: bool, optional
        Determines whether individual contributions from gases, rain, clouds
        and scintillation are returned in addition to the total attenuation
        (True), or just the total atmospheric attenuation (False).
        Default is False
    include_rain: bool, optional
        Determines whether to include the rain contribution in the total
        atmospheric attenuation calculation or not. Default is True
    include_gas: bool, optional
        Determines whether to include the gaseous contribution in the total
        atmospheric attenuation calculation or not. Default is True
    include_scintillation: bool, optional
        Determines whether to include the scintillation contribution in the
        total atmospheric attenuation calculation or not. Default is True
    include_clouds: bool, optional
        Determines whether to include the clouds contribution in the total
        atmospheric attenuation calculation or not. Default is True


    Returns
    -------
    A : Quantity
        Total atmospheric attenuation (dB)

    Ag, Ac, Ar, As, A : tuple
        Gaseous, Cloud, Rain, Scintillation contributions to total attenuation,
        and total attenuation (dB)



    References
    ----------
    [1] Propagation data and prediction methods required for the design of
    Earth-space telecommunication systems:
    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.618-12-201507-I!!PDF-E.pdf


    """
    if np.logical_or(p < 0.001, p > 50).any():
        warnings.warn(
            RuntimeWarning(
                "The method to compute the total "
                "atmospheric attenuation in recommendation ITU-P 618-13 "
                "is only recommended for unavailabilities (p) between "
                "0.001% and 50 %"
            )
        )

    # This takes account of the fact that a large part of the cloud attenuation
    # and gaseous attenuation is already included in the rain attenuation
    # prediction for time percentages below 1%. Eq. 64 and Eq. 65 in
    # Recommendation ITU 618-12
    p_c_g = np.maximum(1, p)

    # Estimate the ground station altitude
    if hs is None:
        hs = topographic_altitude(lat, lon)

    # Surface mean temperature
    if T is None:
        T = surface_mean_temperature(lat, lon)

    # Estimate the surface Pressure
    if P is None:
        P = standard_pressure(hs)

    # Estimate the surface Pressure
    if V_t is None:
        V_t = total_water_vapour_content(lat, lon, p_c_g, hs)

    # Estimate the surface water vapour density
    if rho is None:
        rho = surface_water_vapour_density(lat, lon, p_c_g, hs)

    # Compute the attenuation components
    if include_rain:
        Ar = rain_attenuation(lat, lon, f, el, hs, p, R001, tau, Ls)
    else:
        Ar = 0 * u.dB

    if include_gas:
        Ag = gaseous_attenuation_slant_path(f, el, rho, P, T, V_t, hs, mode)
    else:
        Ag = 0 * u.dB

    if include_clouds:
        Ac = cloud_attenuation(lat, lon, el, f, p_c_g)
    else:
        Ac = 0 * u.dB

    if include_scintillation:
        As = scintillation_attenuation(lat, lon, f, el, p, D, eta, T, H, P, hL)
    else:
        As = 0 * u.dB

    # Compute the total attenuation according to
    A = Ag + np.sqrt((Ar + Ac) ** 2 + As ** 2)

    if return_contributions:
        return Ag, Ac, Ar, As, A
    else:
        return A
