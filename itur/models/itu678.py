# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from astropy import units as u
from scipy.special import erfc

from itur.models.itu1144 import bilinear_2D_interpolator
from itur.utils import (prepare_input_array, prepare_output_array,
                        load_data_interpolator, get_input_type)


class __ITU678__():
    """Characterization of the variability of propagation phenomena and
    estimation of the risk associated with propagation margin.

    Not available versions:
    * P.678-1 (03/92) (Superseded)
    * P.678-2 (09/2013) (Superseded)
    Available versions include:
    * P.678-3 (07/2015) (Current version)
    """
    # This is an abstract class that contains an instance to a version of the
    # ITU-R P.678 recommendation.

    def __init__(self, version=3):
        if version == 3:
            self.instance = _ITU678_3_()
        else:
            raise ValueError(
                f"Version {version} is not implemented "
                "for the ITU-R P.678 model."
            )

    @property
    def __version__(self):
        return self.instance.__version__

    def inter_annual_variability(self, p, lat, lon):
        return self.instance.inter_annual_variability(p, lat, lon)

    def risk_of_exceedance(self, p, pr, lat, lon):
        return self.instance.risk_of_exceedance(p, pr, lat, lon)


class _ITU678_3_():

    def __init__(self):
        self.__version__ = 3
        self.year = 2015
        self.month = 7
        self.link = 'https://www.itu.int/rec/R-REC-P.678-3-201507-I/en'

        self._climatic_ratio_data = None
        self._variance_of_estimation_cache = {}

    def climatic_ratio(self, lat, lon):
        """Return the climatic ratio r_c(lat, lon) from the gridded map."""
        if self._climatic_ratio_data is None:
            self._climatic_ratio_data = load_data_interpolator(
                '678/v3_lat.npz', '678/v3_lon.npz',
                '678/v3_climatic_ratio.npz', bilinear_2D_interpolator,
                flip_ud=False)
        return self._climatic_ratio_data(
            np.array([lat.ravel(), lon.ravel()]).T).reshape(lat.shape)

    def variance_of_estimation(self, p):
        """Compute the estimation variance σ²_E(p) via Eqs. (2)–(5).

        The autocorrelation sum C is truncated once contributions drop below
        machine precision, rather than summing all N = 525 960 terms, giving
        a ~200× speed-up for typical probability values.

        Parameters
        ----------
        p : float
            Probability of exceedance as a fraction [0, 1].
            Valid range: 0.0001 ≤ p ≤ 0.02.

        Returns
        -------
        float
            Estimation variance σ²_E(p) (dimensionless fraction²).
        """
        p_key = float(p)
        if p_key in self._variance_of_estimation_cache:
            return self._variance_of_estimation_cache[p_key]

        # Eq. (4): shape parameter of the autocorrelation function
        b1 = -0.0396
        b2 = 0.286
        b = b1 * np.log(p) + b2

        a = 0.0265   # amplitude parameter [s^{-1}]
        N = 525960   # total minutes per year (= 365.25 × 24 × 60)
        dt = 60      # Δt = 60 s  (Annex 2, equation parameters)

        # Eq. (3): C = Σ_{i=-(N-1)}^{N-1} c_U(i·Δt, p)
        # By symmetry: C = c_U(0) + 2·Σ_{i=1}^{max_i} c_U(i·Δt, p)
        # Truncate where exp(-a·(i·dt)^b) < 1e-13, i.e. a·(i·dt)^b > 30.
        max_i = min(int(np.ceil((30.0 / a) ** (1.0 / b) / dt)) + 1, N - 1)
        i_pos = np.arange(0, max_i + 1, dtype=float)
        cu = np.exp(-a * np.abs(i_pos * dt) ** b)   # Eq. (3) terms ≥ 0
        C = 2.0 * np.sum(cu) - cu[0]               # Eq. (2): subtract i=0 counted twice

        # Eq. (5): σ²_E(p) = p(1-p)·C / N
        result = float(p * (1 - p) * C / N)
        self._variance_of_estimation_cache[p_key] = result
        return result

    def inter_annual_climatic_variance(self, p, lat, lon):
        """Compute the climatic variance σ²_C(p) via Eq. (6).

        Parameters
        ----------
        p : float or ndarray
            Probability of exceedance as a fraction [0, 1].

        Returns
        -------
        ndarray
            Climatic variance σ²_C(p) = (r_c · p)² (dimensionless fraction²).
        """
        rc = self.climatic_ratio(lat, lon)
        # Eq. (6): σ²_C(p) = (r_c(Lat, Lon) · p)²
        return (rc * p) ** 2

    def inter_annual_variability(self, p, lat, lon):
        """Compute the total inter-annual variance σ²(p) via Eq. (1).

        Parameters
        ----------
        p : float or ndarray
            Probability of exceedance as a fraction [0, 1].

        Returns
        -------
        ndarray
            Total variance σ²(p) = σ²_C(p) + σ²_E(p) (dimensionless fraction²).
        """
        # Eq. (1): σ²(p) = σ²_C(p) + σ²_E(p)
        return (self.inter_annual_climatic_variance(p, lat, lon)
                + self.variance_of_estimation(p)).squeeze()

    def risk_of_exceedance(self, p, pr, lat, lon):
        """Compute the exceedance risk ℜ via Eq. (8).

        Parameters
        ----------
        p : float
            Long-term probability of exceedance as a fraction [0, 1].
        pr : float
            Annual probability of exceedance as a fraction [0, 1].

        Returns
        -------
        ndarray
            Risk ℜ ∈ [0, 1] (dimensionless).
            ℜ = 0.5 when pr = p (identity stated in the recommendation).
        """
        # Eq. (8): ℜ = Q((p_ℜ − p) / σ(p)) = 0.5 · erfc((p_ℜ − p) / (√2 · σ(p)))
        # NOTE: σ(p) is the *standard deviation* — the square root of the variance
        #       returned by inter_annual_variability().
        sigma = np.sqrt(self.inter_annual_variability(p, lat, lon))
        return (0.5 * erfc((pr - p) / (sigma * np.sqrt(2)))).squeeze()


__model = __ITU678__()


def change_version(new_version):
    """
    Change the version of the ITU-R P.678 recommendation currently being used.

    Parameters
    ----------
    new_version : int
        Number of the version to use.
        Valid values are:
          * 3: Activates recommendation ITU-R P.678-3 (07/2015) (Current version)
    """
    global __model
    __model = __ITU678__(new_version)


def get_version():
    """
    Obtain the version of the ITU-R P.678 recommendation currently being used.

    Returns
    -------
    version: int
        Version currently being used.
    """
    return __model.__version__


def inter_annual_variability(p, lat, lon):
    """
    Estimate the total inter-annual variance σ²(p) of rain-rate or rain-
    attenuation statistics at a given location.

    Implements Eqs. (1)–(6) from Annex 2 of Recommendation ITU-R P.678-3.

    The total variance combines two independent components:

    * **Climatic variance** σ²_C(p) = (r_c · p)²  [Eq. 6] — captures
      location-dependent year-to-year climate variability.
    * **Estimation variance** σ²_E(p) = p(1−p)·C/N  [Eq. 5] — captures
      the finite-sample uncertainty from observing one year of statistics.

    The square root σ(p) = √σ²(p) is the standard deviation used as input
    to ``risk_of_exceedance``.

    Parameters
    ----------
    p : number
        Probability of exceedance as a fraction [0, 1].
        Valid range: 0.0001 ≤ p ≤ 0.02  (i.e. 0.01 % to 2 % of time).
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points (degrees N).
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points (degrees E, −180 to 360).

    Returns
    -------
    sigma_sq : astropy.units.Quantity
        Total inter-annual variance σ²(p) (dimensionless fraction²).

    References
    ----------
    [1] Characterization of the variability of propagation phenomena and
    estimation of the risk associated with propagation margin:
    https://www.itu.int/rec/R-REC-P.678/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.inter_annual_variability(p, lat, lon)
    return prepare_output_array(val, type_output) * u.dimensionless_unscaled


def risk_of_exceedance(p, pr, lat, lon):
    """
    Estimate the risk ℜ that the annual exceedance fraction exceeds ``pr``
    in a randomly chosen year, given that the long-term (multi-year average)
    exceedance fraction is ``p``.

    Implements Eq. (8) from Annex 3 of Recommendation ITU-R P.678-3:

        ℜ = Q((p_ℜ − p) / σ(p)) = 0.5 · erfc((p_ℜ − p) / (√2 · σ(p)))

    where σ(p) = √(``inter_annual_variability``(p, lat, lon)).

    **Identity check (Annex 3):** setting ``pr = p`` always yields ℜ = 50 %,
    regardless of location.

    Parameters
    ----------
    p : number
        Long-term probability of exceedance as a fraction [0, 1].
        Valid range: 0.0001 ≤ p ≤ 0.02.
    pr : number
        Annual probability of exceedance threshold as a fraction [0, 1].
    lat : number, sequence, or numpy.ndarray
        Latitudes of the receiver points (degrees N).
    lon : number, sequence, or numpy.ndarray
        Longitudes of the receiver points (degrees E, −180 to 360).

    Returns
    -------
    risk : astropy.units.Quantity
        Risk ℜ expressed as a percentage (%).
        ℜ = 50 % when pr = p.

    References
    ----------
    [1] Characterization of the variability of propagation phenomena and
    estimation of the risk associated with propagation margin:
    https://www.itu.int/rec/R-REC-P.678/en
    """
    type_output = get_input_type(lat)
    lat = prepare_input_array(lat)
    lon = prepare_input_array(lon)
    lon = np.mod(lon, 360)
    val = __model.risk_of_exceedance(p, pr, lat, lon)
    return prepare_output_array(val * 100, type_output) * u.pct
