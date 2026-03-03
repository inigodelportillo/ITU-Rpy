# -*- coding: utf-8 -*-
"""
Validation tests for ITU-R P.676-13 (08/2022).

Test data sourced from:
  CG-3M3J-13-ValEx-Rev8.3.0.xlsx
  (ITU-R Study Group 3 validation examples)

Sheets used:
  - P.676-13 SpAtt          : Specific attenuation (Eqs. 1–9, Annex 1)
  - P.676-13 A_Gas_A1_2.2.1a: Slant path exact (Annex 1, Eqs. 13–19b)
  - P.676-13 A_Gas_A2_INST  : Slant path approximate (Annex 2, Eqs. 31–37)

Pressure convention:
  - gamma0_exact / gammaw_exact use DRY air pressure p (hPa)
  - slant_inclined_path_equivalent_height Eq. 31 uses TOTAL surface pressure
    Ps = p + e, where e = rho*T/216.7.  The function accepts p (dry) and
    computes Ps internally.
"""
import unittest
import itur.models.itu676 as itu676


class TestITU676_13_SpAtt(unittest.TestCase):
    """Specific attenuation at standard conditions (Annex 1, Eqs. 1–9).

    Standard conditions from SpAtt sheet:
      p = 1013.25 hPa  (dry air pressure, ITU convention)
      e = 9.97289 hPa  (water vapour partial pressure, computed from rho)
      T = 288.15 K
      rho = 7.5 g/m³
    """

    def setUp(self):
        itu676.change_version(13)
        self.p = 1013.25     # dry air pressure [hPa]
        self.rho = 7.5       # water vapour density [g/m³]
        self.T = 288.15      # temperature [K]

    def test_gamma0_exact_12GHz(self):
        self.assertAlmostEqual(
            itu676.gamma0_exact(12, self.p, self.rho, self.T).value,
            0.008698264068773570, places=5)

    def test_gamma0_exact_20GHz(self):
        self.assertAlmostEqual(
            itu676.gamma0_exact(20, self.p, self.rho, self.T).value,
            0.011883550477807600, places=5)

    def test_gamma0_exact_60GHz(self):
        self.assertAlmostEqual(
            itu676.gamma0_exact(60, self.p, self.rho, self.T).value,
            14.623474796486100, places=4)

    def test_gamma0_exact_90GHz(self):
        self.assertAlmostEqual(
            itu676.gamma0_exact(90, self.p, self.rho, self.T).value,
            0.038869711072423500, places=5)

    def test_gamma0_exact_130GHz(self):
        self.assertAlmostEqual(
            itu676.gamma0_exact(130, self.p, self.rho, self.T).value,
            0.041509083599522800, places=5)

    def test_gammaw_exact_12GHz(self):
        self.assertAlmostEqual(
            itu676.gammaw_exact(12, self.p, self.rho, self.T).value,
            0.009535388220245930, places=5)

    def test_gammaw_exact_20GHz(self):
        self.assertAlmostEqual(
            itu676.gammaw_exact(20, self.p, self.rho, self.T).value,
            0.097047304815111700, places=5)

    def test_gammaw_exact_60GHz(self):
        self.assertAlmostEqual(
            itu676.gammaw_exact(60, self.p, self.rho, self.T).value,
            0.154841841000000, places=5)

    def test_gammaw_exact_90GHz(self):
        self.assertAlmostEqual(
            itu676.gammaw_exact(90, self.p, self.rho, self.T).value,
            0.341973394000000, places=5)

    def test_gammaw_exact_130GHz(self):
        self.assertAlmostEqual(
            itu676.gammaw_exact(130, self.p, self.rho, self.T).value,
            0.751844704000000, places=5)


class TestITU676_13_SlantPathExact(unittest.TestCase):
    """Slant path attenuation using Annex 1 exact layer-by-layer method.

    From A_Gas_A1_2.2.1a sheet, Example 1:
      f = 28 GHz, h1 = 0 km, elevation = 30 deg
      Uses P.835-6 standard atmosphere with surface rho = 7.5 g/m³.
      Ray tracing via Eqs. 13, 14, 15, 17, 18b, 19a.
      Expected A_gas = 0.47081173472870474 dB
    """

    def setUp(self):
        itu676.change_version(13)

    def test_slant_path_exact_28GHz_30deg(self):
        # Surface rho = 7.5 g/m³; p = 1013.25 hPa (dry, standard atmosphere)
        A_gas = itu676.gaseous_attenuation_slant_path(
            28, 30, 7.5, 1013.25, 288.15, mode='exact').value
        self.assertAlmostEqual(A_gas, 0.47081173472870474, places=5)


class TestITU676_13_SlantPathApprox(unittest.TestCase):
    """Slant path attenuation using Annex 2 simplified method.

    From A_Gas_A2_INST sheet (el=45 deg):
      - gamma0 and gammaw computed via exact method at surface conditions
      - h0 from Part 1 coefficient file (Eq. 31): h0 = a0+b0*T+c0*Ps+d0*rho
        where Ps = total surface pressure = p_dry + e
      - hw from Eq. 37 (Table 4 coefficients)
      - A_gas = (gamma0*h0 + gammaw*hw) / sin(theta)

    Input convention: P argument = dry air pressure p = Ps - e
    """

    def setUp(self):
        itu676.change_version(13)
        self.el = 45.0

    def _check(self, f, ps, T, rho, expected, places=5):
        """Helper: compute A_gas and check against expected value."""
        A_gas = itu676.gaseous_attenuation_slant_path(
            f, self.el, rho, ps, T, mode='approx').value
        self.assertAlmostEqual(A_gas, expected, places=places)

    # --- f = 38.5 GHz cases ---

    def test_case01_Ts295_rho14p0(self):
        self._check(38.5, 988.3342860812425, 295.15, 13.998103358274586,
                    0.6724061393008622)

    def test_case02_Ts294_rho14p0(self):
        self._check(38.5, 988.8194616390798, 294.45, 14.04229126442994,
                    0.6783244181497916)

    def test_case03_Ts295_rho14p2(self):
        self._check(38.5, 989.4840337179357, 294.65, 14.205904949340969,
                    0.6833829361095082)

    def test_case04_Ts297_rho14p3(self):
        self._check(38.5, 989.4976770016127, 297.15, 14.295215863202152,
                    0.6721249815284822)

    def test_case05_Ts301_rho13p2(self):
        self._check(38.5, 990.712469428683, 300.85, 13.172371197621427,
                    0.6171528309720985)

    def test_case06_Ts303_rho12p9(self):
        self._check(38.5, 990.6016244371225, 303.25, 12.932952957874935,
                    0.5981110294418397)

    def test_case07_Ts304_rho13p5(self):
        self._check(38.5, 989.3537790809332, 304.05, 13.503193794315951,
                    0.6104970109958014)

    def test_case08_Ts303_rho14p8(self):
        self._check(38.5, 987.9839758399005, 302.65, 14.83285126546697,
                    0.6579285461427736)

    def test_case09_Ts301_rho14p5(self):
        self._check(38.5, 989.3878692172694, 301.35, 14.53449059438436,
                    0.6566089491834671)

    # --- f = 39.5 GHz case ---

    def test_case10_f39p5_Ts298_rho15p9(self):
        self._check(39.5, 988.9020367914704, 297.65, 15.942511766465097,
                    0.7707708981960036)


if __name__ == '__main__':
    unittest.main()
