# -*- coding: utf-8 -*-
import unittest as test
import numpy as np

import itur.models as models
import sys
from astropy import units as u


def suite():
    """ A test suite for the ITU-P Recommendations. Recommendations tested:
    * ITU-P R-676-0
    * ITU-P R-618-12
    * ITU-P R-453-12
    * ITU-P R-837-6
    * ITU-P R-838-3
    * ITU-P R-839-4
    * ITU-P R-840-4
    """
    suite = test.TestSuite()

    # ITU-R P.676 tests (Gaseous attenuation)
    suite.addTest(ITUR676_9TestCase('test_gammaw'))
    suite.addTest(ITUR676_9TestCase('test_gamma0'))
    suite.addTest(ITUR676_9TestCase('test_zenit_water_vapour_attenuation'))

    # ITU-R P.618 tests (Rain attenuation)
    suite.addTest(ITUR618_12TestCase(
        'test_rain_cross_polarization_discrimination'))
    suite.addTest(ITUR618_12TestCase('test_rain_attenuation'))
    suite.addTest(ITUR618_12TestCase('test_scintillation_attenuation'))

    # ITU-R P.83X tests (Tabulated parameters)
    suite.addTest(ITUR453_12TestCase('test_wet_term_radio_refractivity'))
    suite.addTest(ITUR837_6TestCase('test_rainfall_rate'))
    suite.addTest(ITUR838_3TestCase('test_rain_specific_attenuation'))
    suite.addTest(ITUR839_4TestCase('test_rain_height'))

    # ITU-R P.840 tests (Clouds attenuation)
    suite.addTest(ITUR840_4TestCase('test_columnar_content_reduced_liquid'))
    suite.addTest(ITUR840_4TestCase('test_cloud_attenuation'))
    return suite


class ITUR453_12TestCase(test.TestCase):
    def setUp(self):
        models.itu453.change_version(12)

    def test_wet_term_radio_refractivity(self):
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(51.5, 359.86).value,
            45.130667, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(41.9, 12.49).value,
            53.756489, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(33.94, 18.43).value,
            76.349680, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(51.5, 359.86).value,
            45.130667, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(41.9, 12.49).value,
            53.756489, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(33.94, 18.43).value,
            76.349680, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(22.9, 316.77).value,
            87.907733, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(25.78, 279.78).value,
            101.416373, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(22.9, 316.77).value,
            87.907733, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(25.78, 279.78).value,
            101.416373, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(28.717, 77.3).value,
            60.060569, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(3.133, 101.7).value,
            105.920333, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(9.05, 38.7).value,
            50.162000, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(28.717, 77.3).value,
            60.060569, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(3.133, 101.7).value,
            105.920333, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(9.05, 38.7).value,
            50.162000, places=5)


class ITUR676_9TestCase(test.TestCase):

    def setUp(self):
        models.itu676.change_version(9)
        models.itu836.change_version(4)

    def test_gammaw(self):
        # The ITU models are non-sense and believe that the conversion between
        # Kelvin is 273 instead of 273.15
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(12, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.Celsius).value,
            0.00705700000, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(20, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.Celsius).value,
            0.06742720000, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(60, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.Celsius).value,
            0.11538020000, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(90, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.Celsius).value,
            0.25568340000, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(130, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.Celsius).value,
            0.56358380000, places=5)

    def test_gamma0(self):
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        282.724 - 0.15).value,
            0.00941327, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        287.4834667 - 0.15).value,
            0.00898682, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        293.1487022 - 0.15).value,
            0.00851359, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        282.724 - 0.15).value,
            0.00941327, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        287.4834667 - 0.15).value,
            0.00898682, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        293.1487022 - 0.15).value,
            0.00851359, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        282.724 - 0.15).value,
            0.00941327, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        287.4834667 - 0.15).value,
            0.00898682, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        293.1487022 - 0.15).value,
            0.00851359, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        282.724 - 0.15).value,
            0.00941327, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        287.4834667 - 0.15).value,
            0.00898682, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        293.1487022 - 0.15).value,
            0.00851359, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        282.724 - 0.15).value,
            0.02043748, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        287.4834667 - 0.15).value,
            0.01954568, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        293.1487022 - 0.15).value,
            0.01856193, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        282.724 - 0.15).value,
            0.02043748, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        287.4834667 - 0.15).value,
            0.01954568, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        293.1487022 - 0.15).value,
            0.01856193, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        282.724 - 0.15).value,
            0.02043748, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        287.4834667 - 0.15).value,
            0.01954568, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        293.1487022 - 0.15).value,
            0.01856193, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        282.724 - 0.15).value,
            0.02043748, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        287.4834667 - 0.15).value,
            0.01954568, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25,
                                        293.1487022 - 0.15).value,
            0.01856193, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25,
                                        296.602 - 0.15).value,
            0.00824203, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 296.7208533 - 0.15).value,
            0.0082329, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 296.602 - 0.15).value,
            0.00824203, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 296.7208533 - 0.15).value,
            0.0082329, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 296.602 - 0.15).value,
            0.00824203, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 296.7208533 - 0.15).value,
            0.0082329, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 296.602 - 0.15).value,
            0.00824203, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 296.7208533 - 0.15).value,
            0.0082329, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 296.602 - 0.15).value,
            0.01800011, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 296.7208533 - 0.15).value,
            0.01798125, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 296.602 - 0.15).value,
            0.01800011, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 296.7208533 - 0.15).value,
            0.01798125, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 296.602 - 0.15).value,
            0.01800011, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 296.7208533 - 0.15).value,
            0.01798125, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 296.602 - 0.15).value,
            0.01800011, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 296.7208533 - 0.15).value,
            0.01798125, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 299.0966578 - 0.15).value,
            0.00805331, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 297.9322267 - 0.15).value,
            0.00814064, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 287.444 - 0.15).value,
            0.00899025, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 299.0966578 - 0.15).value,
            0.00805331, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 297.9322267 - 0.15).value,
            0.00814064, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 287.444 - 0.15).value,
            0.00899025, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 299.0966578 - 0.15).value,
            0.00805331, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 297.9322267 - 0.15).value,
            0.00814064, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 287.444 - 0.15).value,
            0.00899025, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 299.0966578 - 0.15).value,
            0.00805331, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 297.9322267 - 0.15).value,
            0.00814064, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 287.444 - 0.15).value,
            0.00899025, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 299.0966578 - 0.15).value,
            0.01761077, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 297.9322267 - 0.15).value,
            0.01779083, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 287.444 - 0.15).value,
            0.01955282, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 299.0966578 - 0.15).value,
            0.01761077, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 297.9322267 - 0.15).value,
            0.01779083, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 287.444 - 0.15).value,
            0.01955282, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 299.0966578 - 0.15).value,
            0.01761077, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 297.9322267 - 0.15).value,
            0.01779083, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 287.444 - 0.15).value,
            0.01955282, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 299.0966578 - 0.15).value,
            0.01761077, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 297.9322267 - 0.15).value,
            0.01779083, places=6)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 287.444 - 0.15).value,
            0.01955282, places=6)

    def zenit_water_vapour_attenuation(self, lat, lon, el, p, f, alt):
        gamma_w = models.itu676.zenit_water_vapour_attenuation(lat,
                                                               lon,
                                                               p,
                                                               f,
                                                               None,
                                                               alt=alt).value
        return gamma_w / np.sin(np.deg2rad(el))

    def test_zenit_water_vapour_attenuation(self):
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                51.5, 359.86, 30.87067768, 1, 14.25, 0.06916422),
            0.12789267, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                41.9, 12.49, 40.97052773, 1, 14.25, 0.05670104),
            0.10865204, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                33.94, 18.43, 47.91280491, 1, 14.25, 0),
            0.10205633, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                51.5, 359.86, 30.87067768, 0.1, 14.25, 0.06916422),
            0.15315923, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                41.9, 12.49, 40.97052773, 0.1, 14.25, 0.05670104),
            0.12223686, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                33.94, 18.43, 47.91280491, 0.1, 14.25, 0),
            0.12410189, places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(51.5, 359.86, 30.87067768, 0.01, 14.25, 0.06916422),
#            0.15315923,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(41.9, 12.49, 40.97052773, 0.01, 14.25, 0.05670104),
#            0.12223686,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(33.94, 18.43, 47.91280491, 0.01, 14.25, 0),
#            0.12410189,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(51.5, 359.86, 30.87067768, 0.001, 14.25, 0.06916422),
#            0.15315923,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(41.9, 12.49, 40.97052773, 0.001, 14.25, 0.05670104),
#            0.12223686,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(33.94, 18.43, 47.91280491, 0.001, 14.25, 0),
#            0.12410189,  places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                51.5, 359.86, 30.87067768, 1, 29, 0.06916422),
            0.60896934, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                41.9, 12.49, 40.97052773, 1, 29, 0.05670104),
            0.51690529, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                33.94, 18.43, 47.91280491, 1, 29, 0),
            0.48519817, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                51.5, 359.86, 30.87067768, 0.1, 29, 0.06916422),
            0.72784676, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                41.9, 12.49, 40.97052773, 0.1, 29, 0.05670104),
            0.58076456, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                33.94, 18.43, 47.91280491, 0.1, 29, 0),
            0.58863533, places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(51.5, 359.86, 30.87067768, 0.01, 29, 0.06916422),
#            0.72784676,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(41.9, 12.49, 40.97052773, 0.01, 29, 0.05670104),
#            0.58076456,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(33.94, 18.43, 47.91280491, 0.01, 29, 0),
#            0.58863533,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(51.5, 359.86, 30.87067768, 0.001, 29, 0.06916422),
#            0.72784676,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(41.9, 12.49, 40.97052773, 0.001, 29, 0.05670104),
#            0.58076456,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(33.94, 18.43, 47.91280491, 0.001, 29, 0),
#            0.58863533,  places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                22.9, 316.77, 59.81487174, 1, 14.25, 0),
            0.1181882, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                25.78, 279.78, 49.20900369, 1, 14.25, 0.00007511),
            0.16093386, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                22.9, 316.77, 59.81487174, 0.1, 14.25, 0),
            0.13730617, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                25.78, 279.78, 49.20900369, 0.1, 14.25, 0.00007511),
            0.17798382, places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(22.9, 316.77, 59.81487174, 0.01, 14.25, 0),
#            0.13730617,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(25.78, 279.78, 49.20900369, 0.01, 14.25, 0.00007511),
#            0.17798382,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(22.9, 316.77, 59.81487174, 0.001, 14.25, 0),
#            0.13730617,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(25.78, 279.78, 49.20900369, 0.001, 14.25, 0.00007511),
#            0.17798382,  places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                22.9, 316.77, 59.81487174, 1, 29, 0),
            0.55983815, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                25.78, 279.78, 49.20900369, 1, 29, 0.00007511),
            0.76047761, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                22.9, 316.77, 59.81487174, 0.1, 29, 0),
            0.64906814, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                25.78, 279.78, 49.20900369, 0.1, 29, 0.00007511),
            0.83981774, places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(22.9, 316.77, 59.81487174, 0.01, 29, 0),
#            0.64906814,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(25.78, 279.78, 49.20900369, 0.01, 29, 0.00007511),
#            0.83981774,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(22.9, 316.77, 59.81487174, 0.001, 29, 0),
#            0.64906814,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(25.78, 279.78, 49.20900369, 0.001, 29, 0.00007511),
#            0.83981774,  places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                28.717, 77.3, 55.90591362, 1, 14.25, 0.21755946),
            0.18628614, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                3.133, 101.7, 67.76751981, 1, 14.25, 0.23610446),
            0.13468573, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                9.05, 38.7, 38.14104832, 1, 14.25, 2.45000492),
            0.08369587, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                28.717, 77.3, 55.90591362, 0.1, 14.25, 0.21755946),
            0.20242415, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                3.133, 101.7, 67.76751981, 0.1, 14.25, 0.23610446),
            0.14372476, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                9.05, 38.7, 38.14104832, 0.1, 14.25, 2.45000492),
            0.09153026, places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(28.717, 77.3, 55.90591362, 0.01, 14.25, 0.21755946),
#            0.20242415,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(3.133, 101.7, 67.76751981, 0.01, 14.25, 0.23610446),
#            0.14372476,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(9.05, 38.7, 38.14104832, 0.01, 14.25, 2.45000492),
#            0.09153026,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(28.717, 77.3, 55.90591362, 0.001, 14.25, 0.21755946),
#            0.20242415,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(3.133, 101.7, 67.76751981, 0.001, 14.25, 0.23610446),
#            0.14372476,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(9.05, 38.7, 38.14104832, 0.001, 14.25, 2.45000492),
#            0.09153026,  places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                28.717, 77.3, 55.90591362, 1, 29, 0.21755946),
            0.8771945, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                3.133, 101.7, 67.76751981, 1, 29, 0.23610446),
            0.63623574, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                9.05, 38.7, 38.14104832, 1, 29, 2.45000492),
            0.39942177, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                28.717, 77.3, 55.90591362, 0.1, 29, 0.21755946),
            0.95194476, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                3.133, 101.7, 67.76751981, 0.1, 29, 0.23610446),
            0.67829402, places=2)
        self.assertAlmostEqual(
            self.zenit_water_vapour_attenuation(
                9.05, 38.7, 38.14104832, 0.1, 29, 2.45000492),
            0.43646179, places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(28.717, 77.3, 55.90591362, 0.01, 29, 0.21755946),
#            0.95194476,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(3.133, 101.7, 67.76751981, 0.01, 29, 0.23610446),
#            0.67829402,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(9.05, 38.7, 38.14104832, 0.01, 29, 2.45000492),
#            0.43646179,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(28.717, 77.3, 55.90591362, 0.001, 29, 0.21755946),
#            0.95194476,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(3.133, 101.7, 67.76751981, 0.001, 29, 0.23610446),
#            0.67829402,  places=2)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(9.05, 38.7, 38.14104832, 0.001, 29, 2.45000492),
#            0.43646179,  places=2)


class ITUR838_3TestCase(test.TestCase):

    def setUp(self):
        models.itu838.change_version(3)

    def test_rain_specific_attenuation(self):
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 14.25, 30.87067768, 0).value,
            1.879742, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 14.25, 40.97052773, 0).value,
            3.630988, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 14.25, 47.91280491, 0).value,
            3.503189, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 14.25, 30.87067768, 0).value,
            1.879742, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 14.25, 40.97052773, 0).value,
            3.630988, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 14.25, 47.91280491, 0).value,
            3.503189, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 14.25, 30.87067768, 0).value,
            1.879742, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 14.25, 40.97052773, 0).value,
            3.630988, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 14.25, 47.91280491, 0).value,
            3.503189, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 14.25, 30.87067768, 0).value,
            1.879742, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 14.25, 40.97052773, 0).value,
            3.630988, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 14.25, 47.91280491, 0).value,
            3.503189, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 29.00, 30.87067768, 0).value,
            5.814832, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 29.00, 40.97052773, 0).value,
            10.157375, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 29.00, 47.91280491, 0).value,
            9.846762, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 29.00, 30.87067768, 0).value,
            5.814832, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 29.00, 40.97052773, 0).value,
            10.157375, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 29.00, 47.91280491, 0).value,
            9.846762, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 29.00, 30.87067768, 0).value,
            5.814832, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 29.00, 40.97052773, 0).value,
            10.157375, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 29.00, 47.91280491, 0).value,
            9.846762, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 29.00, 30.87067768, 0).value,
            5.814832, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 29.00, 40.97052773, 0).value,
            10.157375, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 29.00, 47.91280491, 0).value,
            9.846762, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 14.25, 59.81487174, 0).value,
            3.628282, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 14.25, 49.20900369, 0).value,
            5.948478, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 14.25, 59.81487174, 0).value,
            3.628282, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 14.25, 49.20900369, 0).value,
            5.948478, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 14.25, 59.81487174, 0).value,
            3.628282, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 14.25, 49.20900369, 0).value,
            5.948478, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 14.25, 59.81487174, 0).value,
            3.628282, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 14.25, 49.20900369, 0).value,
            5.948478, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 29.00, 59.81487174, 0).value,
            10.132682, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 29.00, 49.20900369, 0).value,
            15.460212, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 29.00, 59.81487174, 0).value,
            10.132682, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 29.00, 49.20900369, 0).value,
            15.460212, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 29.00, 59.81487174, 0).value,
            10.132682, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 29.00, 49.20900369, 0).value,
            15.460212, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 29.00, 59.81487174, 0).value,
            10.132682, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 29.00, 49.20900369, 0).value,
            15.460212, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 14.25, 55.90591362, 0).value,
            3.603569, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 14.25, 67.76751981, 0).value,
            6.06336, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 14.25, 38.14104832, 0).value,
            3.523996, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 14.25, 55.90591362, 0).value,
            3.603569, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 14.25, 67.76751981, 0).value,
            6.06336, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 14.25, 38.14104832, 0).value,
            3.523996, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 14.25, 55.90591362, 0).value,
            3.603569, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 14.25, 67.76751981, 0).value,
            6.06336, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 14.25, 38.14104832, 0).value,
            3.523996, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 14.25, 55.90591362, 0).value,
            3.603569, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 14.25, 67.76751981, 0).value,
            6.06336, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 14.25, 38.14104832, 0).value,
            3.523996, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 29.00, 55.90591362, 0).value,
            10.078266, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 29.00, 67.76751981, 0).value,
            15.712442, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 29.00, 38.14104832, 0).value,
            9.904098, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 29.00, 55.90591362, 0).value,
            10.078266, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 29.00, 67.76751981, 0).value,
            15.712442, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 29.00, 38.14104832, 0).value,
            9.904098, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 29.00, 55.90591362, 0).value,
            10.078266, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 29.00, 67.76751981, 0).value,
            15.712442, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 29.00, 38.14104832, 0).value,
            9.904098, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 29.00, 55.90591362, 0).value,
            10.078266, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 29.00, 67.76751981, 0).value,
            15.712442, places=6)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 29.00, 38.14104832, 0).value,
            9.904098, places=6)


class ITUR837_6TestCase(test.TestCase):

    def setUp(self):
        models.itu837.change_version(6)

    def test_rainfall_rate(self):
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(51.500, 359.86, 0.01).value,
            30.8750240, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(41.900, 12.49, 0.01).value,
            56.3700090, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(33.940, 18.43, 0.01).value,
            55.2316250, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(51.500, 359.86, 0.01).value,
            30.8750240, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(41.900, 12.49, 0.01).value,
            56.3700090, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(33.940, 18.43, 0.01).value,
            55.2316250, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(22.900, 316.77, 0.01).value,
            58.0942160, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(25.780, 279.78, 0.01).value,
            89.1141030, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(22.900, 316.77, 0.01).value,
            58.0942160, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(25.780, 279.78, 0.01).value,
            89.1141030, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(28.717, 77.30, 0.01).value,
            57.3962300, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(3.133, 101.70, 0.01).value,
            93.6070980, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(9.050, 38.70, 0.01).value,
            54.6234110, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(28.717, 77.30, 0.01).value,
            57.3962300, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(3.133, 101.70, 0.01).value,
            93.6070980, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(9.050, 38.70, 0.01).value,
            54.6234110, places=5)


class ITUR839_4TestCase(test.TestCase):

    def setUp(self):
        models.itu839.change_version(4)

    def test_rain_height(self):
        self.assertAlmostEqual(
            models.itu839.rain_height(51.500, 359.86).value,
            2.4527330, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(41.900, 12.49).value,
            3.0474930, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(33.940, 18.43).value,
            2.5633030, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(51.500, 359.86).value,
            2.4527330, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(41.900, 12.49).value,
            3.0474930, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(33.940, 18.43).value,
            2.5633030, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(22.900, 316.77).value,
            4.1587790, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(25.780, 279.78).value,
            4.5694610, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(22.900, 316.77).value,
            4.1587790, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(25.780, 279.78).value,
            4.5694610, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(28.717, 77.30).value,
            5.2582040, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(3.133, 101.70).value,
            4.9579740, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(9.050, 38.70).value,
            4.7839070, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(28.717, 77.30).value,
            5.2582040, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(3.133, 101.70).value,
            4.9579740, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(9.050, 38.70).value,
            4.7839070, places=5)


class ITUR618_12TestCase(test.TestCase):

    def setUp(self):
        models.itu618.change_version(12)
        models.itu838.change_version(3)

    def test_rain_cross_polarization_discrimination(self):
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(16.38308757, 14.25,
                                                                 30.870677680, 0.001, 0).value,
            27.143007980, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(3.89479806, 14.25,
                                                                 40.970527730, 0.1, 0).value,
            37.386086000, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(9.71179484, 14.25,
                                                                 47.912804910, 0.01, 0).value,
            33.812795580, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(71.44613350, 29,
                                                                 40.970527730, 0.001, 0).value,
            21.244470560, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(12.87478397, 29,
                                                                 47.912804910, 0.1, 0).value,
            35.166125690, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(39.07323323, 29,
                                                                 40.970527730, 0.01, 0).value,
            25.180145740, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(23.00197384, 14.25,
                                                                 59.814871740, 0.001, 0).value,
            33.308530550, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(32.74150676, 14.25,
                                                                 49.209003690, 0.001, 0).value,
            25.508227320, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(4.92489694, 14.25,
                                                                 59.814871740, 0.1, 0).value,
            41.798127850, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(6.96606559, 14.25,
                                                                 49.209003690, 0.1, 0).value,
            34.830206060, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(12.76053997, 14.25,
                                                                 59.814871740, 0.01, 0).value,
            36.168649690, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(18.06938866, 14.25,
                                                                 49.209003690, 0.01, 0).value,
            28.803871260, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(23.00197384, 14.25,
                                                                 59.814871740, 0.001, 0).value,
            33.308530550, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(32.74150676, 14.25,
                                                                 49.209003690, 0.001, 0).value,
            25.508227320, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(100.96022257, 29,
                                                                 49.209003690, 0.001, 0).value,
            20.365001500, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(20.43214239, 29,
                                                                 59.814871740, 0.1, 0).value,
            35.581135690, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(27.86774318, 29,
                                                                 49.209003690, 0.1, 0).value,
            28.745547830, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(46.32024457, 29,
                                                                 59.814871740, 0.01, 0).value,
            30.303830010, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(63.46384760, 29,
                                                                 49.209003690, 0.01, 0).value,
            23.046241580, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(73.05533363, 29,
                                                                 59.814871740, 0.001, 0).value,
            28.089155910, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(26.85402570, 14.25,
                                                                 55.905913620, 0.001, 0).value,
            29.993601830, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(4.44533923, 14.25,
                                                                 38.141048320, 0.1, 0).value,
            35.652315760, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(11.06265445, 14.25,
                                                                 38.141048320, 0.01, 0).value,
            30.034285750, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(26.85402570, 14.25,
                                                                 55.905913620, 0.001, 0).value,
            29.993601830, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(21.84602116, 29,
                                                                 55.905913620, 0.1, 0).value,
            33.289964560, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(51.72271818, 29,
                                                                 55.905913620, 0.01, 0).value,
            27.480618010, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(53.61322867, 29,
                                                                 38.141048320, 0.001, 0).value,
            23.354012700, places=5)

    def test_rain_attenuation(self):
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 14.25, 30.87067768,
                                           hs=0.0691640, p=0.01, tau=0.00).value,
            7.5572640, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 14.25, 40.97052773,
                                           hs=0.0567010, p=0.01, tau=0.00).value,
            11.4735460, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 14.25, 47.91280491,
                                           hs=0.0000000, p=0.01, tau=0.00).value,
            9.7117950, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 29.00, 30.87067768,
                                           hs=0.0691640, p=0.01, tau=0.00).value,
            25.7166770, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 29.00, 40.97052773,
                                           hs=0.0567010, p=0.01, tau=0.00).value,
            39.0732330, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 29.00, 47.91280491,
                                           hs=0.0000000, p=0.01, tau=0.00).value,
            33.4169840, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 14.25, 59.81487174,
                                           hs=0.0000000, p=0.01, tau=0.00).value,
            12.7605400, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 14.25, 49.20900369,
                                           hs=0.0000750, p=0.01, tau=0.00).value,
            18.0693890, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 29.00, 59.81487174,
                                           hs=0.0000000, p=0.01, tau=0.00).value,
            46.3202450, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 29.00, 49.20900369,
                                           hs=0.0000750, p=0.01, tau=0.00).value,
            63.4638480, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 14.25, 55.90591362,
                                           hs=0.2175590, p=0.01, tau=0.00).value,
            14.1707990, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 14.25, 67.76751981,
                                           hs=0.2361040, p=0.01, tau=0.00).value,
            19.6617050, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 14.25, 38.14104832,
                                           hs=2.4500050, p=0.01, tau=0.00).value,
            11.0626540, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 29.00, 55.90591362,
                                           hs=0.2175590, p=0.01, tau=0.00).value,
            51.7227180, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 29.00, 67.76751981,
                                           hs=0.2361040, p=0.01, tau=0.00).value,
            70.5396050, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 29.00, 38.14104832,
                                           hs=2.4500050, p=0.01, tau=0.00).value,
            35.1160650, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 14.25, 30.87067768,
                                           hs=0.0691640, p=0.10, tau=0.00).value,
            2.4567600, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 14.25, 40.97052773,
                                           hs=0.0567010, p=0.10, tau=0.00).value,
            3.8947980, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 14.25, 47.91280491,
                                           hs=0.0000000, p=0.10, tau=0.00).value,
            3.2920370, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 29.00, 30.87067768,
                                           hs=0.0691640, p=0.10, tau=0.00).value,
            9.4912070, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 29.00, 40.97052773,
                                           hs=0.0567010, p=0.10, tau=0.00).value,
            15.0594580, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 29.00, 47.91280491,
                                           hs=0.0000000, p=0.10, tau=0.00).value,
            12.8747840, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 14.25, 59.81487174,
                                           hs=0.0000000, p=0.10, tau=0.00).value,
            4.9248970, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 14.25, 49.20900369,
                                           hs=0.0000750, p=0.10, tau=0.00).value,
            6.9660660, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 29.00, 59.81487174,
                                           hs=0.0000000, p=0.10, tau=0.00).value,
            20.4321420, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 29.00, 49.20900369,
                                           hs=0.0000750, p=0.10, tau=0.00).value,
            27.8677430, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 14.25, 55.90591362,
                                           hs=0.2175590, p=0.10, tau=0.00).value,
            5.2338740, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 14.25, 67.76751981,
                                           hs=0.2361040, p=0.10, tau=0.00).value,
            9.6728110, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 14.25, 38.14104832,
                                           hs=2.4500050, p=0.10, tau=0.00).value,
            4.4453390, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 29.00, 55.90591362,
                                           hs=0.2175590, p=0.10, tau=0.00).value,
            21.8460210, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 29.00, 67.76751981,
                                           hs=0.2361040, p=0.10, tau=0.00).value,
            39.6143120, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 29.00, 38.14104832,
                                           hs=2.4500050, p=0.10, tau=0.00).value,
            15.9048720, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 14.25, 30.87067768,
                                           hs=0.0691640, p=1.00, tau=0.00).value,
            0.5628470, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 14.25, 40.97052773,
                                           hs=0.0567010, p=1.00, tau=0.00).value,
            0.9317550, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 14.25, 47.91280491,
                                           hs=0.0000000, p=1.00, tau=0.00).value,
            0.7619040, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 29.00, 30.87067768,
                                           hs=0.0691640, p=1.00, tau=0.00).value,
            2.4686380, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 29.00, 40.97052773,
                                           hs=0.0567010, p=1.00, tau=0.00).value,
            4.0904280, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 29.00, 47.91280491,
                                           hs=0.0000000, p=1.00, tau=0.00).value,
            3.3867490, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 14.25, 59.81487174,
                                           hs=0.0000000, p=1.00, tau=0.00).value,
            1.0593540, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 14.25, 49.20900369,
                                           hs=0.0000750, p=1.00, tau=0.00).value,
            1.6122160, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 29.00, 59.81487174,
                                           hs=0.0000000, p=1.00, tau=0.00).value,
            5.0231130, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 29.00, 49.20900369,
                                           hs=0.0000750, p=1.00, tau=0.00).value,
            7.3463010, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 14.25, 55.90591362,
                                           hs=0.2175590, p=1.00, tau=0.00).value,
            1.2022670, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 14.25, 67.76751981,
                                           hs=0.2361040, p=1.00, tau=0.00).value,
            1.7852610, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 14.25, 38.14104832,
                                           hs=2.4500050, p=1.00, tau=0.00).value,
            0.8916230, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 29.00, 55.90591362,
                                           hs=0.2175590, p=1.00, tau=0.00).value,
            5.7386810, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 29.00, 67.76751981,
                                           hs=0.2361040, p=1.00, tau=0.00).value,
            8.3461990, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 29.00, 38.14104832,
                                           hs=2.4500050, p=1.00, tau=0.00).value,
            3.5957140, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 14.25, 30.87067768,
                                           hs=0.0691640, p=0.01, tau=0.00).value,
            7.5572640, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 14.25, 40.97052773,
                                           hs=0.0567010, p=0.01, tau=0.00).value,
            11.4735460, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 14.25, 47.91280491,
                                           hs=0.0000000, p=0.01, tau=0.00).value,
            9.7117950, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 29.00, 30.87067768,
                                           hs=0.0691640, p=0.01, tau=0.00).value,
            25.7166770, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 29.00, 40.97052773,
                                           hs=0.0567010, p=0.01, tau=0.00).value,
            39.0732330, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 29.00, 47.91280491,
                                           hs=0.0000000, p=0.01, tau=0.00).value,
            33.4169840, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 14.25, 59.81487174,
                                           hs=0.0000000, p=0.01, tau=0.00).value,
            12.7605400, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 14.25, 49.20900369,
                                           hs=0.0000750, p=0.01, tau=0.00).value,
            18.0693890, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 29.00, 59.81487174,
                                           hs=0.0000000, p=0.01, tau=0.00).value,
            46.3202450, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 29.00, 49.20900369,
                                           hs=0.0000750, p=0.01, tau=0.00).value,
            63.4638480, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 14.25, 55.90591362,
                                           hs=0.2175590, p=0.01, tau=0.00).value,
            14.1707990, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 14.25, 67.76751981,
                                           hs=0.2361040, p=0.01, tau=0.00).value,
            19.6617050, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 14.25, 38.14104832,
                                           hs=2.4500050, p=0.01, tau=0.00).value,
            11.0626540, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 29.00, 55.90591362,
                                           hs=0.2175590, p=0.01, tau=0.00).value,
            51.7227180, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 29.00, 67.76751981,
                                           hs=0.2361040, p=0.01, tau=0.00).value,
            70.5396050, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 29.00, 38.14104832,
                                           hs=2.4500050, p=0.01, tau=0.00).value,
            35.1160650, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 14.25, 30.87067768,
                                           hs=0.0691640, p=0.10, tau=0.00).value,
            2.4567600, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 14.25, 40.97052773,
                                           hs=0.0567010, p=0.10, tau=0.00).value,
            3.8947980, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 14.25, 47.91280491,
                                           hs=0.0000000, p=0.10, tau=0.00).value,
            3.2920370, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 29.00, 30.87067768,
                                           hs=0.0691640, p=0.10, tau=0.00).value,
            9.4912070, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 29.00, 40.97052773,
                                           hs=0.0567010, p=0.10, tau=0.00).value,
            15.0594580, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 29.00, 47.91280491,
                                           hs=0.0000000, p=0.10, tau=0.00).value,
            12.8747840, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 14.25, 59.81487174,
                                           hs=0.0000000, p=0.10, tau=0.00).value,
            4.9248970, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 14.25, 49.20900369,
                                           hs=0.0000750, p=0.10, tau=0.00).value,
            6.9660660, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 29.00, 59.81487174,
                                           hs=0.0000000, p=0.10, tau=0.00).value,
            20.4321420, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 29.00, 49.20900369,
                                           hs=0.0000750, p=0.10, tau=0.00).value,
            27.8677430, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 14.25, 55.90591362,
                                           hs=0.2175590, p=0.10, tau=0.00).value,
            5.2338740, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 14.25, 67.76751981,
                                           hs=0.2361040, p=0.10, tau=0.00).value,
            9.6728110, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 14.25, 38.14104832,
                                           hs=2.4500050, p=0.10, tau=0.00).value,
            4.4453390, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 29.00, 55.90591362,
                                           hs=0.2175590, p=0.10, tau=0.00).value,
            21.8460210, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 29.00, 67.76751981,
                                           hs=0.2361040, p=0.10, tau=0.00).value,
            39.6143120, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 29.00, 38.14104832,
                                           hs=2.4500050, p=0.10, tau=0.00).value,
            15.9048720, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 14.25, 30.87067768,
                                           hs=0.0691640, p=1.00, tau=0.00).value,
            0.5628470, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 14.25, 40.97052773,
                                           hs=0.0567010, p=1.00, tau=0.00).value,
            0.9317550, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 14.25, 47.91280491,
                                           hs=0.0000000, p=1.00, tau=0.00).value,
            0.7619040, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(51.500, 359.86, 29.00, 30.87067768,
                                           hs=0.0691640, p=1.00, tau=0.00).value,
            2.4686380, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(41.900, 12.49, 29.00, 40.97052773,
                                           hs=0.0567010, p=1.00, tau=0.00).value,
            4.0904280, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(33.940, 18.43, 29.00, 47.91280491,
                                           hs=0.0000000, p=1.00, tau=0.00).value,
            3.3867490, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 14.25, 59.81487174,
                                           hs=0.0000000, p=1.00, tau=0.00).value,
            1.0593540, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 14.25, 49.20900369,
                                           hs=0.0000750, p=1.00, tau=0.00).value,
            1.6122160, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(22.900, 316.77, 29.00, 59.81487174,
                                           hs=0.0000000, p=1.00, tau=0.00).value,
            5.0231130, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(25.780, 279.78, 29.00, 49.20900369,
                                           hs=0.0000750, p=1.00, tau=0.00).value,
            7.3463010, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 14.25, 55.90591362,
                                           hs=0.2175590, p=1.00, tau=0.00).value,
            1.2022670, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 14.25, 67.76751981,
                                           hs=0.2361040, p=1.00, tau=0.00).value,
            1.7852610, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 14.25, 38.14104832,
                                           hs=2.4500050, p=1.00, tau=0.00).value,
            0.8916230, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(28.717, 77.30, 29.00, 55.90591362,
                                           hs=0.2175590, p=1.00, tau=0.00).value,
            5.7386810, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(3.133, 101.70, 29.00, 67.76751981,
                                           hs=0.2361040, p=1.00, tau=0.00).value,
            8.3461990, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(9.050, 38.70, 29.00, 38.14104832,
                                           hs=2.4500050, p=1.00, tau=0.00).value,
            3.5957140, places=5)

    def test_scintillation_attenuation(self):
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, 359.86, 14.25, 30.87067768, 0.001, 0.9, eta=0.6).value,
            0.866044, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 14.25, 40.97052773, 0.001, 0.9, eta=0.6).value,
            0.710527, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 14.25, 47.91280491, 0.001, 0.9, eta=0.6).value,
            0.764448, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, 359.86, 29, 30.87067768, 0.001, 0.9, eta=0.6).value,
            1.289482, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 29, 40.97052773, 0.001, 0.9, eta=0.6).value,
            1.054611, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 29, 47.91280491, 0.001, 0.9, eta=0.6).value,
            1.132606, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, 316.77, 14.25, 59.81487174, 0.001, 0.9, eta=0.6).value,
            0.699472, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, 279.78, 14.25, 49.20900369, 0.001, 0.9, eta=0.6).value,
            0.912438, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, 316.77, 29, 59.81487174, 0.001, 0.9, eta=0.6).value,
            1.033819, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, 279.78, 29, 49.20900369, 0.001, 0.9, eta=0.6).value,
            1.351457, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 14.25, 55.90591362, 0.001, 0.9, eta=0.6).value,
            0.571530, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 14.25, 67.76751981, 0.001, 0.9, eta=0.6).value,
            0.736636, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 14.25, 38.14104832, 0.001, 0.9, eta=0.6).value,
            0.733740, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 29, 55.90591362, 0.001, 0.9, eta=0.6).value,
            0.845322, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 29, 67.76751981, 0.001, 0.9, eta=0.6).value,
            1.087468, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 29, 38.14104832, 0.001, 0.9, eta=0.6).value,
            1.089954, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, 359.86, 14.25, 30.87067768, 0.1, 0.9, eta=0.6).value,
            0.402326, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 14.25, 40.97052773, 0.1, 0.9, eta=0.6).value,
            0.330080, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 14.25, 47.91280491, 0.1, 0.9, eta=0.6).value,
            0.355129, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, 359.86, 29, 30.87067768, 0.1, 0.9, eta=0.6).value,
            0.599037, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 29, 40.97052773, 0.1, 0.9, eta=0.6).value,
            0.489926, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 29, 47.91280491, 0.1, 0.9, eta=0.6).value,
            0.526159, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, 316.77, 14.25, 59.81487174, 0.1, 0.9, eta=0.6).value,
            0.324944, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, 279.78, 14.25, 49.20900369, 0.1, 0.9, eta=0.6).value,
            0.423879, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, 316.77, 29, 59.81487174, 0.1, 0.9, eta=0.6).value,
            0.480267, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, 279.78, 29, 49.20900369, 0.1, 0.9, eta=0.6).value,
            0.627828, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 14.25, 55.90591362, 0.1, 0.9, eta=0.6).value,
            0.265508, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 14.25, 67.76751981, 0.1, 0.9, eta=0.6).value,
            0.342209, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 14.25, 38.14104832, 0.1, 0.9, eta=0.6).value,
            0.340864, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 29, 55.90591362, 0.1, 0.9, eta=0.6).value,
            0.392700, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 29, 67.76751981, 0.1, 0.9, eta=0.6).value,
            0.505190, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 29, 38.14104832, 0.1, 0.9, eta=0.6).value,
            0.506345, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, 359.86, 14.25, 30.87067768, 1, 0.9, eta=0.6).value,
            0.249221, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 14.25, 40.97052773, 1, 0.9, eta=0.6).value,
            0.204468, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 14.25, 47.91280491, 1, 0.9, eta=0.6).value,
            0.219985, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, 359.86, 29, 30.87067768, 1, 0.9, eta=0.6).value,
            0.371074, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 29, 40.97052773, 1, 0.9, eta=0.6).value,
            0.303485, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 29, 47.91280491, 1, 0.9, eta=0.6).value,
            0.325930, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, 316.77, 14.25, 59.81487174, 1, 0.9, eta=0.6).value,
            0.201287, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, 279.78, 14.25, 49.20900369, 1, 0.9, eta=0.6).value,
            0.262572, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, 316.77, 29, 59.81487174, 1, 0.9, eta=0.6).value,
            0.297502, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, 279.78, 29, 49.20900369, 1, 0.9, eta=0.6).value,
            0.388909, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 14.25, 55.90591362, 1, 0.9, eta=0.6).value,
            0.164469, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 14.25, 67.76751981, 1, 0.9, eta=0.6).value,
            0.211982, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 14.25, 38.14104832, 1, 0.9, eta=0.6).value,
            0.211148, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 29, 55.90591362, 1, 0.9, eta=0.6).value,
            0.243258, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 29, 67.76751981, 1, 0.9, eta=0.6).value,
            0.312940, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 29, 38.14104832, 1, 0.9, eta=0.6).value,
            0.313656, places=5)


class ITUR840_4TestCase(test.TestCase):

    def setUp(self):
        models.itu840.change_version(4)

    def test_columnar_content_reduced_liquid(self):
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                51.5, -0.14, 1.000).value,
            1.26328612, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                41.9, 12.49, 1.000).value,
            0.91467189, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                33.94, 18.43, 1.000).value,
            0.73072098, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                51.5, -0.14, 0.100).value,
            1.90329847, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                41.9, 12.49, 0.100).value,
            1.49845951, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(33.94, 18.43, 0.100).value,
#       1.47628568, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(51.5, -0.14, 0.010).value,
#       1.90329847, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(41.9, 12.49, 0.010).value,
#       1.49845951, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(33.94, 18.43, 0.010).value,
#       1.47628568, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(51.5, -0.14, 0.001).value,
#       1.90329847, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(41.9, 12.49, 0.001).value,
#       1.49845951, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(33.94, 18.43, 0.001).value,
#       1.47628568, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                51.5, -0.14, 1.000).value,
            1.26328612, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                41.9, 12.49, 1.000).value,
            0.91467189, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                33.94, 18.43, 1.000).value,
            0.73072098, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                51.5, -0.14, 0.100).value,
            1.90329847, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                41.9, 12.49, 0.100).value,
            1.49845951, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                33.94, 18.43, 0.100).value,
            1.47628568, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(51.5, -0.14, 0.010).value,
#       1.90329847, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(41.9, 12.49, 0.010).value,
#       1.49845951, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(33.94, 18.43, 0.010).value,
#       1.47628568, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(51.5, -0.14, 0.001).value,
#       1.90329847, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(41.9, 12.49, 0.001).value,
#       1.49845951, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(33.94, 18.43, 0.001).value,
#       1.47628568, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                22.9, -43.23, 1.000).value,
            1.10444871, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                25.78, -80.22, 1.000).value,
            2.27978216, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                22.9, -43.23, 0.100).value,
            2.82993169, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                25.78, -80.22, 0.100).value,
            3.52927516, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(22.9, -43.23, 0.010).value,
#       2.82993169, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(25.78, -80.22, 0.010).value,
#       3.52927516, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(22.9, -43.23, 0.001).value,
#       2.82993169, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(25.78, -80.22, 0.001).value,
#       3.52927516, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                22.9, -43.23, 1.000).value,
            1.10444871, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                25.78, -80.22, 1.000).value,
            2.27978216, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                22.9, -43.23, 0.100).value,
            2.82993169, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                25.78, -80.22, 0.100).value,
            3.52927516, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(22.9, -43.23, 0.010).value,
#       2.82993169, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(25.78, -80.22, 0.010).value,
#       3.52927516, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(22.9, -43.23, 0.001).value,
#       2.82993169, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(25.78, -80.22, 0.001).value,
#       3.52927516, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                28.717, 77.3, 1.000).value,
            2.75109958, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                3.133, 101.7, 1.000).value,
            3.33600769, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                9.05, 38.7, 1.000).value,
            1.21770185, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                28.717, 77.3, 0.100).value,
            4.23072604, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                3.133, 101.7, 0.100).value,
            3.80525123, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                9.05, 38.7, 0.100).value,
            1.49251459, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(28.717, 77.3, 0.010).value,
#       4.23072604, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(3.133, 101.7, 0.010).value,
#       3.80525123, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(9.05, 38.7, 0.010).value,
#       1.49251459, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(28.717, 77.3, 0.001).value,
#       4.23072604, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(3.133, 101.7, 0.001).value,
#       3.80525123, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(9.05, 38.7, 0.001).value,
#       1.49251459, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                28.717, 77.3, 1.000).value,
            2.75109958, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                3.133, 101.7, 1.000).value,
            3.33600769, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                9.05, 38.7, 1.000).value,
            1.21770185, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                28.717, 77.3, 0.100).value,
            4.23072604, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                3.133, 101.7, 0.100).value,
            3.80525123, places=1)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                9.05, 38.7, 0.100).value,
            1.49251459, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(28.717, 77.3, 0.010).value,
#       4.23072604, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(3.133, 101.7, 0.010).value,
#       3.80525123, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(9.05, 38.7, 0.010).value,
#       1.49251459, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(28.717, 77.3, 0.001).value,
#       4.23072604, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(3.133, 101.7, 0.001).value,
#       3.80525123, places=1)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(9.05, 38.7, 0.001).value,
#       1.49251459, places=1)

    def test_cloud_attenuation(self):
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 30.87067768, 14.25, 1.000).value,
            0.45792895, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.97052773, 14.25, 1.000).value,
            0.25946553, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 47.91280491, 14.25, 1.000).value,
            0.18313623, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 30.87067768, 14.25, 0.100).value,
            0.68992722, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.97052773, 14.25, 0.100).value,
            0.42506892, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 47.91280491, 14.25, 0.100).value,
            0.36999265, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(51.5,-0.14,30.87067768,14.25, 0.010).value,
#            0.68992722, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(41.9,12.49,40.97052773,14.25, 0.010).value,
#            0.42506892, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(33.94,18.43,47.91280491,14.25, 0.010).value,
#            0.36999265, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(51.5,-0.14,30.87067768,14.25, 0.001).value,
#            0.68992722, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(41.9,12.49,40.97052773,14.25, 0.001).value,
#            0.42506892, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(33.94,18.43,47.91280491,14.25, 0.001).value,
#            0.36999265, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 30.87067768, 29, 1.000).value,
            1.79599547, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.97052773, 29, 1.000).value,
            1.01762274, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 47.91280491, 29, 1.000).value,
            0.71825953, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 30.87067768, 29, 0.100).value,
            2.70589171, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.97052773, 29, 0.100).value,
            1.66711854, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 47.91280491, 29, 0.100).value,
            1.45110964, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(51.5,-0.14,30.87067768,29, 0.010).value,
#            2.70589171, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(41.9,12.49,40.97052773,29, 0.010).value,
#            1.66711854, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(33.94,18.43,47.91280491,29, 0.010).value,
#            1.45110964, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(51.5,-0.14,30.87067768,29, 0.001).value,
#            2.70589171, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(41.9,12.49,40.97052773,29, 0.001).value,
#            1.66711854, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(33.94,18.43,47.91280491,29, 0.001).value,
#            1.45110964, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 59.81487174, 14.25, 1.000).value,
            0.23764476, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 49.20900369, 14.25, 1.000).value,
            0.56006901, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 59.81487174, 14.25, 0.100).value,
            0.60891776, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 49.20900369, 14.25, 0.100).value,
            0.86702917, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(22.9,-43.23,59.81487174,14.25, 0.010).value,
#            0.60891776, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(25.78,-80.22,49.20900369,14.25, 0.010).value,
#            0.86702917, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(22.9,-43.23,59.81487174,14.25, 0.001).value,
#            0.60891776, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(25.78,-80.22,49.20900369,14.25, 0.001).value,
#            0.86702917, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 59.81487174, 29, 1.000).value,
            0.93204177, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 49.20900369, 29, 1.000).value,
            2.19658834, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 59.81487174, 29, 0.100).value,
            2.38817297, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 49.20900369, 29, 0.100).value,
            3.40048483, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(22.9,-43.23,59.81487174,29, 0.010).value,
#            2.38817297, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(25.78,-80.22,49.20900369,29, 0.010).value,
#            3.40048483, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(22.9,-43.23,59.81487174,29, 0.001).value,
#            2.38817297, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(25.78,-80.22,49.20900369,29, 0.001).value,
#            3.40048483, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.717, 77.3, 55.90591362, 14.25, 1.000).value,
            0.6178942, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.133, 101.7, 67.76751981, 14.25, 1.000).value,
            0.67031269, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 38.14104832, 14.25, 1.000).value,
            0.36671963, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.717, 77.3, 55.90591362, 14.25, 0.100).value,
            0.95021681, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.133, 101.7, 67.76751981, 14.25, 0.100).value,
            0.76459901, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 38.14104832, 14.25, 0.100).value,
            0.44948146, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(28.717,77.3,55.90591362,14.25, 0.010).value,
#            0.95021681, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(3.133,101.7,67.76751981,14.25, 0.010).value,
#            0.76459901, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(9.05,38.7,38.14104832,14.25, 0.010).value,
#            0.44948146, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(28.717,77.3,55.90591362,14.25, 0.001).value,
#            0.95021681, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(3.133,101.7,67.76751981,14.25, 0.001).value,
#            0.76459901, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(9.05,38.7,38.14104832,14.25, 0.001).value,
#            0.44948146, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.717, 77.3, 55.90591362, 29, 1.000).value,
            2.4233785, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.133, 101.7, 67.76751981, 29, 1.000).value,
            2.6289636, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 38.14104832, 29, 1.000).value,
            1.43827289, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.717, 77.3, 55.90591362, 29, 0.100).value,
            3.72674641, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.133, 101.7, 67.76751981, 29, 0.100).value,
            2.99875418, places=1)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 38.14104832, 29, 0.100).value,
            1.76286444, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(28.717,77.3,55.90591362,29, 0.010).value,
#            3.72674641, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(3.133,101.7,67.76751981,29, 0.010).value,
#            2.99875418, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(9.05,38.7,38.14104832,29, 0.010).value,
#            1.76286444, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(28.717,77.3,55.90591362,29, 0.001).value,
#            3.72674641, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(3.133,101.7,67.76751981,29, 0.001).value,
#            2.99875418, places=1)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(9.05,38.7,38.14104832,29, 0.001).value,
#            1.76286444, places=1)


if __name__ == '__main__':
    pass
    suite = suite()
    print('Validation tests for the ITU-R models')
    print('------------------------')
    print(
        'A total of %d test-cases are going to be tested' %
        suite.countTestCases())
    sys.stdout.flush()
    test.TextTestRunner(verbosity=2).run(suite)
