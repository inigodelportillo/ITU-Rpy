# -*- coding: utf-8 -*-
import unittest as test

import itur
import itur.models as models

import sys
from astropy import units as u


def suite():
    """ A test suite for the ITU-P Recommendations. Recommendations tested:
    * ITU-P R-676-9
    * ITU-P R-676-11
    * ITU-P R-618-12
    * ITU-P R-618-13
    * ITU-P R-453-12
    * ITU-P R-837-6
    * ITU-P R-837-7
    * ITU-P R-838-3
    * ITU-P R-839-4
    * ITU-P R-840-4
    * ITU-P R-840-7
    * ITU-P R-1511-1
    """
    suite = test.TestSuite()

    # Ensure models are in the right version
    models.itu453.change_version(13)
    models.itu618.change_version(13)
    models.itu676.change_version(11)
    models.itu836.change_version(6)
    models.itu837.change_version(7)
    models.itu838.change_version(3)
    models.itu839.change_version(4)
    models.itu840.change_version(7)
    models.itu1510.change_version(1)
    models.itu1511.change_version(1)

    # ITU-R P.676 tests (Gaseous attenuation)
    suite.addTest(ITUR676_9TestCase('test_gammaw'))
    suite.addTest(ITUR676_9TestCase('test_gamma0'))
#    suite.addTest(ITUR676_9TestCase('test_zenit_water_vapour_attenuation'))
    suite.addTest(ITUR676_11TestCase('test_gammaw_exact'))
    suite.addTest(ITUR676_11TestCase('test_gamma0_exact'))
    suite.addTest(ITUR676_11TestCase('test_gammaw_approx'))
    suite.addTest(ITUR676_11TestCase('test_gamma0_approx'))
    suite.addTest(ITUR676_11TestCase('test_zenit_water_vapour_attenuation'))

    # ITU-R P.618 tests (Rain attenuation)
    suite.addTest(ITUR618_12TestCase(
        'test_rain_cross_polarization_discrimination'))
    suite.addTest(ITUR618_12TestCase('test_rain_attenuation'))
    suite.addTest(ITUR618_12TestCase('test_scintillation_attenuation'))

    suite.addTest(ITUR618_13TestCase('test_rain_attenuation'))
    suite.addTest(ITUR618_13TestCase('test_probability_of_rain_attenuation'))
#    suite.addTest(ITUR618_13TestCase('test_site_diversity'))
    suite.addTest(ITUR618_13TestCase('test_scintillation_attenuation'))
    suite.addTest(ITUR618_13TestCase(
        'test_rain_cross_polarization_discrimination'))
    suite.addTest(ITUR618_13TestCase('test_total_attenuation'))

    # ITU-R P.453 tests (Wet term radio refractivity)
    suite.addTest(ITUR453_12TestCase('test_wet_term_radio_refractivity'))
    suite.addTest(ITUR453_13TestCase('test_wet_term_radio_refractivity'))

    # ITU-R P.836 tests (Water vapour density)
    suite.addTest(ITUR836_6TestCase('test_surface_water_vapour_density'))
    suite.addTest(ITUR836_6TestCase('test_total_water_vapour_content'))

    # ITU-R P.836 tests (Rainfall rate)
    suite.addTest(ITUR837_6TestCase('test_rainfall_rate'))
    suite.addTest(ITUR837_7TestCase('test_rainfall_rate'))
    suite.addTest(ITUR837_7TestCase('test_rainfall_probability'))
    suite.addTest(ITUR837_7TestCase('test_rainfall_rate_R001'))

    # ITU-R P.836 tests (Rainfall specific attenuation)
    suite.addTest(ITUR838_3TestCase('test_rain_specific_attenuation'))

    # ITU-R P.839 tests (Rain height)
    suite.addTest(ITUR839_4TestCase('test_isoterm_0_deg'))
    suite.addTest(ITUR839_4TestCase('test_rain_height'))

    # ITU-R P.840 tests (Clouds attenuation)
#    suite.addTest(ITUR840_4TestCase('test_columnar_content_reduced_liquid'))
#    suite.addTest(ITUR840_4TestCase('test_cloud_attenuation'))
    suite.addTest(ITUR840_7TestCase('test_columnar_content_reduced_liquid'))
    suite.addTest(ITUR840_7TestCase('test_cloud_attenuation'))

    # ITU-R P.1511 tests (Topographic altitude)
    suite.addTest(ITUR1511_1TestCase('test_topographic_altitude'))
    suite.addTest(ITUR1511_2TestCase('test_topographic_altitude'))

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


class ITUR453_13TestCase(test.TestCase):
    def setUp(self):
        models.itu453.change_version(13)

    def test_wet_term_radio_refractivity(self):
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(
                3.133, 101.7, 50).value,
            128.14080027, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(
                22.9, -43.23, 50).value,
            104.35847467, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(
                23, 30, 50).value,
            36.47166667, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(
                25.78, -80.22, 50).value,
            113.2738672, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(
                28.717, 77.3, 50).value,
            75.66013547, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(
                33.94, 18.43, 50).value,
            80.14015964, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(
                41.9, 12.49, 50).value,
            61.21890044, places=5)
        self.assertAlmostEqual(
            models.itu453.map_wet_term_radio_refractivity(
                51.5, -0.14, 50).value,
            50.38926222, places=5)


class ITUR676_9TestCase(test.TestCase):

    def setUp(self):
        models.itu676.change_version(9)
        models.itu836.change_version(4)

    def test_gammaw(self):
        # The ITU models are non-sense and believe that the conversion between
        # Kelvin is 273 instead of 273.15
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(12, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.deg_C).value,
            0.00705700000, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(20, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.deg_C).value,
            0.06742720000, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(60, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.deg_C).value,
            0.11538020000, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(90, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.deg_C).value,
            0.25568340000, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(130, 1013.25, 4.98154290000,
                                        (5.9435147000 - 0.15) * u.deg_C).value,
            0.56358380000, places=5)

    def test_gamma0(self):
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        282.724 - 0.15).value,
            0.00941327, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        287.4834667 - 0.15).value,
            0.00898682, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        293.1487022 - 0.15).value,
            0.00851359, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        282.724 - 0.15).value,
            0.00941327, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        287.4834667 - 0.15).value,
            0.00898682, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        293.1487022 - 0.15).value,
            0.00851359, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        282.724 - 0.15).value,
            0.00941327, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        287.4834667 - 0.15).value,
            0.00898682, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        293.1487022 - 0.15).value,
            0.00851359, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        282.724 - 0.15).value,
            0.00941327, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        287.4834667 - 0.15).value,
            0.00898682, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        293.1487022 - 0.15).value,
            0.00851359, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        282.724 - 0.15).value,
            0.02043748, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        287.4834667 - 0.15).value,
            0.01954568, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        293.1487022 - 0.15).value,
            0.01856193, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        282.724 - 0.15).value,
            0.02043748, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        287.4834667 - 0.15).value,
            0.01954568, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        293.1487022 - 0.15).value,
            0.01856193, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        282.724 - 0.15).value,
            0.02043748, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        287.4834667 - 0.15).value,
            0.01954568, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        293.1487022 - 0.15).value,
            0.01856193, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        282.724 - 0.15).value,
            0.02043748, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        287.4834667 - 0.15).value,
            0.01954568, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        293.1487022 - 0.15).value,
            0.01856193, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        296.602 - 0.15).value,
            0.00824203, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 7.5, 296.7208533 - 0.15).value,
            0.0082329, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        296.602 - 0.15).value,
            0.00824203, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        296.7208533 - 0.15).value,
            0.0082329, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        296.602 - 0.15).value,
            0.00824203, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        296.7208533 - 0.15).value,
            0.0082329, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        296.602 - 0.15).value,
            0.00824203, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(
                14.25, 1013.25, 7.5, 296.7208533 - 0.15).value,
            0.0082329, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        296.602 - 0.15).value,
            0.01800011, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        296.7208533 - 0.15).value,
            0.01798125, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        296.602 - 0.15).value,
            0.01800011, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        296.7208533 - 0.15).value,
            0.01798125, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        296.602 - 0.15).value,
            0.01800011, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        296.7208533 - 0.15).value,
            0.01798125, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        296.602 - 0.15).value,
            0.01800011, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        296.7208533 - 0.15).value,
            0.01798125, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        299.0966578 - 0.15).value,
            0.00805331, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        297.9322267 - 0.15).value,
            0.00814064, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        287.444 - 0.15).value,
            0.00899025, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        299.0966578 - 0.15).value,
            0.00805331, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        297.9322267 - 0.15).value,
            0.00814064, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        287.444 - 0.15).value,
            0.00899025, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        299.0966578 - 0.15).value,
            0.00805331, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        297.9322267 - 0.15).value,
            0.00814064, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        287.444 - 0.15).value,
            0.00899025, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        299.0966578 - 0.15).value,
            0.00805331, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        297.9322267 - 0.15).value,
            0.00814064, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14.25, 1013.25, 7.5,
                                        287.444 - 0.15).value,
            0.00899025, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        299.0966578 - 0.15).value,
            0.01761077, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        297.9322267 - 0.15).value,
            0.01779083, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        287.444 - 0.15).value,
            0.01955282, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        299.0966578 - 0.15).value,
            0.01761077, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        297.9322267 - 0.15).value,
            0.01779083, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        287.444 - 0.15).value,
            0.01955282, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        299.0966578 - 0.15).value,
            0.01761077, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        297.9322267 - 0.15).value,
            0.01779083, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        287.444 - 0.15).value,
            0.01955282, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        299.0966578 - 0.15).value,
            0.01761077, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        297.9322267 - 0.15).value,
            0.01779083, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5,
                                        287.444 - 0.15).value,
            0.01955282, places=5)

#    def zenit_water_vapour_attenuation(self, lat, lon, el, p, f, alt):
#        gamma_w = models.itu676.zenit_water_vapour_attenuation(lat,
#                                                               lon,
#                                                               p,
#                                                               f,
#                                                               None,
#                                                               alt=alt).value
#        return gamma_w / np.sin(np.deg2rad(el))
#
#    def test_zenit_water_vapour_attenuation(self):
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                51.5, 359.86, 30.87067768, 1, 14.25, 0.06916422),
#            0.12789267, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                41.9, 12.49, 40.97052773, 1, 14.25, 0.05670104),
#            0.10865204, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                33.94, 18.43, 47.91280491, 1, 14.25, 0),
#            0.10205633, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                51.5, 359.86, 30.87067768, 0.1, 14.25, 0.06916422),
#            0.15315923, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                41.9, 12.49, 40.97052773, 0.1, 14.25, 0.05670104),
#            0.12223686, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                33.94, 18.43, 47.91280491, 0.1, 14.25, 0),
#            0.12410189, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                51.5, 359.86, 30.87067768, 0.01, 14.25, 0.06916422),
#            0.15315923,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                41.9, 12.49, 40.97052773, 0.01, 14.25, 0.05670104),
#            0.12223686,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                33.94, 18.43, 47.91280491, 0.01, 14.25, 0),
#            0.12410189,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                51.5, 359.86, 30.87067768, 0.001, 14.25, 0.06916422),
#            0.15315923,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                41.9, 12.49, 40.97052773, 0.001, 14.25, 0.05670104),
#            0.12223686,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                33.94, 18.43, 47.91280491, 0.001, 14.25, 0),
#            0.12410189,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                51.5, 359.86, 30.87067768, 1, 29, 0.06916422),
#            0.60896934, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                41.9, 12.49, 40.97052773, 1, 29, 0.05670104),
#            0.51690529, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                33.94, 18.43, 47.91280491, 1, 29, 0),
#            0.48519817, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                51.5, 359.86, 30.87067768, 0.1, 29, 0.06916422),
#            0.72784676, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                41.9, 12.49, 40.97052773, 0.1, 29, 0.05670104),
#            0.58076456, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                33.94, 18.43, 47.91280491, 0.1, 29, 0),
#            0.58863533, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                51.5, 359.86, 30.87067768, 0.01, 29, 0.06916422),
#            0.72784676,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                41.9, 12.49, 40.97052773, 0.01, 29, 0.05670104),
#            0.58076456,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                33.94, 18.43, 47.91280491, 0.01, 29, 0),
#            0.58863533,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                51.5, 359.86, 30.87067768, 0.001, 29, 0.06916422),
#            0.72784676,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                41.9, 12.49, 40.97052773, 0.001, 29, 0.05670104),
#            0.58076456,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                33.94, 18.43, 47.91280491, 0.001, 29, 0),
#            0.58863533,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                22.9, 316.77, 59.81487174, 1, 14.25, 0),
#            0.1181882, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                25.78, 279.78, 49.20900369, 1, 14.25, 0.00007511),
#            0.16093386, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                22.9, 316.77, 59.81487174, 0.1, 14.25, 0),
#            0.13730617, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                25.78, 279.78, 49.20900369, 0.1, 14.25, 0.00007511),
#            0.17798382, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                22.9, 316.77, 59.81487174, 0.01, 14.25, 0),
#            0.13730617,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                25.78, 279.78, 49.20900369, 0.01, 14.25, 0.00007511),
#            0.17798382,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                22.9, 316.77, 59.81487174, 0.001, 14.25, 0),
#            0.13730617,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                25.78, 279.78, 49.20900369, 0.001, 14.25, 0.00007511),
#            0.17798382,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                22.9, 316.77, 59.81487174, 1, 29, 0),
#            0.55983815, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                25.78, 279.78, 49.20900369, 1, 29, 0.00007511),
#            0.76047761, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                22.9, 316.77, 59.81487174, 0.1, 29, 0),
#            0.64906814, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                25.78, 279.78, 49.20900369, 0.1, 29, 0.00007511),
#            0.83981774, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                22.9, 316.77, 59.81487174, 0.01, 29, 0),
#            0.64906814,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                25.78, 279.78, 49.20900369, 0.01, 29, 0.00007511),
#            0.83981774,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                22.9, 316.77, 59.81487174, 0.001, 29, 0),
#            0.64906814,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                25.78, 279.78, 49.20900369, 0.001, 29, 0.00007511),
#            0.83981774,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                28.717, 77.3, 55.90591362, 1, 14.25, 0.21755946),
#            0.18628614, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                3.133, 101.7, 67.76751981, 1, 14.25, 0.23610446),
#            0.13468573, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                9.05, 38.7, 38.14104832, 1, 14.25, 2.45000492),
#            0.08369587, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                28.717, 77.3, 55.90591362, 0.1, 14.25, 0.21755946),
#            0.20242415, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                3.133, 101.7, 67.76751981, 0.1, 14.25, 0.23610446),
#            0.14372476, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                9.05, 38.7, 38.14104832, 0.1, 14.25, 2.45000492),
#            0.09153026, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                28.717, 77.3, 55.90591362, 0.01, 14.25, 0.21755946),
#            0.20242415,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                3.133, 101.7, 67.76751981, 0.01, 14.25, 0.23610446),
#            0.14372476,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                9.05, 38.7, 38.14104832, 0.01, 14.25, 2.45000492),
#            0.09153026,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                28.717, 77.3, 55.90591362, 0.001, 14.25, 0.21755946),
#            0.20242415,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                3.133, 101.7, 67.76751981, 0.001, 14.25, 0.23610446),
#            0.14372476,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                9.05, 38.7, 38.14104832, 0.001, 14.25, 2.45000492),
#            0.09153026,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                28.717, 77.3, 55.90591362, 1, 29, 0.21755946),
#            0.8771945, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                3.133, 101.7, 67.76751981, 1, 29, 0.23610446),
#            0.63623574, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                9.05, 38.7, 38.14104832, 1, 29, 2.45000492),
#            0.39942177, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                28.717, 77.3, 55.90591362, 0.1, 29, 0.21755946),
#            0.95194476, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                3.133, 101.7, 67.76751981, 0.1, 29, 0.23610446),
#            0.67829402, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                9.05, 38.7, 38.14104832, 0.1, 29, 2.45000492),
#            0.43646179, places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                28.717, 77.3, 55.90591362, 0.01, 29, 0.21755946),
#            0.95194476,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                3.133, 101.7, 67.76751981, 0.01, 29, 0.23610446),
#            0.67829402,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                9.05, 38.7, 38.14104832, 0.01, 29, 2.45000492),
#            0.43646179,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                28.717, 77.3, 55.90591362, 0.001, 29, 0.21755946),
#            0.95194476,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                3.133, 101.7, 67.76751981, 0.001, 29, 0.23610446),
#            0.67829402,  places=5)
#        self.assertAlmostEqual(
#            self.zenit_water_vapour_attenuation(
#                9.05, 38.7, 38.14104832, 0.001, 29, 2.45000492),
#            0.43646179,  places=5)


class ITUR676_11TestCase(test.TestCase):

    def setUp(self):
        models.itu676.change_version(11)
        models.itu836.change_version(6)
        models.itu1511.change_version(1)

    def test_gammaw_exact(self):
        # The ITU models are non-sense and believe that the conversion between
        # Kelvin is 273 instead of 273.15
        self.assertAlmostEqual(
            models.itu676.gammaw_exact(12, 1013.25, 7.5, 288.15).value,
            0.00953539, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_exact(20, 1013.25, 7.5, 288.15).value,
            0.09704730, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_exact(60, 1013.25, 7.5, 288.15).value,
            0.15484184, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_exact(90, 1013.25, 7.5, 288.15).value,
            0.34197339, places=5)
        self.assertAlmostEqual(
            models.itu676.gammaw_exact(130, 1013.25, 7.5, 288.15).value,
            0.75184470, places=5)

    def test_gamma0_exact(self):
        self.assertAlmostEqual(
            models.itu676.gamma0_exact(12, 1013.25, 7.5, 288.15).value,
            0.00869826, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_exact(20, 1013.25, 7.5, 288.15).value,
            0.01188355, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_exact(60, 1013.25, 7.5, 288.15).value,
            14.62347480, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_exact(90, 1013.25, 7.5, 288.15).value,
            0.03886971, places=5)
        self.assertAlmostEqual(
            models.itu676.gamma0_exact(130, 1013.25, 7.5, 288.15).value,
            0.04150908, places=5)

    def test_gammaw_approx(self):
        # The ITU models are non-sense and believe that the conversion between
        # Kelvin is 273 instead of 273.15
        self.assertAlmostEqual(
            models.itu676.gammaw_approx(1, 1013.25, 7.5, 288.15).value,
            5.06e-05, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(2, 1013.25, 7.5, 288.15).value,
            0.000203124, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(3, 1013.25, 7.5, 288.15).value,
            0.000459962, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(4, 1013.25, 7.5, 288.15).value,
            0.000825295, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(5, 1013.25, 7.5, 288.15).value,
            0.001305574, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(6, 1013.25, 7.5, 288.15).value,
            0.001910194, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(7, 1013.25, 7.5, 288.15).value,
            0.00265257, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(8, 1013.25, 7.5, 288.15).value,
            0.00355178, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(9, 1013.25, 7.5, 288.15).value,
            0.00463511, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(10, 1013.25, 7.5, 288.15).value,
            0.005942065, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(11, 1013.25, 7.5, 288.15).value,
            0.007530789, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(12, 1013.25, 7.5, 288.15).value,
            0.009488627, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(13, 1013.25, 7.5, 288.15).value,
            0.01194992, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(14, 1013.25, 7.5, 288.15).value,
            0.015126834, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(15, 1013.25, 7.5, 288.15).value,
            0.019364141, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(16, 1013.25, 7.5, 288.15).value,
            0.025238305, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(17, 1013.25, 7.5, 288.15).value,
            0.033736014, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(18, 1013.25, 7.5, 288.15).value,
            0.04655406, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(19, 1013.25, 7.5, 288.15).value,
            0.066459485, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(20, 1013.25, 7.5, 288.15).value,
            0.096940958, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(21, 1013.25, 7.5, 288.15).value,
            0.137887422, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(22, 1013.25, 7.5, 288.15).value,
            0.17418431, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(23, 1013.25, 7.5, 288.15).value,
            0.180393135, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(24, 1013.25, 7.5, 288.15).value,
            0.15839854, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(25, 1013.25, 7.5, 288.15).value,
            0.130540688, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(26, 1013.25, 7.5, 288.15).value,
            0.108338372, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(27, 1013.25, 7.5, 288.15).value,
            0.092962551, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(28, 1013.25, 7.5, 288.15).value,
            0.082791566, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(29, 1013.25, 7.5, 288.15).value,
            0.076209755, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(30, 1013.25, 7.5, 288.15).value,
            0.072073391, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(31, 1013.25, 7.5, 288.15).value,
            0.069632181, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(32, 1013.25, 7.5, 288.15).value,
            0.06839841, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(33, 1013.25, 7.5, 288.15).value,
            0.068050819, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(34, 1013.25, 7.5, 288.15).value,
            0.068373336, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(35, 1013.25, 7.5, 288.15).value,
            0.069217296, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(36, 1013.25, 7.5, 288.15).value,
            0.070478105, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(37, 1013.25, 7.5, 288.15).value,
            0.072080617, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(38, 1013.25, 7.5, 288.15).value,
            0.073969796, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(39, 1013.25, 7.5, 288.15).value,
            0.076104615, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(40, 1013.25, 7.5, 288.15).value,
            0.078454003, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(41, 1013.25, 7.5, 288.15).value,
            0.080994086, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(42, 1013.25, 7.5, 288.15).value,
            0.08370628, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(43, 1013.25, 7.5, 288.15).value,
            0.086575946, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(44, 1013.25, 7.5, 288.15).value,
            0.089591433, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(45, 1013.25, 7.5, 288.15).value,
            0.092743375, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(46, 1013.25, 7.5, 288.15).value,
            0.096024183, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(47, 1013.25, 7.5, 288.15).value,
            0.099427654, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(48, 1013.25, 7.5, 288.15).value,
            0.102948692, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(49, 1013.25, 7.5, 288.15).value,
            0.106583076, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(50, 1013.25, 7.5, 288.15).value,
            0.110327299, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(51, 1013.25, 7.5, 288.15).value,
            0.11417843, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(52, 1013.25, 7.5, 288.15).value,
            0.118134012, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(53, 1013.25, 7.5, 288.15).value,
            0.122191981, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(54, 1013.25, 7.5, 288.15).value,
            0.126350598, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(55, 1013.25, 7.5, 288.15).value,
            0.130608397, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(56, 1013.25, 7.5, 288.15).value,
            0.134964144, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(57, 1013.25, 7.5, 288.15).value,
            0.139416798, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(58, 1013.25, 7.5, 288.15).value,
            0.143965489, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(59, 1013.25, 7.5, 288.15).value,
            0.148609489, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(60, 1013.25, 7.5, 288.15).value,
            0.153348196, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(61, 1013.25, 7.5, 288.15).value,
            0.158181114, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(62, 1013.25, 7.5, 288.15).value,
            0.163107847, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(63, 1013.25, 7.5, 288.15).value,
            0.168128079, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(64, 1013.25, 7.5, 288.15).value,
            0.173241572, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(65, 1013.25, 7.5, 288.15).value,
            0.178448154, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(66, 1013.25, 7.5, 288.15).value,
            0.183747712, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(67, 1013.25, 7.5, 288.15).value,
            0.18914019, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(68, 1013.25, 7.5, 288.15).value,
            0.194625582, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(69, 1013.25, 7.5, 288.15).value,
            0.200203926, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(70, 1013.25, 7.5, 288.15).value,
            0.205875306, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(71, 1013.25, 7.5, 288.15).value,
            0.211639845, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(72, 1013.25, 7.5, 288.15).value,
            0.217497702, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(73, 1013.25, 7.5, 288.15).value,
            0.223449076, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(74, 1013.25, 7.5, 288.15).value,
            0.229494196, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(75, 1013.25, 7.5, 288.15).value,
            0.235633329, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(76, 1013.25, 7.5, 288.15).value,
            0.241866771, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(77, 1013.25, 7.5, 288.15).value,
            0.248194851, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(78, 1013.25, 7.5, 288.15).value,
            0.254617931, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(79, 1013.25, 7.5, 288.15).value,
            0.261136401, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(80, 1013.25, 7.5, 288.15).value,
            0.267750686, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(81, 1013.25, 7.5, 288.15).value,
            0.274461239, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(82, 1013.25, 7.5, 288.15).value,
            0.281268547, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(83, 1013.25, 7.5, 288.15).value,
            0.28817313, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(84, 1013.25, 7.5, 288.15).value,
            0.295175539, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(85, 1013.25, 7.5, 288.15).value,
            0.302276362, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(86, 1013.25, 7.5, 288.15).value,
            0.309476219, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(87, 1013.25, 7.5, 288.15).value,
            0.316775769, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(88, 1013.25, 7.5, 288.15).value,
            0.324175708, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(89, 1013.25, 7.5, 288.15).value,
            0.331676772, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(90, 1013.25, 7.5, 288.15).value,
            0.339279738, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(91, 1013.25, 7.5, 288.15).value,
            0.346985426, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(92, 1013.25, 7.5, 288.15).value,
            0.354794703, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(93, 1013.25, 7.5, 288.15).value,
            0.362708483, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(94, 1013.25, 7.5, 288.15).value,
            0.370727732, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(95, 1013.25, 7.5, 288.15).value,
            0.378853468, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(96, 1013.25, 7.5, 288.15).value,
            0.387086768, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(97, 1013.25, 7.5, 288.15).value,
            0.395428769, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(98, 1013.25, 7.5, 288.15).value,
            0.403880673, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(99, 1013.25, 7.5, 288.15).value,
            0.412443748, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(100, 1013.25, 7.5, 288.15).value,
            0.421119341, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(101, 1013.25, 7.5, 288.15).value,
            0.429908872, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(102, 1013.25, 7.5, 288.15).value,
            0.438813848, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(103, 1013.25, 7.5, 288.15).value,
            0.447835866, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(104, 1013.25, 7.5, 288.15).value,
            0.456976619, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(105, 1013.25, 7.5, 288.15).value,
            0.466237905, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(106, 1013.25, 7.5, 288.15).value,
            0.475621633, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(107, 1013.25, 7.5, 288.15).value,
            0.485129833, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(108, 1013.25, 7.5, 288.15).value,
            0.494764666, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(109, 1013.25, 7.5, 288.15).value,
            0.504528432, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(110, 1013.25, 7.5, 288.15).value,
            0.514423584, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(111, 1013.25, 7.5, 288.15).value,
            0.524452741, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(112, 1013.25, 7.5, 288.15).value,
            0.5346187, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(113, 1013.25, 7.5, 288.15).value,
            0.54492445, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(114, 1013.25, 7.5, 288.15).value,
            0.555373195, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(115, 1013.25, 7.5, 288.15).value,
            0.565968366, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(116, 1013.25, 7.5, 288.15).value,
            0.576713646, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(117, 1013.25, 7.5, 288.15).value,
            0.58761299, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(118, 1013.25, 7.5, 288.15).value,
            0.598670654, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(119, 1013.25, 7.5, 288.15).value,
            0.609891221, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(120, 1013.25, 7.5, 288.15).value,
            0.621279631, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(121, 1013.25, 7.5, 288.15).value,
            0.63284122, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(122, 1013.25, 7.5, 288.15).value,
            0.644581758, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(123, 1013.25, 7.5, 288.15).value,
            0.656507491, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(124, 1013.25, 7.5, 288.15).value,
            0.668625191, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(125, 1013.25, 7.5, 288.15).value,
            0.680942215, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(126, 1013.25, 7.5, 288.15).value,
            0.69346656, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(127, 1013.25, 7.5, 288.15).value,
            0.70620694, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(128, 1013.25, 7.5, 288.15).value,
            0.719172861, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(129, 1013.25, 7.5, 288.15).value,
            0.73237471, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(130, 1013.25, 7.5, 288.15).value,
            0.745823861, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(131, 1013.25, 7.5, 288.15).value,
            0.759532783, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(132, 1013.25, 7.5, 288.15).value,
            0.773515178, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(133, 1013.25, 7.5, 288.15).value,
            0.787786128, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(134, 1013.25, 7.5, 288.15).value,
            0.802362262, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(135, 1013.25, 7.5, 288.15).value,
            0.817261961, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(136, 1013.25, 7.5, 288.15).value,
            0.832505575, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(137, 1013.25, 7.5, 288.15).value,
            0.848115693, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(138, 1013.25, 7.5, 288.15).value,
            0.864117433, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(139, 1013.25, 7.5, 288.15).value,
            0.880538802, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(140, 1013.25, 7.5, 288.15).value,
            0.897411097, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(141, 1013.25, 7.5, 288.15).value,
            0.914769381, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(142, 1013.25, 7.5, 288.15).value,
            0.93265304, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(143, 1013.25, 7.5, 288.15).value,
            0.951106434, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(144, 1013.25, 7.5, 288.15).value,
            0.970179674, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(145, 1013.25, 7.5, 288.15).value,
            0.989929528, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(146, 1013.25, 7.5, 288.15).value,
            1.010420514, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(147, 1013.25, 7.5, 288.15).value,
            1.03172619, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(148, 1013.25, 7.5, 288.15).value,
            1.053930717, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(149, 1013.25, 7.5, 288.15).value,
            1.077130727, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(150, 1013.25, 7.5, 288.15).value,
            1.101437596, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(151, 1013.25, 7.5, 288.15).value,
            1.126980206, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(152, 1013.25, 7.5, 288.15).value,
            1.153908327, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(153, 1013.25, 7.5, 288.15).value,
            1.182396776, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(154, 1013.25, 7.5, 288.15).value,
            1.212650574, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(155, 1013.25, 7.5, 288.15).value,
            1.244911365, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(156, 1013.25, 7.5, 288.15).value,
            1.279465482, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(157, 1013.25, 7.5, 288.15).value,
            1.316654129, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(158, 1013.25, 7.5, 288.15).value,
            1.356886363, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(159, 1013.25, 7.5, 288.15).value,
            1.400655759, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(160, 1013.25, 7.5, 288.15).value,
            1.448562004, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(161, 1013.25, 7.5, 288.15).value,
            1.501339131, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(162, 1013.25, 7.5, 288.15).value,
            1.559892824, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(163, 1013.25, 7.5, 288.15).value,
            1.625350216, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(164, 1013.25, 7.5, 288.15).value,
            1.699127159, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(165, 1013.25, 7.5, 288.15).value,
            1.783020212, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(166, 1013.25, 7.5, 288.15).value,
            1.87933414, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(167, 1013.25, 7.5, 288.15).value,
            1.991061177, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(168, 1013.25, 7.5, 288.15).value,
            2.122137016, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(169, 1013.25, 7.5, 288.15).value,
            2.277812508, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(170, 1013.25, 7.5, 288.15).value,
            2.465203115, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(171, 1013.25, 7.5, 288.15).value,
            2.694116678, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(172, 1013.25, 7.5, 288.15).value,
            2.978325696, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(173, 1013.25, 7.5, 288.15).value,
            3.337563176, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(174, 1013.25, 7.5, 288.15).value,
            3.80071648, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(175, 1013.25, 7.5, 288.15).value,
            4.411026238, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(176, 1013.25, 7.5, 288.15).value,
            5.23462829, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(177, 1013.25, 7.5, 288.15).value,
            6.374446918, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(178, 1013.25, 7.5, 288.15).value,
            7.991434174, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(179, 1013.25, 7.5, 288.15).value,
            10.33006475, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(180, 1013.25, 7.5, 288.15).value,
            13.71659631, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(181, 1013.25, 7.5, 288.15).value,
            18.39188186, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(182, 1013.25, 7.5, 288.15).value,
            23.83194406, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(183, 1013.25, 7.5, 288.15).value,
            27.67449812, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(184, 1013.25, 7.5, 288.15).value,
            27.03213321, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(185, 1013.25, 7.5, 288.15).value,
            22.60135009, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(186, 1013.25, 7.5, 288.15).value,
            17.47071693, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(187, 1013.25, 7.5, 288.15).value,
            13.32603388, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(188, 1013.25, 7.5, 288.15).value,
            10.35715037, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(189, 1013.25, 7.5, 288.15).value,
            8.290514155, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(190, 1013.25, 7.5, 288.15).value,
            6.842342894, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(191, 1013.25, 7.5, 288.15).value,
            5.808188543, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(192, 1013.25, 7.5, 288.15).value,
            5.053346248, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(193, 1013.25, 7.5, 288.15).value,
            4.490555513, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(194, 1013.25, 7.5, 288.15).value,
            4.062799398, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(195, 1013.25, 7.5, 288.15).value,
            3.73214099, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(196, 1013.25, 7.5, 288.15).value,
            3.472798768, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(197, 1013.25, 7.5, 288.15).value,
            3.266872193, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(198, 1013.25, 7.5, 288.15).value,
            3.101674612, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(199, 1013.25, 7.5, 288.15).value,
            2.968040373, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(200, 1013.25, 7.5, 288.15).value,
            2.859229328, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(201, 1013.25, 7.5, 288.15).value,
            2.77020383, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(202, 1013.25, 7.5, 288.15).value,
            2.697142323, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(203, 1013.25, 7.5, 288.15).value,
            2.637106098, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(204, 1013.25, 7.5, 288.15).value,
            2.587807034, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(205, 1013.25, 7.5, 288.15).value,
            2.547443114, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(206, 1013.25, 7.5, 288.15).value,
            2.514580199, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(207, 1013.25, 7.5, 288.15).value,
            2.488065853, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(208, 1013.25, 7.5, 288.15).value,
            2.466965732, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(209, 1013.25, 7.5, 288.15).value,
            2.450516051, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(210, 1013.25, 7.5, 288.15).value,
            2.43808768, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(211, 1013.25, 7.5, 288.15).value,
            2.429158726, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(212, 1013.25, 7.5, 288.15).value,
            2.423293411, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(213, 1013.25, 7.5, 288.15).value,
            2.420125627, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(214, 1013.25, 7.5, 288.15).value,
            2.419346051, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(215, 1013.25, 7.5, 288.15).value,
            2.420691947, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(216, 1013.25, 7.5, 288.15).value,
            2.423939044, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(217, 1013.25, 7.5, 288.15).value,
            2.428895021, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(218, 1013.25, 7.5, 288.15).value,
            2.435394244, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(219, 1013.25, 7.5, 288.15).value,
            2.443293492, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(220, 1013.25, 7.5, 288.15).value,
            2.452468459, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(221, 1013.25, 7.5, 288.15).value,
            2.462810881, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(222, 1013.25, 7.5, 288.15).value,
            2.474226165, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(223, 1013.25, 7.5, 288.15).value,
            2.486631411, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(224, 1013.25, 7.5, 288.15).value,
            2.499953772, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(225, 1013.25, 7.5, 288.15).value,
            2.514129072, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(226, 1013.25, 7.5, 288.15).value,
            2.529100646, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(227, 1013.25, 7.5, 288.15).value,
            2.544818361, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(228, 1013.25, 7.5, 288.15).value,
            2.561237781, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(229, 1013.25, 7.5, 288.15).value,
            2.578319459, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(230, 1013.25, 7.5, 288.15).value,
            2.596028333, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(231, 1013.25, 7.5, 288.15).value,
            2.614333207, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(232, 1013.25, 7.5, 288.15).value,
            2.633206307, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(233, 1013.25, 7.5, 288.15).value,
            2.652622892, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(234, 1013.25, 7.5, 288.15).value,
            2.672560927, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(235, 1013.25, 7.5, 288.15).value,
            2.693000794, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(236, 1013.25, 7.5, 288.15).value,
            2.713925039, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(237, 1013.25, 7.5, 288.15).value,
            2.735318161, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(238, 1013.25, 7.5, 288.15).value,
            2.757166418, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(239, 1013.25, 7.5, 288.15).value,
            2.779457666, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(240, 1013.25, 7.5, 288.15).value,
            2.80218121, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(241, 1013.25, 7.5, 288.15).value,
            2.825327683, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(242, 1013.25, 7.5, 288.15).value,
            2.848888936, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(243, 1013.25, 7.5, 288.15).value,
            2.87285794, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(244, 1013.25, 7.5, 288.15).value,
            2.897228704, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(245, 1013.25, 7.5, 288.15).value,
            2.921996202, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(246, 1013.25, 7.5, 288.15).value,
            2.947156311, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(247, 1013.25, 7.5, 288.15).value,
            2.972705756, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(248, 1013.25, 7.5, 288.15).value,
            2.998642066, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(249, 1013.25, 7.5, 288.15).value,
            3.024963531, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(250, 1013.25, 7.5, 288.15).value,
            3.051669175, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(251, 1013.25, 7.5, 288.15).value,
            3.078758722, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(252, 1013.25, 7.5, 288.15).value,
            3.106232585, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(253, 1013.25, 7.5, 288.15).value,
            3.13409184, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(254, 1013.25, 7.5, 288.15).value,
            3.162338224, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(255, 1013.25, 7.5, 288.15).value,
            3.190974123, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(256, 1013.25, 7.5, 288.15).value,
            3.220002576, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(257, 1013.25, 7.5, 288.15).value,
            3.249427276, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(258, 1013.25, 7.5, 288.15).value,
            3.27925258, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(259, 1013.25, 7.5, 288.15).value,
            3.309483522, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(260, 1013.25, 7.5, 288.15).value,
            3.340125831, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(261, 1013.25, 7.5, 288.15).value,
            3.371185956, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(262, 1013.25, 7.5, 288.15).value,
            3.402671096, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(263, 1013.25, 7.5, 288.15).value,
            3.434589233, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(264, 1013.25, 7.5, 288.15).value,
            3.466949175, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(265, 1013.25, 7.5, 288.15).value,
            3.499760606, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(266, 1013.25, 7.5, 288.15).value,
            3.533034141, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(267, 1013.25, 7.5, 288.15).value,
            3.566781392, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(268, 1013.25, 7.5, 288.15).value,
            3.601015043, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(269, 1013.25, 7.5, 288.15).value,
            3.635748934, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(270, 1013.25, 7.5, 288.15).value,
            3.670998161, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(271, 1013.25, 7.5, 288.15).value,
            3.706779184, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(272, 1013.25, 7.5, 288.15).value,
            3.743109957, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(273, 1013.25, 7.5, 288.15).value,
            3.780010072, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(274, 1013.25, 7.5, 288.15).value,
            3.817500924, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(275, 1013.25, 7.5, 288.15).value,
            3.855605898, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(276, 1013.25, 7.5, 288.15).value,
            3.894350591, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(277, 1013.25, 7.5, 288.15).value,
            3.933763053, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(278, 1013.25, 7.5, 288.15).value,
            3.97387408, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(279, 1013.25, 7.5, 288.15).value,
            4.014717535, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(280, 1013.25, 7.5, 288.15).value,
            4.056330735, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(281, 1013.25, 7.5, 288.15).value,
            4.098754887, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(282, 1013.25, 7.5, 288.15).value,
            4.142035602, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(283, 1013.25, 7.5, 288.15).value,
            4.186223487, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(284, 1013.25, 7.5, 288.15).value,
            4.231374849, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(285, 1013.25, 7.5, 288.15).value,
            4.277552506, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(286, 1013.25, 7.5, 288.15).value,
            4.324826757, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(287, 1013.25, 7.5, 288.15).value,
            4.373276518, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(288, 1013.25, 7.5, 288.15).value,
            4.422990681, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(289, 1013.25, 7.5, 288.15).value,
            4.474069728, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(290, 1013.25, 7.5, 288.15).value,
            4.526627666, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(291, 1013.25, 7.5, 288.15).value,
            4.580794366, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(292, 1013.25, 7.5, 288.15).value,
            4.63671838, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(293, 1013.25, 7.5, 288.15).value,
            4.694570386, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(294, 1013.25, 7.5, 288.15).value,
            4.754547391, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(295, 1013.25, 7.5, 288.15).value,
            4.816877916, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(296, 1013.25, 7.5, 288.15).value,
            4.881828419, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(297, 1013.25, 7.5, 288.15).value,
            4.949711312, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(298, 1013.25, 7.5, 288.15).value,
            5.020895036, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(299, 1013.25, 7.5, 288.15).value,
            5.095816817, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(300, 1013.25, 7.5, 288.15).value,
            5.174998967, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(301, 1013.25, 7.5, 288.15).value,
            5.259069863, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(302, 1013.25, 7.5, 288.15).value,
            5.348791238, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(303, 1013.25, 7.5, 288.15).value,
            5.445094008, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(304, 1013.25, 7.5, 288.15).value,
            5.549125837, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(305, 1013.25, 7.5, 288.15).value,
            5.662315008, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(306, 1013.25, 7.5, 288.15).value,
            5.786457272, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(307, 1013.25, 7.5, 288.15).value,
            5.923835584, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(308, 1013.25, 7.5, 288.15).value,
            6.077387588, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(309, 1013.25, 7.5, 288.15).value,
            6.250943769, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(310, 1013.25, 7.5, 288.15).value,
            6.449572043, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(311, 1013.25, 7.5, 288.15).value,
            6.6800861, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(312, 1013.25, 7.5, 288.15).value,
            6.951811263, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(313, 1013.25, 7.5, 288.15).value,
            7.277764947, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(314, 1013.25, 7.5, 288.15).value,
            7.676520947, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(315, 1013.25, 7.5, 288.15).value,
            8.175226571, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(316, 1013.25, 7.5, 288.15).value,
            8.814587905, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(317, 1013.25, 7.5, 288.15).value,
            9.657150573, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(318, 1013.25, 7.5, 288.15).value,
            10.80040832, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(319, 1013.25, 7.5, 288.15).value,
            12.39287203, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(320, 1013.25, 7.5, 288.15).value,
            14.63291434, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(321, 1013.25, 7.5, 288.15).value,
            17.71848802, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(322, 1013.25, 7.5, 288.15).value,
            21.89833011, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(323, 1013.25, 7.5, 288.15).value,
            27.52921207, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(324, 1013.25, 7.5, 288.15).value,
            33.93584273, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(325, 1013.25, 7.5, 288.15).value,
            37.82487596, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(326, 1013.25, 7.5, 288.15).value,
            35.8615979, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(327, 1013.25, 7.5, 288.15).value,
            29.89188489, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(328, 1013.25, 7.5, 288.15).value,
            23.80724266, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(329, 1013.25, 7.5, 288.15).value,
            19.19466647, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(330, 1013.25, 7.5, 288.15).value,
            16.01196137, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(331, 1013.25, 7.5, 288.15).value,
            13.85529573, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(332, 1013.25, 7.5, 288.15).value,
            12.38126427, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(333, 1013.25, 7.5, 288.15).value,
            11.35803945, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(334, 1013.25, 7.5, 288.15).value,
            10.63752623, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(335, 1013.25, 7.5, 288.15).value,
            10.12544, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(336, 1013.25, 7.5, 288.15).value,
            9.760866789, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(337, 1013.25, 7.5, 288.15).value,
            9.503637932, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(338, 1013.25, 7.5, 288.15).value,
            9.326698304, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(339, 1013.25, 7.5, 288.15).value,
            9.211467317, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(340, 1013.25, 7.5, 288.15).value,
            9.144973785, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(341, 1013.25, 7.5, 288.15).value,
            9.118051719, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(342, 1013.25, 7.5, 288.15).value,
            9.124181675, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(343, 1013.25, 7.5, 288.15).value,
            9.158733184, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(344, 1013.25, 7.5, 288.15).value,
            9.218462001, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(345, 1013.25, 7.5, 288.15).value,
            9.301173194, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(346, 1013.25, 7.5, 288.15).value,
            9.405494946, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(347, 1013.25, 7.5, 288.15).value,
            9.53072851, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(348, 1013.25, 7.5, 288.15).value,
            9.676752463, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(349, 1013.25, 7.5, 288.15).value,
            9.843967537, places=5)

        self.assertAlmostEqual(
            models.itu676.gammaw_approx(350, 1013.25, 7.5, 288.15).value,
            10.03327368, places=5)

    def test_gamma0_approx(self):
        self.assertAlmostEqual(
            models.itu676.gamma0_approx(1, 1013.25, 7.5, 288.15).value,
            0.005388658, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(2, 1013.25, 7.5, 288.15).value,
            0.006716038, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(3, 1013.25, 7.5, 288.15).value,
            0.00707596, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(4, 1013.25, 7.5, 288.15).value,
            0.007258969, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(5, 1013.25, 7.5, 288.15).value,
            0.007400426, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(6, 1013.25, 7.5, 288.15).value,
            0.007537212, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(7, 1013.25, 7.5, 288.15).value,
            0.007682905, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(8, 1013.25, 7.5, 288.15).value,
            0.007843794, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(9, 1013.25, 7.5, 288.15).value,
            0.008023466, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(10, 1013.25, 7.5, 288.15).value,
            0.008224416, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(11, 1013.25, 7.5, 288.15).value,
            0.008448705, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(12, 1013.25, 7.5, 288.15).value,
            0.008698263, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(13, 1013.25, 7.5, 288.15).value,
            0.008975056, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(14, 1013.25, 7.5, 288.15).value,
            0.009281177, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(15, 1013.25, 7.5, 288.15).value,
            0.009618923, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(16, 1013.25, 7.5, 288.15).value,
            0.009990845, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(17, 1013.25, 7.5, 288.15).value,
            0.010399811, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(18, 1013.25, 7.5, 288.15).value,
            0.010849054, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(19, 1013.25, 7.5, 288.15).value,
            0.011342243, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(20, 1013.25, 7.5, 288.15).value,
            0.011883547, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(21, 1013.25, 7.5, 288.15).value,
            0.012477725, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(22, 1013.25, 7.5, 288.15).value,
            0.013130219, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(23, 1013.25, 7.5, 288.15).value,
            0.013847273, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(24, 1013.25, 7.5, 288.15).value,
            0.014636078, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(25, 1013.25, 7.5, 288.15).value,
            0.015504937, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(26, 1013.25, 7.5, 288.15).value,
            0.016463481, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(27, 1013.25, 7.5, 288.15).value,
            0.017522921, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(28, 1013.25, 7.5, 288.15).value,
            0.018696367, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(29, 1013.25, 7.5, 288.15).value,
            0.019999221, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(30, 1013.25, 7.5, 288.15).value,
            0.021449673, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(31, 1013.25, 7.5, 288.15).value,
            0.023069328, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(32, 1013.25, 7.5, 288.15).value,
            0.024883993, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(33, 1013.25, 7.5, 288.15).value,
            0.026924712, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(34, 1013.25, 7.5, 288.15).value,
            0.029229084, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(35, 1013.25, 7.5, 288.15).value,
            0.031843013, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(36, 1013.25, 7.5, 288.15).value,
            0.034823023, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(37, 1013.25, 7.5, 288.15).value,
            0.038239374, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(38, 1013.25, 7.5, 288.15).value,
            0.042180317, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(39, 1013.25, 7.5, 288.15).value,
            0.046757999, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(40, 1013.25, 7.5, 288.15).value,
            0.052116797, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(41, 1013.25, 7.5, 288.15).value,
            0.058445339, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(42, 1013.25, 7.5, 288.15).value,
            0.065994232, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(43, 1013.25, 7.5, 288.15).value,
            0.075102941, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(44, 1013.25, 7.5, 288.15).value,
            0.086241846, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(45, 1013.25, 7.5, 288.15).value,
            0.100080659, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(46, 1013.25, 7.5, 288.15).value,
            0.117605188, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(47, 1013.25, 7.5, 288.15).value,
            0.140329657, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(48, 1013.25, 7.5, 288.15).value,
            0.170719546, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(49, 1013.25, 7.5, 288.15).value,
            0.213175063, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(50, 1013.25, 7.5, 288.15).value,
            0.277268297, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(51, 1013.25, 7.5, 288.15).value,
            0.389670239, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(52, 1013.25, 7.5, 288.15).value,
            0.618429331, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(53, 1013.25, 7.5, 288.15).value,
            1.126611463, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(54, 1013.25, 7.5, 288.15).value,
            2.211541194, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(55, 1013.25, 7.5, 288.15).value,
            4.193281287, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(56, 1013.25, 7.5, 288.15).value,
            7.055044748, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(57, 1013.25, 7.5, 288.15).value,
            10.0652395, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(58, 1013.25, 7.5, 288.15).value,
            12.35314971, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(59, 1013.25, 7.5, 288.15).value,
            13.63529754, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(60, 1013.25, 7.5, 288.15).value,
            14.62347701, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(61, 1013.25, 7.5, 288.15).value,
            15.00716194, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(62, 1013.25, 7.5, 288.15).value,
            13.99621411, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(63, 1013.25, 7.5, 288.15).value,
            10.83108919, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(64, 1013.25, 7.5, 288.15).value,
            6.844588337, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(65, 1013.25, 7.5, 288.15).value,
            3.80880229, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(66, 1013.25, 7.5, 288.15).value,
            1.966616477, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(67, 1013.25, 7.5, 288.15).value,
            1.033387448, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(68, 1013.25, 7.5, 288.15).value,
            0.60546544, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(69, 1013.25, 7.5, 288.15).value,
            0.406984877, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(70, 1013.25, 7.5, 288.15).value,
            0.304104518, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(71, 1013.25, 7.5, 288.15).value,
            0.24160024, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(72, 1013.25, 7.5, 288.15).value,
            0.198531458, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(73, 1013.25, 7.5, 288.15).value,
            0.167045465, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(74, 1013.25, 7.5, 288.15).value,
            0.143141978, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(75, 1013.25, 7.5, 288.15).value,
            0.124484922, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(76, 1013.25, 7.5, 288.15).value,
            0.109604088, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(77, 1013.25, 7.5, 288.15).value,
            0.09752563, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(78, 1013.25, 7.5, 288.15).value,
            0.087579095, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(79, 1013.25, 7.5, 288.15).value,
            0.079288509, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(80, 1013.25, 7.5, 288.15).value,
            0.072307337, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(81, 1013.25, 7.5, 288.15).value,
            0.066377906, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(82, 1013.25, 7.5, 288.15).value,
            0.061305161, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(83, 1013.25, 7.5, 288.15).value,
            0.05693918, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(84, 1013.25, 7.5, 288.15).value,
            0.053163238, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(85, 1013.25, 7.5, 288.15).value,
            0.049885471, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(86, 1013.25, 7.5, 288.15).value,
            0.047032946, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(87, 1013.25, 7.5, 288.15).value,
            0.044547391, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(88, 1013.25, 7.5, 288.15).value,
            0.042382069, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(89, 1013.25, 7.5, 288.15).value,
            0.040499471, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(90, 1013.25, 7.5, 288.15).value,
            0.038869622, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(91, 1013.25, 7.5, 288.15).value,
            0.037468818, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(92, 1013.25, 7.5, 288.15).value,
            0.036278727, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(93, 1013.25, 7.5, 288.15).value,
            0.035285753, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(94, 1013.25, 7.5, 288.15).value,
            0.034480646, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(95, 1013.25, 7.5, 288.15).value,
            0.033858315, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(96, 1013.25, 7.5, 288.15).value,
            0.033417859, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(97, 1013.25, 7.5, 288.15).value,
            0.033162815, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(98, 1013.25, 7.5, 288.15).value,
            0.033101677, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(99, 1013.25, 7.5, 288.15).value,
            0.033248738, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(100, 1013.25, 7.5, 288.15).value,
            0.033625377, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(101, 1013.25, 7.5, 288.15).value,
            0.034261951, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(102, 1013.25, 7.5, 288.15).value,
            0.03520058, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(103, 1013.25, 7.5, 288.15).value,
            0.036499225, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(104, 1013.25, 7.5, 288.15).value,
            0.038237778, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(105, 1013.25, 7.5, 288.15).value,
            0.040527282, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(106, 1013.25, 7.5, 288.15).value,
            0.043524209, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(107, 1013.25, 7.5, 288.15).value,
            0.047453183, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(108, 1013.25, 7.5, 288.15).value,
            0.05264422, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(109, 1013.25, 7.5, 288.15).value,
            0.059596011, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(110, 1013.25, 7.5, 288.15).value,
            0.069087844, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(111, 1013.25, 7.5, 288.15).value,
            0.082387054, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(112, 1013.25, 7.5, 288.15).value,
            0.101654574, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(113, 1013.25, 7.5, 288.15).value,
            0.130786962, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(114, 1013.25, 7.5, 288.15).value,
            0.177281273, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(115, 1013.25, 7.5, 288.15).value,
            0.25660834, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(116, 1013.25, 7.5, 288.15).value,
            0.402453591, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(117, 1013.25, 7.5, 288.15).value,
            0.683016431, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(118, 1013.25, 7.5, 288.15).value,
            1.134866447, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(119, 1013.25, 7.5, 288.15).value,
            1.306379447, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(120, 1013.25, 7.5, 288.15).value,
            0.886108944, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(121, 1013.25, 7.5, 288.15).value,
            0.509171816, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(122, 1013.25, 7.5, 288.15).value,
            0.307768488, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(123, 1013.25, 7.5, 288.15).value,
            0.202100995, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(124, 1013.25, 7.5, 288.15).value,
            0.142570138, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(125, 1013.25, 7.5, 288.15).value,
            0.106445548, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(126, 1013.25, 7.5, 288.15).value,
            0.083103974, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(127, 1013.25, 7.5, 288.15).value,
            0.067232038, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(128, 1013.25, 7.5, 288.15).value,
            0.055982184, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(129, 1013.25, 7.5, 288.15).value,
            0.047732332, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(130, 1013.25, 7.5, 288.15).value,
            0.041509034, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(131, 1013.25, 7.5, 288.15).value,
            0.036701642, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(132, 1013.25, 7.5, 288.15).value,
            0.032912343, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(133, 1013.25, 7.5, 288.15).value,
            0.029873422, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(134, 1013.25, 7.5, 288.15).value,
            0.027399581, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(135, 1013.25, 7.5, 288.15).value,
            0.025359375, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(136, 1013.25, 7.5, 288.15).value,
            0.02365753, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(137, 1013.25, 7.5, 288.15).value,
            0.022223652, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(138, 1013.25, 7.5, 288.15).value,
            0.021004834, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(139, 1013.25, 7.5, 288.15).value,
            0.019960701, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(140, 1013.25, 7.5, 288.15).value,
            0.019060012, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(141, 1013.25, 7.5, 288.15).value,
            0.018278288, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(142, 1013.25, 7.5, 288.15).value,
            0.017596128, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(143, 1013.25, 7.5, 288.15).value,
            0.016997997, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(144, 1013.25, 7.5, 288.15).value,
            0.016471333, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(145, 1013.25, 7.5, 288.15).value,
            0.016005887, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(146, 1013.25, 7.5, 288.15).value,
            0.015593233, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(147, 1013.25, 7.5, 288.15).value,
            0.015226386, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(148, 1013.25, 7.5, 288.15).value,
            0.014899516, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(149, 1013.25, 7.5, 288.15).value,
            0.014607727, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(150, 1013.25, 7.5, 288.15).value,
            0.01434688, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(151, 1013.25, 7.5, 288.15).value,
            0.014113455, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(152, 1013.25, 7.5, 288.15).value,
            0.013904444, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(153, 1013.25, 7.5, 288.15).value,
            0.013717263, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(154, 1013.25, 7.5, 288.15).value,
            0.013549678, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(155, 1013.25, 7.5, 288.15).value,
            0.013399755, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(156, 1013.25, 7.5, 288.15).value,
            0.013265808, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(157, 1013.25, 7.5, 288.15).value,
            0.013146362, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(158, 1013.25, 7.5, 288.15).value,
            0.013040122, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(159, 1013.25, 7.5, 288.15).value,
            0.012945947, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(160, 1013.25, 7.5, 288.15).value,
            0.012862828, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(161, 1013.25, 7.5, 288.15).value,
            0.012789869, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(162, 1013.25, 7.5, 288.15).value,
            0.012726273, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(163, 1013.25, 7.5, 288.15).value,
            0.012671328, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(164, 1013.25, 7.5, 288.15).value,
            0.012624397, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(165, 1013.25, 7.5, 288.15).value,
            0.012584907, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(166, 1013.25, 7.5, 288.15).value,
            0.012552343, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(167, 1013.25, 7.5, 288.15).value,
            0.012526238, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(168, 1013.25, 7.5, 288.15).value,
            0.012506174, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(169, 1013.25, 7.5, 288.15).value,
            0.012491766, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(170, 1013.25, 7.5, 288.15).value,
            0.012482668, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(171, 1013.25, 7.5, 288.15).value,
            0.012478563, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(172, 1013.25, 7.5, 288.15).value,
            0.012479162, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(173, 1013.25, 7.5, 288.15).value,
            0.012484201, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(174, 1013.25, 7.5, 288.15).value,
            0.012493438, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(175, 1013.25, 7.5, 288.15).value,
            0.01250665, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(176, 1013.25, 7.5, 288.15).value,
            0.012523632, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(177, 1013.25, 7.5, 288.15).value,
            0.012544196, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(178, 1013.25, 7.5, 288.15).value,
            0.012568166, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(179, 1013.25, 7.5, 288.15).value,
            0.012595383, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(180, 1013.25, 7.5, 288.15).value,
            0.012625696, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(181, 1013.25, 7.5, 288.15).value,
            0.012658968, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(182, 1013.25, 7.5, 288.15).value,
            0.01269507, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(183, 1013.25, 7.5, 288.15).value,
            0.012733882, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(184, 1013.25, 7.5, 288.15).value,
            0.012775293, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(185, 1013.25, 7.5, 288.15).value,
            0.012819198, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(186, 1013.25, 7.5, 288.15).value,
            0.012865502, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(187, 1013.25, 7.5, 288.15).value,
            0.012914112, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(188, 1013.25, 7.5, 288.15).value,
            0.012964945, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(189, 1013.25, 7.5, 288.15).value,
            0.01301792, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(190, 1013.25, 7.5, 288.15).value,
            0.013072963, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(191, 1013.25, 7.5, 288.15).value,
            0.013130004, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(192, 1013.25, 7.5, 288.15).value,
            0.013188976, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(193, 1013.25, 7.5, 288.15).value,
            0.013249818, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(194, 1013.25, 7.5, 288.15).value,
            0.013312471, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(195, 1013.25, 7.5, 288.15).value,
            0.01337688, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(196, 1013.25, 7.5, 288.15).value,
            0.013442993, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(197, 1013.25, 7.5, 288.15).value,
            0.01351076, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(198, 1013.25, 7.5, 288.15).value,
            0.013580135, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(199, 1013.25, 7.5, 288.15).value,
            0.013651074, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(200, 1013.25, 7.5, 288.15).value,
            0.013723536, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(201, 1013.25, 7.5, 288.15).value,
            0.01379748, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(202, 1013.25, 7.5, 288.15).value,
            0.013872869, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(203, 1013.25, 7.5, 288.15).value,
            0.013949668, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(204, 1013.25, 7.5, 288.15).value,
            0.014027843, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(205, 1013.25, 7.5, 288.15).value,
            0.014107361, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(206, 1013.25, 7.5, 288.15).value,
            0.014188192, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(207, 1013.25, 7.5, 288.15).value,
            0.014270307, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(208, 1013.25, 7.5, 288.15).value,
            0.014353678, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(209, 1013.25, 7.5, 288.15).value,
            0.014438278, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(210, 1013.25, 7.5, 288.15).value,
            0.014524083, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(211, 1013.25, 7.5, 288.15).value,
            0.014611069, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(212, 1013.25, 7.5, 288.15).value,
            0.014699212, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(213, 1013.25, 7.5, 288.15).value,
            0.01478849, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(214, 1013.25, 7.5, 288.15).value,
            0.014878883, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(215, 1013.25, 7.5, 288.15).value,
            0.01497037, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(216, 1013.25, 7.5, 288.15).value,
            0.015062932, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(217, 1013.25, 7.5, 288.15).value,
            0.015156551, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(218, 1013.25, 7.5, 288.15).value,
            0.01525121, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(219, 1013.25, 7.5, 288.15).value,
            0.015346891, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(220, 1013.25, 7.5, 288.15).value,
            0.015443579, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(221, 1013.25, 7.5, 288.15).value,
            0.015541258, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(222, 1013.25, 7.5, 288.15).value,
            0.015639912, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(223, 1013.25, 7.5, 288.15).value,
            0.015739529, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(224, 1013.25, 7.5, 288.15).value,
            0.015840094, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(225, 1013.25, 7.5, 288.15).value,
            0.015941595, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(226, 1013.25, 7.5, 288.15).value,
            0.016044018, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(227, 1013.25, 7.5, 288.15).value,
            0.016147352, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(228, 1013.25, 7.5, 288.15).value,
            0.016251585, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(229, 1013.25, 7.5, 288.15).value,
            0.016356706, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(230, 1013.25, 7.5, 288.15).value,
            0.016462705, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(231, 1013.25, 7.5, 288.15).value,
            0.016569571, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(232, 1013.25, 7.5, 288.15).value,
            0.016677295, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(233, 1013.25, 7.5, 288.15).value,
            0.016785866, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(234, 1013.25, 7.5, 288.15).value,
            0.016895277, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(235, 1013.25, 7.5, 288.15).value,
            0.017005518, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(236, 1013.25, 7.5, 288.15).value,
            0.017116581, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(237, 1013.25, 7.5, 288.15).value,
            0.017228459, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(238, 1013.25, 7.5, 288.15).value,
            0.017341142, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(239, 1013.25, 7.5, 288.15).value,
            0.017454625, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(240, 1013.25, 7.5, 288.15).value,
            0.0175689, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(241, 1013.25, 7.5, 288.15).value,
            0.01768396, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(242, 1013.25, 7.5, 288.15).value,
            0.017799799, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(243, 1013.25, 7.5, 288.15).value,
            0.017916411, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(244, 1013.25, 7.5, 288.15).value,
            0.018033789, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(245, 1013.25, 7.5, 288.15).value,
            0.018151929, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(246, 1013.25, 7.5, 288.15).value,
            0.018270824, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(247, 1013.25, 7.5, 288.15).value,
            0.01839047, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(248, 1013.25, 7.5, 288.15).value,
            0.018510862, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(249, 1013.25, 7.5, 288.15).value,
            0.018631995, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(250, 1013.25, 7.5, 288.15).value,
            0.018753865, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(251, 1013.25, 7.5, 288.15).value,
            0.018876467, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(252, 1013.25, 7.5, 288.15).value,
            0.018999798, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(253, 1013.25, 7.5, 288.15).value,
            0.019123854, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(254, 1013.25, 7.5, 288.15).value,
            0.019248631, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(255, 1013.25, 7.5, 288.15).value,
            0.019374127, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(256, 1013.25, 7.5, 288.15).value,
            0.019500338, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(257, 1013.25, 7.5, 288.15).value,
            0.019627261, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(258, 1013.25, 7.5, 288.15).value,
            0.019754894, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(259, 1013.25, 7.5, 288.15).value,
            0.019883235, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(260, 1013.25, 7.5, 288.15).value,
            0.02001228, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(261, 1013.25, 7.5, 288.15).value,
            0.020142029, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(262, 1013.25, 7.5, 288.15).value,
            0.02027248, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(263, 1013.25, 7.5, 288.15).value,
            0.020403631, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(264, 1013.25, 7.5, 288.15).value,
            0.020535481, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(265, 1013.25, 7.5, 288.15).value,
            0.020668028, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(266, 1013.25, 7.5, 288.15).value,
            0.020801273, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(267, 1013.25, 7.5, 288.15).value,
            0.020935214, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(268, 1013.25, 7.5, 288.15).value,
            0.021069851, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(269, 1013.25, 7.5, 288.15).value,
            0.021205184, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(270, 1013.25, 7.5, 288.15).value,
            0.021341213, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(271, 1013.25, 7.5, 288.15).value,
            0.021477939, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(272, 1013.25, 7.5, 288.15).value,
            0.021615362, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(273, 1013.25, 7.5, 288.15).value,
            0.021753483, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(274, 1013.25, 7.5, 288.15).value,
            0.021892304, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(275, 1013.25, 7.5, 288.15).value,
            0.022031825, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(276, 1013.25, 7.5, 288.15).value,
            0.02217205, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(277, 1013.25, 7.5, 288.15).value,
            0.022312979, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(278, 1013.25, 7.5, 288.15).value,
            0.022454615, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(279, 1013.25, 7.5, 288.15).value,
            0.022596961, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(280, 1013.25, 7.5, 288.15).value,
            0.022740019, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(281, 1013.25, 7.5, 288.15).value,
            0.022883795, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(282, 1013.25, 7.5, 288.15).value,
            0.02302829, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(283, 1013.25, 7.5, 288.15).value,
            0.02317351, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(284, 1013.25, 7.5, 288.15).value,
            0.023319459, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(285, 1013.25, 7.5, 288.15).value,
            0.023466142, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(286, 1013.25, 7.5, 288.15).value,
            0.023613565, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(287, 1013.25, 7.5, 288.15).value,
            0.023761733, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(288, 1013.25, 7.5, 288.15).value,
            0.023910653, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(289, 1013.25, 7.5, 288.15).value,
            0.024060332, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(290, 1013.25, 7.5, 288.15).value,
            0.024210778, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(291, 1013.25, 7.5, 288.15).value,
            0.024361999, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(292, 1013.25, 7.5, 288.15).value,
            0.024514003, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(293, 1013.25, 7.5, 288.15).value,
            0.024666801, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(294, 1013.25, 7.5, 288.15).value,
            0.024820402, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(295, 1013.25, 7.5, 288.15).value,
            0.024974817, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(296, 1013.25, 7.5, 288.15).value,
            0.025130058, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(297, 1013.25, 7.5, 288.15).value,
            0.025286138, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(298, 1013.25, 7.5, 288.15).value,
            0.025443071, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(299, 1013.25, 7.5, 288.15).value,
            0.025600871, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(300, 1013.25, 7.5, 288.15).value,
            0.025759555, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(301, 1013.25, 7.5, 288.15).value,
            0.025919138, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(302, 1013.25, 7.5, 288.15).value,
            0.026079639, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(303, 1013.25, 7.5, 288.15).value,
            0.026241079, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(304, 1013.25, 7.5, 288.15).value,
            0.026403477, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(305, 1013.25, 7.5, 288.15).value,
            0.026566857, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(306, 1013.25, 7.5, 288.15).value,
            0.026731244, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(307, 1013.25, 7.5, 288.15).value,
            0.026896663, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(308, 1013.25, 7.5, 288.15).value,
            0.027063143, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(309, 1013.25, 7.5, 288.15).value,
            0.027230715, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(310, 1013.25, 7.5, 288.15).value,
            0.027399412, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(311, 1013.25, 7.5, 288.15).value,
            0.02756927, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(312, 1013.25, 7.5, 288.15).value,
            0.027740328, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(313, 1013.25, 7.5, 288.15).value,
            0.027912629, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(314, 1013.25, 7.5, 288.15).value,
            0.028086218, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(315, 1013.25, 7.5, 288.15).value,
            0.028261145, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(316, 1013.25, 7.5, 288.15).value,
            0.028437464, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(317, 1013.25, 7.5, 288.15).value,
            0.028615235, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(318, 1013.25, 7.5, 288.15).value,
            0.028794523, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(319, 1013.25, 7.5, 288.15).value,
            0.028975399, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(320, 1013.25, 7.5, 288.15).value,
            0.029157939, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(321, 1013.25, 7.5, 288.15).value,
            0.029342231, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(322, 1013.25, 7.5, 288.15).value,
            0.029528367, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(323, 1013.25, 7.5, 288.15).value,
            0.029716451, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(324, 1013.25, 7.5, 288.15).value,
            0.029906599, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(325, 1013.25, 7.5, 288.15).value,
            0.030098937, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(326, 1013.25, 7.5, 288.15).value,
            0.030293607, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(327, 1013.25, 7.5, 288.15).value,
            0.030490765, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(328, 1013.25, 7.5, 288.15).value,
            0.030690588, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(329, 1013.25, 7.5, 288.15).value,
            0.030893273, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(330, 1013.25, 7.5, 288.15).value,
            0.031099041, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(331, 1013.25, 7.5, 288.15).value,
            0.031308143, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(332, 1013.25, 7.5, 288.15).value,
            0.031520859, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(333, 1013.25, 7.5, 288.15).value,
            0.031737512, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(334, 1013.25, 7.5, 288.15).value,
            0.031958467, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(335, 1013.25, 7.5, 288.15).value,
            0.032184144, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(336, 1013.25, 7.5, 288.15).value,
            0.032415022, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(337, 1013.25, 7.5, 288.15).value,
            0.032651659, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(338, 1013.25, 7.5, 288.15).value,
            0.032894698, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(339, 1013.25, 7.5, 288.15).value,
            0.033144893, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(340, 1013.25, 7.5, 288.15).value,
            0.033403124, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(341, 1013.25, 7.5, 288.15).value,
            0.033670433, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(342, 1013.25, 7.5, 288.15).value,
            0.033948053, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(343, 1013.25, 7.5, 288.15).value,
            0.034237459, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(344, 1013.25, 7.5, 288.15).value,
            0.034540422, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(345, 1013.25, 7.5, 288.15).value,
            0.034859093, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(346, 1013.25, 7.5, 288.15).value,
            0.035196097, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(347, 1013.25, 7.5, 288.15).value,
            0.03555467, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(348, 1013.25, 7.5, 288.15).value,
            0.03593884, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(349, 1013.25, 7.5, 288.15).value,
            0.036353672, places=5)

        self.assertAlmostEqual(
            models.itu676.gamma0_approx(350, 1013.25, 7.5, 288.15).value,
            0.036805605, places=5)

    def test_zenit_water_vapour_attenuation(self):
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                51.5, -0.14, 1.0, 14.25).value,
            0.064981043, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                41.9, 12.49, 1.0, 14.25).value,
            0.070360091, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                33.94, 18.43, 1.0, 14.25).value,
            0.074660262, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                51.5, -0.14, 0.5, 14.25).value,
            0.06911297, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                41.9, 12.49, 0.5, 14.25).value,
            0.073434531, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                33.94, 18.43, 0.5, 14.25).value,
            0.080098077, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                51.5, -0.14, 0.3, 14.25).value,
            0.072394726, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                41.9, 12.49, 0.3, 14.25).value,
            0.075162715, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                33.94, 18.43, 0.3, 14.25).value,
            0.083750389, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                51.5, -0.14, 0.2, 14.25).value,
            0.074394064, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                41.9, 12.49, 0.2, 14.25).value,
            0.076695287, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                33.94, 18.43, 0.2, 14.25).value,
            0.086350752, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                51.5, -0.14, 1.0, 29).value,
            0.305636526, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                41.9, 12.49, 1.0, 29).value,
            0.331425898, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                33.94, 18.43, 1.0, 29).value,
            0.355205229, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                51.5, -0.14, 0.5, 29).value,
            0.324977228, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                41.9, 12.49, 0.5, 29).value,
            0.345830132, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                33.94, 18.43, 0.5, 29).value,
            0.38091961, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                51.5, -0.14, 0.3, 29).value,
            0.340327583, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                41.9, 12.49, 0.3, 29).value,
            0.353923317, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                33.94, 18.43, 0.3, 29).value,
            0.398176611, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                51.5, -0.14, 0.2, 29).value,
            0.349674822, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                41.9, 12.49, 0.2, 29).value,
            0.361098289, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                33.94, 18.43, 0.2, 29).value,
            0.410456469, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                22.9, -43.23, 1.0, 14.25).value,
            0.099820608, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                25.78, -80.22, 1.0, 14.25).value,
            0.118484695, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                22.9, -43.23, 0.5, 14.25).value,
            0.105446054, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                25.78, -80.22, 0.5, 14.25).value,
            0.12252307, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                22.9, -43.23, 0.3, 14.25).value,
            0.108812058, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                25.78, -80.22, 0.3, 14.25).value,
            0.125093339, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                22.9, -43.23, 0.2, 14.25).value,
            0.111441086, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                25.78, -80.22, 0.2, 14.25).value,
            0.127090376, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                22.9, -43.23, 1.0, 29).value,
            0.473979935, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                25.78, -80.22, 1.0, 29).value,
            0.561753331, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                22.9, -43.23, 0.5, 29).value,
            0.500468518, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                25.78, -80.22, 0.5, 29).value,
            0.580717641, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                22.9, -43.23, 0.3, 29).value,
            0.516307047, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                25.78, -80.22, 0.3, 29).value,
            0.592782098, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                22.9, -43.23, 0.2, 29).value,
            0.528672179, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                25.78, -80.22, 0.2, 29).value,
            0.602152942, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                28.717, 77.3, 1.0, 14.25).value,
            0.149156898, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                3.133, 101.7, 1.0, 14.25).value,
            0.121165007, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                9.05, 38.7, 1.0, 14.25).value,
            0.051589359, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                28.717, 77.3, 0.5, 14.25).value,
            0.153859398, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                3.133, 101.7, 0.5, 14.25).value,
            0.123550552, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                9.05, 38.7, 0.5, 14.25).value,
            0.052996133, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                28.717, 77.3, 0.3, 14.25).value,
            0.156616572, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                3.133, 101.7, 0.3, 14.25).value,
            0.125325192, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                9.05, 38.7, 0.3, 14.25).value,
            0.053871006, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                28.717, 77.3, 0.2, 14.25).value,
            0.158958354, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                3.133, 101.7, 0.2, 14.25).value,
            0.126766365, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                9.05, 38.7, 0.2, 14.25).value,
            0.054721343, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                28.717, 77.3, 1.0, 29).value,
            0.683528163, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                3.133, 101.7, 1.0, 29).value,
            0.555168022, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                9.05, 38.7, 1.0, 29).value,
            0.188559832, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                28.717, 77.3, 0.5, 29).value,
            0.704836196, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                3.133, 101.7, 0.5, 29).value,
            0.565993797, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                9.05, 38.7, 0.5, 29).value,
            0.193687836, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                28.717, 77.3, 0.3, 29).value,
            0.717323975, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                3.133, 101.7, 0.3, 29).value,
            0.574044911, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                9.05, 38.7, 0.3, 29).value,
            0.19687619, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                28.717, 77.3, 0.2, 29).value,
            0.727927181, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                3.133, 101.7, 0.2, 29).value,
            0.580581723, places=5)
        self.assertAlmostEqual(
            models.itu676.zenit_water_vapour_attenuation(
                9.05, 38.7, 0.2, 29).value,
            0.19997458, places=5)


class ITUR836_6TestCase(test.TestCase):

    def setUp(self):
        models.itu836.change_version(6)

    def test_surface_water_vapour_density(self):
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                3.133, 101.7, 0.1, 0.236104459).value,
            22.93756598, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                3.133, 101.7, 0.15, 0.236104459).value,
            22.80534575, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                3.133, 101.7, 0.3, 0.236104459).value,
            22.55507955, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                3.133, 101.7, 0.35, 0.236104459).value,
            22.49361957, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                22.9, -43.23, 0.1, 0).value,
            21.59164912, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                22.9, -43.23, 0.15, 0).value,
            21.46164369, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                22.9, -43.23, 0.3, 0).value,
            21.24753319, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                22.9, -43.23, 0.35, 0).value,
            21.18676013, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                23, 30, 0.1, 0.247).value,
            11.88170822, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                23, 30, 0.15, 0.247).value,
            11.61777268, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                23, 30, 0.3, 0.247).value,
            11.12235912, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                23, 30, 0.35, 0.247).value,
            11.00877052, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                25.78, -80.22, 0.1, 7.51071e-05).value,
            23.50748104, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                25.78, -80.22, 0.15, 7.51071e-05).value,
            23.34324475, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                25.78, -80.22, 0.3, 7.51071e-05).value,
            23.06574222, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                25.78, -80.22, 0.35, 7.51071e-05).value,
            23.00327243, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                28.717, 77.3, 0.1, 0.217559455).value,
            25.95287453, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                28.717, 77.3, 0.15, 0.217559455).value,
            25.71217873, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                28.717, 77.3, 0.3, 0.217559455).value,
            25.34018758, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                28.717, 77.3, 0.35, 0.217559455).value,
            25.2557054, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                33.94, 18.43, 0.1, 0).value,
            24.00156532, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                33.94, 18.43, 0.15, 0).value,
            23.85987554, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                33.94, 18.43, 0.3, 0).value,
            23.51464505, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                33.94, 18.43, 0.35, 0).value,
            23.41954477, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                41.9, 12.49, 0.1, 0.056701045).value,
            19.78501126, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                41.9, 12.49, 0.15, 0.056701045).value,
            19.48948848, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                41.9, 12.49, 0.3, 0.056701045).value,
            19.02450953, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                41.9, 12.49, 0.35, 0.056701045).value,
            18.92055161, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                51.5, -0.14, 0.1, 0.069164224).value,
            15.21351315, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                51.5, -0.14, 0.15, 0.069164224).value,
            15.0172773, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                51.5, -0.14, 0.3, 0.069164224).value,
            14.6189506, places=5)
        self.assertAlmostEqual(
            models.itu836.surface_water_vapour_density(
                51.5, -0.14, 0.35, 0.069164224).value,
            14.50640729, places=5)

    def test_total_water_vapour_content(self):
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                3.133, 101.7, 0.1, 0.23610446).value,
            62.16532093, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                3.133, 101.7, 0.15, 0.23610446).value,
            61.59527521, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                3.133, 101.7, 0.3, 0.23610446).value,
            60.58285243, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                3.133, 101.7, 0.35, 0.23610446).value,
            60.35619302, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                22.9, -43.23, 0.1, 0.0).value,
            56.38788554, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                22.9, -43.23, 0.15, 0.0).value,
            55.36064664, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                22.9, -43.23, 0.3, 0.0).value,
            53.4851113, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                22.9, -43.23, 0.35, 0.0).value,
            53.03918259, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                23, 30, 0.1, 0.247).value,
            38.47288189, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                23, 30, 0.15, 0.247).value,
            37.21449337, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                23, 30, 0.3, 0.247).value,
            34.63093178, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                23, 30, 0.35, 0.247).value,
            34.06569649, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                25.78, -80.22, 0.1, 7.511e-05).value,
            62.84315177, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                25.78, -80.22, 0.15, 7.511e-05).value,
            61.95641322, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                25.78, -80.22, 0.3, 7.511e-05).value,
            60.48487688, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                25.78, -80.22, 0.35, 7.511e-05).value,
            60.1561742, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                28.717, 77.3, 0.1, 0.21755946).value,
            75.44891006, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                28.717, 77.3, 0.15, 0.21755946).value,
            74.79639702, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                28.717, 77.3, 0.3, 0.21755946).value,
            73.40408393, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                28.717, 77.3, 0.35, 0.21755946).value,
            73.07234727, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                33.94, 18.43, 0.1, 0.0).value,
            45.19895208, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                33.94, 18.43, 0.15, 0.0).value,
            44.15275162, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                33.94, 18.43, 0.3, 0.0).value,
            42.21022387, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                33.94, 18.43, 0.35, 0.0).value,
            41.69772633, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                41.9, 12.49, 0.1, 0.05670104).value,
            39.93693588, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                41.9, 12.49, 0.15, 0.05670104).value,
            39.33984158, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                41.9, 12.49, 0.3, 0.05670104).value,
            38.19321515, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                41.9, 12.49, 0.35, 0.05670104).value,
            37.94621912, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                51.5, -0.14, 0.1, 0.06916422).value,
            39.23803432, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                51.5, -0.14, 0.15, 0.06916422).value,
            38.41414987, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                51.5, -0.14, 0.3, 0.06916422).value,
            36.88058222, places=5)
        self.assertAlmostEqual(
            models.itu836.total_water_vapour_content(
                51.5, -0.14, 0.35, 0.06916422).value,
            36.4074561, places=5)


class ITUR838_3TestCase(test.TestCase):

    def setUp(self):
        models.itu838.change_version(3)

    def test_rain_specific_attenuation(self):
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 14.25, 30.87067768, 0).value,
            1.879742, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 14.25, 40.97052773, 0).value,
            3.630988, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 14.25, 47.91280491, 0).value,
            3.503189, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 14.25, 30.87067768, 0).value,
            1.879742, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 14.25, 40.97052773, 0).value,
            3.630988, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 14.25, 47.91280491, 0).value,
            3.503189, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 14.25, 30.87067768, 0).value,
            1.879742, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 14.25, 40.97052773, 0).value,
            3.630988, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 14.25, 47.91280491, 0).value,
            3.503189, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 14.25, 30.87067768, 0).value,
            1.879742, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 14.25, 40.97052773, 0).value,
            3.630988, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 14.25, 47.91280491, 0).value,
            3.503189, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 29.00, 30.87067768, 0).value,
            5.814832, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 29.00, 40.97052773, 0).value,
            10.157375, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 29.00, 47.91280491, 0).value,
            9.846762, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 29.00, 30.87067768, 0).value,
            5.814832, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 29.00, 40.97052773, 0).value,
            10.157375, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 29.00, 47.91280491, 0).value,
            9.846762, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 29.00, 30.87067768, 0).value,
            5.814832, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 29.00, 40.97052773, 0).value,
            10.157375, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 29.00, 47.91280491, 0).value,
            9.846762, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                30.875024, 29.00, 30.87067768, 0).value,
            5.814832, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                56.370009, 29.00, 40.97052773, 0).value,
            10.157375, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                55.231625, 29.00, 47.91280491, 0).value,
            9.846762, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 14.25, 59.81487174, 0).value,
            3.628282, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 14.25, 49.20900369, 0).value,
            5.948478, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 14.25, 59.81487174, 0).value,
            3.628282, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 14.25, 49.20900369, 0).value,
            5.948478, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 14.25, 59.81487174, 0).value,
            3.628282, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 14.25, 49.20900369, 0).value,
            5.948478, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 14.25, 59.81487174, 0).value,
            3.628282, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 14.25, 49.20900369, 0).value,
            5.948478, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 29.00, 59.81487174, 0).value,
            10.132682, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 29.00, 49.20900369, 0).value,
            15.460212, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 29.00, 59.81487174, 0).value,
            10.132682, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 29.00, 49.20900369, 0).value,
            15.460212, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 29.00, 59.81487174, 0).value,
            10.132682, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 29.00, 49.20900369, 0).value,
            15.460212, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                58.094216, 29.00, 59.81487174, 0).value,
            10.132682, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                89.114103, 29.00, 49.20900369, 0).value,
            15.460212, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 14.25, 55.90591362, 0).value,
            3.603569, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 14.25, 67.76751981, 0).value,
            6.06336, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 14.25, 38.14104832, 0).value,
            3.523996, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 14.25, 55.90591362, 0).value,
            3.603569, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 14.25, 67.76751981, 0).value,
            6.06336, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 14.25, 38.14104832, 0).value,
            3.523996, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 14.25, 55.90591362, 0).value,
            3.603569, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 14.25, 67.76751981, 0).value,
            6.06336, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 14.25, 38.14104832, 0).value,
            3.523996, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 14.25, 55.90591362, 0).value,
            3.603569, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 14.25, 67.76751981, 0).value,
            6.06336, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 14.25, 38.14104832, 0).value,
            3.523996, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 29.00, 55.90591362, 0).value,
            10.078266, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 29.00, 67.76751981, 0).value,
            15.712442, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 29.00, 38.14104832, 0).value,
            9.904098, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 29.00, 55.90591362, 0).value,
            10.078266, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 29.00, 67.76751981, 0).value,
            15.712442, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 29.00, 38.14104832, 0).value,
            9.904098, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 29.00, 55.90591362, 0).value,
            10.078266, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 29.00, 67.76751981, 0).value,
            15.712442, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 29.00, 38.14104832, 0).value,
            9.904098, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                57.39623, 29.00, 55.90591362, 0).value,
            10.078266, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                93.607098, 29.00, 67.76751981, 0).value,
            15.712442, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                54.623411, 29.00, 38.14104832, 0).value,
            9.904098, places=5)
        # New values in validation 4
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                26.48052, 14.25, 31.07694309, 0).value,
            1.581308489, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                33.936232, 14.25, 40.23202374, 0).value,
            2.06173217, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                27.13586832, 14.25, 46.35969261, 0).value,
            1.592084199, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                26.48052, 14.25, 31.07694309, 0).value,
            1.581308489, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                33.936232, 14.25, 40.23202374, 0).value,
            2.06173217, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                27.13586832, 14.25, 46.35969261, 0).value,
            1.592084199, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                26.48052, 14.25, 31.07694309, 0).value,
            1.581308489, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                33.936232, 14.25, 40.23202374, 0).value,
            2.06173217, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                27.13586832, 14.25, 46.35969261, 0).value,
            1.592084199, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                26.48052, 14.25, 31.07694309, 0).value,
            1.581308489, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                33.936232, 14.25, 40.23202374, 0).value,
            2.06173217, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                27.13586832, 14.25, 46.35969261, 0).value,
            1.592084199, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                26.48052, 29, 31.07694309, 0).value,
            5.021802196, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                33.936232, 29, 40.23202374, 0).value,
            6.278460355, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                27.13586832, 29, 46.35969261, 0).value,
            5.031354793, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                26.48052, 29, 31.07694309, 0).value,
            5.021802196, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                33.936232, 29, 40.23202374, 0).value,
            6.278460355, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                27.13586832, 29, 46.35969261, 0).value,
            5.031354793, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                26.48052, 29, 31.07694309, 0).value,
            5.021802196, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                33.936232, 29, 40.23202374, 0).value,
            6.278460355, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                27.13586832, 29, 46.35969261, 0).value,
            5.031354793, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                26.48052, 29, 31.07694309, 0).value,
            5.021802196, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                33.936232, 29, 40.23202374, 0).value,
            6.278460355, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                27.13586832, 29, 46.35969261, 0).value,
            5.031354793, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                50.639304, 14.25, 22.27833468, 0).value,
            3.321396378, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                78.2994993, 14.25, 52.6789929, 0).value,
            5.11503455, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                50.639304, 14.25, 22.27833468, 0).value,
            3.321396378, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                78.2994993, 14.25, 52.6789929, 0).value,
            5.11503455, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                50.639304, 14.25, 22.27833468, 0).value,
            3.321396378, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                78.2994993, 14.25, 52.6789929, 0).value,
            5.11503455, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                50.639304, 14.25, 22.27833468, 0).value,
            3.321396378, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                78.2994993, 14.25, 52.6789929, 0).value,
            5.11503455, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                50.639304, 29, 22.27833468, 0).value,
            9.424302438, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                78.2994993, 29, 52.6789929, 0).value,
            13.59290067, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                50.639304, 29, 22.27833468, 0).value,
            9.424302438, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                78.2994993, 29, 52.6789929, 0).value,
            13.59290067, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                50.639304, 29, 22.27833468, 0).value,
            9.424302438, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                78.2994993, 29, 52.6789929, 0).value,
            13.59290067, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                50.639304, 29, 22.27833468, 0).value,
            9.424302438, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                78.2994993, 29, 52.6789929, 0).value,
            13.59290067, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                63.62668149, 14.25, 48.23861222, 90).value,
            3.72899602, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                99.13558978, 14.25, 85.80767474, 90).value,
            6.340652096, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                42.91007183, 14.25, 20.14348033, 90).value,
            2.350323497, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                63.62668149, 14.25, 48.23861222, 90).value,
            3.72899602, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                99.13558978, 14.25, 85.80767474, 90).value,
            6.340652096, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                42.91007183, 14.25, 20.14348033, 90).value,
            2.350323497, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                63.62668149, 14.25, 48.23861222, 90).value,
            3.72899602, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                99.13558978, 14.25, 85.80767474, 90).value,
            6.340652096, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                42.91007183, 14.25, 20.14348033, 90).value,
            2.350323497, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                63.62668149, 14.25, 48.23861222, 90).value,
            3.72899602, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                99.13558978, 14.25, 85.80767474, 90).value,
            6.340652096, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                42.91007183, 14.25, 20.14348033, 90).value,
            2.350323497, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                63.62668149, 29, 48.23861222, 90).value,
            10.28694456, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                99.13558978, 29, 85.80767474, 90).value,
            16.31838263, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                42.91007183, 29, 20.14348033, 90).value,
            6.833646475, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                63.62668149, 29, 48.23861222, 90).value,
            10.28694456, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                99.13558978, 29, 85.80767474, 90).value,
            16.31838263, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                42.91007183, 29, 20.14348033, 90).value,
            6.833646475, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                63.62668149, 29, 48.23861222, 90).value,
            10.28694456, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                99.13558978, 29, 85.80767474, 90).value,
            16.31838263, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                42.91007183, 29, 20.14348033, 90).value,
            6.833646475, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                63.62668149, 29, 48.23861222, 90).value,
            10.28694456, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                99.13558978, 29, 85.80767474, 90).value,
            16.31838263, places=5)
        self.assertAlmostEqual(
            models.itu838.rain_specific_attenuation(
                42.91007183, 29, 20.14348033, 90).value,
            6.833646475, places=5)


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


class ITUR837_7TestCase(test.TestCase):

    def setUp(self):
        models.itu837.change_version(7)

    def test_rainfall_rate(self):
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(3.133, 101.7, 0.1).value,
            34.64798123, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(3.133, 101.7, 0.15).value,
            27.7636201, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(3.133, 101.7, 0.3).value,
            18.26254364, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(3.133, 101.7, 0.35).value,
            16.49493229, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(22.9, -43.23, 0.1).value,
            14.58963041, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(22.9, -43.23, 0.15).value,
            11.00510082, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(22.9, -43.23, 0.3).value,
            6.23796236, places=2)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(22.9, -43.23, 0.35).value,
            5.38239642, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(23, 30, 0.1).value,
            0.0, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(23, 30, 0.15).value,
            0.0, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(23, 30, 0.3).value,
            0.0, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(23, 30, 0.35).value,
            0.0, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(25.78, -80.22, 0.1).value,
            25.33888119, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(25.78, -80.22, 0.15).value,
            19.86683577, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(25.78, -80.22, 0.3).value,
            12.43676554, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(25.78, -80.22, 0.35).value,
            11.07566126, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(28.717, 77.3, 0.1).value,
            16.53857378, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(28.717, 77.3, 0.15).value,
            12.04651363, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(28.717, 77.3, 0.3).value,
            6.21600589, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(28.717, 77.3, 0.35).value,
            5.19609765, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(33.94, 18.43, 0.1).value,
            7.43193175, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(33.94, 18.43, 0.15).value,
            5.53031864, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(33.94, 18.43, 0.3).value,
            3.03506603, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(33.94, 18.43, 0.35).value,
            2.59276061, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(41.9, 12.49, 0.1).value,
            11.19798305, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(41.9, 12.49, 0.15).value,
            8.88472572, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(41.9, 12.49, 0.3).value,
            5.75356253, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(41.9, 12.49, 0.35).value,
            5.18058827, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(51.5, -0.14, 0.1).value,
            8.9924712, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(51.5, -0.14, 0.15).value,
            7.17369312, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(51.5, -0.14, 0.3).value,
            4.69033625, places=3)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(51.5, -0.14, 0.35).value,
            4.23258601, places=3)

    def test_rainfall_probability(self):

        self.assertAlmostEqual(
            models.itu837.rainfall_probability(3.133, 101.7).value,
            4.53654368, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_probability(22.9, -43.23).value,
            1.41773353, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_probability(23, 30).value,
            0.00051911, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_probability(25.78, -80.22).value,
            2.90785192, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_probability(28.717, 77.3).value,
            1.07089363, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_probability(33.94, 18.43).value,
            1.27567391, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_probability(41.9, 12.49).value,
            5.26971907, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_probability(51.5, -0.14).value,
            5.3615096, places=5)

    def test_rainfall_rate_R001(self):

        self.assertAlmostEqual(
            models.itu837.rainfall_rate(3.133, 101.7, 0.01).value,
            99.1481136, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(22.9, -43.23, 0.01).value,
            50.639304, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(23.0, 30.0, 0.01).value,
            0.0, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(25.78, -80.22, 0.01).value,
            78.2982928, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(28.717, 77.3, 0.01).value,
            63.5972464, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(33.94, 18.43, 0.01).value,
            27.1349664, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(41.9, 12.49, 0.01).value,
            33.936232, places=5)
        self.assertAlmostEqual(
            models.itu837.rainfall_rate(51.5, -0.14, 0.01).value,
            26.48052, places=5)


class ITUR839_4TestCase(test.TestCase):

    def setUp(self):
        models.itu839.change_version(4)

    def test_isoterm_0_deg(self):
        self.assertAlmostEqual(
            models.itu839.isoterm_0(3.133, 101.7).value,
            4.5979744, places=5)
        self.assertAlmostEqual(
            models.itu839.isoterm_0(22.9, -43.23).value,
            3.79877867, places=5)
        self.assertAlmostEqual(
            models.itu839.isoterm_0(23, 30).value,
            4.168, places=5)
        self.assertAlmostEqual(
            models.itu839.isoterm_0(25.78, -80.22).value,
            4.20946133, places=5)
        self.assertAlmostEqual(
            models.itu839.isoterm_0(28.717, 77.3).value,
            4.89820404, places=5)
        self.assertAlmostEqual(
            models.itu839.isoterm_0(33.94, 18.43).value,
            2.20330276, places=5)
        self.assertAlmostEqual(
            models.itu839.isoterm_0(41.9, 12.49).value,
            2.68749333, places=5)
        self.assertAlmostEqual(
            models.itu839.isoterm_0(51.5, -0.14).value,
            2.09273333, places=5)

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
        self.assertAlmostEqual(
            models.itu839.rain_height(3.133, 101.7).value,
            4.9579744, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(22.9, -43.23).value,
            4.15877867, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(23, 30).value,
            4.528, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(25.78, -80.22).value,
            4.56946133, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(28.717, 77.3).value,
            5.25820404, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(33.94, 18.43).value,
            2.56330276, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(41.9, 12.49).value,
            3.04749333, places=5)
        self.assertAlmostEqual(
            models.itu839.rain_height(51.5, -0.14).value,
            2.45273333, places=5)


class ITUR618_12TestCase(test.TestCase):

    def setUp(self):
        models.itu618.change_version(12)
        models.itu453.change_version(12)
        models.itu838.change_version(3)
        models.itu836.change_version(5)
        models.itu837.change_version(6)
        models.itu840.change_version(6)

    def test_rain_cross_polarization_discrimination(self):
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                16.38308757, 14.25, 30.870677680, 0.001, 0).value,
            27.143007980, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                3.89479806, 14.25, 40.970527730, 0.1, 0).value,
            37.386086000, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                9.71179484, 14.25, 47.912804910, 0.01, 0).value,
            33.812795580, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                71.44613350, 29, 40.970527730, 0.001, 0).value,
            21.244470560, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                12.87478397, 29, 47.912804910, 0.1, 0).value,
            35.166125690, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                39.07323323, 29, 40.970527730, 0.01, 0).value,
            25.180145740, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                23.00197384, 14.25, 59.814871740, 0.001, 0).value,
            33.308530550, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                32.74150676, 14.25, 49.209003690, 0.001, 0).value,
            25.508227320, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                4.92489694, 14.25, 59.814871740, 0.1, 0).value,
            41.798127850, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                6.96606559, 14.25, 49.209003690, 0.1, 0).value,
            34.830206060, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                12.76053997, 14.25, 59.814871740, 0.01, 0).value,
            36.168649690, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                18.06938866, 14.25, 49.209003690, 0.01, 0).value,
            28.803871260, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                23.00197384, 14.25, 59.814871740, 0.001, 0).value,
            33.308530550, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                32.74150676, 14.25, 49.209003690, 0.001, 0).value,
            25.508227320, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                100.96022257, 29, 49.209003690, 0.001, 0).value,
            20.365001500, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                20.43214239, 29, 59.814871740, 0.1, 0).value,
            35.581135690, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                27.86774318, 29, 49.209003690, 0.1, 0).value,
            28.745547830, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                46.32024457, 29, 59.814871740, 0.01, 0).value,
            30.303830010, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                63.46384760, 29, 49.209003690, 0.01, 0).value,
            23.046241580, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                73.05533363, 29, 59.814871740, 0.001, 0).value,
            28.089155910, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                26.85402570, 14.25, 55.905913620, 0.001, 0).value,
            29.993601830, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                4.44533923, 14.25, 38.141048320, 0.1, 0).value,
            35.652315760, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                11.06265445, 14.25, 38.141048320, 0.01, 0).value,
            30.034285750, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                26.85402570, 14.25, 55.905913620, 0.001, 0).value,
            29.993601830, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                21.84602116, 29, 55.905913620, 0.1, 0).value,
            33.289964560, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                51.72271818, 29, 55.905913620, 0.01, 0).value,
            27.480618010, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                53.61322867, 29, 38.141048320, 0.001, 0).value,
            23.354012700, places=5)

    def test_rain_attenuation(self):
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 14.25, 30.87067768,
                hs=0.0691640, p=0.01, tau=0.00).value,
            7.5572640, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 14.25, 40.97052773,
                hs=0.0567010, p=0.01, tau=0.00).value,
            11.4735460, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 14.25, 47.91280491,
                hs=0.0000000, p=0.01, tau=0.00).value,
            9.7117950, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 29.00, 30.87067768,
                hs=0.0691640, p=0.01, tau=0.00).value,
            25.7166770, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 29.00, 40.97052773,
                hs=0.0567010, p=0.01, tau=0.00).value,
            39.0732330, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 29.00, 47.91280491,
                hs=0.0000000, p=0.01, tau=0.00).value,
            33.4169840, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 14.25, 59.81487174,
                hs=0.0000000, p=0.01, tau=0.00).value,
            12.7605400, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 14.25, 49.20900369,
                hs=0.0000750, p=0.01, tau=0.00).value,
            18.0693890, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 29.00, 59.81487174,
                hs=0.0000000, p=0.01, tau=0.00).value,
            46.3202450, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 29.00, 49.20900369,
                hs=0.0000750, p=0.01, tau=0.00).value,
            63.4638480, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 14.25, 55.90591362,
                hs=0.2175590, p=0.01, tau=0.00).value,
            14.1707990, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 14.25, 67.76751981,
                hs=0.2361040, p=0.01, tau=0.00).value,
            19.6617050, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 14.25, 38.14104832,
                hs=2.4500050, p=0.01, tau=0.00).value,
            11.0626540, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 29.00, 55.90591362,
                hs=0.2175590, p=0.01, tau=0.00).value,
            51.7227180, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 29.00, 67.76751981,
                hs=0.2361040, p=0.01, tau=0.00).value,
            70.5396050, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 29.00, 38.14104832,
                hs=2.4500050, p=0.01, tau=0.00).value,
            35.1160650, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 14.25, 30.87067768,
                hs=0.0691640, p=0.10, tau=0.00).value,
            2.4567600, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 14.25, 40.97052773,
                hs=0.0567010, p=0.10, tau=0.00).value,
            3.8947980, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 14.25, 47.91280491,
                hs=0.0000000, p=0.10, tau=0.00).value,
            3.2920370, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 29.00, 30.87067768,
                hs=0.0691640, p=0.10, tau=0.00).value,
            9.4912070, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 29.00, 40.97052773,
                hs=0.0567010, p=0.10, tau=0.00).value,
            15.0594580, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 29.00, 47.91280491,
                hs=0.0000000, p=0.10, tau=0.00).value,
            12.8747840, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 14.25, 59.81487174,
                hs=0.0000000, p=0.10, tau=0.00).value,
            4.9248970, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 14.25, 49.20900369,
                hs=0.0000750, p=0.10, tau=0.00).value,
            6.9660660, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 29.00, 59.81487174,
                hs=0.0000000, p=0.10, tau=0.00).value,
            20.4321420, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 29.00, 49.20900369,
                hs=0.0000750, p=0.10, tau=0.00).value,
            27.8677430, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 14.25, 55.90591362,
                hs=0.2175590, p=0.10, tau=0.00).value,
            5.2338740, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 14.25, 67.76751981,
                hs=0.2361040, p=0.10, tau=0.00).value,
            9.6728110, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 14.25, 38.14104832,
                hs=2.4500050, p=0.10, tau=0.00).value,
            4.4453390, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 29.00, 55.90591362,
                hs=0.2175590, p=0.10, tau=0.00).value,
            21.8460210, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 29.00, 67.76751981,
                hs=0.2361040, p=0.10, tau=0.00).value,
            39.6143120, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 29.00, 38.14104832,
                hs=2.4500050, p=0.10, tau=0.00).value,
            15.9048720, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 14.25, 30.87067768,
                hs=0.0691640, p=1.00, tau=0.00).value,
            0.5628470, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 14.25, 40.97052773,
                hs=0.0567010, p=1.00, tau=0.00).value,
            0.9317550, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 14.25, 47.91280491,
                hs=0.0000000, p=1.00, tau=0.00).value,
            0.7619040, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 29.00, 30.87067768,
                hs=0.0691640, p=1.00, tau=0.00).value,
            2.4686380, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 29.00, 40.97052773,
                hs=0.0567010, p=1.00, tau=0.00).value,
            4.0904280, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 29.00, 47.91280491,
                hs=0.0000000, p=1.00, tau=0.00).value,
            3.3867490, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 14.25, 59.81487174,
                hs=0.0000000, p=1.00, tau=0.00).value,
            1.0593540, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 14.25, 49.20900369,
                hs=0.0000750, p=1.00, tau=0.00).value,
            1.6122160, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 29.00, 59.81487174,
                hs=0.0000000, p=1.00, tau=0.00).value,
            5.0231130, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 29.00, 49.20900369,
                hs=0.0000750, p=1.00, tau=0.00).value,
            7.3463010, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 14.25, 55.90591362,
                hs=0.2175590, p=1.00, tau=0.00).value,
            1.2022670, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 14.25, 67.76751981,
                hs=0.2361040, p=1.00, tau=0.00).value,
            1.7852610, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 14.25, 38.14104832,
                hs=2.4500050, p=1.00, tau=0.00).value,
            0.8916230, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 29.00, 55.90591362,
                hs=0.2175590, p=1.00, tau=0.00).value,
            5.7386810, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 29.00, 67.76751981,
                hs=0.2361040, p=1.00, tau=0.00).value,
            8.3461990, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 29.00, 38.14104832,
                hs=2.4500050, p=1.00, tau=0.00).value,
            3.5957140, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 14.25, 30.87067768,
                hs=0.0691640, p=0.01, tau=0.00).value,
            7.5572640, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 14.25, 40.97052773,
                hs=0.0567010, p=0.01, tau=0.00).value,
            11.4735460, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 14.25, 47.91280491,
                hs=0.0000000, p=0.01, tau=0.00).value,
            9.7117950, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 29.00, 30.87067768,
                hs=0.0691640, p=0.01, tau=0.00).value,
            25.7166770, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 29.00, 40.97052773,
                hs=0.0567010, p=0.01, tau=0.00).value,
            39.0732330, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 29.00, 47.91280491,
                hs=0.0000000, p=0.01, tau=0.00).value,
            33.4169840, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 14.25, 59.81487174,
                hs=0.0000000, p=0.01, tau=0.00).value,
            12.7605400, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 14.25, 49.20900369,
                hs=0.0000750, p=0.01, tau=0.00).value,
            18.0693890, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 29.00, 59.81487174,
                hs=0.0000000, p=0.01, tau=0.00).value,
            46.3202450, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 29.00, 49.20900369,
                hs=0.0000750, p=0.01, tau=0.00).value,
            63.4638480, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 14.25, 55.90591362,
                hs=0.2175590, p=0.01, tau=0.00).value,
            14.1707990, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 14.25, 67.76751981,
                hs=0.2361040, p=0.01, tau=0.00).value,
            19.6617050, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 14.25, 38.14104832,
                hs=2.4500050, p=0.01, tau=0.00).value,
            11.0626540, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 29.00, 55.90591362,
                hs=0.2175590, p=0.01, tau=0.00).value,
            51.7227180, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 29.00, 67.76751981,
                hs=0.2361040, p=0.01, tau=0.00).value,
            70.5396050, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 29.00, 38.14104832,
                hs=2.4500050, p=0.01, tau=0.00).value,
            35.1160650, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 14.25, 30.87067768,
                hs=0.0691640, p=0.10, tau=0.00).value,
            2.4567600, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 14.25, 40.97052773,
                hs=0.0567010, p=0.10, tau=0.00).value,
            3.8947980, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 14.25, 47.91280491,
                hs=0.0000000, p=0.10, tau=0.00).value,
            3.2920370, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 29.00, 30.87067768,
                hs=0.0691640, p=0.10, tau=0.00).value,
            9.4912070, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 29.00, 40.97052773,
                hs=0.0567010, p=0.10, tau=0.00).value,
            15.0594580, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 29.00, 47.91280491,
                hs=0.0000000, p=0.10, tau=0.00).value,
            12.8747840, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 14.25, 59.81487174,
                hs=0.0000000, p=0.10, tau=0.00).value,
            4.9248970, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 14.25, 49.20900369,
                hs=0.0000750, p=0.10, tau=0.00).value,
            6.9660660, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 29.00, 59.81487174,
                hs=0.0000000, p=0.10, tau=0.00).value,
            20.4321420, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 29.00, 49.20900369,
                hs=0.0000750, p=0.10, tau=0.00).value,
            27.8677430, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 14.25, 55.90591362,
                hs=0.2175590, p=0.10, tau=0.00).value,
            5.2338740, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 14.25, 67.76751981,
                hs=0.2361040, p=0.10, tau=0.00).value,
            9.6728110, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 14.25, 38.14104832,
                hs=2.4500050, p=0.10, tau=0.00).value,
            4.4453390, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 29.00, 55.90591362,
                hs=0.2175590, p=0.10, tau=0.00).value,
            21.8460210, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 29.00, 67.76751981,
                hs=0.2361040, p=0.10, tau=0.00).value,
            39.6143120, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 29.00, 38.14104832,
                hs=2.4500050, p=0.10, tau=0.00).value,
            15.9048720, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 14.25, 30.87067768,
                hs=0.0691640, p=1.00, tau=0.00).value,
            0.5628470, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 14.25, 40.97052773,
                hs=0.0567010, p=1.00, tau=0.00).value,
            0.9317550, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 14.25, 47.91280491,
                hs=0.0000000, p=1.00, tau=0.00).value,
            0.7619040, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.500, 359.86, 29.00, 30.87067768,
                hs=0.0691640, p=1.00, tau=0.00).value,
            2.4686380, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.900, 12.49, 29.00, 40.97052773,
                hs=0.0567010, p=1.00, tau=0.00).value,
            4.0904280, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.940, 18.43, 29.00, 47.91280491,
                hs=0.0000000, p=1.00, tau=0.00).value,
            3.3867490, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 14.25, 59.81487174,
                hs=0.0000000, p=1.00, tau=0.00).value,
            1.0593540, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 14.25, 49.20900369,
                hs=0.0000750, p=1.00, tau=0.00).value,
            1.6122160, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.900, 316.77, 29.00, 59.81487174,
                hs=0.0000000, p=1.00, tau=0.00).value,
            5.0231130, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.780, 279.78, 29.00, 49.20900369,
                hs=0.0000750, p=1.00, tau=0.00).value,
            7.3463010, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 14.25, 55.90591362,
                hs=0.2175590, p=1.00, tau=0.00).value,
            1.2022670, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 14.25, 67.76751981,
                hs=0.2361040, p=1.00, tau=0.00).value,
            1.7852610, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 14.25, 38.14104832,
                hs=2.4500050, p=1.00, tau=0.00).value,
            0.8916230, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.30, 29.00, 55.90591362,
                hs=0.2175590, p=1.00, tau=0.00).value,
            5.7386810, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.70, 29.00, 67.76751981,
                hs=0.2361040, p=1.00, tau=0.00).value,
            8.3461990, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.050, 38.70, 29.00, 38.14104832,
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


class ITUR618_13TestCase(test.TestCase):

    def setUp(self):
        models.itu453.change_version(13)
        models.itu618.change_version(13)
        models.itu676.change_version(11)
        models.itu836.change_version(6)
        models.itu837.change_version(7)
        models.itu838.change_version(3)
        models.itu839.change_version(4)
        models.itu840.change_version(7)
        models.itu1510.change_version(1)
        models.itu1511.change_version(1)

    def test_rain_attenuation(self):
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.5, -0.14, 14.25, 31.07694309, p=1.0,
                tau=0, R001=26.48052).value,
            0.4891464, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.9, 12.49, 14.25, 40.23202374, p=1.0,
                tau=0, R001=33.936232).value,
            0.62159245, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.94, 18.43, 14.25, 46.35969261, p=1.0,
                tau=0, R001=27.13586832).value,
            0.42101702, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.5, -0.14, 14.25, 31.07694309, p=0.1,
                tau=0, R001=26.48052).value,
            2.16093996, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.9, 12.49, 14.25, 40.23202374, p=0.1,
                tau=0, R001=33.936232).value,
            2.69015654, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.94, 18.43, 14.25, 46.35969261, p=0.1,
                tau=0, R001=27.13586832).value,
            1.91338757, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.5, -0.14, 14.25, 31.07694309, p=0.01,
                tau=0, R001=26.48052).value,
            6.72784425, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.9, 12.49, 14.25, 40.23202374, p=0.01,
                tau=0, R001=33.936232).value,
            8.20500328, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.94, 18.43, 14.25, 46.35969261, p=0.01,
                tau=0, R001=27.13586832).value,
            5.9418061, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.5, -0.14, 14.25, 31.07694309, p=0.001,
                tau=0, R001=26.48052).value,
            14.76177358, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.9, 12.49, 14.25, 40.23202374, p=0.001,
                tau=0, R001=33.936232).value,
            17.636376, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.94, 18.43, 14.25, 46.35969261, p=0.001,
                tau=0, R001=27.13586832).value,
            12.98151687, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.5, -0.14, 29.0, 31.07694309, p=1.0,
                tau=0, R001=26.48052).value,
            2.17898357, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.9, 12.49, 29.0, 40.23202374, p=1.0,
                tau=0, R001=33.936232).value,
            2.81537632, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.94, 18.43, 29.0, 46.35969261, p=1.0,
                tau=0, R001=27.13586832).value,
            1.96063611, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.5, -0.14, 29.0, 31.07694309, p=0.1,
                tau=0, R001=26.48052).value,
            8.46779316, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.9, 12.49, 29.0, 40.23202374, p=0.1,
                tau=0, R001=33.936232).value,
            10.70289842, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.94, 18.43, 29.0, 46.35969261, p=0.1,
                tau=0, R001=27.13586832).value,
            7.80832251, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.5, -0.14, 29.0, 31.07694309, p=0.01,
                tau=0, R001=26.48052).value,
            23.1908096, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.9, 12.49, 29.0, 40.23202374, p=0.01,
                tau=0, R001=33.936232).value,
            28.67449232, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.94, 18.43, 29.0, 46.35969261, p=0.01,
                tau=0, R001=27.13586832).value,
            21.24861968, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                51.5, -0.14, 29.0, 31.07694309, p=0.001,
                tau=0, R001=26.48052).value,
            44.76009125, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                41.9, 12.49, 29.0, 40.23202374, p=0.001,
                tau=0, R001=33.936232).value,
            54.14015005, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                33.94, 18.43, 29.0, 46.35969261, p=0.001,
                tau=0, R001=27.13586832).value,
            40.68133015, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.9, -43.23, 14.25, 22.27833468, p=1.0,
                tau=0, R001=50.639304).value,
            1.70690128, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.78, -80.22, 14.25, 52.6789929, p=1.0,
                tau=0, R001=78.2994993).value,
            1.43904149, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.9, -43.23, 14.25, 22.27833468, p=0.1,
                tau=0, R001=50.639304).value,
            8.27164744, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.78, -80.22, 14.25, 52.6789929, p=0.1,
                tau=0, R001=78.2994993).value,
            6.30417186, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.9, -43.23, 14.25, 22.27833468, p=0.01,
                tau=0, R001=50.639304).value,
            18.94410356, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.78, -80.22, 14.25, 52.6789929, p=0.01,
                tau=0, R001=78.2994993).value,
            16.44617644, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.9, -43.23, 14.25, 22.27833468, p=0.001,
                tau=0, R001=50.639304).value,
            29.91171296, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.78, -80.22, 14.25, 52.6789929, p=0.001,
                tau=0, R001=78.2994993).value,
            29.95767701, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.9, -43.23, 29.0, 22.27833468, p=1.0,
                tau=0, R001=50.639304).value,
            6.81336808, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.78, -80.22, 29.0, 52.6789929, p=1.0,
                tau=0, R001=78.2994993).value,
            6.66385625, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.9, -43.23, 29.0, 22.27833468, p=0.1,
                tau=0, R001=50.639304).value,
            29.31896844, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.78, -80.22, 29.0, 52.6789929, p=0.1,
                tau=0, R001=78.2994993).value,
            25.59455941, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.9, -43.23, 29.0, 22.27833468, p=0.01,
                tau=0, R001=50.639304).value,
            59.62576355, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.78, -80.22, 29.0, 52.6789929, p=0.01,
                tau=0, R001=78.2994993).value,
            58.53988572, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                22.9, -43.23, 29.0, 22.27833468, p=0.001,
                tau=0, R001=50.639304).value,
            83.5996391, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                25.78, -80.22, 29.0, 52.6789929, p=0.001,
                tau=0, R001=78.2994993).value,
            93.48939944, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.3, 14.25, 48.24116215, p=1.0,
                tau=90, R001=63.61888808).value,
            1.2731081, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.7, 14.25, 85.80457401, p=1.0,
                tau=90, R001=99.15117186).value,
            1.93713255, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.05, 38.7, 14.25, 20.14348033, p=1.0,
                tau=90, R001=42.91007183).value,
            1.04440572, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.3, 14.25, 48.24116215, p=0.1,
                tau=90, R001=63.61888808).value,
            5.48101228, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.7, 14.25, 85.80457401, p=0.1,
                tau=90, R001=99.15117186).value,
            10.67987642, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.05, 38.7, 14.25, 20.14348033, p=0.1,
                tau=90, R001=42.91007183).value,
            6.0510347, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.3, 14.25, 48.24116215, p=0.01,
                tau=90, R001=63.61888808).value,
            14.85903351, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.7, 14.25, 85.80457401, p=0.01,
                tau=90, R001=99.15117186).value,
            21.03740448, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.05, 38.7, 14.25, 20.14348033, p=0.01,
                tau=90, R001=42.91007183).value,
            12.61120361, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.3, 14.25, 48.24116215, p=0.001,
                tau=90, R001=63.61888808).value,
            28.21372983, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.7, 14.25, 85.80457401, p=0.001,
                tau=90, R001=99.15117186).value,
            28.13337932, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.05, 38.7, 14.25, 20.14348033, p=0.001,
                tau=90, R001=42.91007183).value,
            17.85045772, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.3, 29.0, 48.24116215, p=1.0,
                tau=90, R001=63.61888808).value,
            5.88085649, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.7, 29.0, 85.80457401, p=1.0,
                tau=90, R001=99.15117186).value,
            9.84052929, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.05, 38.7, 29.0, 20.14348033, p=1.0,
                tau=90, R001=42.91007183).value,
            3.8213237, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.3, 29.0, 48.24116215, p=0.1,
                tau=90, R001=63.61888808).value,
            22.20219047, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.7, 29.0, 85.80457401, p=0.1,
                tau=90, R001=99.15117186).value,
            47.18910296, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.05, 38.7, 29.0, 20.14348033, p=0.1,
                tau=90, R001=42.91007183).value,
            19.80717661, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.3, 29.0, 48.24116215, p=0.01,
                tau=90, R001=63.61888808).value,
            52.7819415, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.7, 29.0, 85.80457401, p=0.01,
                tau=90, R001=99.15117186).value,
            80.85074503, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.05, 38.7, 29.0, 20.14348033, p=0.01,
                tau=90, R001=42.91007183).value,
            36.93157357, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                28.717, 77.3, 29.0, 48.24116215, p=0.001,
                tau=90, R001=63.61888808).value,
            87.88505965, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                3.133, 101.7, 29.0, 85.80457401, p=0.001,
                tau=90, R001=99.15117186).value,
            94.0437949, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation(
                9.05, 38.7, 29.0, 20.14348033, p=0.001,
                tau=90, R001=42.91007183).value,
            46.76694226, places=5)

    def test_probability_of_rain_attenuation(self):
        self.assertAlmostEqual(
            models.itu618.rain_attenuation_probability(
                51.5, -0.14, 31.07694309).value,
            7.32466089, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation_probability(
                41.9, 12.49, 40.23202374).value,
            7.08992377, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation_probability(
                33.94, 18.43, 46.35969261).value,
            1.74467895, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation_probability(
                22.9, -43.23, 22.27833468).value,
            2.5828985, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation_probability(
                25.78, -80.22, 52.6789929).value,
            4.0392312, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation_probability(
                28.717, 77.3, 48.24116215).value,
            1.64420965, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation_probability(
                3.133, 101.7, 85.80457401).value,
            5.00075505, places=4)
        self.assertAlmostEqual(
            models.itu618.rain_attenuation_probability(
                9.05, 38.7, 20.14348033).value,
            7.0357202, places=5)

#    def test_site_diversity(self):
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.768, -80.205, 9, 52.40999326,
#                    25.463, -80.486, 9, 52.48526958, 14.5, tau=0).value,
#            0.00098637, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.768, -80.205, 9, 52.40999326,
#                    25.463, -80.486, 3, 52.48526958, 14.5, tau=0).value,
#            0.0049444, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.768, -80.205, 3, 52.40999326,
#                    25.463, -80.486, 9, 52.48526958, 14.5, tau=0).value,
#            0.00503721, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.768, -80.205, 9, 52.40999326,
#                    25.463, -80.486, 9, 52.48526958, 18, tau=0).value,
#            0.00513052, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.768, -80.205, 9, 52.40999326,
#                    25.463, -80.486, 3, 52.48526958, 18, tau=0).value,
#            0.01982845, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.768, -80.205, 3, 52.40999326,
#                    25.463, -80.486, 9, 52.48526958, 18, tau=0).value,
#            0.02027952, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.768, -80.205, 9, 52.40999326,
#                    25.463, -80.486, 9, 52.48526958, 29, tau=0).value,
#            0.07543135, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.768, -80.205, 9, 52.40999326,
#                    25.463, -80.486, 3, 52.48526958, 29, tau=0).value,
#            0.16564191, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.768, -80.205, 3, 52.40999326,
#                    25.463, -80.486, 9, 52.48526958, 29, tau=0).value,
#            0.17005653, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.796, -80.287, 9, 52.33141826,
#                    25.889, -80.278, 9, 52.25682688, 29, tau=0).value,
#            0.25228844, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.796, -80.287, 9, 52.33141826,
#                    25.889, -80.278, 3, 52.25682688, 29, tau=0).value,
#            0.40360211, places=5)
#        self.assertAlmostEqual(
#            models.itu618.site_diversity_rain_outage_probability(
#                    25.796, -80.287, 3, 52.33141826,
#                    25.889, -80.278, 9, 52.25682688, 29, tau=0).value,
#            0.39740505, places=5)

    def test_scintillation_attenuation(self):
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, -0.14, 14.25, 31.07694309, 1, 1, 0.65).value,
            0.26193234, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 14.25, 40.23202374, 1, 1, 0.65).value,
            0.22405226, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 14.25, 46.35969261, 1, 1, 0.65).value,
            0.23279942, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, -0.14, 14.25, 31.07694309, 0.1, 1, 0.65).value,
            0.4228461, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 14.25, 40.23202374, 0.1, 1, 0.65).value,
            0.36169504, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 14.25, 46.35969261, 0.1, 1, 0.65).value,
            0.37581586, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, -0.14, 14.25, 31.07694309, 0.01, 1, 0.65).value,
            0.62828836, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 14.25, 40.23202374, 0.01, 1, 0.65).value,
            0.5374267, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 14.25, 46.35969261, 0.01, 1, 0.65).value,
            0.55840821, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, -0.14, 14.25, 31.07694309, 0.001, 1, 0.65).value,
            0.91021486, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 14.25, 40.23202374, 0.001, 1, 0.65).value,
            0.77858162, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 14.25, 46.35969261, 0.001, 1, 0.65).value,
            0.80897798, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, -0.14, 29, 31.07694309, 1, 1, 0.65).value,
            0.38849319, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 29, 40.23202374, 1, 1, 0.65).value,
            0.33115269, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 29, 46.35969261, 1, 1, 0.65).value,
            0.34339899, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, -0.14, 29, 31.07694309, 0.1, 1, 0.65).value,
            0.62715751, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 29, 40.23202374, 0.1, 1, 0.65).value,
            0.53459083, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 29, 46.35969261, 0.1, 1, 0.65).value,
            0.55436043, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, -0.14, 29, 31.07694309, 0.01, 1, 0.65).value,
            0.93186567, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 29, 40.23202374, 0.01, 1, 0.65).value,
            0.79432493, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 29, 46.35969261, 0.01, 1, 0.65).value,
            0.82369971, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                51.5, -0.14, 29, 31.07694309, 0.001, 1, 0.65).value,
            1.35001384, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                41.9, 12.49, 29, 40.23202374, 0.001, 1, 0.65).value,
            1.15075561, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                33.94, 18.43, 29, 46.35969261, 0.001, 1, 0.65).value,
            1.19331148, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, -43.23, 14.25, 22.27833468, 1, 1, 0.65).value,
            0.62009744, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, -80.22, 14.25, 52.6789929, 1, 1, 0.65).value,
            0.2664749, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, -43.23, 14.25, 22.27833468, 0.1, 1, 0.65).value,
            1.00104396, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, -80.22, 14.25, 52.6789929, 0.1, 1, 0.65).value,
            0.43017931, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, -43.23, 14.25, 22.27833468, 0.01, 1, 0.65).value,
            1.48740705, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, -80.22, 14.25, 52.6789929, 0.01, 1, 0.65).value,
            0.63918446, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, -43.23, 14.25, 22.27833468, 0.001, 1, 0.65).value,
            2.15483859, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, -80.22, 14.25, 52.6789929, 0.001, 1, 0.65).value,
            0.92600027, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, -43.23, 29, 22.27833468, 1, 1, 0.65).value,
            0.92341029, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, -80.22, 29, 52.6789929, 1, 1, 0.65).value,
            0.39237999, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, -43.23, 29, 22.27833468, 0.1, 1, 0.65).value,
            1.49069201, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, -80.22, 29, 52.6789929, 0.1, 1, 0.65).value,
            0.63343209, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, -43.23, 29, 22.27833468, 0.01, 1, 0.65).value,
            2.21495349, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, -80.22, 29, 52.6789929, 0.01, 1, 0.65).value,
            0.9411888, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                22.9, -43.23, 29, 22.27833468, 0.001, 1, 0.65).value,
            3.20885076, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                25.78, -80.22, 29, 52.6789929, 0.001, 1, 0.65).value,
            1.36352046, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 14.25, 48.24116215, 1, 1, 0.65).value,
            0.2156413, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 14.25, 85.80457401, 1, 1, 0.65).value,
            0.22167129, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 14.25, 20.14348033, 1, 1, 0.65).value,
            0.48533645, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 14.25, 48.24116215, 0.1, 1, 0.65).value,
            0.34811693, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 14.25, 85.80457401, 0.1, 1, 0.65).value,
            0.35785136, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 14.25, 20.14348033, 0.1, 1, 0.65).value,
            0.78349481, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 14.25, 48.24116215, 0.01, 1, 0.65).value,
            0.51725159, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 14.25, 85.80457401, 0.01, 1, 0.65).value,
            0.53171554, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 14.25, 20.14348033, 0.01, 1, 0.65).value,
            1.16416037, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 14.25, 48.24116215, 0.001, 1, 0.65).value,
            0.7493535, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 14.25, 85.80457401, 0.001, 1, 0.65).value,
            0.77030774, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 14.25, 20.14348033, 0.001, 1, 0.65).value,
            1.68654418, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 29, 48.24116215, 1, 1, 0.65).value,
            0.31791278, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 29, 85.80457401, 1, 1, 0.65).value,
            0.32486881, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 29, 20.14348033, 1, 1, 0.65).value,
            0.72351614, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 29, 48.24116215, 0.1, 1, 0.65).value,
            0.5132172, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 29, 85.80457401, 0.1, 1, 0.65).value,
            0.52444655, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 29, 20.14348033, 0.1, 1, 0.65).value,
            1.16799623, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 29, 48.24116215, 0.01, 1, 0.65).value,
            0.76256679, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 29, 85.80457401, 0.01, 1, 0.65).value,
            0.77925198, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 29, 20.14348033, 0.01, 1, 0.65).value,
            1.73547406, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                28.717, 77.3, 29, 48.24116215, 0.001, 1, 0.65).value,
            1.10474691, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                3.133, 101.7, 29, 85.80457401, 0.001, 1, 0.65).value,
            1.12891911, places=5)
        self.assertAlmostEqual(
            models.itu618.scintillation_attenuation(
                9.05, 38.7, 29, 20.14348033, 0.001, 1, 0.65).value,
            2.5142186, places=5)

    def test_rain_cross_polarization_discrimination(self):
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                0.4891464, 14.25, 31.07694309, 1.0, 0).value,
            49.57582307, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                0.62159245, 14.25, 40.23202374, 1.0, 0).value,
            49.3981550, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                0.42101702, 14.25, 46.35969261, 1.0, 0).value,
            53.93857057, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                2.16093996, 14.25, 31.07694309, 0.1, 0).value,
            40.29800396, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                2.69015654, 14.25, 40.23202374, 0.1, 0).value,
            40.28034662, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                1.91338757, 14.25, 46.35969261, 0.1, 0).value,
            44.68265675, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                6.72784425, 14.25, 31.07694309, 0.01, 0).value,
            32.97842659, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                8.20500328, 14.25, 40.23202374, 0.01, 0).value,
            33.13972017, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                5.9418061, 14.25, 46.35969261, 0.01, 0).value,
            37.62918682, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                14.76177358, 14.25, 31.07694309, 0.001, 0).value,
            28.14021762, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                17.636376, 14.25, 40.23202374, 0.001, 0).value,
            28.49940232, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                12.98151687, 14.25, 46.35969261, 0.001, 0).value,
            33.07510332, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                2.17898357, 29.0, 31.07694309, 1.0, 0).value,
            44.30006506, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                2.81537632, 29.0, 40.23202374, 1.0, 0).value,
            43.8603725, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                1.96063611, 29.0, 46.35969261, 1.0, 0).value,
            48.36964892, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                8.46779316, 29.0, 31.07694309, 0.1, 0).value,
            35.03444, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                10.70289842, 29.0, 40.23202374, 0.1, 0).value,
            34.76315732, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                7.80832251, 29.0, 46.35969261, 0.1, 0).value,
            39.12690283, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                23.1908096, 29.0, 31.07694309, 0.01, 0).value,
            27.96431726, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                28.67449232, 29.0, 40.23202374, 0.01, 0).value,
            27.8830305, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                21.24861968, 29.0, 46.35969261, 0.01, 0).value,
            32.34366876, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                44.76009125, 29.0, 31.07694309, 0.001, 0).value,
            23.64462724, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                54.14015005, 29.0, 40.23202374, 0.001, 0).value,
            23.7749224, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                40.68133015, 29.0, 46.35969261, 0.001, 0).value,
            28.33381119, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                1.70690128, 14.25, 22.27833468, 1.0, 0).value,
            38.65072987, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                1.43904149, 14.25, 52.6789929, 1.0, 0).value,
            46.23051298, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                8.27164744, 14.25, 22.27833468, 0.1, 0).value,
            27.9634536, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                6.30417186, 14.25, 52.6789929, 0.1, 0).value,
            36.82555192, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                18.94410356, 14.25, 22.27833468, 0.01, 0).value,
            22.64492814, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                16.44617644, 14.25, 52.6789929, 0.01, 0).value,
            30.86009092, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                29.91171296, 14.25, 22.27833468, 0.001, 0).value,
            20.29292318, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                29.95767701, 14.25, 52.6789929, 0.001, 0).value,
            27.62415271, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                6.81336808, 29.0, 22.27833468, 1.0, 0).value,
            33.64688473, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                6.66385625, 29.0, 52.6789929, 1.0, 0).value,
            40.0755612, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                29.31896844, 29.0, 22.27833468, 0.1, 0).value,
            22.85413903, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                25.59455941, 29.0, 52.6789929, 0.1, 0).value,
            30.6650529, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                59.62576355, 29.0, 22.27833468, 0.01, 0).value,
            17.88255372, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                58.53988572, 29.0, 52.6789929, 0.01, 0).value,
            25.03203051, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                83.5996391, 29.0, 22.27833468, 0.001, 0).value,
            16.16922861, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                93.48939944, 29.0, 52.6789929, 0.001, 0).value,
            22.41718851, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                1.2731081, 14.25, 48.24116215, 1.0, 90).value,
            45.80237934, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                1.93713255, 14.25, 85.80457401, 1.0, 90).value,
            75.12972446, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                1.04440572, 14.25, 20.14348033, 1.0, 90).value,
            42.28242577, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                5.48101228, 14.25, 48.24116215, 0.1, 90).value,
            36.51649699, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                10.67987642, 14.25, 85.80457401, 0.1, 90).value,
            65.51910372, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                6.0510347, 14.25, 20.14348033, 0.1, 90).value,
            30.32827626, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                14.85903351, 14.25, 48.24116215, 0.01, 90).value,
            30.19759496, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                21.03740448, 14.25, 85.80457401, 0.01, 90).value,
            63.60558975, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                12.61120361, 14.25, 20.14348033, 0.01, 90).value,
            25.96615447, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                28.21372983, 14.25, 48.24116215, 0.001, 90).value,
            26.54453432, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                28.13337932, 14.25, 85.80457401, 0.001, 90).value,
            64.9390721, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                17.85045772, 14.25, 20.14348033, 0.001, 90).value,
            24.79563761, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                5.88085649, 29.0, 48.24116215, 1.0, 90).value,
            39.73121592, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                9.84052929, 29.0, 85.80457401, 1.0, 90).value,
            68.04931733, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                3.8213237, 29.0, 20.14348033, 1.0, 90).value,
            38.25788658, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                22.20219047, 29.0, 48.24116215, 0.1, 90).value,
            30.45232663, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                47.18910296, 29.0, 85.80457401, 0.1, 90).value,
            58.32352317, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                19.80717661, 29.0, 20.14348033, 0.1, 90).value,
            26.0924586, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                52.7819415, 29.0, 48.24116215, 0.01, 90).value,
            24.4471052, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                80.85074503, 29.0, 85.80457401, 0.01, 90).value,
            56.92074963, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                36.93157357, 29.0, 20.14348033, 0.01, 90).value,
            22.11041426, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                87.88505965, 29.0, 48.24116215, 0.001, 90).value,
            21.39198393, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                94.0437949, 29.0, 85.80457401, 0.001, 90).value,
            59.09547625, places=5)
        self.assertAlmostEqual(
            models.itu618.rain_cross_polarization_discrimination(
                46.76694226, 29.0, 20.14348033, 0.001, 90).value,
            21.61918999, places=5)

    def total_attenuation_fcn(self, lat, lon, f, el, p, D, eta, tau,
                              val_g, val_c, val_r, val_s, val_t):

        # The validation function uses the exact method to compute the rainfall
        # rate exceeded for 0.01% of the time
        R001 = models.itu837.rainfall_rate(lat, lon, 0.01000000001)
        A_g, A_c, A_r, A_s, A = itur.atmospheric_attenuation_slant_path(
            lat, lon, f, el, p, D, eta=eta, tau=tau, R001=R001,
            return_contributions=True)

        self.assertAlmostEqual(A_g.value, val_g, places=5)
        self.assertAlmostEqual(A_c.value, val_c, places=5)
        self.assertAlmostEqual(A_r.value, val_r, places=5)
        self.assertAlmostEqual(A_s.value, val_s, places=5)
        self.assertAlmostEqual(A.value, val_t, places=5)

    def test_total_attenuation(self):
        self.total_attenuation_fcn(
            51.5, -0.14, 14.25, 31.07694309, 1, 1, 0.65, 0,
            0.223693782, 0.45517046, 0.48914539, 0.26193234, 1.203663661)
        self.total_attenuation_fcn(
            41.9, 12.49, 14.25, 40.23202374, 1, 1, 0.65, 0,
            0.184499507, 0.26338517, 0.62159459, 0.22405226, 1.097400703)
        self.total_attenuation_fcn(
            33.94, 18.43, 14.25, 46.35969261, 1, 1, 0.65, 0,
            0.168635988, 0.18779409, 0.42101546, 0.23279942, 0.820437057)
        self.total_attenuation_fcn(
            51.5, -0.14, 14.25, 31.07694309, 0.1, 1, 0.65, 0,
            0.223693782, 0.45517046, 2.16093588, 0.4228461, 2.873752501)
        self.total_attenuation_fcn(
            41.9, 12.49, 14.25, 40.23202374, 0.1, 1, 0.65, 0,
            0.184499507, 0.26338517, 2.69016502, 0.36169504, 3.16011407)
        self.total_attenuation_fcn(
            33.94, 18.43, 14.25, 46.35969261, 0.1, 1, 0.65, 0,
            0.168635988, 0.18779409, 1.91338106, 0.37581586, 2.303155743)
        self.total_attenuation_fcn(
            51.5, -0.14, 14.25, 31.07694309, 0.01, 1, 0.65, 0,
            0.223693782, 0.45517046, 6.72783273, 0.62828836, 7.434122415)
        self.total_attenuation_fcn(
            41.9, 12.49, 14.25, 40.23202374, 0.01, 1, 0.65, 0,
            0.184499507, 0.26338517, 8.20502671, 0.5374267, 8.669947478)
        self.total_attenuation_fcn(
            33.94, 18.43, 14.25, 46.35969261, 0.01, 1, 0.65, 0,
            0.168635988, 0.18779409, 5.94178778, 0.55840821, 6.323600947)
        self.total_attenuation_fcn(
            51.5, -0.14, 14.25, 31.07694309, 0.001, 1, 0.65, 0,
            0.223693782, 0.45517046, 14.76175093, 0.91021486, 15.46781355)
        self.total_attenuation_fcn(
            41.9, 12.49, 14.25, 40.23202374, 0.001, 1, 0.65, 0,
            0.184499507, 0.26338517, 17.63642115, 0.77858162, 18.10123067)
        self.total_attenuation_fcn(
            33.94, 18.43, 14.25, 46.35969261, 0.001, 1, 0.65, 0,
            0.168635988, 0.18779409, 12.981481, 0.80897798, 13.36273512)
        self.total_attenuation_fcn(
            51.5, -0.14, 29, 31.07694309, 1, 1, 0.65, 0,
            0.799999368, 1.77247154, 2.17897957, 0.38849319, 4.770502219)
        self.total_attenuation_fcn(
            41.9, 12.49, 29, 40.23202374, 1, 1, 0.65, 0,
            0.673619867, 1.0256437, 2.81538514, 0.33115269, 4.528897381)
        self.total_attenuation_fcn(
            33.94, 18.43, 29, 46.35969261, 1, 1, 0.65, 0,
            0.62972417, 0.73128577, 1.96062953, 0.34339899, 3.343454224)
        self.total_attenuation_fcn(
            51.5, -0.14, 29, 31.07694309, 0.1, 1, 0.65, 0,
            0.799999368, 1.77247154, 8.46777895, 0.62715751, 11.05943682)
        self.total_attenuation_fcn(
            41.9, 12.49, 29, 40.23202374, 0.1, 1, 0.65, 0,
            0.673619867, 1.0256437, 10.70292908, 0.53459083, 12.41436971)
        self.total_attenuation_fcn(
            33.94, 18.43, 29, 46.35969261, 0.1, 1, 0.65, 0,
            0.62972417, 0.73128577, 7.80829852, 0.55436043, 9.18728313)
        self.total_attenuation_fcn(
            51.5, -0.14, 29, 31.07694309, 0.01, 1, 0.65, 0,
            0.799999368, 1.77247154, 23.19077435, 0.93186567, 25.78063225)
        self.total_attenuation_fcn(
            41.9, 12.49, 29, 40.23202374, 0.01, 1, 0.65, 0,
            0.673619867, 1.0256437, 28.67456675, 0.79432493, 30.38445044)
        self.total_attenuation_fcn(
            33.94, 18.43, 29, 46.35969261, 0.01, 1, 0.65, 0,
            0.62972417, 0.73128577, 21.24856054, 0.82369971, 22.62499923)
        self.total_attenuation_fcn(
            51.5, -0.14, 29, 31.07694309, 0.001, 1, 0.65, 0,
            0.799999368, 1.77247154, 44.76003026, 1.35001384, 47.35208054)
        self.total_attenuation_fcn(
            41.9, 12.49, 29, 40.23202374, 0.001, 1, 0.65, 0,
            0.673619867, 1.0256437, 54.14027603, 1.15075561, 55.85154062)
        self.total_attenuation_fcn(
            33.94, 18.43, 29, 46.35969261, 0.001, 1, 0.65, 0,
            0.62972417, 0.73128577, 40.68122866, 1.19331148, 42.05942781)
        self.total_attenuation_fcn(
            22.9, -43.23, 14.25, 22.27833468, 1, 1, 0.65, 0,
            0.383178724, 0.54183293, 1.70690691, 0.62009744, 2.715849229)
        self.total_attenuation_fcn(
            25.78, -80.22, 14.25, 52.6789929, 1, 1, 0.65, 0,
            0.206227197, 0.53317506, 1.43904233, 0.2664749, 2.196365451)
        self.total_attenuation_fcn(
            22.9, -43.23, 14.25, 22.27833468, 0.1, 1, 0.65, 0,
            0.383178724, 0.54183293, 8.27167236, 1.00104396, 9.253351467)
        self.total_attenuation_fcn(
            25.78, -80.22, 14.25, 52.6789929, 0.1, 1, 0.65, 0,
            0.206227197, 0.53317506, 6.30417519, 0.43017931, 7.057096675)
        self.total_attenuation_fcn(
            22.9, -43.23, 14.25, 22.27833468, 0.01, 1, 0.65, 0,
            0.383178724, 0.54183293, 18.94415527, 1.48740705, 19.92585295)
        self.total_attenuation_fcn(
            25.78, -80.22, 14.25, 52.6789929, 0.01, 1, 0.65, 0,
            0.206227197, 0.53317506, 16.44618432, 0.63918446, 17.1976133)
        self.total_attenuation_fcn(
            22.9, -43.23, 14.25, 22.27833468, 0.001, 1, 0.65, 0,
            0.383178724, 0.54183293, 29.91178614, 2.15483859, 30.91293869)
        self.total_attenuation_fcn(
            25.78, -80.22, 14.25, 52.6789929, 0.001, 1, 0.65, 0,
            0.206227197, 0.53317506, 29.95768987, 0.92600027, 30.71115009)
        self.total_attenuation_fcn(
            22.9, -43.23, 29, 22.27833468, 1, 1, 0.65, 0,
            1.504259763, 2.1099424, 6.81338837, 0.92341029, 10.4752418)
        self.total_attenuation_fcn(
            25.78, -80.22, 29, 52.6789929, 1, 1, 0.65, 0,
            0.827675954, 2.07622792, 6.66385994, 0.39237999, 9.576567189)
        self.total_attenuation_fcn(
            22.9, -43.23, 29, 22.27833468, 0.1, 1, 0.65, 0,
            1.504259763, 2.1099424, 29.31904828, 1.49069201, 32.9685827)
        self.total_attenuation_fcn(
            25.78, -80.22, 29, 52.6789929, 0.1, 1, 0.65, 0,
            0.827675954, 2.07622792, 25.59457239, 0.63343209, 28.50572549)
        self.total_attenuation_fcn(
            22.9, -43.23, 29, 22.27833468, 0.01, 1, 0.65, 0,
            1.504259763, 2.1099424, 59.62591067, 2.21495349, 63.27983401)
        self.total_attenuation_fcn(
            25.78, -80.22, 29, 52.6789929, 0.01, 1, 0.65, 0,
            0.827675954, 2.07622792, 58.53991262, 0.9411888, 61.45112298)
        self.total_attenuation_fcn(
            22.9, -43.23, 29, 22.27833468, 0.001, 1, 0.65, 0,
            1.504259763, 2.1099424, 83.59982398, 3.20885076, 87.2740725)
        self.total_attenuation_fcn(
            25.78, -80.22, 29, 52.6789929, 0.001, 1, 0.65, 0,
            0.827675954, 2.07622792, 93.48943794, 1.36352046, 96.4030686)
        self.total_attenuation_fcn(
            28.717, 77.3, 14.25, 48.24116215, 1, 1, 0.65, 90,
            0.257653026, 0.68592197, 1.27311232, 0.2156413, 2.228519972)
        self.total_attenuation_fcn(
            3.133, 101.7, 14.25, 85.80457401, 1, 1, 0.65, 90,
            0.163655312, 0.62211863, 1.93712821, 0.22167129, 2.732484342)
        self.total_attenuation_fcn(
            9.05, 38.7, 14.25, 20.14348033, 1, 1, 0.65, 90,
            0.22310495, 0.65764822, 1.04440674, 0.48533645, 1.993003982)
        self.total_attenuation_fcn(
            28.717, 77.3, 14.25, 48.24116215, 0.1, 1, 0.65, 90,
            0.257653026, 0.68592197, 5.48102886, 0.34811693, 6.434421439)
        self.total_attenuation_fcn(
            3.133, 101.7, 14.25, 85.80457401, 0.1, 1, 0.65, 90,
            0.163655312, 0.62211863, 10.67985456, 0.35785136, 11.47129236)
        self.total_attenuation_fcn(
            9.05, 38.7, 14.25, 20.14348033, 0.1, 1, 0.65, 90,
            0.22310495, 0.65764822, 6.05104013, 0.78349481, 6.977389778)
        self.total_attenuation_fcn(
            28.717, 77.3, 14.25, 48.24116215, 0.01, 1, 0.65, 90,
            0.257653026, 0.68592197, 14.85907425, 0.51725159, 15.81125251)
        self.total_attenuation_fcn(
            3.133, 101.7, 14.25, 85.80457401, 0.01, 1, 0.65, 90,
            0.163655312, 0.62211863, 21.03736546, 0.53171554, 21.82966493)
        self.total_attenuation_fcn(
            9.05, 38.7, 14.25, 20.14348033, 0.01, 1, 0.65, 90,
            0.22310495, 0.65764822, 12.61121387, 1.16416037, 13.54293868)
        self.total_attenuation_fcn(
            28.717, 77.3, 14.25, 48.24116215, 0.001, 1, 0.65, 90,
            0.257653026, 0.68592197, 28.21379917, 0.7493535, 29.16708769)
        self.total_attenuation_fcn(
            3.133, 101.7, 14.25, 85.80457401, 0.001, 1, 0.65, 90,
            0.163655312, 0.62211863, 28.13333254, 0.77030774, 28.92942223)
        self.total_attenuation_fcn(
            9.05, 38.7, 14.25, 20.14348033, 0.001, 1, 0.65, 90,
            0.22310495, 0.65764822, 17.85047073, 1.68654418, 18.80790784)
        self.total_attenuation_fcn(
            28.717, 77.3, 29, 48.24116215, 1, 1, 0.65, 90,
            1.038585522, 2.67103709, 5.88087518, 0.31791278, 9.596404871)
        self.total_attenuation_fcn(
            3.133, 101.7, 29, 85.80457401, 1, 1, 0.65, 90,
            0.645959831, 2.42258159, 9.84050759, 0.32486881, 12.9133514)
        self.total_attenuation_fcn(
            9.05, 38.7, 29, 20.14348033, 1, 1, 0.65, 90,
            0.703128217, 2.56093676, 3.82132703, 0.72351614, 7.126271267)
        self.total_attenuation_fcn(
            28.717, 77.3, 29, 48.24116215, 0.1, 1, 0.65, 90,
            1.038585522, 2.67103709, 22.20225497, 0.5132172, 25.9171717)
        self.total_attenuation_fcn(
            3.133, 101.7, 29, 85.80457401, 0.1, 1, 0.65, 90,
            0.645959831, 2.42258159, 47.18900784, 0.52444655, 50.26032116)
        self.total_attenuation_fcn(
            9.05, 38.7, 29, 20.14348033, 0.1, 1, 0.65, 90,
            0.703128217, 2.56093676, 19.80719237, 1.16799623, 23.10173121)
        self.total_attenuation_fcn(
            28.717, 77.3, 29, 48.24116215, 0.01, 1, 0.65, 90,
            1.038585522, 2.67103709, 52.78208044, 0.76256679, 56.49694605)
        self.total_attenuation_fcn(
            3.133, 101.7, 29, 85.80457401, 0.01, 1, 0.65, 90,
            0.645959831, 2.42258159, 80.85059735, 0.77925198, 83.92278473)
        self.total_attenuation_fcn(
            9.05, 38.7, 29, 20.14348033, 0.01, 1, 0.65, 90,
            0.703128217, 2.56093676, 36.9316002, 1.73547406, 40.23377893)
        self.total_attenuation_fcn(
            28.717, 77.3, 29, 48.24116215, 0.001, 1, 0.65, 90,
            1.038585522, 2.67103709, 87.88526702, 1.10474691, 91.6016281)
        self.total_attenuation_fcn(
            3.133, 101.7, 29, 85.80457401, 0.001, 1, 0.65, 90,
            0.645959831, 2.42258159, 94.04364093, 1.12891911, 97.11878785)
        self.total_attenuation_fcn(
            9.05, 38.7, 29, 20.14348033, 0.001, 1, 0.65, 90,
            0.703128217, 2.56093676, 46.76697248, 2.5142186, 50.09507012)


# class ITUR840_4TestCase(test.TestCase):
#
#    def setUp(self):
#        models.itu840.change_version(4)
#
#    def test_columnar_content_reduced_liquid(self):
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                51.5, -0.14, 1.000).value,
#            1.26328612, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                41.9, 12.49, 1.000).value,
#            0.91467189, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                33.94, 18.43, 1.000).value,
#            0.73072098, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                51.5, -0.14, 0.100).value,
#            1.90329847, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                41.9, 12.49, 0.100).value,
#            1.49845951, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    33.94, 18.43, 0.100).value,
#            1.47628568, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    51.5, -0.14, 0.010).value,
#            1.90329847, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    41.9, 12.49, 0.010).value,
#            1.49845951, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    33.94, 18.43, 0.010).value,
#            1.47628568, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    51.5, -0.14, 0.001).value,
#            1.90329847, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    41.9, 12.49, 0.001).value,
#            1.49845951, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    33.94, 18.43, 0.001).value,
#            1.47628568, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                51.5, -0.14, 1.000).value,
#            1.26328612, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                41.9, 12.49, 1.000).value,
#            0.91467189, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                33.94, 18.43, 1.000).value,
#            0.73072098, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                51.5, -0.14, 0.100).value,
#            1.90329847, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                41.9, 12.49, 0.100).value,
#            1.49845951, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                33.94, 18.43, 0.100).value,
#            1.47628568, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    51.5, -0.14, 0.010).value,
#            1.90329847, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    41.9, 12.49, 0.010).value,
#            1.49845951, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    33.94, 18.43, 0.010).value,
#            1.47628568, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    51.5, -0.14, 0.001).value,
#            1.90329847, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    41.9, 12.49, 0.001).value,
#            1.49845951, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    33.94, 18.43, 0.001).value,
#            1.47628568, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                22.9, -43.23, 1.000).value,
#            1.10444871, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                25.78, -80.22, 1.000).value,
#            2.27978216, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                22.9, -43.23, 0.100).value,
#            2.82993169, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                25.78, -80.22, 0.100).value,
#            3.52927516, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    22.9, -43.23, 0.010).value,
#            2.82993169, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    25.78, -80.22, 0.010).value,
#            3.52927516, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    22.9, -43.23, 0.001).value,
#            2.82993169, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    25.78, -80.22, 0.001).value,
#            3.52927516, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                22.9, -43.23, 1.000).value,
#            1.10444871, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                25.78, -80.22, 1.000).value,
#            2.27978216, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                22.9, -43.23, 0.100).value,
#            2.82993169, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                25.78, -80.22, 0.100).value,
#            3.52927516, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    22.9, -43.23, 0.010).value,
#            2.82993169, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    25.78, -80.22, 0.010).value,
#            3.52927516, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    22.9, -43.23, 0.001).value,
#            2.82993169, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                    25.78, -80.22, 0.001).value,
#            3.52927516, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                28.717, 77.3, 1.000).value,
#            2.75109958, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                3.133, 101.7, 1.000).value,
#            3.33600769, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                9.05, 38.7, 1.000).value,
#            1.21770185, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                28.717, 77.3, 0.100).value,
#            4.23072604, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                3.133, 101.7, 0.100).value,
#            3.80525123, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                9.05, 38.7, 0.100).value,
#            1.49251459, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                28.717, 77.3, 0.010).value,
#            4.23072604, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                3.133, 101.7, 0.010).value,
#            3.80525123, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                9.05, 38.7, 0.010).value,
#            1.49251459, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                28.717, 77.3, 0.001).value,
#            4.23072604, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                3.133, 101.7, 0.001).value,
#            3.80525123, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                9.05, 38.7, 0.001).value,
#            1.49251459, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                28.717, 77.3, 1.000).value,
#            2.75109958, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                3.133, 101.7, 1.000).value,
#            3.33600769, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                9.05, 38.7, 1.000).value,
#            1.21770185, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                28.717, 77.3, 0.100).value,
#            4.23072604, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                3.133, 101.7, 0.100).value,
#            3.80525123, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                9.05, 38.7, 0.100).value,
#            1.49251459, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                28.717, 77.3, 0.010).value,
#            4.23072604, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                3.133, 101.7, 0.010).value,
#            3.80525123, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                9.05, 38.7, 0.010).value,
#            1.49251459, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                28.717, 77.3, 0.001).value,
#            4.23072604, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                3.133, 101.7, 0.001).value,
#            3.80525123, places=5)
#        self.assertAlmostEqual(
#            models.itu840.columnar_content_reduced_liquid(
#                9.05, 38.7, 0.001).value,
#            1.49251459, places=5)
#
#    def test_cloud_attenuation(self):
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                51.5, -0.14, 30.87067768, 14.25, 1.000).value,
#            0.45792895, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                41.9, 12.49, 40.97052773, 14.25, 1.000).value,
#            0.25946553, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                33.94, 18.43, 47.91280491, 14.25, 1.000).value,
#            0.18313623, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                51.5, -0.14, 30.87067768, 14.25, 0.100).value,
#            0.68992722, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                41.9, 12.49, 40.97052773, 14.25, 0.100).value,
#            0.42506892, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                33.94, 18.43, 47.91280491, 14.25, 0.100).value,
#            0.36999265, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                51.5, -0.14, 30.87067768, 14.25, 0.010).value,
#            0.68992722, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                41.9, 12.49, 40.97052773, 14.25, 0.010).value,
#            0.42506892, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                33.94, 18.43, 47.91280491, 14.25, 0.010).value,
#            0.36999265, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                51.5, -0.14, 30.87067768, 14.25, 0.001).value,
#            0.68992722, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                41.9, 12.49, 40.97052773, 14.25, 0.001).value,
#            0.42506892, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                33.94, 18.43, 47.91280491, 14.25, 0.001).value,
#            0.36999265, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                51.5, -0.14, 30.87067768, 29, 1.000).value,
#            1.79599547, places=2)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                41.9, 12.49, 40.97052773, 29, 1.000).value,
#            1.01762274, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                33.94, 18.43, 47.91280491, 29, 1.000).value,
#            0.71825953, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                51.5, -0.14, 30.87067768, 29, 0.100).value,
#            2.70589171, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                41.9, 12.49, 40.97052773, 29, 0.100).value,
#            1.66711854, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                33.94, 18.43, 47.91280491, 29, 0.100).value,
#            1.45110964, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                51.5, -0.14, 30.87067768, 29, 0.010).value,
#            2.70589171, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                41.9, 12.49, 40.97052773, 29, 0.010).value,
#            1.66711854, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                33.94, 18.43, 47.91280491, 29, 0.010).value,
#            1.45110964, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                51.5, -0.14, 30.87067768, 29, 0.001).value,
#            2.70589171, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                41.9, 12.49, 40.97052773, 29, 0.001).value,
#            1.66711854, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                33.94, 18.43, 47.91280491, 29, 0.001).value,
#            1.45110964, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                22.9, -43.23, 59.81487174, 14.25, 1.000).value,
#            0.23764476, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                25.78, -80.22, 49.20900369, 14.25, 1.000).value,
#            0.56006901, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                22.9, -43.23, 59.81487174, 14.25, 0.100).value,
#            0.60891776, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                25.78, -80.22, 49.20900369, 14.25, 0.100).value,
#            0.86702917, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                22.9, -43.23, 59.81487174, 14.25, 0.010).value,
#            0.60891776, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                25.78, -80.22, 49.20900369, 14.25, 0.010).value,
#            0.86702917, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                22.9, -43.23, 59.81487174, 14.25, 0.001).value,
#            0.60891776, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                25.78, -80.22, 49.20900369, 14.25, 0.001).value,
#            0.86702917, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                22.9, -43.23, 59.81487174, 29, 1.000).value,
#            0.93204177, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                25.78, -80.22, 49.20900369, 29, 1.000).value,
#            2.19658834, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                22.9, -43.23, 59.81487174, 29, 0.100).value,
#            2.38817297, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                25.78, -80.22, 49.20900369, 29, 0.100).value,
#            3.40048483, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                22.9, -43.23, 59.81487174, 29, 0.010).value,
#            2.38817297, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                25.78, -80.22, 49.20900369, 29, 0.010).value,
#            3.40048483, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                22.9, -43.23, 59.81487174, 29, 0.001).value,
#            2.38817297, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                25.78, -80.22, 49.20900369, 29, 0.001).value,
#            3.40048483, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                28.717, 77.3, 55.90591362, 14.25, 1.000).value,
#            0.6178942, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                3.133, 101.7, 67.76751981, 14.25, 1.000).value,
#            0.67031269, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                9.05, 38.7, 38.14104832, 14.25, 1.000).value,
#            0.36671963, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                28.717, 77.3, 55.90591362, 14.25, 0.100).value,
#            0.95021681, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                3.133, 101.7, 67.76751981, 14.25, 0.100).value,
#            0.76459901, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                9.05, 38.7, 38.14104832, 14.25, 0.100).value,
#            0.44948146, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                28.717, 77.3, 55.90591362, 14.25, 0.010).value,
#            0.95021681, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                3.133, 101.7, 67.76751981, 14.25, 0.010).value,
#            0.76459901, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                9.05, 38.7, 38.14104832, 14.25, 0.010).value,
#            0.44948146, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                28.717, 77.3, 55.90591362, 14.25, 0.001).value,
#            0.95021681, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                3.133, 101.7, 67.76751981, 14.25, 0.001).value,
#            0.76459901, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                9.05, 38.7, 38.14104832, 14.25, 0.001).value,
#            0.44948146, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                28.717, 77.3, 55.90591362, 29, 1.000).value,
#            2.4233785, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                3.133, 101.7, 67.76751981, 29, 1.000).value,
#            2.6289636, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                9.05, 38.7, 38.14104832, 29, 1.000).value,
#            1.43827289, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                28.717, 77.3, 55.90591362, 29, 0.100).value,
#            3.72674641, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                3.133, 101.7, 67.76751981, 29, 0.100).value,
#            2.99875418, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                9.05, 38.7, 38.14104832, 29, 0.100).value,
#            1.76286444, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                28.717, 77.3, 55.90591362, 29, 0.010).value,
#            3.72674641, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                3.133, 101.7, 67.76751981, 29, 0.010).value,
#            2.99875418, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                9.05, 38.7, 38.14104832, 29, 0.010).value,
#            1.76286444, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                28.717, 77.3, 55.90591362, 29, 0.001).value,
#            3.72674641, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                3.133, 101.7, 67.76751981, 29, 0.001).value,
#            2.99875418, places=3)
#        self.assertAlmostEqual(
#            models.itu840.cloud_attenuation(
#                9.05, 38.7, 38.14104832, 29, 0.001).value,
#            1.76286444, places=3)


class ITUR840_7TestCase(test.TestCase):

    def setUp(self):
        models.itu840.change_version(7)

    def test_columnar_content_reduced_liquid(self):
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                3.133, 101.7, 0.1).value,
            3.805251208, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                3.133, 101.7, 0.15).value,
            3.744512329, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                3.133, 101.7, 0.3).value,
            3.630957766, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                3.133, 101.7, 0.35).value,
            3.594946111, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                22.9, -43.23, 0.1).value,
            2.829931669, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                22.9, -43.23, 0.15).value,
            2.615428331, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                22.9, -43.23, 0.3).value,
            2.152560931, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                22.9, -43.23, 0.35).value,
            2.030424796, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                23, 30, 0.1).value,
            0.443821013, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                23, 30, 0.15).value,
            0.367758574, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                23, 30, 0.3).value,
            0.25249597, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                23, 30, 0.35).value,
            0.230476914, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                25.78, -80.22, 0.1).value,
            3.52927514, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                25.78, -80.22, 0.15).value,
            3.368053109, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                25.78, -80.22, 0.3).value,
            3.090031167, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                25.78, -80.22, 0.35).value,
            2.98280226, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                28.717, 77.3, 0.1).value,
            4.230726014, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                28.717, 77.3, 0.15).value,
            4.004951665, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                28.717, 77.3, 0.3).value,
            3.641943304, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                28.717, 77.3, 0.35).value,
            3.550068054, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                33.94, 18.43, 0.1).value,
            1.476285677, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                33.94, 18.43, 0.15).value,
            1.342662497, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                33.94, 18.43, 0.3).value,
            1.117630129, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                33.94, 18.43, 0.35).value,
            1.061278891, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                41.9, 12.49, 0.1).value,
            1.498459518, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                41.9, 12.49, 0.15).value,
            1.411411719, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                41.9, 12.49, 0.3).value,
            1.254176128, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                41.9, 12.49, 0.35).value,
            1.214239524, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                51.5, -0.14, 0.1).value,
            1.903298487, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                51.5, -0.14, 0.15).value,
            1.803803604, places=5)
        self.assertAlmostEqual(
            models.itu840.columnar_content_reduced_liquid(
                51.5, -0.14, 0.3).value,
            1.641289077, places=5)

    def test_cloud_attenuation(self):
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 31.07694309, 14.25, 1.0).value,
            0.45517046, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.23202374, 14.25, 1.0).value,
            0.26338517, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 46.35969261, 14.25, 1.0).value,
            0.18779409, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 31.07694309, 14.25, 0.5).value,
            0.53457216, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.23202374, 14.25, 0.5).value,
            0.3230387, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 46.35969261, 14.25, 0.5).value,
            0.23923797, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 31.07694309, 14.25, 0.3).value,
            0.59136745, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.23202374, 14.25, 0.3).value,
            0.36114741, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 46.35969261, 14.25, 0.3).value,
            0.2872291, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 31.07694309, 14.25, 0.2).value,
            0.62448748, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.23202374, 14.25, 0.2).value,
            0.38863977, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 46.35969261, 14.25, 0.2).value,
            0.32069677, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 31.07694309, 29, 1.0).value,
            1.77247154, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.23202374, 29, 1.0).value,
            1.0256437, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 46.35969261, 29, 1.0).value,
            0.73128577, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 31.07694309, 29, 0.5).value,
            2.08166837, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.23202374, 29, 0.5).value,
            1.2579395, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 46.35969261, 29, 0.5).value,
            0.9316125, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 31.07694309, 29, 0.3).value,
            2.30283391, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.23202374, 29, 0.3).value,
            1.40633801, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 46.35969261, 29, 0.3).value,
            1.11849396, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                51.5, -0.14, 31.07694309, 29, 0.2).value,
            2.43180607, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                41.9, 12.49, 40.23202374, 29, 0.2).value,
            1.51339553, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                33.94, 18.43, 46.35969261, 29, 0.2).value,
            1.24881983, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 22.27833468, 14.25, 1.0).value,
            0.54183293, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 52.6789929, 14.25, 1.0).value,
            0.53317506, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 22.27833468, 14.25, 0.5).value,
            0.85746792, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 52.6789929, 14.25, 0.5).value,
            0.63956606, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 22.27833468, 14.25, 0.3).value,
            1.05602769, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 52.6789929, 14.25, 0.3).value,
            0.72266885, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 22.27833468, 14.25, 0.2).value,
            1.20844208, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 52.6789929, 14.25, 0.2).value,
            0.76093789, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 22.27833468, 29, 1.0).value,
            2.1099424, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 52.6789929, 29, 1.0).value,
            2.07622792, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 22.27833468, 29, 0.5).value,
            3.33905126, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 52.6789929, 29, 0.5).value,
            2.49052334, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 22.27833468, 29, 0.3).value,
            4.11225948, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 52.6789929, 29, 0.3).value,
            2.81413248, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                22.9, -43.23, 22.27833468, 29, 0.2).value,
            4.70577375, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                25.78, -80.22, 52.6789929, 29, 0.2).value,
            2.96315532, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.72, 77.3, 48.24116215, 14.25, 1.0).value,
            0.68560078, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.13, 101.7, 85.80457401, 14.25, 1.0).value,
            0.62214817, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 20.14348033, 14.25, 1.0).value,
            0.65764822, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.72, 77.3, 48.24116215, 14.25, 0.5).value,
            0.83179446, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.13, 101.7, 85.80457401, 14.25, 0.5).value,
            0.65489922, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 20.14348033, 14.25, 0.5).value,
            0.7181604, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.72, 77.3, 48.24116215, 14.25, 0.3).value,
            0.90773089, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.13, 101.7, 85.80457401, 14.25, 0.3).value,
            0.6771593, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 20.14348033, 14.25, 0.3).value,
            0.75244454, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.72, 77.3, 48.24116215, 14.25, 0.2).value,
            0.95830261, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.13, 101.7, 85.80457401, 14.25, 0.2).value,
            0.69030616, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 20.14348033, 14.25, 0.2).value,
            0.77111549, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.72, 77.3, 48.24116215, 29, 1.0).value,
            2.66978635, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.13, 101.7, 85.80457401, 29, 1.0).value,
            2.42269662, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 20.14348033, 29, 1.0).value,
            2.56093676, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.72, 77.3, 48.24116215, 29, 0.5).value,
            3.23907665, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.13, 101.7, 85.80457401, 29, 0.5).value,
            2.55023192, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 20.14348033, 29, 0.5).value,
            2.79657622, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.72, 77.3, 48.24116215, 29, 0.3).value,
            3.53477943, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.13, 101.7, 85.80457401, 29, 0.3).value,
            2.63691452, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 20.14348033, 29, 0.3).value,
            2.93008149, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                28.72, 77.3, 48.24116215, 29, 0.2).value,
            3.73170991, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                3.13, 101.7, 85.80457401, 29, 0.2).value,
            2.68810948, places=5)
        self.assertAlmostEqual(
            models.itu840.cloud_attenuation(
                9.05, 38.7, 20.14348033, 29, 0.2).value,
            3.00278773, places=5)


class ITUR1511_1TestCase(test.TestCase):

    def setUp(self):
        models.itu1511.change_version(1)

    def test_topographic_altitude(self):
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(3.133, 101.7).value,
            0.23610446, places=5)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(22.9, -43.23).value,
            0.0, places=5)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(23.0, 30.0).value,
            0.247, places=5)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(25.78, -80.22).value,
            7.511e-05, places=5)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(28.717, 77.3).value,
            0.21755946, places=5)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(33.94, 18.43).value,
            0.0, places=5)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(41.9, 12.49).value,
            0.05670104, places=5)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(51.5, -0.14).value,
            0.06916422, places=5)


class ITUR1511_2TestCase(test.TestCase):

    def setUp(self):
        models.itu1511.change_version(2)

    def test_topographic_altitude(self):
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(51.5, -0.14).value,
            0.031382983999999, places=4)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(41.9, 12.49).value,
            0.0461229880100015, places=4)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(33.94, 18.43).value,
            0, places=5)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(22.9, -43.23).value,
            0, places=5)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(25.78, -80.22).value,
            0.00861727999508758, places=4)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(28.717, 77.3).value,
            0.209383698952704, places=4)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(3.133, 101.7).value,
            0.0512514559528945, places=4)
        self.assertAlmostEqual(
            models.itu1511.topographic_altitude(9.05, 38.7).value,
            2.5398618775, places=4)


if __name__ == '__main__':
    suite = suite()
    print('Validation tests for the ITU-R models')
    print('------------------------')
    print(
        'A total of %d test-cases are going to be tested' %
        suite.countTestCases())
    sys.stdout.flush()
    test.TextTestRunner(verbosity=2).run(suite)
