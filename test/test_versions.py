# -*- coding: utf-8 -*-
import unittest as test

import itur
from itur import utils
import itur.models as models

import sys
from astropy import units as u


def suite():
    """ A test suite for the ITU-P Recommendations. Recommendations tested:
    """
    suite = test.TestSuite()

    # Test valid versions
    suite.addTest(TestVersions('change_versions'))

    # For each version test all functions for vectorization and for
    suite.addTest(TestFunctionsRecommendation453('test_453'))
    suite.addTest(TestFunctionsRecommendation618('test_618'))
    suite.addTest(TestFunctionsRecommendation676('test_676'))
    suite.addTest(TestFunctionsRecommendation835('test_835'))
    suite.addTest(TestFunctionsRecommendation836('test_836'))
    suite.addTest(TestFunctionsRecommendation837('test_837'))
    suite.addTest(TestFunctionsRecommendation838('test_838'))
    suite.addTest(TestFunctionsRecommendation839('test_839'))
    suite.addTest(TestFunctionsRecommendation840('test_840'))
    suite.addTest(TestFunctionsRecommendation1510('test_1510'))
    suite.addTest(TestFunctionsRecommendation1511('test_1511'))

    return suite


class TestVersions(test.TestCase):

    def change_versions(self):

        # For
        for i in range(1, 12):
            self.assertRaises(ValueError,
                              models.itu453.change_version, i)

        for i in range(1, 16):
            self.assertRaises(ValueError,
                              models.itu530.change_version, i)

        for i in range(1, 12):
            self.assertRaises(ValueError,
                              models.itu618.change_version, i)

        for i in range(1, 2):
            self.assertRaises(ValueError,
                              models.itu839.change_version, i)

        for i in range(1, 9):
            self.assertRaises(ValueError,
                              models.itu676.change_version, i)

        for i in range(1, 5):
            self.assertRaises(ValueError,
                              models.itu835.change_version, i)

        for i in range(1, 12):
            self.assertRaises(ValueError,
                              models.itu453.change_version, i)

        for i in range(1, 5):
            self.assertRaises(ValueError,
                              models.itu835.change_version, i)

        for i in range(1, 4):
            self.assertRaises(ValueError,
                              models.itu836.change_version, i)

        for i in range(1, 6):
            self.assertRaises(ValueError,
                              models.itu837.change_version, i)

        for i in range(1, 4):
            self.assertRaises(ValueError,
                              models.itu840.change_version, i)


class TestFunctionsRecommendation453(test.TestCase):
    def setUp(self):
        self.versions = [12, 13]

    def test_all_functions_453(self):

        T = 15 * itur.u.deg_C
        e = (T.value * 7.5 / 216.7) * u.hPa
        Pd = 1013.15 * itur.u.hPa
        P = 1013.15 * itur.u.hPa
        H = 60 * itur.u.percent
        lat = 51
        lon = -53
        p = 0.51
        p_exact = 0.5

        models.itu453.wet_term_radio_refractivity(e, T)
        models.itu453.wet_term_radio_refractivity([e, e], [T, T])
        models.itu453.wet_term_radio_refractivity(e, [T, T])

        models.itu453.radio_refractive_index(P, e, T)
        models.itu453.radio_refractive_index([P, P], [e, e], [T, T])
        models.itu453.radio_refractive_index([P, P], [e, e], T)
        models.itu453.radio_refractive_index([P, P], e, [T, T])
        models.itu453.radio_refractive_index(P, [e, e], [T, T])

        models.itu453.dry_term_radio_refractivity(Pd, T)
        models.itu453.dry_term_radio_refractivity([Pd, Pd], [T, T])
        models.itu453.dry_term_radio_refractivity([Pd, Pd], [T, T])

        models.itu453.water_vapour_pressure(T, P, H, type_hydrometeor='water')
        models.itu453.water_vapour_pressure(T, P, H, type_hydrometeor='ice')
        models.itu453.water_vapour_pressure([T, T], [P, P], [H, H],
                                            type_hydrometeor='water')
        models.itu453.water_vapour_pressure([T, T], [P, P], [H, H],
                                            type_hydrometeor='ice')
        models.itu453.water_vapour_pressure([T, T], P, [H, H],
                                            type_hydrometeor='water')
        models.itu453.water_vapour_pressure([T, T], P, [H, H],
                                            type_hydrometeor='ice')
        models.itu453.water_vapour_pressure([T, T], [P, P], H,
                                            type_hydrometeor='water')
        models.itu453.water_vapour_pressure([T, T], [P, P], H,
                                            type_hydrometeor='ice')

        models.itu453.saturation_vapour_pressure(T, P,
                                                 type_hydrometeor='water')
        models.itu453.saturation_vapour_pressure(T, P,
                                                 type_hydrometeor='ice')
        models.itu453.saturation_vapour_pressure([T, T], [P, P],
                                                 type_hydrometeor='water')
        models.itu453.saturation_vapour_pressure([T, T], [P, P],
                                                 type_hydrometeor='ice')
        models.itu453.saturation_vapour_pressure([T, T], P,
                                                 type_hydrometeor='water')
        models.itu453.saturation_vapour_pressure([T, T], P,
                                                 type_hydrometeor='ice')
        models.itu453.saturation_vapour_pressure(T, [P, P],
                                                 type_hydrometeor='water')
        models.itu453.saturation_vapour_pressure(T, [P, P],
                                                 type_hydrometeor='ice')

        models.itu453.map_wet_term_radio_refractivity(lat, lon, p).value
        models.itu453.map_wet_term_radio_refractivity(
                [lat, lat], [lon, lon], [p, p])
        models.itu453.map_wet_term_radio_refractivity(lat, lon, [p, p])

        models.itu453.DN65(lat, lon, p_exact)
        models.itu453.DN65([lat, lat], [lon, lon], [p_exact, p_exact])
        models.itu453.DN65(lat, lon, [p_exact, p_exact])

        models.itu453.DN1(lat, lon, p_exact)
        models.itu453.DN1([lat, lat], [lon, lon], [p_exact, p_exact])
        models.itu453.DN1(lat, lon, [p_exact, p_exact])

    def test_453(self):

        for version in self.versions:
            utils.memory.clear()
            models.itu453.change_version(version)
            self.test_all_functions_453()
            self.assertEqual(models.itu453.get_version(), version)


class TestFunctionsRecommendation618(test.TestCase):
    def setUp(self):
        self.versions = [12, 13]

    def test_all_functions_618(self):

        f = 29 * itur.u.GHz
        el = 31
        tau = 45
        Ls = 3 * itur.u.km
        hs = 0.05 * itur.u.km
        D = 1.2 * itur.u.m
        R001 = 34 * itur.u.mm / itur.u.hr
        lat = 51
        lon = -53
        p = 0.51

        a1 = 10
        a2 = 12
        lat1 = 51
        lon1 = -53
        lat2 = 52
        lon2 = -53
        el1 = 30
        el2 = 54

        Ap = 10

        P_k = P0 = 0.5

        models.itu618.rain_attenuation(lat, lon, f, el, p=p)
        models.itu618.rain_attenuation([lat, lat], [lon, lon], f, el, p=p)
        models.itu618.rain_attenuation([lat, lat], [lon, lon], [f, f],
                                       [el, el], p=p)
        models.itu618.rain_attenuation([lat, lat], [lon, lon], [f, f],
                                       [el, el], p=[p, p])

        models.itu618.rain_attenuation(lat, lon, f, el, hs=hs, p=p, R001=R001,
                                       tau=tau, Ls=Ls)
        models.itu618.rain_attenuation([lat, lat], [lon, lon], f, el, hs=hs,
                                       p=p, R001=R001, tau=tau, Ls=Ls)
        models.itu618.rain_attenuation([lat, lat], [lon, lon], [f, f],
                                       [el, el], hs=hs, p=p, R001=R001,
                                       tau=tau, Ls=Ls)
        models.itu618.rain_attenuation([lat, lat], [lon, lon], [f, f],
                                       [el, el], hs=[hs, hs], p=[p, p],
                                       R001=[R001, R001], Ls=[Ls, Ls],
                                       tau=[tau, tau])

        models.itu618.rain_attenuation_probability(lat, lon, el)
        models.itu618.rain_attenuation_probability([lat, lat], [lon, lon], el)
        models.itu618.rain_attenuation_probability(
            [lat, lat], [lon, lon], [el, el])

        models.itu618.rain_attenuation_probability(
            lat, lon, el, hs=hs, Ls=Ls, P0=P0)
        models.itu618.rain_attenuation_probability([lat, lat], [lon, lon], el,
                                                   hs=hs, Ls=Ls, P0=P0)
        models.itu618.rain_attenuation_probability(
            [lat, lat], [lon, lon], [el, el], hs=[hs, hs], Ls=[Ls, Ls],
            P0=[P0, P0])

        models.itu618.site_diversity_rain_outage_probability(
            lat1, lon1, a1, el1, lat2, lon2, a2, el2, f,
            tau=45, hs1=None, hs2=None)

        models.itu618.scintillation_attenuation(lat, lon, f, el, p, D)
        models.itu618.scintillation_attenuation([lat, lat], [lon, lon], [f, f],
                                                [el, el], p, D)
        models.itu618.scintillation_attenuation([lat, lat], [lon, lon], [f, f],
                                                [el, el], [p, p], [D, D])

        models.itu618.rain_cross_polarization_discrimination(
            Ap, f, el, p, tau=45)
        models.itu618.rain_cross_polarization_discrimination(
            [Ap, Ap], [f, f], [el, el], [p, p], tau=45)
        models.itu618.rain_cross_polarization_discrimination(
            [Ap, Ap], [f, f], [el, el], [p, p], tau=tau)
        models.itu618.rain_cross_polarization_discrimination(
            [Ap, Ap], [f, f], [el, el], [p, p], tau=[tau, tau])

        models.itu618.fit_rain_attenuation_to_lognormal(
            lat, lon, f, el, hs, P_k, tau)
        models.itu618.fit_rain_attenuation_to_lognormal(
            [lat, lat], [lon, lon], [f, f], [el, el], hs, P_k, tau)
        models.itu618.fit_rain_attenuation_to_lognormal(
            [lat, lat], [lon, lon], [f, f], [el, el], [hs, hs], [P_k, P_k],
            [tau, tau])

    def test_618(self):

        for version in self.versions:
            utils.memory.clear()
            models.itu618.change_version(version)
            self.test_all_functions_618()
            self.assertEqual(models.itu618.get_version(), version)


class TestFunctionsRecommendation676(test.TestCase):
    def setUp(self):
        self.versions = [9, 10, 11]

    def test_all_functions_676(self):

        r = 5 * itur.u.km
        f = 29 * itur.u.GHz
        el = 31
        rho = 7.5
        P = 1013 * itur.u.hPa
        T = 15 * itur.u.deg_C
        V_t = 20 * itur.u.kg / itur.u.m**2
        h = 0.05 * itur.u.km
        lat = 51
        lon = -53
        p = 0.51

        models.itu676.gaseous_attenuation_terrestrial_path(
            r, f, el, rho, P, T, 'approx')
        models.itu676.gaseous_attenuation_terrestrial_path(
            r, f, el, rho, P, T, 'exact')
        models.itu676.gaseous_attenuation_terrestrial_path(
            r, [f, f], [el, el], [rho, rho], [P, P], [T, T], 'approx')
        models.itu676.gaseous_attenuation_terrestrial_path(
            r, [f, f], [el, el], [rho, rho], [P, P], [T, T], 'exact')

        models.itu676.gaseous_attenuation_slant_path(
            f, el, rho, P, T, V_t=None, h=None, mode='approx')
        models.itu676.gaseous_attenuation_slant_path(
            f, el, rho, P, T, V_t=None, h=None, mode='exact')
        models.itu676.gaseous_attenuation_slant_path(
            [f, f], [el, el], [rho, rho], [P, P], [T, T],
            V_t=None, h=None, mode='approx')
        models.itu676.gaseous_attenuation_slant_path(
            [f, f], [el, el], [rho, rho], [P, P], [T, T],
            V_t=None, h=None, mode='exact')
        models.itu676.gaseous_attenuation_slant_path(
            f, el, rho, P, T, V_t=V_t, h=h, mode='approx')
        models.itu676.gaseous_attenuation_slant_path(
            f, el, rho, P, T, V_t=V_t, h=h, mode='exact')
        models.itu676.gaseous_attenuation_slant_path(
            [f, f], [el, el], [rho, rho], [P, P], [T, T],
            V_t=[V_t, V_t], h=[h, h], mode='approx')
        models.itu676.gaseous_attenuation_slant_path(
            [f, f], [el, el], [rho, rho], [P, P], [T, T],
            V_t=[V_t, V_t], h=[h, h], mode='exact')

        models.itu676.zenit_water_vapour_attenuation(
            lat, lon, p, f, V_t=None, h=None)
        models.itu676.zenit_water_vapour_attenuation(
            [lat, lat], [lon, lon], [p, p], [f, f], V_t=None, h=None)
        models.itu676.zenit_water_vapour_attenuation(
            lat, lon, p, f, V_t=V_t, h=h)
        models.itu676.zenit_water_vapour_attenuation(
            [lat, lat], [lon, lon], [p, p], [f, f], V_t=[V_t, V_t], h=[h, h])

        models.itu676.gammaw_approx(f, P, rho, T)
        models.itu676.gammaw_approx([f, f], [P, P], [rho, rho], [T, T])

        models.itu676.gamma0_approx(f, P, rho, T)
        models.itu676.gamma0_approx([f, f], [P, P], [rho, rho], [T, T])

        if models.itu676.get_version() == 11:
            models.itu676.gamma0_exact(f, P, rho, T)
            models.itu676.gamma0_exact([f, f], [P, P], [rho, rho], [T, T])

            models.itu676.gammaw_exact(f, P, rho, T)
            models.itu676.gammaw_exact([f, f], [P, P], [rho, rho], [T, T])

    def test_676(self):

        for version in self.versions:
            utils.memory.clear()
            models.itu676.change_version(version)
            self.test_all_functions_676()
            self.assertEqual(models.itu676.get_version(), version)


class TestFunctionsRecommendation835(test.TestCase):
    def setUp(self):
        self.versions = [5, 6]

    def test_all_functions_835(self):

        T_0 = 15 * itur.u.deg_C
        h_0 = 2 * itur.u.km
        P_0 = 1013.25 * itur.u.hPa
        rho_0 = 7.5 * itur.u.g / itur.u.m**3
        h = 0.05 * itur.u.km
        lat = 51

        models.itu835.standard_water_vapour_pressure(h, h_0, rho_0)
        models.itu835.standard_water_vapour_density([h, h], h_0, rho_0)

        models.itu835.standard_pressure(h, T_0, P_0)
        models.itu835.standard_pressure([h, h], T_0, P_0)

        models.itu835.standard_temperature(h, T_0)
        models.itu835.standard_temperature([h, h], T_0)

        models.itu835.water_vapour_density(lat, h, season='summer')
        models.itu835.water_vapour_density(lat, h, season='winter')
        models.itu835.water_vapour_density([lat, lat], [h, h], season='summer')
        models.itu835.water_vapour_density([lat, lat], [h, h], season='winter')

        models.itu835.pressure(lat, h, season='summer')
        models.itu835.pressure(lat, h, season='winter')
        models.itu835.pressure([lat, lat], [h, h], season='summer')
        models.itu835.pressure([lat, lat], [h, h], season='winter')

        models.itu835.temperature(lat, h, season='summer')
        models.itu835.temperature(lat, h, season='winter')
        models.itu835.temperature([lat, lat], [h, h], season='summer')
        models.itu835.temperature([lat, lat], [h, h], season='winter')

    def test_835(self):

        for version in self.versions:
            utils.memory.clear()
            models.itu835.change_version(version)
            self.test_all_functions_835()
            self.assertEqual(models.itu835.get_version(), version)


class TestFunctionsRecommendation836(test.TestCase):
    def setUp(self):
        self.versions = [4, 5, 6]

    def test_all_functions_836(self):

        lat = 51
        lon = -63
        p = 0.51
        alt = 0.5

        models.itu836.surface_water_vapour_density(lat, lon, p)
        models.itu836.surface_water_vapour_density(
                [lat, lat], [lon, lon], p)
        models.itu836.surface_water_vapour_density(
                [lat, lat], [lon, lon], [p, p])

        models.itu836.surface_water_vapour_density(lat, lon, p, alt)
        models.itu836.surface_water_vapour_density(
                [lat, lat], [lon, lon], p, [alt, alt])
        models.itu836.surface_water_vapour_density(
                [lat, lat], [lon, lon], [p, p], [alt, alt])

        models.itu836.total_water_vapour_content(lat, lon, p)
        models.itu836.total_water_vapour_content(
                [lat, lat], [lon, lon], p)
        models.itu836.total_water_vapour_content(
                [lat, lat], [lon, lon], [p, p])

        models.itu836.total_water_vapour_content(lat, lon, p, alt)
        models.itu836.total_water_vapour_content(
                [lat, lat], [lon, lon], p, [alt, alt])
        models.itu836.total_water_vapour_content(
                [lat, lat], [lon, lon], [p, p], [alt, alt])

    def test_836(self):

        for version in self.versions:
            utils.memory.clear()
            models.itu836.change_version(version)
            self.test_all_functions_836()
            self.assertEqual(models.itu836.get_version(), version)


class TestFunctionsRecommendation837(test.TestCase):
    def setUp(self):
        self.versions = [6, 7]

    def test_all_functions_837(self):
        lat = 51
        lon = -63
        p = 0.51

        models.itu837.rainfall_probability(lat, lon)
        models.itu837.rainfall_probability([lat, lat], [lon, lon])

        models.itu837.rainfall_rate(lat, lon, p)
        models.itu837.rainfall_rate([lat, lat], [lon, lon], p)
        models.itu837.rainfall_rate([lat, lat], [lon, lon], [p, p])

    def test_837(self):

        for version in self.versions:
            utils.memory.clear()
            models.itu837.change_version(version)
            self.test_all_functions_837()
            self.assertEqual(models.itu837.get_version(), version)


class TestFunctionsRecommendation838(test.TestCase):
    def setUp(self):
        self.versions = [0, 1, 2, 3]

    def test_all_functions_838(self):

        f = 29 * itur.u.GHz
        el = 32
        tau = 45
        R = 10 * itur.u.mm / itur.u.hr

        models.itu838.rain_specific_attenuation_coefficients(f, el, tau)
        models.itu838.rain_specific_attenuation_coefficients(
                [f, f], [el, el], [tau, tau])
        models.itu838.rain_specific_attenuation_coefficients(
                f, [el, el], [tau, tau])

        models.itu838.rain_specific_attenuation(R, f, el, tau)
        models.itu838.rain_specific_attenuation(R, [f, f], el, [tau, tau])
        models.itu838.rain_specific_attenuation(
                R, [f, f], [el, el], [tau, tau])
        models.itu838.rain_specific_attenuation(
                [R, R], [f, f], [el, el], [tau, tau])

    def test_838(self):
        for version in self.versions:
            utils.memory.clear()
            models.itu838.change_version(version)
            self.test_all_functions_838()
            self.assertEqual(models.itu838.get_version(), version)


class TestFunctionsRecommendation839(test.TestCase):
    def setUp(self):
        self.versions = [2, 3, 4]

    def test_all_functions_839(self):
        lat = 51
        lon = -63

        models.itu839.isoterm_0(lat, lon)
        models.itu839.isoterm_0([lat, lat], [lon, lon])

        models.itu839.rain_height(lat, lon)
        models.itu839.rain_height([lat, lat], [lon, lon])

    def test_839(self):

        for version in self.versions:
            utils.memory.clear()
            models.itu839.change_version(version)
            self.test_all_functions_839()
            self.assertEqual(models.itu839.get_version(), version)


class TestFunctionsRecommendation840(test.TestCase):
    def setUp(self):
        self.versions = [4, 5, 6, 7]

    def test_all_functions_840(self):

        T = 15 * itur.u.deg_C
        f = 29 * itur.u.GHz
        el = 32
        p = 0.51
        lat = 51
        lon = -54

        models.itu840.specific_attenuation_coefficients(f, T)
        models.itu840.specific_attenuation_coefficients([f, f], [T, T])
        models.itu840.specific_attenuation_coefficients(f, [T, T])

        models.itu840.columnar_content_reduced_liquid(lat, lon, p)
        models.itu840.columnar_content_reduced_liquid(
                [lat, lat], [lon, lon], p)
        models.itu840.columnar_content_reduced_liquid(
                [lat, lat], [lon, lon], [p, p])

        models.itu840.cloud_attenuation(lat, lon, el, f, p)
        models.itu840.cloud_attenuation(
                [lat, lat], [lon, lon], [el, el], [f, f], [p, p])
        models.itu840.cloud_attenuation(
                [lat, lat], [lon, lon], el, [f, f], [p, p])
        models.itu840.cloud_attenuation(
                [lat, lat], [lon, lon], el, f, [p, p])

        models.itu840.lognormal_approximation_coefficient(lat, lon)
        models.itu840.lognormal_approximation_coefficient(
                [lat, lat], [lon, lon])

    def test_840(self):

        for version in self.versions:
            utils.memory.clear()
            models.itu840.change_version(version)
            self.test_all_functions_840()
            self.assertEqual(models.itu840.get_version(), version)


class TestFunctionsRecommendation1510(test.TestCase):
    def setUp(self):
        self.versions = [0, 1]

    def test_all_functions_1510(self):
        lat = 51
        lon = -63

        models.itu1510.surface_mean_temperature(lat, lon)
        models.itu1510.surface_mean_temperature([lat, lat], [lon, lon])

        if models.itu1510.get_version() == 1:
            for m in range(1, 13):
                models.itu1510.surface_month_mean_temperature(lat, lon, m)
                models.itu1510.surface_month_mean_temperature(
                        [lat, lat], [lon, lon], m)
                models.itu1510.surface_month_mean_temperature(
                        [lat, lat], [lon, lon], [m, m])

    def test_1510(self):
        for version in self.versions:
            utils.memory.clear()
            models.itu1510.change_version(version)
            self.test_all_functions_1510()
            self.assertEqual(models.itu1510.get_version(), version)


class TestFunctionsRecommendation1511(test.TestCase):
    def setUp(self):
        self.versions = [0, 1]

    def test_all_functions_1511(self):
        lat = 51
        lon = -63

        models.itu1511.topographic_altitude(lat, lon)
        models.itu1511.topographic_altitude([lat, lat], [lon, lon])

    def test_1511(self):

        for version in self.versions:
            utils.memory.clear()
            models.itu1511.change_version(version)
            self.test_all_functions_1511()
            self.assertEqual(models.itu1511.get_version(), version)

if __name__ == '__main__':
    pass
    suite = suite()
    print('Test versioning of the code')
    print('------------------------')
    print(
        'A total of %d test-cases are going to be tested' %
        suite.countTestCases())
    sys.stdout.flush()
    test.TextTestRunner(verbosity=2).run(suite)
