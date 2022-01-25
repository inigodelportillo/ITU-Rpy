# -*- coding: utf-8 -*-
import os
import sys
import warnings
import numpy as np
import unittest as test

import itur
import itur.models as models

from astropy import units as u

basepath = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(basepath, 'test_data')


def suite():
    """ A test suite for the ITU-P Recommendations. Recommendations tested:
    """
    suite = test.TestSuite()

    # Test valid versions
    suite.addTest(TestVersions('test_change_to_not_implemented_versions'))

    # For each version test all functions for vectorization and for
    suite.addTest(TestFunctionsRecommendation453('test_453'))
    suite.addTest(TestFunctionsRecommendation530('test_530'))
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
    suite.addTest(TestFunctionsRecommendation1623('test_1623'))
    suite.addTest(TestFunctionsRecommendation1853('test_1853'))

    # Basic import module functionality
    suite.addTest(TestImportModules('test_import_itur'))
    suite.addTest(TestImportModules('test_import_itur_utils'))
    suite.addTest(TestImportModules('test_import_itur_plotting'))
    suite.addTest(TestImportModules('test_import_itur_models'))
    suite.addTest(TestImportModules('test_import_itur_models_1853'))
    suite.addTest(TestImportModules('test_import_itur_models_1623'))
    suite.addTest(TestImportModules('test_import_itur_models_1511'))
    suite.addTest(TestImportModules('test_import_itur_models_1510'))
    suite.addTest(TestImportModules('test_import_itur_models_453'))
    suite.addTest(TestImportModules('test_import_itur_models_530'))
    suite.addTest(TestImportModules('test_import_itur_models_618'))
    suite.addTest(TestImportModules('test_import_itur_models_676'))
    suite.addTest(TestImportModules('test_import_itur_models_835'))
    suite.addTest(TestImportModules('test_import_itur_models_836'))
    suite.addTest(TestImportModules('test_import_itur_models_837'))
    suite.addTest(TestImportModules('test_import_itur_models_838'))
    suite.addTest(TestImportModules('test_import_itur_models_839'))
    suite.addTest(TestImportModules('test_import_itur_models_840'))

    # Test version of itur
    suite.addTest(TestImportModules('test_version'))

    # Test slant_path_attenuation calls
    suite.addTest(TestIturMainFunctions('test_slant_path_attenuation'))
    suite.addTest(TestIturMainFunctions('test_slant_path_attenuation_p_below'))
    suite.addTest(TestIturMainFunctions('test_slant_path_attenuation_p_above'))
    suite.addTest(TestIturMainFunctions(
        'test_slant_path_attenuation_without_rain'))
    suite.addTest(TestIturMainFunctions(
        'test_slant_path_attenuation_without_gas'))
    suite.addTest(TestIturMainFunctions(
        'test_slant_path_attenuation_without_clouds'))
    suite.addTest(TestIturMainFunctions(
        'test_slant_path_attenuation_without_scintillation'))

    # Test utils library
    suite.addTest(TestIturUtils('test_read_file_npz'))
    suite.addTest(TestIturUtils('test_read_file_npy'))
    suite.addTest(TestIturUtils('test_read_file_txt'))
    suite.addTest(TestIturUtils('test_read_file_txt_as_txt'))
    suite.addTest(TestIturUtils('test_distance_wsg84'))
    suite.addTest(TestIturUtils('test_distance_haversine'))
    suite.addTest(TestIturUtils('test_prepare_quantity'))
    suite.addTest(TestIturUtils('test_prepare_output_array'))
    suite.addTest(TestIturUtils('test_regular_lat_lon_grid'))

    return suite


class TestVersions(test.TestCase):

    def test_change_to_not_implemented_versions(self):

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


class TestImportModules(test.TestCase):
    """ Tests that all submodules are imporable.
    """

    @staticmethod
    def test_version():
        import itur as itu
        print(itu.__version__)

    @staticmethod
    def test_import_itur():
        import itur

    @staticmethod
    def test_import_itur_utils():
        import itur.utils

    @staticmethod
    def test_import_itur_plotting():
        import itur.plotting

    @staticmethod
    def test_import_itur_models():
        import itur.models

    @staticmethod
    def test_import_itur_models_1510():
        import itur.models.itu1510

    @staticmethod
    def test_import_itur_models_1511():
        import itur.models.itu1511

    @staticmethod
    def test_import_itur_models_1623():
        import itur.models.itu1853

    @staticmethod
    def test_import_itur_models_1853():
        import itur.models.itu1853

    @staticmethod
    def test_import_itur_models_453():
        import itur.models.itu453

    @staticmethod
    def test_import_itur_models_530():
        import itur.models.itu530

    @staticmethod
    def test_import_itur_models_618():
        import itur.models.itu618

    @staticmethod
    def test_import_itur_models_676():
        import itur.models.itu676

    @staticmethod
    def test_import_itur_models_835():
        import itur.models.itu835

    @staticmethod
    def test_import_itur_models_836():
        import itur.models.itu836

    @staticmethod
    def test_import_itur_models_837():
        import itur.models.itu837

    @staticmethod
    def test_import_itur_models_838():
        import itur.models.itu838

    @staticmethod
    def test_import_itur_models_839():
        import itur.models.itu839

    @staticmethod
    def test_import_itur_models_840():
        import itur.models.itu840


class TestIturMainFunctions(test.TestCase):

    def setUp(self):
        self.lat = 0
        self.lon = 0
        self.f = 22 * u.GHz
        self.el = 45
        self.p = 0.01
        self.D = 1

    def test_slant_path_attenuation(self):
        itur.atmospheric_attenuation_slant_path(
            lat=self.lat, lon=self.lon, f=self.f, el=self.el, p=self.p,
            D=self.D)

    @test.skipIf(sys.version_info[0] < 3, "Only supported in Python 3+")
    def test_slant_path_attenuation_p_below(self):
        with self.assertWarns(RuntimeWarning):
            itur.atmospheric_attenuation_slant_path(
                lat=self.lat, lon=self.lon, f=self.f, el=self.el, p=1e-4,
                D=self.D)

    @test.skipIf(sys.version_info[0] < 3, "Only supported in Python 3+")
    def test_slant_path_attenuation_p_above(self):
        with self.assertWarns(RuntimeWarning):
            itur.atmospheric_attenuation_slant_path(
                lat=self.lat, lon=self.lon, f=self.f, el=self.el, p=90,
                D=self.D)

    def test_slant_path_attenuation_without_rain(self):
        itur.atmospheric_attenuation_slant_path(
            lat=self.lat, lon=self.lon, f=self.f, el=self.el, p=self.p,
            D=self.D, include_rain=False)

    def test_slant_path_attenuation_without_clouds(self):
        itur.atmospheric_attenuation_slant_path(
            lat=self.lat, lon=self.lon, f=self.f, el=self.el, p=self.p,
            D=self.D, include_clouds=False)

    def test_slant_path_attenuation_without_gas(self):
        itur.atmospheric_attenuation_slant_path(
            lat=self.lat, lon=self.lon, f=self.f, el=self.el, p=self.p,
            D=self.D, include_gas=False)

    def test_slant_path_attenuation_without_scintillation(self):
        itur.atmospheric_attenuation_slant_path(
            lat=self.lat, lon=self.lon, f=self.f, el=self.el, p=self.p,
            D=self.D, include_scintillation=False)


class TestIturUtils(test.TestCase):

    def test_read_file_npy(self):
        path = os.path.join(test_data, 'v3_esalat.npy')
        itur.utils.load_data(path)

    def test_read_file_npz(self):
        path = os.path.join(test_data, 'v3_esalat.npz')
        itur.utils.load_data(path)

    def test_read_file_txt(self):
        path = os.path.join(test_data, 'v12_lines_oxygen.txt')
        itur.utils.load_data(path, skip_header=1)

    def test_read_file_txt_as_txt(self):
        path = os.path.join(test_data, 'v12_lines_oxygen.txt')
        itur.utils.load_data(path, is_text=True)

    def test_distance_haversine(self):
        val = itur.utils.compute_distance_earth_to_earth_haversine(
            lat_p=0, lon_p=0, lat_grid=10, lon_grid=10)
        np.testing.assert_allclose(val, 1568.5205567985759)

        val = itur.utils.compute_distance_earth_to_earth_haversine(
            lat_p=0, lon_p=0, lat_grid=np.array([10, 20]),
            lon_grid=np.array([10, 20]))
        np.testing.assert_allclose(val, np.array([1568.5205567985759,
                                                  3112.445040079722]))

        val = itur.utils.compute_distance_earth_to_earth_haversine(
            lat_p=0, lon_p=0, lat_grid=np.array([[10], [20]]),
            lon_grid=np.array([[10], [20]]))
        np.testing.assert_allclose(val, np.array([[1568.5205567985759],
                                                  [3112.445040079722]]))

    def test_distance_wsg84(self):
        val = itur.utils.compute_distance_earth_to_earth(
            lat_p=0, lon_p=0, lat_grid=10, lon_grid=10,
            method='WGS84')
        self.assertEqual(val, 1565.10909921789)

        val = itur.utils.compute_distance_earth_to_earth(
            lat_p=0, lon_p=0, lat_grid=np.array([10, 20]),
            lon_grid=np.array([10, 20]),
            method='WGS84')
        np.testing.assert_allclose(val, np.array([1565.10909922,
                                                  3106.12677679]))

        val = itur.utils.compute_distance_earth_to_earth(
            lat_p=0, lon_p=0, lat_grid=np.array([[10], [20]]),
            lon_grid=np.array([[10], [20]]),
            method='WGS84')
        np.testing.assert_allclose(val, np.array([[1565.10909922],
                                                  [3106.12677679]]))

    def test_prepare_quantity(self):
        # Test temperature conversion
        val = itur.utils.prepare_quantity((273.15, 373.15) * itur.u.K,
                                          units=itur.u.Celsius)
        np.testing.assert_array_equal(val, np.array([0, 100]))

        # Test individual numbers
        val = itur.utils.prepare_quantity(1, units=itur.u.m)
        self.assertEqual(val, 1)

        val = itur.utils.prepare_quantity(1 * itur.u.km, units=itur.u.m)
        self.assertEqual(val, 1000)

        val = itur.utils.prepare_quantity(None, units=itur.u.m)
        self.assertEqual(val, None)

        # Test tuples of values
        val = itur.utils.prepare_quantity((1, 2), units=itur.u.m)
        np.testing.assert_array_equal(val, np.array([1, 2]))

        val = itur.utils.prepare_quantity((1, 2) * itur.u.km, units=itur.u.m)
        np.testing.assert_array_equal(val, np.array([1000, 2000]))

        # Test numpy arrays
        val = itur.utils.prepare_quantity(np.array([1, 2]), units=itur.u.m)
        np.testing.assert_array_equal(val, np.array([1, 2]))

        val = itur.utils.prepare_quantity(np.array([1, 2]) * itur.u.km,
                                          units=itur.u.m)
        np.testing.assert_array_equal(val, np.array([1000, 2000]))

        # Test lists of values
        val = itur.utils.prepare_quantity([1, 2], units=itur.u.m)
        np.testing.assert_array_equal(val, np.array([1, 2]))

        val = itur.utils.prepare_quantity([1, 2] * itur.u.km, units=itur.u.m)
        np.testing.assert_array_equal(val, np.array([1000, 2000]))

        # Check that invalid formats raise an exception
        with self.assertRaises(ValueError):
            itur.utils.prepare_quantity({}, units=itur.u.m)

    def test_prepare_output_array(self):
        out_array = np.array([[1, 2], [3, 4]])

        # Check values without units
        val = itur.utils.prepare_output_array(out_array,
                                              type_input=list)
        self.assertEqual(val, out_array.tolist())

        val = itur.utils.prepare_output_array(out_array,
                                              type_input=np.ndarray)
        np.testing.assert_array_equal(val, out_array)

        val = itur.utils.prepare_output_array(5, type_input=float)
        self.assertEqual(val, 5.0)

        val = itur.utils.prepare_output_array([5, 10], type_input=list)
        self.assertEqual(val, [5, 10])

        # Check values with units
        val = itur.utils.prepare_output_array(5 * itur.u.m, type_input=float)
        self.assertEqual(val, 5.0 * itur.u.m)

        val = itur.utils.prepare_output_array([5, 10] * itur.u.m,
                                              type_input=list)
        np.testing.assert_array_equal(val, [5, 10] * itur.u.m)

        val = itur.utils.prepare_output_array(out_array * itur.u.m,
                                              type_input=list)
        np.testing.assert_array_equal(val, out_array.tolist() * itur.u.m)

        val = itur.utils.prepare_output_array(out_array * itur.u.m,
                                              type_input=np.ndarray)
        np.testing.assert_array_equal(val, out_array.tolist() * itur.u.m)

    def test_regular_lat_lon_grid(self):
        itur.utils.regular_lat_lon_grid(lon_start_0=True)
        itur.utils.regular_lat_lon_grid(lon_start_0=False)


class TestFunctionsRecommendation453(test.TestCase):
    def setUp(self):
        self.versions = [12, 13]

    @staticmethod
    def test_all_functions_453():

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
            models.itu453.change_version(version)
            self.test_all_functions_453()
            self.assertEqual(models.itu453.get_version(), version)


class TestFunctionsRecommendation530(test.TestCase):
    def setUp(self):
        self.versions = [16, 17]

    @staticmethod
    def test_all_functions_530():

        d1 = d2 = d = 10 * itur.u.km
        f = 29 * itur.u.GHz
        h = h_e = h_r = 100 * itur.u.m
        A = Ap = 10 * itur.u.dB
        el = 45
        XPD_g = 20 * itur.u.dB
        C0_I = 20 * itur.u.dB
        lat = 51
        lon = -53
        p = 0.05

        models.itu530.fresnel_ellipse_radius(d1, d2, f)
        models.itu530.diffraction_loss(d1, d2, h, f)
        models.itu530.multipath_loss_for_A(lat, lon, h_e, h_r, d, f, A)
        models.itu530.multipath_loss(lat, lon, h_e, h_r, d, f, A)
        models.itu530.rain_attenuation(lat, lon, d, f, el, p)
        models.itu530.inverse_rain_attenuation(lat, lon, d, f, el, Ap)

        models.itu530.rain_event_count(lat, lon, d, f, el, A)
        models.itu530.XPD_outage_clear_air(
            lat, lon, h_e, h_r, d, f, XPD_g, C0_I)
        models.itu530.XPD_outage_precipitation(lat, lon, d, f, el, C0_I)

    def test_530(self):

        for version in self.versions:
            models.itu530.change_version(version)
            self.test_all_functions_530()
            self.assertEqual(models.itu530.get_version(), version)


class TestFunctionsRecommendation618(test.TestCase):
    def setUp(self):
        self.versions = [12, 13]

    @staticmethod
    def test_all_functions_618():

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
            models.itu618.change_version(version)
            self.test_all_functions_618()
            self.assertEqual(models.itu618.get_version(), version)


class TestFunctionsRecommendation676(test.TestCase):
    def setUp(self):
        self.versions = [9, 10, 11, 12]

    @staticmethod
    def test_all_functions_676():

        r = 5 * itur.u.km
        f = 29 * itur.u.GHz
        el = 71
        el_low = 4
        rho = 7.5
        P = 1013 * itur.u.hPa
        T = 15 * itur.u.deg_C
        V_t = 20 * itur.u.kg / itur.u.m**2
        h = 0.05 * itur.u.km
        h1 = 0.05 * itur.u.km
        h2 = 0.15 * itur.u.km
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

        models.itu676.gaseous_attenuation_inclined_path(
            f, el, rho, P, T, h1=h1, h2=h2, mode='approx')
        models.itu676.gaseous_attenuation_inclined_path(
            f, el, rho, P, T, h1=h1, h2=h2, mode='exact')
        models.itu676.gaseous_attenuation_inclined_path(
            [f, f], [el, el], [rho, rho], [P, P], [T, T],
            h1=h1, h2=h2, mode='approx')
        models.itu676.gaseous_attenuation_inclined_path(
            [f, f], [el, el], [rho, rho], [P, P], [T, T],
            h1=h1, h2=h2, mode='exact')
        models.itu676.gaseous_attenuation_inclined_path(
            f, el, rho, P, T, h1=h1, h2=h2, mode='approx')
        models.itu676.gaseous_attenuation_inclined_path(
            f, el, rho, P, T, h1=h1, h2=h2, mode='exact')

        with warnings.catch_warnings(record=True) as w:
            models.itu676.gaseous_attenuation_inclined_path(
                f, el_low, rho, P, T, h1=h1, h2=h2, mode='approx')
            models.itu676.gaseous_attenuation_inclined_path(
                f, el_low, rho, P, T, h1=h1, h2=h2, mode='exact')
            models.itu676.gaseous_attenuation_inclined_path(
                [f, f], [el_low, el_low], [rho, rho], [P, P], [T, T],
                h1=h1, h2=h2, mode='approx')
            models.itu676.gaseous_attenuation_inclined_path(
                [f, f], [el_low, el_low], [rho, rho], [P, P], [T, T],
                h1=h1, h2=h2, mode='exact')
            models.itu676.gaseous_attenuation_inclined_path(
                f, el_low, rho, P, T, h1=h1, h2=h2, mode='approx')
            models.itu676.gaseous_attenuation_inclined_path(
                f, el_low, rho, P, T, h1=h1, h2=h2, mode='exact')

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
            models.itu676.change_version(version)
            self.test_all_functions_676()
            self.assertEqual(models.itu676.get_version(), version)


class TestFunctionsRecommendation835(test.TestCase):
    def setUp(self):
        self.versions = [5, 6]

    @staticmethod
    def test_all_functions_835():

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
            models.itu835.change_version(version)
            self.test_all_functions_835()
            self.assertEqual(models.itu835.get_version(), version)


class TestFunctionsRecommendation836(test.TestCase):
    def setUp(self):
        self.versions = [4, 5, 6]

    @staticmethod
    def test_all_functions_836():

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
            models.itu836.change_version(version)
            self.test_all_functions_836()
            self.assertEqual(models.itu836.get_version(), version)


class TestFunctionsRecommendation837(test.TestCase):
    def setUp(self):
        self.versions = [6, 7]

    @staticmethod
    def test_all_functions_837():
        lat = 51
        lon = -63
        p = 0.51
        R = 10

        models.itu837.rainfall_probability(lat, lon)
        models.itu837.rainfall_probability([lat, lat], [lon, lon])

        models.itu837.rainfall_rate(lat, lon, p)
        models.itu837.rainfall_rate([lat, lat], [lon, lon], p)
        models.itu837.rainfall_rate([lat, lat], [lon, lon], [p, p])

        models.itu837.unavailability_from_rainfall_rate(lat, lon, R)

    def test_837(self):

        for version in self.versions:
            models.itu837.change_version(version)
            self.test_all_functions_837()
            self.assertEqual(models.itu837.get_version(), version)


class TestFunctionsRecommendation838(test.TestCase):
    def setUp(self):
        self.versions = [0, 1, 2, 3]

    @staticmethod
    def test_all_functions_838():

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
            models.itu838.change_version(version)
            self.test_all_functions_838()
            self.assertEqual(models.itu838.get_version(), version)


class TestFunctionsRecommendation839(test.TestCase):
    def setUp(self):
        self.versions = [2, 3, 4]

    @staticmethod
    def test_all_functions_839():
        lat = 51
        lon = -63

        models.itu839.isoterm_0(lat, lon)
        models.itu839.isoterm_0([lat, lat], [lon, lon])

        models.itu839.rain_height(lat, lon)
        models.itu839.rain_height([lat, lat], [lon, lon])

    def test_839(self):

        for version in self.versions:
            models.itu839.change_version(version)
            self.test_all_functions_839()
            self.assertEqual(models.itu839.get_version(), version)


class TestFunctionsRecommendation840(test.TestCase):
    def setUp(self):
        self.versions = [4, 5, 6, 7, 8]

    @staticmethod
    def test_all_functions_840():

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
            models.itu840.change_version(version)
            self.test_all_functions_840()
            self.assertEqual(models.itu840.get_version(), version)


class TestFunctionsRecommendation1510(test.TestCase):
    def setUp(self):
        self.versions = [0, 1]

    @staticmethod
    def test_all_functions_1510():
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
            models.itu1510.change_version(version)
            self.test_all_functions_1510()
            self.assertEqual(models.itu1510.get_version(), version)


class TestFunctionsRecommendation1511(test.TestCase):
    def setUp(self):
        self.versions = [0, 1, 2]

    @staticmethod
    def test_all_functions_1511():
        lat = 51
        lon = -63

        models.itu1511.topographic_altitude(lat, lon)
        models.itu1511.topographic_altitude([lat, lat], [lon, lon])

    def test_1511(self):

        for version in self.versions:
            models.itu1511.change_version(version)
            self.test_all_functions_1511()
            self.assertEqual(models.itu1511.get_version(), version)


class TestFunctionsRecommendation1623(test.TestCase):
    def setUp(self):
        self.versions = [0, 1]

    @staticmethod
    def test_all_functions_1623():
        T_tot = 0.224 * 365.25 * 24 * 3600
        D = np.array([1, 10, 30, 60, 120, 180, 300, 600, 900, 1200, 1500,
                      1800, 2400, 3600])

        z = np.linspace(-2, 2, 100)
        A = 10
        f_B = 0.02
        delta_t = 1
        N_target = 25
        D_target = 60
        PofA = np.array([50, 30, 20, 10, 5, 3, 2, 1, .5, .3, .2, .1, .05, .03,
                         .02, .01, .005, .003, .002, .001])
        A_arr = np.array([0.4, 0.6, 0.8, 1.8, 2.70, 3.5, 4.20, 5.7, 7.4, 9,
                          10.60, 14, 18.3, 22.3, 25.8, 32.6, 40.1, 46.1, 50.8,
                          58.8])
        el = 38.5
        f = 28

        models.itu1623.fade_duration(D, A, el, f, T_tot)
        models.itu1623.fade_slope(z, A, f_B, delta_t)
        models.itu1623.fade_depth(N_target, D_target, A_arr, PofA, el, f)

    def test_1623(self):
        for version in self.versions:
            models.itu1623.change_version(version)
            self.test_all_functions_1623()
            self.assertEqual(models.itu1623.get_version(), version)


class TestFunctionsRecommendation1853(test.TestCase):

    def setUp(self):
        self.versions = [0, 1]

    @staticmethod
    def test_all_functions_1853():
        lat = 51
        lon = -63
        f = 29 * itur.u.GHz
        el = 32
        p = 0.05
        D = 1 * itur.u.m
        hs = 100 * itur.u.m
        Ns = 60 * 60 * 24

        models.itu1853.set_seed(42)
        models.itu1853.scintillation_attenuation_synthesis(Ns=Ns)
        models.itu1853.rain_attenuation_synthesis(
            lat=lat, lon=lon, f=f, el=el, hs=hs, Ns=Ns)

        if models.itu1853.get_version() > 0:
            models.itu1853.cloud_liquid_water_synthesis(lat, lon, Ns)
            models.itu1853.integrated_water_vapour_synthesis(lat, lon, Ns)
            models.itu1853.total_attenuation_synthesis(
                lat, lon, f, el, p, D, Ns)
            models.itu1853.total_attenuation_synthesis(
                lat, lon, f, el, p, D, Ns, return_contributions=True)

    def test_1853(self):

        for version in self.versions:
            models.itu1853.change_version(version)
            self.test_all_functions_1853()
            self.assertEqual(models.itu1853.get_version(), version)


if __name__ == '__main__':

    # Create the test suite and run all the tests
    suite = suite()
    print('Test versioning of the code')
    print('------------------------')
    print(
        'A total of %d test-cases are going to be tested' %
        suite.countTestCases())
    sys.stdout.flush()
    test.TextTestRunner(verbosity=2).run(suite)
