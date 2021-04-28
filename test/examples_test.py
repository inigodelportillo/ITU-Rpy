# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('agg')

import itur.models.itu835 as itu835
import itur.models.itu676 as itu676
import itur
import unittest as test
import numpy as np
import sys
import matplotlib.pyplot as plt


def suite():
    """A test suite for the examples includes in itur. """
    suite = test.TestSuite()

    # Test valid versions
    suite.addTest(TestMapAfrica('test_map_africa'))
    suite.addTest(TestGaseousAttenuation('test_gaseous_attenuation'))
    suite.addTest(TestMultipleLocations('test_multiple_locations'))
    suite.addTest(TestSingleLocation('test_single_location'))
    suite.addTest(TestSingleLocationVsFrequency('test_single_location_vs_f'))
    suite.addTest(TestSingleLocationVsUnavailability(
        'test_single_location_vs_p'))

    return suite


class TestMapAfrica(test.TestCase):

    @staticmethod
    def test_map_africa():
        # Generate a regular grid of latitude and longitudes with 0.1
        # degree resolution for the region of interest.
        lat, lon = itur.utils.regular_lat_lon_grid(lat_max=60,
                                                   lat_min=-60,
                                                   lon_max=65,
                                                   lon_min=-35,
                                                   resolution_lon=1,
                                                   resolution_lat=1)

        # Satellite coordinates (GEO, 4 E)
        lat_sat = 0
        lon_sat = 4
        h_sat = 35786 * itur.u.km

        # Compute the elevation angle between satellite and ground stations
        el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat, lat, lon)

        # Set the link parameters
        f = 22.5 * itur.u.GHz    # Link frequency
        D = 1.2 * itur.u.m       # Antenna diameters
        p = 0.1                  # Unavailability (Vals exceeded 0.1% of time)

        # Compute the atmospheric attenuation
        Att = itur.atmospheric_attenuation_slant_path(lat, lon, f, el, p, D)

        # Now we show the surface mean temperature distribution
        T = itur.surface_mean_temperature(lat, lon)\
            .to(itur.u.Celsius, equivalencies=itur.u.temperature())

        # Plot the results
        try:
            m = itur.plotting.plot_in_map(Att, lat, lon,
                                          cbar_text='Atmospheric attenuation [dB]',
                                          cmap='magma')

            # Plot the satellite location
            m.scatter(lon_sat, lat_sat, c='white', s=20)

            m = itur.plotting.plot_in_map(
                T, lat, lon, cbar_text='Surface mean temperature [C]',
                cmap='RdBu_r')
        except RuntimeError as e:
            print(e)


class TestMultipleLocations(test.TestCase):

    @staticmethod
    def test_multiple_locations():
        # Obtain the coordinates of the different cities
        cities = {'Boston': (42.36, -71.06),
                  'New York': (40.71, -74.01),
                  'Los Angeles': (34.05, -118.24),
                  'Denver': (39.74, -104.99),
                  'Las Vegas': (36.20, -115.14),
                  'Seattle': (47.61, -122.33),
                  'Washington DC': (38.91, -77.04)}

        lat = [coords[0] for coords in cities.values()]
        lon = [coords[1] for coords in cities.values()]

        # Satellite coordinates (GEO, 4 E)
        lat_sat = 0
        lon_sat = -77
        h_sat = 35786 * itur.u.km

        # Compute the elevation angle between satellite and ground stations
        el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat, lat, lon)

        # Set the link parameters
        f = 22.5 * itur.u.GHz    # Link frequency
        D = 1.2 * itur.u.m       # Antenna diameters
        p = 0.1                  # Unavailability (Vals exceeded 0.1% of time)

        # Compute the atmospheric attenuation
        Ag, Ac, Ar, As, Att = itur.atmospheric_attenuation_slant_path(
            lat, lon, f, el, p, D, return_contributions=True)

        # Plot the results
        city_idx = np.arange(len(cities))
        width = 0.15

        fig, ax = plt.subplots(1, 1)
        ax.bar(city_idx, Att.value, 0.6, label='Total atmospheric Attenuation')
        ax.bar(city_idx - 1.5 * width, Ar.value, width,
               label='Rain attenuation')
        ax.bar(city_idx - 0.5 * width, Ag.value, width,
               label='Gaseous attenuation')
        ax.bar(city_idx + 0.5 * width, Ac.value, width,
               label='Clouds attenuation')
        ax.bar(city_idx + 1.5 * width, As.value, width,
               label='Scintillation attenuation')

        # Set the labels
        ticks = ax.set_xticklabels([''] + list(cities.keys()))
        for t in ticks:
            t.set_rotation(45)
        ax.set_ylabel('Atmospheric attenuation exceeded for 0.1% [dB]')

        # Format image
        ax.yaxis.grid(which='both', linestyle=':')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=2)
        plt.tight_layout(rect=(0, 0, 1, 0.85))


class TestSingleLocation(test.TestCase):

    @staticmethod
    def test_single_location():

        # Location of the receiver ground stations
        lat = 41.39
        lon = -71.05

        # Link parameters
        el = 60                # Elevation angle equal to 60 degrees
        f = 22.5 * itur.u.GHz  # Frequency equal to 22.5 GHz
        D = 1 * itur.u.m       # Receiver antenna diameter of 1 m
        p = 0.1                # We compute values exceeded during 0.1 % of
        # the average year

        # Compute atmospheric parameters
        hs = itur.topographic_altitude(lat, lon)
        T = itur.surface_mean_temperature(lat, lon)
        P = itur.models.itu835.pressure(lat, hs)
        rho_p = itur.surface_water_vapour_density(lat, lon, p, hs)
        itur.models.itu835.water_vapour_density(lat, hs)
        itur.models.itu835.temperature(lat, hs)
        itur.models.itu836.total_water_vapour_content(lat, lon, p, hs)

        # Compute rain and cloud-related parameters
        itur.models.itu618.rain_attenuation_probability(
            lat, lon, el, hs)
        itur.models.itu837.rainfall_probability(lat, lon)
        itur.models.itu837.rainfall_rate(lat, lon, p)
        itur.models.itu839.isoterm_0(lat, lon)
        itur.models.itu839.rain_height(lat, lon)
        itur.models.itu840.columnar_content_reduced_liquid(
            lat, lon, p)
        itur.models.itu676.zenit_water_vapour_attenuation(
            lat, lon, p, f, h=hs)

        # Compute attenuation values
        itur.gaseous_attenuation_slant_path(f, el, rho_p, P, T)
        itur.rain_attenuation(lat, lon, f, el, hs=hs, p=p)
        itur.cloud_attenuation(lat, lon, el, f, p)
        itur.scintillation_attenuation(lat, lon, f, el, p, D)
        itur.atmospheric_attenuation_slant_path(lat, lon, f, el, p, D)


class TestSingleLocationVsFrequency(test.TestCase):

    @staticmethod
    def test_single_location_vs_f():
        # Ground station coordinates (Boston)
        lat_GS = 42.3601
        lon_GS = -71.0942

        ################################################
        # First case: Attenuation vs. frequency        #
        ################################################

        # Satellite coordinates (GEO, 77 W)
        lat_sat = 0
        lon_sat = -77
        h_sat = 35786 * itur.u.km

        # Compute the elevation angle between satellite and ground station
        el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat,
                                        lat_GS, lon_GS)

        f = 22.5 * itur.u.GHz    # Link frequency
        D = 1.2 * itur.u.m       # Antenna diameters
        p = 1

        f = np.logspace(-0.2, 2, 100) * itur.u.GHz

        Ag, Ac, Ar, As, A =\
            itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f,
                                                    el, p, D,
                                                    return_contributions=True)

        # Plot the results
        fig, ax = plt.subplots(1, 1)
        ax.loglog(f, Ag, label='Gaseous attenuation')
        ax.loglog(f, Ac, label='Cloud attenuation')
        ax.loglog(f, Ar, label='Rain attenuation')
        ax.loglog(f, As, label='Scintillation attenuation')
        ax.loglog(f, A, label='Total atmospheric attenuation')

        ax.set_xlabel('Frequency [GHz]')
        ax.set_ylabel('Atmospheric attenuation [dB]')
        ax.grid(which='both', linestyle=':')
        plt.legend()

        ################################################
        # Second case: Attenuation vs. elevation angle #
        ################################################

        f = 22.5 * itur.u.GHz
        el = np.linspace(5, 90, 100)

        Ag, Ac, Ar, As, A =\
            itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS,
                                                    f, el, p, D,
                                                    return_contributions=True)

        # Plot the results
        fig, ax = plt.subplots(1, 1)
        ax.plot(el, Ag, label='Gaseous attenuation')
        ax.plot(el, Ac, label='Cloud attenuation')
        ax.plot(el, Ar, label='Rain attenuation')
        ax.plot(el, As, label='Scintillation attenuation')
        ax.plot(el, A, label='Total atmospheric attenuation')

        ax.set_xlabel('Elevation angle [deg]')
        ax.set_ylabel('Atmospheric attenuation [dB]')
        ax.grid(which='both', linestyle=':')
        plt.legend()


class TestSingleLocationVsUnavailability(test.TestCase):

    @staticmethod
    def test_single_location_vs_p():
        # Ground station coordinates (Boston)
        lat_GS = 42.3601
        lon_GS = -71.0942

        # Satellite coordinates (GEO, 77 W)
        lat_sat = 0
        lon_sat = -77
        h_sat = 35786 * itur.u.km

        # Compute the elevation angle between satellite and ground station
        el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat,
                                        lat_GS, lon_GS)

        f = 22.5 * itur.u.GHz    # Link frequency
        D = 1.2 * itur.u.m       # Antenna diameters

        # Define unavailabilities vector in logarithmic scale
        p = np.logspace(-1.5, 1.5, 100)

        A_g, A_c, A_r, A_s, A_t = \
            itur.atmospheric_attenuation_slant_path(
                lat_GS, lon_GS, f, el, p, D, return_contributions=True)

        # Plot the results using matplotlib
        f, ax = plt.subplots(1, 1)
        ax.semilogx(p, A_g.value, label='Gaseous attenuation')
        ax.semilogx(p, A_c.value, label='Cloud attenuation')
        ax.semilogx(p, A_r.value, label='Rain attenuation')
        ax.semilogx(p, A_s.value, label='Scintillation attenuation')
        ax.semilogx(p, A_t.value, label='Total atmospheric attenuation')

        ax.set_xlabel('Percentage of time attenuation value is exceeded [%]')
        ax.set_ylabel('Attenuation [dB]')
        ax.grid(which='both', linestyle=':')
        plt.legend()


class TestGaseousAttenuation(test.TestCase):

    @staticmethod
    def test_gaseous_attenuation():
        # Define atmospheric parameters
        rho_wet = 7.5 * itur.u.g / itur.u.m**3
        rho_dry = 0 * itur.u.g / itur.u.m**3
        P = 1013.25 * itur.u.hPa
        T = 15 * itur.u.deg_C

        # Define frequency logspace parameters
        N_freq = 1000
        fs = np.linspace(0, 1000, N_freq)

        # Compute the attenuation values
        att_wet = itu676.gamma_exact(fs, P, rho_wet, T)
        att_dry = itu676.gamma_exact(fs, P, rho_dry, T)

        # Plot the results
        plt.figure()
        plt.plot(fs, att_wet.value, 'b--', label='Wet atmosphere')
        plt.plot(fs, att_dry.value, 'r', label='Dry atmosphere')
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Specific attenuation [dB/km]')
        plt.yscale('log')
        plt.xscale('linear')
        plt.xlim(0, 1000)
        plt.ylim(1e-3, 1e5)
        plt.legend()
        plt.grid(which='both', linestyle=':', color='gray',
                 linewidth=0.3, alpha=0.5)
        plt.grid(which='major', linestyle=':', color='black')
        plt.title('FIGURE 1. - Specific attenuation due to atmospheric gases,'
                  '\ncalculated at 1 GHz intervals, including line centres')
        plt.tight_layout()

        #######################################################################
        #               Specific attenuation at different altitudes           #
        #######################################################################

        # Define atmospheric parameters
        hs = np.array([0, 5, 10, 15, 20]) * itur.u.km

        # Define frequency logspace parameters
        N_freq = 2001
        fs = np.linspace(50, 70, N_freq)

        # Plot the results
        plt.figure()

        # Loop over heights and compute values
        for h in hs:
            rho = itu835.standard_water_vapour_density(h)
            P = itu835.standard_pressure(h)
            T = itu835.standard_temperature(h)
            atts = itu676.gamma_exact(fs * itur.u.GHz, P, rho, T)
            plt.plot(fs, atts.value, label='Altitude {0} km'.format(h.value))

        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Specific attenuation [dB/km]')
        plt.yscale('log')
        plt.xscale('linear')
        plt.xlim(50, 70)
        plt.ylim(1e-3, 1e2)
        plt.legend()
        plt.grid(which='both', linestyle=':', color='gray',
                 linewidth=0.3, alpha=0.5)
        plt.grid(which='major', linestyle=':', color='black')
        plt.title('FIGURE 2. - Specific attenuation in the range 50-70 GHz'
                  ' at the\n altitudes indicated, calculated at intervals of'
                  ' 10 MHz\nincluding line centers (0, 5, 10 15, 20) km')
        plt.tight_layout()

        #######################################################################
        #           Comparison of line-by-line and approximate method         #
        #######################################################################
        # Define atmospheric parameters
        el = 90
        rho = 7.5 * itur.u.g / itur.u.m**3
        P = 1013.25 * itur.u.hPa
        T = 15 * itur.u.deg_C

        # Define frequency logspace parameters
        N_freq = 350
        fs = np.linspace(0, 350, N_freq)

        # Initialize result vectors
        atts_approx = []
        atts_exact = []

        # Loop over frequencies and compute values
        atts_approx = itu676.gaseous_attenuation_slant_path(
            fs, el, rho, P, T, mode='approx')

        atts_exact = itu676.gaseous_attenuation_slant_path(
            fs, el, rho, P, T, mode='exact')

        # Plot the results
        plt.figure()
        plt.plot(fs, atts_approx.value, 'b--',
                 label='Approximate method Annex 2')
        plt.plot(fs, atts_exact.value, 'r', label='Exact line-by-line method')
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Attenuation [dB]')
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        plt.grid(which='both', linestyle=':', color='gray',
                 linewidth=0.3, alpha=0.5)
        plt.grid(which='major', linestyle=':', color='black')
        plt.title('Comparison of line-by-line method to approximate method')
        plt.tight_layout()


if __name__ == '__main__':
    suite = suite()
    print('Test examples of the code')
    print('------------------------')
    print(
        'A total of %d test-cases are going to be tested' %
        suite.countTestCases())
    sys.stdout.flush()
    test.TextTestRunner(verbosity=2).run(suite)
