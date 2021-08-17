# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import os.path as path
from collections import defaultdict, OrderedDict

import unittest as test

import itur.models as models
from itur import atmospheric_attenuation_slant_path


if sys.version[:3] <= '3.5':
    pd.set_option('display.max_colwidth', -1)
else:
    pd.set_option('display.max_colwidth', None)

basepath = path.dirname(path.realpath(__file__))
test_data = path.join(basepath, 'test_data')
html_path = path.join(basepath, '../docs/validation')


desc_validation =\
    """This page contains the validation examples for Recommendation {0}: {1}.

All test cases were extracted from the
`ITU Validation examples file (rev 5.1) <https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev5_1.xlsx>`_.

.. contents:: Functions tested
    :depth: 2

"""

desc_test_case =\
    """The table below contains the results of testing function ``{0}``.
The test cases were extracted from spreadsheet ``{1}`` from the
`ITU Validation examples file (rev 5.1) <https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev5_1.xlsx>`_.
In addition to the input-arguments, expected result (``ITU Validation``), and
ITU-Rpy computed result (``ITUR-py Result``), the absolute and relative errors
are shown. Each test case is color-coded depending on the magnitude of the
errors (green = pass, errors are negligible, red = fail, relative error is
above 0.01%).

In addition, the code snippet below shows an example of how to generate the
first row of the results in the table:

.. code-block:: python

{2}

"""

html_header = """
       <html>
         <head>
             <title>Validation results {0}</title>
             <style>
                   table {{
                       border-collapse: collapse;
                       font-size: 10px;
                       width: 100%;
                       padding: 2px;
                     }}

                     th {{
                         background-color: black;
                         color: white;
                     }}

                     th, td {{
                       text-align: center;
                       padding: 2px;
                       width: 1%;
                       font-family: Arial, Helvetica, sans-serif;
                       white-space: nowrap;
                     }}

                     tr:nth-child(even) {{background-color: #f2f2f2;}}
                     tr:hover {{background-color: khaki;}}

               </style>
         </head>
         <body>
   """

html_footer = """
            <br><br>
          </body>
        </html>
        """


def create_ITU_suite():
    """ A test suite for the ITU-P Recommendations. Recommendations tested:
    * ITU-P R-676-11
    * ITU-P R-618-13
    * ITU-P R-453-14
    * ITU-P R-837-7
    * ITU-P R-838-3
    * ITU-P R-839-4
    * ITU-P R-840-8
    * ITU-P R-1511-1
    * ITU-P R-1511-2
    """
    suite = ITU_Suite()

    # ITU-R P.453 tests (Gaseous attenuation)
    suite.add_test(ITUR453_14TestCase('test_wet_term_radio_refractivity'))

    # ITU-R P.618
    suite.add_test(ITUR618_13TestCase('test_rain_attenuation'))
    suite.add_test(ITUR618_13TestCase('test_rain_probability'))
    suite.add_test(ITUR618_13TestCase('test_scintillation_attenuation'))
    suite.add_test(ITUR618_13TestCase('test_total_attenuation'))
    suite.add_test(ITUR618_13TestCase(
        'test_cross_polarization_discrimination'))

    # ITU-R P.676
    suite.add_test(ITUR676_12TestCase('test_gamma0'))
    suite.add_test(ITUR676_12TestCase('test_gammaw'))
    suite.add_test(ITUR676_12TestCase('test_gamma'))
    suite.add_test(ITUR676_12TestCase('test_zenith_attenuation'))
    suite.add_test(ITUR676_12TestCase('test_attenuation_gas'))

    # ITU-R P.836
    suite.add_test(ITUR836_6TestCase(
        'test_surface_water_vapour_density_annual'))
    suite.add_test(ITUR836_6TestCase('test_total_water_vapour_content_annual'))

    # ITU-R P.837
    suite.add_test(ITUR837_7TestCase('test_rainfall_rate'))
    suite.add_test(ITUR837_7TestCase('test_rainfall_rate_probability'))
    suite.add_test(ITUR837_7TestCase('test_rainfall_rate_R001'))

    # ITU-R P.838
    suite.add_test(ITUR838_3TestCase('test_rain_specific_attenuation'))

    # ITU-R P.839
    suite.add_test(ITUR839_4TestCase('test_isoterm_0_deg'))
    suite.add_test(ITUR839_4TestCase('test_rain_height'))

    # ITU-R P.840
    suite.add_test(ITUR840_8TestCase('test_columnar_content_reduced_liquid'))
    suite.add_test(ITUR840_8TestCase('test_cloud_attenuation'))

    # ITU-R P.1510
    suite.add_test(ITUR1510_1TestCase('test_surface_mean_temperature'))

    # ITU-R P.1511
    suite.add_test(ITUR1511_1TestCase('test_topographic_altitude'))
    suite.add_test(ITUR1511_2TestCase('test_topographic_altitude'))

    # ITU-R P.1623
    suite.add_test(ITUR1623_1TestCase(
        'test_fade_duration_cummulative_probability'))
    suite.add_test(ITUR1623_1TestCase('test_fade_duration_number_fades'))
    suite.add_test(ITUR1623_1TestCase('test_fade_duration_probability'))
    suite.add_test(ITUR1623_1TestCase(
        'test_fade_duration_total_exceedance_time'))

    return suite


# Format HTML code
def formatter_fcn(s):
    return '\t\t\t<td style="text-align:left">' + str(s)


def formatter_rel_error_cell(s):
    if np.isnan(float(s)) or np.isinf(float(s)):
        return '\t\t\t<td bgcolor="cornflowerblue">{0:.3f}'.format(s)
    elif abs(float(s)) < 0.01:
        return '\t\t\t<td bgcolor="lightgreen">{0:.3f}'.format(s)
    else:
        return '\t\t\t<td bgcolor="salmon">{0:.3f}'.format(s)


def formatter_error(s):
    if np.isnan(float(s)):
        return '\t\t\t<td bgcolor="cornflowerblue">{0:.2e}'.format(s)
    elif abs(float(s)) < 0.1:
        return '\t\t\t<td bgcolor="lightgreen">{0:.2e}'.format(s)
    else:
        return '\t\t\t<td bgcolor="salmon">{0:.3e}'.format(s)


def format_table(table):
    # Fix cells with a cell within
    table = table.replace('<td><td', '<td')

    # Format headers
    table = table.replace('res_val', 'ITU Validation')
    table = table.replace('res_fcn', 'ITU-Rpy Result')
    table = table.replace('error_rel', 'Relative Error')
    table = table.replace('error', 'Absolute Error')
    table = table.replace('fcn', 'ITU-Rpy Function')

    return table


def formatter_digits(fcn, val):
    ret = []
    COL_STR = '<span style="color: {1}">{0}</span>'
    for f, v in zip(fcn, val):
        i_equal = 0

        # Convert numbers to strings
        s_f = '{:0.6f}'.format(f)
        s_v = '{:0.6f}'.format(v)

        # Determine how many numbers are equal
        for c_f, c_v in zip(s_f, s_v):
            if c_f == c_v:
                i_equal += 1
            else:
                break

        # Format the digits by coloring equal and different sections
        if i_equal > 0:
            s = COL_STR.format(s_f[:i_equal], 'darkgreen')
            if i_equal < len(s_f):
                s += COL_STR.format(s_f[i_equal:], 'darkred')
        else:
            s = COL_STR.format(s_f, 'darkred')

        ret.append(s)

    return ret


class ITU_Suite(test.TestSuite):

    def __init__(self):
        test.TestSuite.__init__(self)
        self.test_cases = OrderedDict({})

    def add_test(self, test_case):
        self.test_cases[test_case.__class__.__name__] = test_case
        self.addTest(test_case)

    def rst_reports(self, path_report=None):
        for test_name, test_case in self.test_cases.items():
            ret = test_case.produce_rst_report()
            if path_report:
                fpath = path.join(
                    path_report,
                    test_name.lower().replace('testcase', '') + '.rst')
                with open(fpath, 'w', encoding='utf-8') as fd:
                    fd.write(ret)


class ITU_TestCase(test.TestCase):
    report = defaultdict(dict)

    def read_csv(self, path_name, columns):
        self.path_name = path_name
        df = pd.read_csv(path_name, sep=',', skiprows=range(1, 2), encoding='cp1252')
        units = pd.read_csv(path_name, sep=',', nrows=2, encoding='cp1252')
        self.units = dict(units[columns].iloc[0])
        return df[columns]

    def setUp(self):
        self.tests = []

    def __run__(self, test_name, test_fcn, df, attributes, result_value,
                n_places=5):

        test_fcn_name = test_fcn
        test_fcn = eval(test_fcn)

        # Add units for the result value and function
        self.units['res_val'] = self.units[result_value]
        self.units['res_fcn'] = self.units[result_value]
        # self.units['error'] = self.units[result_value]
        # self.units['error_rel'] = '(%)'

        # Evaluate all the functions
        res = []
        for i, row in df.iterrows():
            args = {a: row[a] for a in attributes}
            # Evaluate function
            res_fcn = test_fcn(**args)
            res_val = row[result_value]

            # Format dictionary to be added to report
            line = dict(args)
            line['fcn'] = test_fcn_name
            line['res_fcn'] = res_fcn.value
            line['res_val'] = res_val
            line['error'] = (res_val - res_fcn.value)
            if res_val == 0 and abs(line['error']) < 1e-6:
                line['error_rel'] = 0
            else:
                line['error_rel'] = round(
                    (res_val - res_fcn.value) / res_val * 100, 3)

            res.append(line)

        # Create data frame with the report
        order = ['fcn'] + attributes + ['res_val', 'res_fcn', 'error',
                                        'error_rel']
        df = pd.DataFrame(res)
        self.report[self.__class__.__name__][test_name] = {
            'df': df,
            'class_name': self.__class__.__name__,
            'test_name': test_name,
            'units': self.units,
            'path_csv': self.path_name,
            'test_fcn': test_fcn_name,
            'attributes': attributes,
            'n_places': n_places,
            'report_html': df[order]
        }

        # Do the assert equal for all the tests
        for ret in res:
            res_val = ret['res_val']
            res_fcn = ret['res_fcn']
            try:
                self.assertAlmostEqual(res_val, res_fcn, places=n_places)
            except AssertionError as e:
                print(e)

    def generate_code_example(self, report):
        """Generate a code example of the call used for this function

        Parameters
        ----------
        report : dict
          The dictionary containing the parameters used in a test case.

        Returns
        -------
        str
          A string with an example of code used to call this test.
        """
        ret = ["    import itur", "", "    # Define input attributes"]
        attributes = report['attributes']
        test_fcn = report['test_fcn']
        test_fcn_name = report['test_fcn'].split('.')[-1]

        row_1 = report['df'].iloc[0]
        units = report['units']

        # Write the attributesf
        for attr_name, attr_val in zip(attributes, row_1[attributes]):
            ret.append('    {0} = {1}  # {2}'.format(
                attr_name, attr_val, units[attr_name]))

        # Add call to test-function
        ret.extend([
            "", "    # Make call to test-function {0}".format(test_fcn_name),
            '    itur_val = itur.{0}({1})'.format(
                test_fcn, ", ".join([att + '=' + att for att in attributes]))])

        # Compute errors
        ret.extend([
            "",
            "    # Compute error with respect to value in ITU example file",
            '    ITU_example_val = {0}  # {1}'.format(
                row_1['res_val'], units['res_val']),
            '    error = ITU_example_val - itur_val.value',
            '    error_rel = error / ITU_example_val * 100  # (%)'])

        return "\n".join(ret)

    def produce_rst_report(self):

        ret = []
        title = "Validation results {0}".format(self.itu_name)
        ret.append(title)
        ret.append('=' * len(title))
        ret.append('')
        ret.append(desc_validation.format(self.itu_name, self.itu_description))

        for test_name in self.report[self.__class__.__name__]:

            # Create HTML table for this test
            report = self.report[self.__class__.__name__][test_name]
            table = self.create_html_table(report, include_header=True)
            html_file = path.join(html_path, test_name + '_table.html')
            with open(html_file, 'w', encoding='utf-8') as fd:
                fd.write(table)

            # Create HTML table for this test
            test_fcn_name = report['test_fcn'].split('.')[-1]
            test_case_name = "Function {0}".format(test_fcn_name)
            ret.append(test_case_name)
            ret.append("-" * len(test_case_name))
            ret.append("")
            ret.append(desc_test_case.format(
                test_fcn_name, path.basename(report['path_csv']),
                self.generate_code_example(report)))

            ret.extend(['.. raw:: html',
                        '    :file: {0}'.format(test_name + '_table.html'),
                        '', ''])

        return "\n".join(ret)

    def create_html_table(self, report, include_header=False):
        fmtrs = {'error_rel': formatter_rel_error_cell,
                 'error': formatter_error,
                 'fcn': formatter_fcn}

        if include_header:
            table = html_header.format(report['test_fcn'])

        df = report['report_html']
        df['res_fcn'] = formatter_digits(df['res_fcn'], df['res_val'])

        # Add units to header attributes
        col_dict = {}
        for col in df.columns:
            if col in report['units']:
                col_dict[col] = '{0} {1}'.format(col, report['units'][col])
            else:
                col_dict[col] = col

        df.rename(columns=col_dict, inplace=True)
        table += df.to_html(bold_rows=True, index=False, justify='center',
                            table_id=report['test_name'].lower(), escape=False,
                            formatters=fmtrs)

        table = format_table(table)
        if include_header:
            table += html_footer
        return table


class ITUR453_14TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.453-14'
    itu_description = 'TBD'

    def test_wet_term_radio_refractivity(self):
        # Set the version to the
        models.itu453.change_version(13)

        # Read the test data
        df = self.read_csv(path.join(test_data, '453/ITURP453-14_Nwet.csv'),
                           columns=['lat', 'lon', 'p', 'Nwet'])

        # Run test and generate the report
        self.__run__('test_wet_term_radio_refractivity',
                     test_fcn='models.itu453.map_wet_term_radio_refractivity',
                     df=df, attributes=['lat', 'lon', 'p'],
                     result_value='Nwet',
                     n_places=5)


class ITUR618_13TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.618-13'
    itu_description = 'Propagation data and prediction methods required for' +\
                      ' the design of Earth-space telecommunication systems'

    def test_rain_attenuation(self):
        # Set the version to the
        models.itu618.change_version(13)

        # Read the test data
        df = self.read_csv(path.join(test_data, '618/ITURP618-13_A_rain.csv'),
                           columns=['lat', 'lon', 'hs', 'el', 'f', 'tau', 'p',
                                    'R001', 'A_rain'])

        # Run test and generate the report
        self.__run__('test_rain_attenuation',
                     test_fcn='models.itu618.rain_attenuation',
                     df=df, attributes=['lat', 'lon', 'hs', 'el', 'f',
                                        'tau', 'p', 'R001'],
                     result_value='A_rain',
                     n_places=5)

    def test_rain_probability(self):
        # Set the version to the
        models.itu618.change_version(13)

        # Read the test data
        df = self.read_csv(path.join(test_data, '618/ITURP618-13_A_rain.csv'),
                           columns=['lat', 'lon', 'hs', 'el', 'Ls', 'P0',
                                    'P_rain'])

        # Run test and generate the report
        self.__run__('test_rain_probability',
                     test_fcn='models.itu618.rain_attenuation_probability',
                     df=df, attributes=['lat', 'lon', 'hs', 'el', 'Ls', 'P0'],
                     result_value='P_rain',
                     n_places=5)

    def test_scintillation_attenuation(self):
        # Set the version to the
        models.itu618.change_version(13)

        # Read the test data
        df = self.read_csv(path.join(test_data, '618/ITURP618-13_A_sci.csv'),
                           columns=['lat', 'lon', 'f', 'el', 'p', 'D', 'eta',
                                    'A_scin'])

        # Run test and generate the report
        self.__run__('test_scintillation_attenuation',
                     test_fcn='models.itu618.scintillation_attenuation',
                     df=df, attributes=['lat', 'lon', 'f', 'el', 'p', 'D',
                                        'eta'],
                     result_value='A_scin',
                     n_places=5)

    def test_cross_polarization_discrimination(self):
        # Set the version to the
        models.itu618.change_version(13)

        # Read the test data
        df = self.read_csv(path.join(test_data, '618/ITURP618-13_A_xpd.csv'),
                           columns=['f', 'el', 'p', 'tau', 'Ap', 'XPD'])

        # Run test and generate the report
        self.__run__('test_cross_polarization_discrimination',
                     test_fcn='models.itu618.rain_cross_polarization_discrimination',
                     df=df, attributes=['f', 'el', 'p', 'tau', 'Ap'],
                     result_value='XPD',
                     n_places=5)

    def test_total_attenuation(self):
        # Set the version to the
        models.itu618.change_version(13)

        # Read the test data
        df = self.read_csv(path.join(test_data, '618/ITURP618-13_A_total.csv'),
                           columns=['lat', 'lon', 'f', 'el', 'p', 'D', 'eta',
                                    'tau', 'hs', 'A_total'])

        # Run test and generate the report
        self.__run__('test_total_attenuation',
                     test_fcn='atmospheric_attenuation_slant_path',
                     df=df, attributes=['lat', 'lon', 'f', 'el', 'p', 'D',
                                        'eta', 'tau', 'hs'],
                     result_value='A_total',
                     n_places=4)


class ITUR676_12TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.676-12'
    itu_description = 'Attenuation by atmospheric gases and related effects'

    def test_gamma0(self):
        # Set the version to the
        models.itu676.change_version(12)

        path_file = '676/ITURP676-12_gamma.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['f', 'P', 'rho', 'T', 'gamma0'])

        # Run test and generate the report
        self.__run__('test_gamma0',
                     test_fcn='models.itu676.gamma0_exact',
                     df=df, attributes=['f', 'P', 'rho', 'T'],
                     result_value='gamma0',
                     n_places=5)

    def test_gammaw(self):
        # Set the version to the
        models.itu676.change_version(12)

        path_file = '676/ITURP676-12_gamma.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['f', 'P', 'rho', 'T', 'gammaw'])

        # Run test and generate the report
        self.__run__('test_gammaw',
                     test_fcn='models.itu676.gammaw_exact',
                     df=df, attributes=['f', 'P', 'rho', 'T'],
                     result_value='gammaw',
                     n_places=5)

    def test_gamma(self):
        # Set the version to the
        models.itu676.change_version(12)

        path_file = '676/ITURP676-12_gamma.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['f', 'P', 'rho', 'T', 'gamma'])

        # Run test and generate the report
        self.__run__('test_gamma',
                     test_fcn='models.itu676.gamma_exact',
                     df=df, attributes=['f', 'P', 'rho', 'T'],
                     result_value='gamma',
                     n_places=5)

    def test_attenuation_gas(self):
        # Set the version to the
        models.itu676.change_version(12)

        path_file = '676/ITURP676-12_A_gas.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['f', 'el', 'P', 'rho', 'T', 'h', 'V_t',
                                    'A_gas'])

        # Run test and generate the report
        self.__run__('test_attenuation_gas',
                     test_fcn='models.itu676.gaseous_attenuation_slant_path',
                     df=df, attributes=['f', 'el', 'rho', 'P', 'T', 'h',
                                        'V_t'],
                     result_value='A_gas',
                     n_places=5)

    def test_zenith_attenuation(self):
        # Set the version to the
        models.itu676.change_version(12)

        path_file = '676/ITURP676-12_zenith_attenuation.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'p', 'f', 'h', 'V_t', 'Aw'])

        # Run test and generate the report
        self.__run__('test_zenith_attenuation',
                     test_fcn='models.itu676.zenit_water_vapour_attenuation',
                     df=df, attributes=['lat', 'lon', 'p', 'f', 'h', 'V_t'],
                     result_value='Aw',
                     n_places=5)


class ITUR836_6TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.836-6'
    itu_description = 'Water vapour: surface density and total columnar content'

    def test_surface_water_vapour_density_annual(self):
        # Set the version to the
        models.itu836.change_version(6)

        path_file = '836/ITURP836-6_surface_water_vapour_density_annual.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'alt', 'p', 'rho'])

        # Run test and generate the report
        self.__run__('test_surface_water_vapour_density_annual',
                     test_fcn='models.itu836.surface_water_vapour_density',
                     df=df, attributes=['lat', 'lon', 'alt', 'p'],
                     result_value='rho',
                     n_places=5)

    def test_total_water_vapour_content_annual(self):
        # Set the version to the
        models.itu836.change_version(6)

        path_file = '836/ITURP836-6_total_water_vapour_content_annual.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'alt', 'p', 'V'])

        # Run test and generate the report
        self.__run__('test_total_water_vapour_content_annual',
                     test_fcn='models.itu836.total_water_vapour_content',
                     df=df, attributes=['lat', 'lon', 'alt', 'p'],
                     result_value='V',
                     n_places=5)


class ITUR837_7TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.837-7'
    itu_description = 'Characteristics of precipitation for propagation ' +\
                      'modelling'

    def test_rainfall_rate(self):
        # Set the version to the
        models.itu837.change_version(7)

        path_file = '837/ITURP837-7_rainfall_rate.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'p', 'Rp'])

        # Run test and generate the report
        self.__run__('test_rainfall_rate',
                     test_fcn='models.itu837.rainfall_rate',
                     df=df, attributes=['lat', 'lon', 'p'],
                     result_value='Rp',
                     n_places=3)

    def test_rainfall_rate_R001(self):
        # Set the version to the
        models.itu837.change_version(7)

        path_file = '837/ITURP837-7_rainfall_rate_R001.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'p', 'Rp'])

        # Run test and generate the report
        self.__run__('test_rainfall_rate_R001',
                     test_fcn='models.itu837.rainfall_rate',
                     df=df, attributes=['lat', 'lon', 'p'],
                     result_value='Rp',
                     n_places=5)

    def test_rainfall_rate_probability(self):
        # Set the version to the
        models.itu837.change_version(7)

        path_file = '837/ITURP837-7_rainfall_rate_probability.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'p'])

        # Run test and generate the report
        self.__run__('test_rainfall_rate_probability',
                     test_fcn='models.itu837.rainfall_probability',
                     df=df, attributes=['lat', 'lon'],
                     result_value='p',
                     n_places=5)


class ITUR838_3TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.838-3'
    itu_description = 'Specific attenuation model for rain for use in ' +\
                      'prediction methods'

    def test_rain_specific_attenuation(self):
        # Set the version to the
        models.itu838.change_version(3)

        path_file = '838/ITURP838-3_rain_specific_attenuation.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['el', 'f', 'R', 'tau', 'gamma_r'])

        # Run test and generate the report
        self.__run__('test_rain_specific_attenuation',
                     test_fcn='models.itu838.rain_specific_attenuation',
                     df=df, attributes=['el', 'f', 'R', 'tau'],
                     result_value='gamma_r',
                     n_places=5)


class ITUR839_4TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.839-4'
    itu_description = 'Rain height model for prediction methods'

    def test_isoterm_0_deg(self):
        # Set the version to the
        models.itu839.change_version(4)

        path_file = '839/ITURP839-4_rain_height.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'h0'])

        # Run test and generate the report
        self.__run__('test_isoterm_0_deg',
                     test_fcn='models.itu839.isoterm_0',
                     df=df, attributes=['lat', 'lon'],
                     result_value='h0',
                     n_places=5)

    def test_rain_height(self):
        # Set the version to the
        models.itu839.change_version(4)

        path_file = '839/ITURP839-4_rain_height.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'hr'])

        # Run test and generate the report
        self.__run__('test_rain_height',
                     test_fcn='models.itu839.rain_height',
                     df=df, attributes=['lat', 'lon'],
                     result_value='hr',
                     n_places=5)


class ITUR840_8TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.840-8'
    itu_description = 'Attenuation due to clouds and fog'

    def test_columnar_content_reduced_liquid(self):
        # Set the version to the
        models.itu840.change_version(8)

        path_file = '840/ITURP840-8_columnar_content_reduced_liquid.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'p', 'Lred'])

        # Run test and generate the report
        self.__run__('test_columnar_content_reduced_liquid',
                     test_fcn='models.itu840.columnar_content_reduced_liquid',
                     df=df, attributes=['lat', 'lon', 'p'],
                     result_value='Lred',
                     n_places=5)

    def test_cloud_attenuation(self):
        # Set the version to the
        models.itu840.change_version(8)

        path_file = '840/ITURP840-8_cloud_attenuation.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'f', 'el', 'p', 'Ac'])

        # Run test and generate the report
        self.__run__('test_cloud_attenuation',
                     test_fcn='models.itu840.cloud_attenuation',
                     df=df, attributes=['lat', 'lon', 'f', 'el', 'p'],
                     result_value='Ac',
                     n_places=5)


class ITUR1510_1TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.1510-1'
    itu_description = 'Mean surface temperature'

    def test_surface_mean_temperature(self):
        # Set the version to the
        models.itu1510.change_version(1)

        path_file = '1510/ITURP1510-1_temperature.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'T'])

        # Run test and generate the report
        self.__run__('test_surface_mean_temperature',
                     test_fcn='models.itu1510.surface_mean_temperature',
                     df=df, attributes=['lat', 'lon'],
                     result_value='T',
                     n_places=5)


class ITUR1511_1TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.1511-1'
    itu_description = 'Topography for Earth-to-space propagation modelling'

    def test_topographic_altitude(self):
        # Set the version to the
        models.itu1511.change_version(1)

        path_file = '1511/ITURP1511-1_topographic_altitude.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'hs'])

        # Run test and generate the report
        self.__run__('test_topographic_altitude',
                     test_fcn='models.itu1511.topographic_altitude',
                     df=df, attributes=['lat', 'lon'],
                     result_value='hs',
                     n_places=5)


class ITUR1511_2TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.1511-2'
    itu_description = 'Topography for Earth-to-space propagation modelling'

    def test_topographic_altitude(self):
        # Set the version to the
        models.itu1511.change_version(2)

        path_file = '1511/ITURP1511-2_topographic_altitude.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['lat', 'lon', 'hs'])

        # Run test and generate the report
        self.__run__('test_topographic_altitude',
                     test_fcn='models.itu1511.topographic_altitude',
                     df=df, attributes=['lat', 'lon'],
                     result_value='hs',
                     n_places=5)


class ITUR1623_1TestCase(ITU_TestCase):

    itu_name = 'ITU-R P.1623-1'
    itu_description = 'Prediction method of fade dynamics on Earth-space paths'

    def test_fade_duration_probability(self):
        # Set the version to the
        models.itu1623.change_version(1)

        path_file = '1623/ITURP1623-1_fade_duration_params.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['D', 'A', 'el', 'f', 'P'])

        # Run test and generate the report
        self.__run__('test_fade_duration_probability',
                     test_fcn='models.itu1623.fade_duration_probability',
                     df=df, attributes=['D', 'A', 'el', 'f'],
                     result_value='P',
                     n_places=5)

    def test_fade_duration_cummulative_probability(self):
        # Set the version to the
        models.itu1623.change_version(1)

        path_file = '1623/ITURP1623-1_fade_duration_params.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['D', 'A', 'el', 'f', 'F'])

        # Run test and generate the report
        self.__run__('test_fade_duration_cummulative_probability',
                     test_fcn='models.itu1623.fade_duration_cummulative_probability',
                     df=df, attributes=['D', 'A', 'el', 'f'],
                     result_value='F',
                     n_places=5)

    def test_fade_duration_total_exceedance_time(self):
        # Set the version to the
        models.itu1623.change_version(1)

        path_file = '1623/ITURP1623-1_fade_duration_params.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['D', 'A', 'el', 'f', 'T_tot', 'T'])

        # Run test and generate the report
        self.__run__('test_fade_duration_total_exceedance_time',
                     test_fcn='models.itu1623.fade_duration_total_exceedance_time',
                     df=df, attributes=['D', 'A', 'el', 'f', 'T_tot'],
                     result_value='T',
                     n_places=5)

    def test_fade_duration_number_fades(self):
        # Set the version to the
        models.itu1623.change_version(1)

        path_file = '1623/ITURP1623-1_number_of_fades.csv'
        # Read the test data
        df = self.read_csv(path.join(test_data, path_file),
                           columns=['D', 'A', 'el', 'f', 'T_tot', 'N'])

        # Run test and generate the report
        self.__run__('test_fade_duration_number_fades',
                     test_fcn='models.itu1623.fade_duration_number_fades',
                     df=df, attributes=['D', 'A', 'el', 'f', 'T_tot'],
                     result_value='N',
                     n_places=5)


class TestValidationReports(test.TestCase):

    @test.skipIf(sys.version_info[0] < 3,
                 "Only supported in Python 3+ (open does not have encoding)")
    def test_validation_reports(self):
        # Test valid versions
        suite = create_ITU_suite()
        test.TextTestRunner(verbosity=0).run(suite)
        suite.rst_reports(html_path)


if __name__ == '__main__':
    suite = create_ITU_suite()
    print('Validation tests for the ITU-R models')
    print('------------------------')
    print('A total of %d test-cases are going to be tested' %
          suite.countTestCases())
    sys.stdout.flush()
    test.TextTestRunner(verbosity=2).run(suite)
    # suite.html_reports(html_path)
    suite.rst_reports(html_path)
