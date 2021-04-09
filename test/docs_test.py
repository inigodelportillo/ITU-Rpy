# -*- coding: utf-8 -*-
import sys
import unittest as test

from ITU_validation_report_test import create_ITU_suite, html_path


def suite():
    """A test suite for the documentation of ITUR-py. """
    suite = test.TestSuite()

    suite.addTest(TestValidationReports('test_validation_reports'))
    return suite


class TestValidationReports(test.TestCase):

    def test_validation_reports(self):
        # Test valid versions
        suite = create_ITU_suite()
        test.TextTestRunner(verbosity=0).run(suite)
        suite.rst_reports(html_path)


if __name__ == '__main__':

    # Create the test suite and run all the tests
    suite = suite()
    print('Test dcoumentation works')
    print('------------------------')
    print(
        'A total of %d test-cases are going to be tested' %
        suite.countTestCases())
    sys.stdout.flush()
    test.TextTestRunner(verbosity=2).run(suite)
