Validation
==========

ITU-Rpy has been validated using the `ITU Validation examples (rev 8.3) <https://www.itu.int/en/ITU-R/study-groups/rsg3/rwp3m/Validation%20Example/CG-3M3J-13-ValEx-Rev8.3.0.xlsx>`_ , which provides test cases for parts of Recommendations ITU-R P.453-14, P.618-14, P.676-13, P.836-6, P.837-7, P.838-3, P.839-4, P.840-9, P.1511-3, P.1623-1, P.1815-1, P.2145-0, and end to end tests.

For each part of the recommendation that has a counterpart function in ITU-Rpy, the results of the executing the ITUR-py function are compared against the values provided in the ITU Validation examples. The absolute and relative errors between these two quantities are computed, and each test case is color-coded (green = pass, errors are negligible, red = fail, errors are above 0.01%).

The links below show the validation results for each of the recommendations:

.. toctree::

    validation/itur453_14
    validation/itur618_13
    validation/itur618_14
    validation/itur676_12
    validation/itur836_6
    validation/itur837_7
    validation/itur838_3
    validation/itur839_4
    validation/itur840_8
    validation/itur1510_1
    validation/itur1511_1
    validation/itur1511_2
    validation/itur1511_3
    validation/itur1623_1
    validation/itur2145_0