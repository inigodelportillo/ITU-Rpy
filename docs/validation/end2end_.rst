Validation results End_to_end_attenuation_components
====================================================

This page contains the validation examples for Recommendation End_to_end_attenuation_components: Tests for end to end atmospheric attenuation with components.

All test cases were extracted from the
`ITU Validation examples file (rev 5.1) <https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev5_1.xlsx>`_.

.. contents:: Functions tested
    :depth: 2


Function atmospheric_attenuation_slant_path
-------------------------------------------

The table below contains the results of testing function ``atmospheric_attenuation_slant_path``.
The test cases were extracted from spreadsheet ``e2e_A_total.csv`` from the
`ITU Validation examples file (rev 5.1) <https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev5_1.xlsx>`_.
In addition to the input-arguments, expected result (``ITU Validation``), and
ITU-Rpy computed result (``ITUR-py Result``), the absolute and relative errors
are shown. Each test case is color-coded depending on the magnitude of the
errors (green = pass, errors are negligible, red = fail, relative error is
above 0.01%).

In addition, the code snippet below shows an example of how to generate the
first row of the results in the table:

.. code-block:: python

    import itur

    # Define input attributes
    lat = 51.5  #  (Â°N)
    lon = -0.14  # (Â°E)
    f = 13.75  # (GHz)
    el = 31.07699124  # (Â°)
    p = 1.0  # (%)
    D = 1  # (m)
    eta = 0.65  # (0-1)
    tau = 0  # (Â°)
    hs = 0.031382984  #  (km)
    rain_rate_mode = exact  # nan

    # Make call to test-function atmospheric_attenuation_slant_path
    itur_val = itur.atmospheric_attenuation_slant_path(lat=lat, lon=lon, f=f, el=el, p=p, D=D, eta=eta, tau=tau, hs=hs, rain_rate_mode=rain_rate_mode)

    # Compute error with respect to value in ITU example file
    ITU_example_val = 0.822003497  # (dB)
    error = ITU_example_val - itur_val.value
    error_rel = error / ITU_example_val * 100  # (%)


.. raw:: html
    :file: end2end_total_attenuation_table.html


Function atmospheric_attenuation_slant_path
-------------------------------------------

The table below contains the results of testing function ``atmospheric_attenuation_slant_path``.
The test cases were extracted from spreadsheet ``e2e_A_cloud.csv`` from the
`ITU Validation examples file (rev 5.1) <https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev5_1.xlsx>`_.
In addition to the input-arguments, expected result (``ITU Validation``), and
ITU-Rpy computed result (``ITUR-py Result``), the absolute and relative errors
are shown. Each test case is color-coded depending on the magnitude of the
errors (green = pass, errors are negligible, red = fail, relative error is
above 0.01%).

In addition, the code snippet below shows an example of how to generate the
first row of the results in the table:

.. code-block:: python

    import itur

    # Define input attributes
    lat = 51.5  #  (??N)
    lon = -0.14  # (??E)
    f = 13.75  # (GHz)
    el = 31.07699124  # (??)
    p = 1.0  # (%)
    D = 1.0  # (m)
    eta = 0.65  # (0-1)
    tau = 0.0  # (??)
    hs = 0.031382984  #  (km)
    include_rain = 0.0  # (0-1)
    include_gas = 0.0  # (0-1)
    include_scintillation = 0.0  # (0-1)
    include_clouds = 1.0  # (0-1)

    # Make call to test-function atmospheric_attenuation_slant_path
    itur_val = itur.atmospheric_attenuation_slant_path(lat=lat, lon=lon, f=f, el=el, p=p, D=D, eta=eta, tau=tau, hs=hs, include_rain=include_rain, include_gas=include_gas, include_scintillation=include_scintillation, include_clouds=include_clouds)

    # Compute error with respect to value in ITU example file
    ITU_example_val = 0.12277338  # (dB)
    error = ITU_example_val - itur_val.value
    error_rel = error / ITU_example_val * 100  # (%)


.. raw:: html
    :file: end2end_cloud_attenuation_table.html


Function atmospheric_attenuation_slant_path
-------------------------------------------

The table below contains the results of testing function ``atmospheric_attenuation_slant_path``.
The test cases were extracted from spreadsheet ``e2e_A_gaseous.csv`` from the
`ITU Validation examples file (rev 5.1) <https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev5_1.xlsx>`_.
In addition to the input-arguments, expected result (``ITU Validation``), and
ITU-Rpy computed result (``ITUR-py Result``), the absolute and relative errors
are shown. Each test case is color-coded depending on the magnitude of the
errors (green = pass, errors are negligible, red = fail, relative error is
above 0.01%).

In addition, the code snippet below shows an example of how to generate the
first row of the results in the table:

.. code-block:: python

    import itur

    # Define input attributes
    lat = 51.5  #  (??N)
    lon = -0.14  # (??E)
    f = 13.75  # (GHz)
    el = 31.07699124  # (??)
    p = 1.0  # (%)
    D = 1.0  # (m)
    eta = 0.65  # (0-1)
    tau = 0.0  # (??)
    hs = 0.031382984  #  (km)
    include_rain = 0.0  # (0-1)
    include_gas = 1.0  # (0-1)
    include_scintillation = 0.0  # (0-1)
    include_clouds = 0.0  # (0-1)

    # Make call to test-function atmospheric_attenuation_slant_path
    itur_val = itur.atmospheric_attenuation_slant_path(lat=lat, lon=lon, f=f, el=el, p=p, D=D, eta=eta, tau=tau, hs=hs, include_rain=include_rain, include_gas=include_gas, include_scintillation=include_scintillation, include_clouds=include_clouds)

    # Compute error with respect to value in ITU example file
    ITU_example_val = 0.19084749  # (dB)
    error = ITU_example_val - itur_val.value
    error_rel = error / ITU_example_val * 100  # (%)


.. raw:: html
    :file: end2end_gaseous_attenuation_table.html


Function atmospheric_attenuation_slant_path
-------------------------------------------

The table below contains the results of testing function ``atmospheric_attenuation_slant_path``.
The test cases were extracted from spreadsheet ``e2e_A_rain.csv`` from the
`ITU Validation examples file (rev 5.1) <https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev5_1.xlsx>`_.
In addition to the input-arguments, expected result (``ITU Validation``), and
ITU-Rpy computed result (``ITUR-py Result``), the absolute and relative errors
are shown. Each test case is color-coded depending on the magnitude of the
errors (green = pass, errors are negligible, red = fail, relative error is
above 0.01%).

In addition, the code snippet below shows an example of how to generate the
first row of the results in the table:

.. code-block:: python

    import itur

    # Define input attributes
    lat = 51.5  #  (??N)
    lon = -0.14  # (??E)
    f = 13.75  # (GHz)
    el = 31.07699124  # (??)
    p = 1.0  # (%)
    D = 1  # (m)
    eta = 0.65  # (0-1)
    tau = 0  # (??)
    hs = 0.031382984  #  (km)
    rain_rate_mode = exact  # nan
    include_rain = 1  # (0-1)
    include_gas = 0  # (0-1)
    include_scintillation = 0  # (0-1)
    include_clouds = 0  # (0-1)

    # Make call to test-function atmospheric_attenuation_slant_path
    itur_val = itur.atmospheric_attenuation_slant_path(lat=lat, lon=lon, f=f, el=el, p=p, D=D, eta=eta, tau=tau, hs=hs, rain_rate_mode=rain_rate_mode, include_rain=include_rain, include_gas=include_gas, include_scintillation=include_scintillation, include_clouds=include_clouds)

    # Compute error with respect to value in ITU example file
    ITU_example_val = 0.45380962  # (dB)
    error = ITU_example_val - itur_val.value
    error_rel = error / ITU_example_val * 100  # (%)


.. raw:: html
    :file: end2end_rain_attenuation_table.html


Function atmospheric_attenuation_slant_path
-------------------------------------------

The table below contains the results of testing function ``atmospheric_attenuation_slant_path``.
The test cases were extracted from spreadsheet ``e2e_A_scintillation.csv`` from the
`ITU Validation examples file (rev 5.1) <https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev5_1.xlsx>`_.
In addition to the input-arguments, expected result (``ITU Validation``), and
ITU-Rpy computed result (``ITUR-py Result``), the absolute and relative errors
are shown. Each test case is color-coded depending on the magnitude of the
errors (green = pass, errors are negligible, red = fail, relative error is
above 0.01%).

In addition, the code snippet below shows an example of how to generate the
first row of the results in the table:

.. code-block:: python

    import itur

    # Define input attributes
    lat = 51.5  #  (??N)
    lon = -0.14  # (??E)
    f = 13.75  # (GHz)
    el = 31.07699124  # (??)
    p = 1.0  # (%)
    D = 1.0  # (m)
    eta = 0.65  # (0-1)
    tau = 0.0  # (??)
    hs = 0.031382984  #  (km)
    include_rain = 0.0  # (0-1)
    include_gas = 0.0  # (0-1)
    include_scintillation = 1.0  # (0-1)
    include_clouds = 0.0  # (0-1)

    # Make call to test-function atmospheric_attenuation_slant_path
    itur_val = itur.atmospheric_attenuation_slant_path(lat=lat, lon=lon, f=f, el=el, p=p, D=D, eta=eta, tau=tau, hs=hs, include_rain=include_rain, include_gas=include_gas, include_scintillation=include_scintillation, include_clouds=include_clouds)

    # Compute error with respect to value in ITU example file
    ITU_example_val = 0.25672934  # (dB)
    error = ITU_example_val - itur_val.value
    error_rel = error / ITU_example_val * 100  # (%)


.. raw:: html
    :file: end2end_scintillation_attenuation_table.html

