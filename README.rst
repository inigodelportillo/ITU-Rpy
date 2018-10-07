ITU-Rpy
=======

|GitHub license| |Build Status| |PyPI version| |Coverage Status|

A python implementation of the ITU-R P. Recommendations to compute
atmospheric attenuation in slant and horizontal paths.

The propagation loss on an Earth-space path and a horizontal-path,
relative to the free-space loss, is the sum of different contributions,
namely: attenuation by atmospheric gases; attenuation by rain, other
precipitation and clouds; scintillation and multipath effects;
attenuation by sand and dust storms. Each of these contributions has its
own characteristics as a function of frequency, geographic location and
elevation angle. ITU-Rpy allows for fast, vectorial computation of the
different contributions to the atmospheric attenuation.

Documentation
-------------

The documentation can be found at `ITU-Rpy
documentation <http://itu-rpy.readthedocs.io/en/latest/index.html>`__ in
Read the docs.

Examples of use cases can be found in the `examples
folder <https://github.com/iportillo/ITU-Rpy/tree/master/examples>`__.

Installation
------------

ITU-Rpy has the followind dependencies: ``numpy``, ``scipy``,
``joblib``, ``pyproj``, and ``astropy``. Installation of ``basemap`` and
``matplotlib`` is recommended to display results in a map.

Using pip, you can install all of them by running:

::

    pip install itur

More information about the installation process can be found on the
`documentation <https://github.com/iportillo/ITU-Rpy/blob/master/docs/installation.rst>`__.

ITU-R Recommendations implemented:
----------------------------------

The following ITU-R Recommendations are implemented in ITU-Rpy \*
**ITU-R P.453-13:** The radio refractive index: its formula and
refractivity data \* **ITU-R P.618-13:** Propagation data and prediction
methods required for the design of Earth-space telecommunication systems
\* **ITU-R P.676-11:** Attenuation by atmospheric gases \* **ITU-R
P.835-6:** Reference Standard Atmospheres \* **ITU-R P.836-6:** Water
vapour: surface density and total columnar content \* **ITU-R P.837-7:**
Characteristics of precipitation for propagation modelling \* **ITU-R
P.838-3:** Specific attenuation model for rain for use in prediction
methods \* **ITU-R P.839-4:** Rain height model for prediction methods.
\* **ITU-R P.840-7:** Attenuation due to clouds and fog \* **ITU-R
P.1144-7:** Interpolation methods for the geophysical properties used to
compute propagation effects \* **ITU-R P.1511-1:** Topography for
Earth-to-space propagation modelling \* **ITU-R P.1853-1:** Tropospheric
attenuation time series synthesis

The individual models can be accessed using the ``itur.models`` package.

Usage
-----

The following code example shows the usage of ITU-Rpy. More examples can
be found in the `examples
folder <https://github.com/iportillo/ITU-Rpy/tree/master/examples>`__.

.. code:: python

    import itur

    f = 22.5 * itur.u.GHz    # Link frequency
    D = 1 * itur.u.m         # Size of the receiver antenna
    el = 60                  # Elevation angle constant of 60 degrees
    p = 3                    # Percentage of time that attenuation values are exceeded.
        
    # Generate a regular grid latitude and longitude points with 1 degrees resolution   
    lat, lon = itur.utils.regular_lat_lon_grid() 

    # Comute the atmospheric attenuation
    Att = itur.atmospheric_attenuation_slant_path(lat, lon, f, el, p, D) 
    itur.utils.plot_in_map(Att.value, lat, lon, 
                           cbar_text='Atmospheric attenuation [dB]')

which produces: |Attenuation worldmap|

Citation
--------

If you use ITU-Rpy in one of your research projects, please cite it as:

::

    @misc{iturpy-2017,
          title={ITU-Rpy: A python implementation of the ITU-R P. Recommendations to compute atmospheric
             attenuation in slant and horizontal paths.},
          author={Inigo del Portillo},
          year={2017},
          publisher={GitHub},
          howpublished={\url{https://github.com/iportillo/ITU-Rpy/}}
    }

.. |GitHub license| image:: https://img.shields.io/badge/license-MIT-lightgrey.svg
   :target: https://raw.githubusercontent.com/Carthage/Carthage/master/LICENSE.md
.. |Build Status| image:: https://travis-ci.org/iportillo/ITU-Rpy.svg?branch=master
   :target: https://travis-ci.org/iportillo/ITU-Rpy
.. |PyPI version| image:: https://badge.fury.io/py/itur.svg
   :target: https://badge.fury.io/py/itur
.. |Coverage Status| image:: https://coveralls.io/repos/github/iportillo/ITU-Rpy/badge.svg?branch=master
   :target: https://coveralls.io/github/iportillo/ITU-Rpy?branch=master
.. |Attenuation worldmap| image:: https://raw.githubusercontent.com/iportillo/ITU-Rpy/master/docs/images/att_world.png

