ITU-Rpy |GitHub license|
========================

A python implementation of the ITU-P Recommendations

Dependencies
------------

-  numpy
-  scipy

Installation
------------

Using pip, simply run:

::

    pip install iturpy

Using conda:

::

    conda install iturpy

ITU-P Recommendations implemented:
----------------------------------

The following ITU-P Recommendations are implemented in ITU-Rpy \*
**ITU-P R.453-12:** The radio refractive index: its formula and
refractivity data \* **ITU-P R.618-12:** Propagation data and prediction
methods required for the design of Earth-space telecommunication systems
\* **ITU-P R.676-11:** Attenuation by atmospheric gases \* **ITU-P
R.835-12:** Reference Standard Atmospheres \* **ITU-P R.836-5:** Water
vapour: surface density and total columnar content \* **ITU-P R.837-6:**
Characteristics of precipitation for propagation modelling \* **ITU-P
R.838-3:** Specific attenuation model for rain for use in prediction
methods \* **ITU-P R.839-4:** Rain height model for prediction methods.
\* **ITU-P R.840-6:** Attenuation due to clouds and fog \* **ITU-P
R.1144-7:** Interpolation methods for the geophysical properties used to
compute propagation effects \* **ITU-P R.1511-1:** Topography for
Earth-to-space propagation modelling \* **ITU-P R.1853-1:** Tropospheric
attenuation time series synthesis

Usage
-----

.. code:: python

    import itur
    from astropy import units as u
    f = 86 * u.GHz
    lat, lon = regular_lat_lon_grid() # Produces grids of latitude and longitudes points with 0.1 degrees in resolution

    Att = itur.compute_total_attenuation(lat, lon, f, el) 

.. |GitHub license| image:: https://img.shields.io/badge/license-MIT-lightgrey.svg
   :target: https://raw.githubusercontent.com/Carthage/Carthage/master/LICENSE.md
