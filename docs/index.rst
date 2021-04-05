
ITU-Rpy documentation
=====================

ITU-Rpy is a python implementation of the ITU-R P. Recommendations to compute atmospheric attenuation in slant and horizontal paths.

A complete overview of the contents of this documentation can be found in the :ref:`Table of Contents` at the bottom of this page.

Instructions on how to install ITU-Rpy can be found at the :ref:`Installation` page.

Results of validating ITU-Rpy against the validation examples provided by the ITU (where available) are at the :ref:`Validation` page.

Citation
--------

If you use ITU-Rpy in one of your research projects, please cite it as:

::

    @misc{iturpy-2017,
          title={ITU-Rpy: A python implementation of the ITU-R P. Recommendations to compute 
             atmospheric attenuation in slant and horizontal paths.},
          author={Inigo del Portillo},
          year={2017},
          publisher={GitHub},
          howpublished={\url{https://github.com/inigodelportillo/ITU-Rpy/}}
    }

Usage and examples
------------------

The :ref:`Quick Start` guide provides different examples on how to use ITUR-py. 

Additional examples can be found in the `examples folder <https://github.com/inigodelportillo/ITU-Rpy/tree/master/examples>`_, and the snippet of code below.

.. code-block:: python
     
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
						   
which produces 

.. figure:: https://raw.githubusercontent.com/inigodelportillo/ITU-Rpy/master/docs/images/att_world.png
   :align: center   
   :figclass: align-center   
   :alt: attenuation to single ground station

   Atmospheric attenuation worldmap @ 22.5 GHz.
   
Table of Contents
-----------------

.. toctree::
   :titlesonly:
   :caption: First steps
   :maxdepth: 1

   installation
   quickstart

.. toctree::
   :titlesonly:
   :caption: API documentation
   :maxdepth: 3
      
   api
   validation
   faq
   
.. toctree::
   :titlesonly:
   :caption: Miscellaneous
   :maxdepth: 1
      
   
   contributing
   license
   contact

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Other
-----

ITU-Rpy is mainly written in Python 3 and continuously tested with Python 3.5-3.9.

ITU-Rpy has the following dependencies: `numpy`, `scipy`, `joblib`, `pyproj`, and `astropy`. Installing `basemap` and `matplotlib` is recommended to display results in a map.


