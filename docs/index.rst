
ITU-Rpy documentation
=====================
ITU-RPy is a python implementation of the ITU-P R Recommendations to compute 
atmospheric attenuation for Earth-to-space and horizontal paths,
for frequencies in the GHz range.

The propagation loss on an Earth-space path and a horizontal-path, relative to the free-space loss, is the sum of different contributions, namely: attenuation by atmospheric gases; attenuation by rain, other precipitation and clouds; scintillation and multipath effects; attenuation by sand and dust storms. Each of these contributions has its own characteristics as a function of frequency, geographic location and elevation angle. ITU-Rpy allows for fast, vectorial computation of the different contributions to the atmospheric attenuation.

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
      
   license
   contact

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


ITU-Rpy is mainly written in Python 3 and continuously tested with Python 3.5-3.9.

ITU-Rpy has the followind dependencies: `numpy`, `scipy`, `joblib`, `pyproj`, and `astropy`. Installation of `basemap` and `matplotlib` is recommended to display results in a map.


