
ITU-Rpy documentation
=================================
ITU-RPy  is a  python implementation of the ITU-P R Recommendations to compute 
atmospheric attenuation for Earth-to-space and horizontal paths,
for frequencies in the Ghz range.


The propagation loss on an Earth-space path and a horizontal-path, relative to the free-space loss, is the sum of different contributions, namely: attenuation by atmospheric gases; attenuation by rain, other precipitation and clouds; scintillation and multipath effects; attenuation by sand and dust storms. Each of these contributions has its own characteristics as a function of frequency, geographic location and elevation angle. ITUR-py allows for fast, vectorial computation of the different contributions to the atmospheric attenuation.

Contents:
---------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   installation
   quickstart
   api
   faq
   contact
   license

ITU-Rpy is mainly written in Python 3 and continuously tested with Python 3.5-3.6.
It depends on `numpy`, `scipy`, `joblib`, and `astropy`.


.. Indices and tables
.. ------------------
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
