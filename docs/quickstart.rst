Quick Start
===========

If you have not installed *ITU-Rpy* yet take a look at the `installation instructions <installation.html>`_ and make sure that all the requirements are fulfilled.
Examples to illustrate the usage of *ITU-Rpy* are provided in the examples-folder. In order to understand the models used, check the `models section <apidoc/models.html>`_ .

To get started, we will walk you through a few examples that show the most illustrative case studies for 

* First, we explain the basic usage of *ITU-Rpy* by computing the attenuation at a single location.
* Second, we explain the usage of *ITU-Rpy* using vectorized operations.
* Finally, we summarize other useful atmospheric functions in the library. The whole API description can be accessed through the `API section <api.html>`_ .

Single location attenuation
---------------------------

Here we will compute the link attenuation vs. at a frequency of 22.5 GHz for a link between a satellite in GEO at the orbital slot 77 W and a ground station in Boston.

In addition, we will show how to compute other parameters of interest such as 0 degree isotherm, rainfall intensity exceeded during 0.01 % of the time, total columnar content liquid water or temperature.

First, let's define the coordinates of our ground station and compute the elevation angle of the link

.. code-block:: python

	import itur
	import astropy.units as u
	
	# Ground station coordinates (Boston)
	lat_GS = 42.3601
	lon_GS = -71.0942
	
	# Satellite coordinates (GEO, 77 W)
	lat_sat = 0
	lon_sat = -77 
	h_sat = 35786 * u.km
	
	# Compute the elevation angle between satellite and ground station
	el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat, lat_GS, lon_GS)
	

Next, we define the link parameters

.. code-block:: python

	f = 22.5 * u.GHz    # Link frequency
	D = 1.2 * u.m       # Antenna diameters

Finally, we compute the total atmospheric attenuation as well as the different contributions for a set of unavailability values and plot the results. Note the flag `return_contributions = True` when calling function `itur.atmospheric_attenuation_slant_path`.

.. code-block:: python

	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.ticker import ScalarFormatter

	# Define unavailabilities vector
	unavailabilities = np.logspace(-1.5, 1.5, 100)
	
	# Compute the 
	A_g, A_c, A_r, A_s, A_t = [], [], [], [], []
	for p in unavailabilities:
		a_g, a_c, a_r, a_s, a_t = itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, 
                                                                                  f, el, p, D,
                                                                                  return_contributions=True)
		A_g.append(a_g.value)
		A_c.append(a_c.value)
		A_r.append(a_r.value)
		A_s.append(a_s.value)
		A_t.append(a_t.value)
		
	# Plot the results
	ax = plt.subplot(1,1,1)
	ax.semilogx(unavailabilities, A_g, label='Gaseous attenuation')
	ax.semilogx(unavailabilities, A_c, label='Cloud attenuation')
	ax.semilogx(unavailabilities, A_r, label='Rain attenuation')
	ax.semilogx(unavailabilities, A_s, label='Scintillation attenuation')
	ax.semilogx(unavailabilities, A_t, label='Total atmospheric attenuation')

	ax.xaxis.set_major_formatter(ScalarFormatter())
	ax.set_xlabel('Percentage of time attenuation value is exceeded [%]')
	ax.set_ylabel('Attenuation [dB]')
	ax.grid(which='both', linestyle=':')
	plt.legend()
	
which results in the following plot image:	
	
.. figure:: images/att_single_location.png
   :scale: 80 %
   :align: center   
   :figclass: align-center     
   :alt: attenuation to single ground station

   Atmospheric attenuation at Boston for a link to GEO - 77 W.
   
Note the by default, *ITU-Rpy* returns Quantity type objects, which are based on `astropy.units` module. Quantity objects are special objects that contain a `value` and `unit` attributes. Conversion among units is possible using the `.to()` method.

Atmospheric parameters such as temperature, pressure, or water-vapor density can be passed to function `itur.atmospheric_attenuation_slant_path` manually if known, otherwise *ITU-Rpy* will compute them automatically using the appropriate ITU Recommendation models. Similarly, if the ground station height above mean sea level is known, it can also be introduced manually.


Vectorial operations
--------------------

One of the main characteristics of *ITU-Rpy* is that it allows for broadcasting of operations when using vectors. Let's say that we want to compute the attenuation over Africa of a new Ka-band satellite located in GEO at slot 4 E.


.. code-block:: python

	import itur
	import astropy.units as u
	
	# Generate regular grid of latitude and longitudes with 1 degree resolution
	lat, lon = itur.utils.regular_lat_lon_grid(lat_max=60, 
                                                   lat_min=-60,
                                                   lon_max=65,
                                                   lon_min=-35)
	
	# Satellite coordinates (GEO, 77 W)
	lat_sat = 0
	lon_sat = 4 
	h_sat = 35786 * u.km
	
	# Compute the elevation angle between satellite and ground station
	el = itur.utils.elevation_angle(h_sat, lat_sat, lon_sat, lat, lon)

	f = 22.5 * u.GHz    # Link frequency
	D = 1.2 * u.m       # Antenna diameters
	p = 1               # Unavailability (Values exceeded 1% of time)
	Att = itur.atmospheric_attenuation_slant_path(lat, lon, f, el, p, D)

If you have installed Basemap (see `installation instructions <installation.html>`_ ), you can use function `itur.utils.plot_in_map()` to display the results as an image:

.. code-block:: python

	# Plot the results
	m = itur.utils.plot_in_map(Att.value, lat, lon, 
                                   cbar_text='Atmospheric attenuation [dB]',
                                   cmap='magma')
	
	# Plot the satellite location
	m.scatter(lon_sat, lat_sat, c='white', s=20)
						   

.. figure:: images/att_africa.png
   :scale: 133 %
   :align: center   
   :figclass: align-center   
   :alt: attenuation to single ground station

   Atmospheric attenuation over Africa @ 22.5 GHz.
   
   
Vectorization works in a lot of different ways. For example, let's say that we want to compute the dependency of the attenuation with the elevation angle for a the example presented in the `Single location attenuation example <>`_. We can just do it using the code below

.. code-block:: python

	import itur
	import astropy.units as u
	
	itur.deactivate_output_quantity()
	
	# Ground station coordinates (Boston)
	lat_GS = 42.3601
	lon_GS = -71.0942
	
	# Vectorize the elevation angle
	el = np.linspace(10, 90, 50)
	
	f = 22.5 * u.GHz    # Link frequency
	D = 1.2 * u.m       # Antenna diameters
	p = 1               # Unavailability (Values exceeded 1% of time)
	Att = itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D)
	
	plt.figure()
	plt.plot(el, Att.value)
	plt.xlabel('Elevation angle [deg]')
	plt.ylabel('Attenuation [dB]')
	plt.grid(which='major', linestyle=':')
	
.. figure:: images/elevation_angle.png
   :scale: 70 %
   :align: center   
   :figclass: align-center   
   :alt: attenuation to single ground station

   Atmospheric attenuation at Boston vs. elevation angle.
