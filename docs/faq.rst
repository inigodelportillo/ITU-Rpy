F.A.Q.
======

I cannot install Basemap
^^^^^^^^^^^^^^^^^^^^^^^^

This happens most likely because you are using python  version > 3.X. You can try to install from conda-forge  ``conda install -c conda-forge basemap`` or, if you are using Windows, using the appropriate pre-compiled wheels file from `this webpage <https://www.lfd.uci.edu/~gohlke/pythonlibs/#basemap>`_. Once you download the .whl file you can install it using ``pip install name_of_whl_file.whl``.

The first time I run *ITU-Rpy* is considerable slower
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*ITU-Rpy* loads in memory several datasets upon first execution. This process might take up to 30 seconds. Once that datasets are loaded into memory *ITU-Rpy* uses cached versions to reduce execution time.

I cannot operate with the values returned by *ITU-Rpy*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*ITU-Rpy* returns Quantity objects, which consist of a value and a unit. Only quantities with compatible dimensions can be added / subtracted.

.. code-block:: python
	
	import itur
	d1 = 300 * itur.u.m
	d2 = 0.2 * itur.u.km
	
	print(d1 + d2)     # prints 500.0 m
	
	p1 = 1013 * itur.u.hPa
	print(d1 + p1)     # Generates an error.


The user can transform between compatible units using the ``.to()`` method.

.. code-block:: python
	
	print(d1.to(itur.u.km))    # prints 0.3 km
	

One can access to the values and units using the ``.value`` and ``.unit`` methods respectively. Some matplotlib functions accept Quantities as inputs (``plt.plot``, ``plt.scatter``), whereas others require plain values (`plt.bar`).
 

I discovered a bug/have criticism or ideas on *ITU-Rpy*. Where should I report to?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*ITU-Rpy* uses the `GitHub issue-tracker <https://github.com/iportillo/ITU-Rpy/issues>`_ to take care of bugs and questions. If you experience problems with *ITU-Rpy*, try to provide a full error report with all the typical information (OS, version, console-output, minimum working example, ...). This makes it a lot easier to reproduce the error and locate the problem.
