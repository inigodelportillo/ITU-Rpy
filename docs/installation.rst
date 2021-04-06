Installation
============

.. _installation_pypi:

Installation from pypi
----------------------
To install ITU-Rpy from pypi, please use the following command on the command
line:

.. code-block:: bash

    pip install itur
    
.. _manual_installation:

Manual Installation
-------------------
To install the development version of ITU-Rpy, please type the following commands on the command line:

.. code-block:: bash

    git clone https://github.com/inigodelportillo/ITU-Rpy
    cd ITU-Rpy
    pip install -U -r requirements.txt
    python setup.py install
	
	
Installing Cartopy
------------------

`Cartopy <https://scitools.org.uk/cartopy/docs/latest/installing.html#installing>`_ can be used to plot results in maps. Installation of Cartopy is optional, and ITU-Rpy will still work without it. However, some plotting capabilities will be deactivated. A quick overview of `Cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ if provided below:

`Cartopy <https://scitools.org.uk/cartopy/docs/latest/installing.html#installing>`_ is a Python package designed for geospatial data processing in order to produce maps and other geospatial data analyses. Cartopy has the ability to transform points, lines, vectors, polygons and images between different projections, and it can be combined with  Matplotlib to plot contours, images, vectors, lines or points in the transformed coordinates. 

To install Cartopy from pypi, please use the following command on the command line:

.. code-block:: bash

    pip install cartopy

	
If that does not work, you can try to download it using conda:

.. code-block:: bash

    conda -c conda-forge install cartopy

If you are using Windows, you can also install cartopy using the appropriate pre-compiled wheels file from `this webpage <https://www.lfd.uci.edu/~gohlke/pythonlibs/#cartopy>`_. 
After downloading the .whl file, cartopy can be installed running:


.. code-block:: bash

	pip install name_of_whl_file.whl