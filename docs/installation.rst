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
To install ITU-Rpy from command line, please type the following commands on the
command line:

.. code-block:: bash

    git clone https://github.com/iportillo/ITU-Rpy
    cd ITU-Rpy
    cat requirements.txt | xargs -n 1 -L 1 pip install
    python setup.py install
	
	
Installing Basemap
------------------

`Basemap <https://matplotlib.org/basemap/users/intro.html>`_ can be used to plot results in maps. Installation of Basemap is optional, and ITU-Rpy will still work without it. However, some plotting capabilities will be deactivated. A quick overview of `Basemap <https://matplotlib.org/basemap/users/intro.html>`_ if provided below:

	The `matplotlib basemap toolkit <https://matplotlib.org/basemap/users/intro.html>`_ is a library for plotting 2D data on maps in Python. Basemap does not do any plotting on it’s own, but provides the facilities to transform coordinates to one of 25 different map projections. Matplotlib is then used to plot contours, images, vectors, lines or points in the transformed coordinates. 

To install Basemap from pypi, please use the following command on the command
line:

.. code-block:: bash

    pip install basemap

	
If that does not work, you can try to download it using conda:

.. code-block:: bash

    conda -c conda-forge install basemap

