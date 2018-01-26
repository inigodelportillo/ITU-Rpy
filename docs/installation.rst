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