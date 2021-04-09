itur.models package
===================

The ``itur.models`` package contains the implementation of the models described in the different ITU-R P. recommendations.

Individual modules can be imported from ``itur.models`` using ``itu<recomendation_number>``, as shown below:

.. code-block:: python

    import itur.models.itu618

Each module contains two functions, ``get_version()`` and ``change_version()``, that allow to obtain or change the version currently being used.
The script below provides an example for the module that implements Recommendation ITU-R P.618.

.. code-block:: python

    import itur.models.itu618 as itu618
    
    print('Current version of ITU-R P.618: ', itu618.get_version())
    itu618.change_version(12)  # Change to version ITU-R P.618-12
    
By default, each module contains a ``__model`` object that acts as a singleton instance for the recommendation model. For most use cases, it is recommended 
that the developers interact with the model using the publicly documented functions. This way, it is ensured that all functions called belong to the 
same recommendation version, and that the parameters passed to the functions have the right units and format.

However, if a developer wants to instantiate a new ``model`` object (i.e., to compare results from different version, for development purposes), a new object
can be instantiated as shown in the example below:

.. code-block:: python

    import itur.models.itu618 as itu618
    
    model_618_12 = itu618._ITU618(12)
    model_618_13 = itu618._ITU618(13)
    
Package contents
----------------

.. automodule:: itur.models
    :members:
    :undoc-members:
    :show-inheritance:

.. toctree::
   :maxdepth: 1
    
   itu453
   itu530
   itu618
   itu676
   itu835
   itu836
   itu837
   itu838
   itu839
   itu840
   itu1144
   itu1510
   itu1511
   itu1623
   itu1853