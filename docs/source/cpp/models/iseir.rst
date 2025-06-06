IDE SEIR model
==============

This model is based on integro-differential equations. The four compartments:

- **Susceptible**: may become exposed at any time.
- **Exposed**: becomes infected after some time.
- **Infected**: will recover after some time.
- **Recovered**

are used to simulate the spread of the disease.

An example can be found in the
`examples/ide_seir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ide_seir.cpp>`_.


Overview of the ``iseir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::iseir