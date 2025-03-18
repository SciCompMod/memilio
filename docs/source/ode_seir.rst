ODE SEIR Compartment Model
==========================

This model is a very simple ODE model with only four compartments and few parameters, mostly for demonstration of the MEmilio framework:

- **Susceptible**: may become exposed at any time.
- **Exposed**: becomes infected after some time.
- **Infected**: will recover after some time.
- **Recovered**

Below is an overview of the model architecture and its compartments.

.. image:: https://github.com/SciCompMod/memilio/assets/69154294/80a36be5-57d9-4012-9b5f-25eb08ec8837
   :alt: SEIR_model

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Mathematical variable
     - C++ variable name
     - Description
   * - :math:`\phi`
     - ``ContactPatterns``
     - Daily contact rate / Number of daily contacts.
   * - :math:`\rho`
     - ``TransmissionProbabilityOnContact``
     - Transmission risk for people located in the Susceptible compartment.
   * - :math:`N`
     - ``populations.get_total()``
     - Total population.
   * - :math:`T_{E}`
     - ``TimeExposed``
     - Time in days an individual stays in the Exposed compartment.
   * - :math:`T_{I}`
     - ``TimeInfected``
     - Time in days an individual stays in the Infected compartment.

An example can be found in the
`examples/ode_seir.cpp <https://github.com/SciCompMod/memilio/tree/main/cpp/examples/ode_seir.cpp>`_.


Overview of the ``oseir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::oseir