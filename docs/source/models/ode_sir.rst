ODE SIR Compartment Model
=========================

This model is a very simple ODE model with only three compartments and few parameters, mostly for demonstration of the MEmilio framework:

- **Susceptible**: may become infected at any time.
- **Infected**: will be recovered after some time.
- **Recovered**: recovered from the infectious process (dead or recovered).

We assume simulations over short periods of time, so that the population size can be considered constant and birth as well as (natural) mortality rates can be ignored.

Below is an overview of the model architecture and its compartments.

.. image:: https://github.com/SciCompMod/memilio/assets/69154294/01c9a2ae-2f5c-4bad-b7f0-34de651f2c73
   :alt: SIR_model

+-----------------------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------+
| Mathematical variable       | C++ variable name                             | Description                                                                                   |
+=============================+===============================================+===============================================================================================+
| :math:`\phi`               | ``ContactPatterns``                           | Daily contact rate / Number of daily contacts.                                                |
+-----------------------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------+
| :math:`\rho`               | ``TransmissionProbabilityOnContact``        | Transmission risk for people located in the Susceptible compartment.                          |
+-----------------------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------+
| :math:`N`                  | ``populations.get_total()``                   | Total population.                                                                             |
+-----------------------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------+
| :math:`T_{I}`              | ``TimeInfected``                              | Time in days an individual stays in the Infected compartment.                                |
+-----------------------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------+

An example can be found in the
`examples/ode_sir.cpp <https://github.com/SciCompMod/memilio/tree/main/cpp/examples/ode_sir.cpp>`_.


Overview of the ``osir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::osir