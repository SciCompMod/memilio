SDE SIRS Compartment Model
==========================

This model is a simple stochastic SDE SIRS model with only three compartments and few parameters, addressing waning immunity and allowing reinfection after some time of immunity.

- **Susceptible**: may become infected at any time.
- **Infected**: will be recovered after some time.
- **Recovered**: recovered from previous infection, temporarily immune.

We assume simulations over short periods of time, so that the population size can be considered constant and birth as well as (natural) mortality rates can be ignored.

Important note on the solution: The numerical integration method is an Euler-Maruyama which reuses the EulerIntegratorCore of MEmilio accordingly. The (re)use of other implemented integrators for ODEs is neither considered nor suggested at the moment.

Below is an overview of the model architecture and its compartments.

.. image:: https://github.com/SciCompMod/memilio/assets/69154294/bdde5ca3-cd7d-4695-8704-03a712c92ff7
   :alt: SIR_model

+-------------------------------+-----------------------------------------------+--------------------------------------------------------------------------------------------------+
| Mathematical variable         | C++ variable name                             | Description                                                                                      |
+===============================+===============================================+==================================================================================================+
| :math:`\phi`                  | ``ContactPatterns``                           | Daily contact rate / Number of daily contacts.                                                   |
+-------------------------------+-----------------------------------------------+--------------------------------------------------------------------------------------------------+
| :math:`\rho`                  | ``TransmissionProbabilityOnContact``          | Transmission risk for people located in the Susceptible compartment.                             |
+-------------------------------+-----------------------------------------------+--------------------------------------------------------------------------------------------------+
| :math:`N`                     | ``populations.get_total()``                   | Total population.                                                                                |
+-------------------------------+-----------------------------------------------+--------------------------------------------------------------------------------------------------+
| :math:`T_{I}`                 | ``TimeInfected``                              | Time in days an individual stays in the Infected compartment.                                    |
+-------------------------------+-----------------------------------------------+--------------------------------------------------------------------------------------------------+
| :math:`T_{R}`                 | ``TimeImmune``                                | Time in days an individual stays in the Recovered compartment before becoming Susceptible again. |
+-------------------------------+-----------------------------------------------+--------------------------------------------------------------------------------------------------+

An example can be found in the
`examples/ode_sir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_sir.cpp>`_.


Overview of the ``ssirs`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::ssirs
