SDE-based SEIR-type model with two variants
===========================================

This module models and simulates the epidemic using an SDE-based SEIR-type model approach with two variants. After a recovered infection from the first variant (considered the wild type) you can still get infected with the second variant. Infection with the second variant immunizes you against the first variant. Unlike the agent based model that uses particular agents, this model simulates the spread of a communicable disease in a population with subpopulations being in different compartments. The model also introduces stochasticity compared to an ODE-based model. The compartments are as follows: ``Susceptible``,  ``ExposedV1``, ``ExposedV2``, ``ExposedV1V2``,  ``InfectedV1``, ``InfectedV2``, ``InfectedV1V2``, ``RecoveredV1``, ``RecoveredV2`` and ``RecoveredV1V2``. The compartment names are appended by the relevant variants. The addendum `V1` means that the individual is (or has been) infected with the first variant, the addendum `V2` means that the individual is (or has been) infected with the second variant with no prior infections and the addendum `V1V2` means that the individual is (or has been) infected with the second variant after a successful recovery from the first variant. Only individuals in the infected compartments are infectious.

Below is an overview of the model architecture.

Simulation
------------

The simulation runs in discrete time steps using an Euler-Maruyama integration scheme. The Simulation class handles the parameters and the numerical integrator. It also stores the result.

.. image:: https://github.com/user-attachments/assets/55258e5d-05f5-4b16-93b0-f089f8f70782
   :alt: SEIRVV_model

.. list-table::
   :header-rows: 1

   * - Mathematical variable
     - C++ variable name
     - Description
   * - :math:`\phi`
     - ``ContactPatterns``
     - Daily contact rate / Number of daily contacts.
   * - :math:`\rho_1`
     - ``TransmissionProbabilityOnContactV1``
     - Transmission risk for people located in the Susceptible compartment (susceptible to infection with variant 1).
   * - :math:`\rho_2`
     - ``TransmissionProbabilityOnContactV2``
     - Transmission risk for people located in the Susceptible compartment or in the RecoveredV1 compartment (susceptible to infection with variant 2).
   * - :math:`N`
     - ``populations.get_total()``
     - Total population.
   * - :math:`T_{E,1}`
     - ``TimeExposed1``
     - Average time in days an individual stays in the ExposedV1 compartment.
   * - :math:`T_{E,2}`
     - ``TimeExposed2``
     - Average time in days an individual stays in the ExposedV2 or in the ExposedV1V2 compartment.
   * - :math:`T_{I,1}`
     - ``TimeInfected1``
     - Average time in days an individual stays in the InfectedV1 compartment.
   * - :math:`T_{I,2}`
     - ``TimeInfected2``
     - Average time in days an individual stays in the InfectedV2 or in the InfectedV1V2 compartment.

An example can be found in `examples/sde_seirvv.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/sde_seirvv.cpp>`_


Overview of the ``sseirvv`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::sseirvv