
SDE SIR compartment model
============================

This model is a simple stochastic SDE model with only three compartments and few parameters, mostly for demonstration of the MEmilio framework:

* Susceptible, may become infected at any time.
* Infected, will be recovered after some time.
* Recovered, recovered from previous infection, permanently immune. (Sometimes also Removed for Dead or Recovered together)

We assume simulations over short periods of time, so that the population size can be considered constant and birth as well as (natural) mortality rates can be ignored. 

Important note on the solution: The numerical integration method is a Euler-Maruyama which reuses the EulerIntegratorCore of MEmilio accordingly. The (re)use of other implemented integrators for ODEs is neither considered nor suggested at the moment.

Below is an overview of the model architecture and its compartments.

.. image:: https://martinkuehn.eu/research/images/sir.pdf
    :alt: SIR_model

=====================   ====================================      =====================================================================
Mathematical variable   C++ variable name                         Description
=====================   ====================================      =====================================================================
:math:`\phi`            ``ContactPatterns``                       Daily contact rate / Number of daily contacts.
:math:`\rho`            ``TransmissionProbabilityOnContact``      Transmission risk for people located in the Susceptible compartment.
:math:`N`               ``populations.get_total()``               Total population.
:math:`T_{I}`           ``TimeInfected``                          Time in days an individual stays in the Infected compartment.
=====================   ====================================      =====================================================================



An example can be found in `sde_sir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/sde_sir.cpp>`_

Overview of the ``SSIR`` namespace:
------------------------------------

.. doxygennamespace:: mio::ssir
