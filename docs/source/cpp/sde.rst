SDE-based models
================

These models are a type of aggregated or compartmental model, which is described by a system of initial value problems
(IVP) given by stochastic differential equations (SDE). In MEmilio, they are implemented as an ODE-based model with an
additional function to compute the random noise, as can be seen :doc:`here <sde_creation>`. Hence, for the most part,
SDE models are used exactly like ODE-based models. They mostly differ in how they are simulated, see the :ref:`Simulation SDE`
section below. For everything else, check out the page on :doc:`ODE-based model usage <ode>`.

.. _Simulation SDE:
Simulation
----------

Once the model is set up, one can run a simple simulation from time ``t0`` to ``tmax`` with an initial step size ``dt``
using the  ``mio::simulate_stochastic()`` function. This will run a simulation of type **StochasticSimulation** that
saves the sizes of each compartment over time.
The simulation uses an Euler-Maruyama scheme by default, so the step size does not change over time.

Flow information cannot be obtained even when the **StochasticModel** is defined using a **FlowModel**, as the
integrator may need to rescale results with respect to compartments to avoid negative values.

List of models
--------------
.. toctree::
    :titlesonly:
    
    models/ssir    
    models/ssirs
    models/sseir
