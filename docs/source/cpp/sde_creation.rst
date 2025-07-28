SDE model creation
==================

The mathematical model
----------------------

We assume the following form for our stochastic models:

.. math::

    \begin{aligned}
        \mathrm{d}X_t = a(X_t, t) \mathrm{d}t + b(X_t, t)\mathrm{d}W_t
    \end{aligned}

with initial values :math:`X_0 = x_0` where :math:`W_t` denotes the Wiener process.

Analogously to the ODE models, :math:`X_t` describes the sizes of compartments at time point :math:`t`.
The first summand of the equation describes the deterministic part of the model and the second summand represents a
noise term.

How to define an SDE model
-------------------------

To define an SDE model in MEmilio, you define a **StochasticModel**. It is derived from either a **CompartmentalModel**
or **FlowModel** which differ in the methods to define the right hand side of the deterministic part of the mathematical
model. For both choices you need to define a list of **InfectionState**\s, **Parameter**\s and initial conditions via a 
**Population**. We refer to the corresponding :doc:`ODE model <ode_creation>` documentation for more details.
The stochastic part of the model is determined in a :code:`get_noise()` function of the **StochasticModel**, which 
represents the result of the deterministic noise matrix multiplied by a white noise vector. 

Noise contributions
~~~~~~~~~~~~~~~~~~~

In infectious disease models, the noise matrix typically assigns noise contributions from transitions to the respective
compartments. In that case, the white noise vector usually has the same size as the number of transitions. The resulting
noise vector must have the same dimension as the state vector, since the integration of stochastic models is performed 
on the compartments. This is necessary, because the applied noise can occasionally push compartments into negative 
values. While this can be mitigated by removing negative values and rescaling the population (see `map_to_nonnegative`),
this mitigation can not (in general) be applied to flows. However, the noise is applied per transition, with inflows and
outflows of a compartment using the same random value to conserve the total population.
