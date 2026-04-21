Aggregated models
=================

There are different equation-based or compartmental models implemented in MEmilio that consider an aggregated population as shown in the following figure.

.. image:: http://martinkuehn.eu/research/images/aggregated.png
   :alt: Overview on MEmilio's models using aggregated populations.
   :width: 100%

A key distinction between the model types lies in their assumptions about the dwelling time distribution of individuals within a compartment. ODE-based models implicitly assume exponentially distributed dwelling times, leading to memoryless transitions between compartments. The Linear Chain Trick (LCT) and its generalization (GLCT) allow Erlang or Gamma distributions, providing more realistic, non-exponential dwelling times. IDE-based models offer the greatest flexibility, supporting arbitrary transition distributions. SDE models extend ODE models by introducing stochastic fluctuations.

.. toctree::
    :maxdepth: 1

    Ordinary differential equations (ODE) <ode>
    Linear Chain Trick (LCT) <lct>
    Generalized Linear Chain Trick (GLCT) <glct>
    Integro-differential equations (IDE) <ide>    
    Stochastic-differential equations (SDE) <sde>