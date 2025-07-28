Model Creation
==============

While MEmilio already preimplements many different models, it is possible to create new models. This section describes how to create a new model in MEmilio. All of MEmilio's models have been designed to share a maximum of structure and functionality, however, the creation of a new aggregated model differs fundamentally from the creation of a new inidividual-based model. New metapopulation models are generally created by implementing an aggregated model in a graph structure.

.. toctree::
   :maxdepth: 1

   Common structure <structure>
   Ordinary differential equations <ode_creation>
   Linear Chain Trick <lct_creation>
   Stochastic-differential equations <sde_creation>
   Integro-differential equations <ide_creation>
    