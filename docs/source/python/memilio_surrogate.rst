MEmilio Surrogate Model Package
===============================

MEmilio Surrogate Model contains machine learning based surrogate models that make predictions based on the MEmilio simulation models. 
Currently there are only surrogate models for ODE models. These simulations of these equation-based models are used for data generation. 
The goal is to create a powerful tool that predicts the dynamics faster than a simulation of an expert model, 
e.g., a metapopulation or agent-based model while still having acceptable errors with respect to the original simulations.
 
The package is contained inside the folder `pycode/memilio-surrogatemodel <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-surrogatemodel>`_.

Dependencies
------------

Required python packages:


Since we are running simulations to generate the data, the MEmilio :doc:`memilio-simulation <memilio_simulation>` package
also needs to be installed.

Usage
-----

The package currently provides the following modules:

- `models`: models for different specific tasks
   Currently we have the following models: 
   - `ode_secir_simple`: A simple model allowing for asymptomatic as well as symptomatic infection not stratified by age groups.
   - `ode_secir_groups`: A model allowing for asymptomatic as well as symptomatic infection stratified by age groups and including one damping.

     Each model folder contains the following files: 
     - `data_generation`: data generated from expert model simulation.
     - `model`: training and evaluation of the model. 
     - `network_architectures`: multiple network architectures are saved in this file.


- `tests`: this file contains all tests 

Testing
-------

The package provides a test suite in `memilio/surrogatemodel_test <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-surrogatemodel/memilio/surrogatemodel_test>`_.
To run the tests, simply use the following command.

.. code-block:: console 
    
    python -m unittest

