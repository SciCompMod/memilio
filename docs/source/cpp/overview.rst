Overview
============


The MEmilio C++ project is organized as follows:

Main Directory Structure
---------------------------

The main directory structure in the ``cpp`` directory includes:

- **memilio/**: Contains the core framework for developing epidemiological models
  
  - **ad/**: Algorithmic differentiation framework
  - **compartments/**: Classes for compartment models and simulations
  - **data/**: Data analysis functions
  - **epidemiology/**: Base classes for epidemiological modeling
  - **geography/**: Geographic region data (e.g. school vacations) and functions
  - **io/**: Input/output utilities for various formats (JSON, HDF5...)
  - **math/**: Mathematical utilities such as integrators (Euler, RK)
  - **mobility/**: Different Metapopulation mobility approaches
  - **utils/**: General helper functions (logging, etc.)

- **models/**: Concrete implementation of epidemiological models

- **simulations/**: Applications for scientific publications

- **examples/**: Example applications demonstrating the use of the framework

- **tests/**: Unit tests for framework and models

- **thirdparty/**: Configuration of external dependencies

- **benchmarks/**: Analyzing runtime performance

Model Structure
-----------------

The MEmilio library uses a modular organization of models, where generic implementations are inherited by specific implementations:

.. image:: http://martinkuehn.eu/research/images/overview.png
   :alt: Model Hierarchy
   :width: 100%

**CompartmentalModel**: The base class for all compartment-based models in MEmilio. It defines the fundamental structure for epidemiological models with compartments (e.g., SEIR, SECIR) and provides methods like ``eval_right_hand_side`` and ``get_initial_values`` required for ODE solvers.

**FlowModel**: Inherits from CompartmentalModel and extends it with the concept of flows between compartments. Instead of directly defining derivatives, it specifies the flows between compartments.

**Specific Model Implementations**:

- **ODE Model** (Ordinary Differential Equations): Deterministic models for continuous populations described by ordinary differential equations.
  
- **IDE Model** (Integro-Differential Equations): Extends the ODE model integration terms.
  
- **SDE Model** (Stochastic Differential Equations): Adds stochastic components to model uncertainties and random effects.

**Individual-based Model**: Stands separate from the compartmental hierarchy and models each individual explicitly with its own properties and interactions. This enables more detailed simulations.


Build System
-------------

The project uses CMake as a build system with various configuration options. 
For an explanation, refer to :doc:`installation`.
