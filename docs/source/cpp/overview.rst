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

Build System
-------------

The project uses CMake as a build system with various configuration options such as:

- ``MEMILIO_BUILD_TESTS``: Enables building unit tests
- ``MEMILIO_BUILD_EXAMPLES``: Enables building example applications
- ``MEMILIO_BUILD_MODELS``: Enables building model libraries
- ``MEMILIO_ENABLE_OPENMP``: Enables multithreading with OpenMP

For more details on the build system, refer to TBD
