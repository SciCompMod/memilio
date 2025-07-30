MEmilio Simulation
===================

MEmilio Simulation is a Python interface to the MEmilio C++ library. Using python bindings, it allows to specify and run original C++ routines directly from Python. It can be used with basic Python skills and without prior knowledge of C++.

The package is contained inside the folder `pycode/memilio-simulation <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation>`_.

Dependencies
------------

Required Python packages:

* scikit-build

For a successful build, the development libraries for Python need to be installed, i.e. python3.x-dev. 
Additionally, as this package builds upon the MEmilio C++ library, 
all dependencies of the main library need to be met. Read more about
the C++ dependencies at :doc:`cpp  <../getting_started>`. 

Main directory structure
------------------------

The main directory structure in the ``pycode/memilio-simulation`` directory includes:

- **memilio/**: Folder containing the source code
  
    - **simulation/**: Contains the core framework for developing epidemiological models

        - **bindings/**: Binding code for creating the MEmilio Simulation package

            - **compartments/**: Classes for compartment models and simulations
            - **epidemiology/**: Base classes for epidemiological modeling
            - **io/**: Input/output utilities for various formats
            - **math/**: Mathematical utilities such as integrators (Euler, RK)
            - **mobility/**: Different Metapopulation mobility approaches
            - **models/**: Defines submodules for model specific bindings
            - **utils/**: General helper functions

    - **simulation_test/**: Unit tests for framework and models

- **tools/**: Additional scripts for memilio-simulation

Package structure
-----------------

The package ``memilio.simulation`` and its submodules aim to mirror the namespaces of the C++ library. 
The main module contains general code for modelling infectious disease, while each submodule contains a specific model.
The overall package structure reflects the directory layout of the Python files, starting from the memilio-simulation <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation>_ folder.
These Python files may include native Python extensions, while the binded C++ code is imported as binary modules into these.

Usage
-----

For understanding the package you can follow the introductions on model usage and creation.
Additionally, multiple `examples <https://github.com/SciCompMod/memilio/blob/main/pycode/examples/simulation>`_ for the 
different features are provided. Lastly, this documentation provides tutorials on applications for the MEmilio Python interface,
coupling it with other libraries, showing the advantages of in language simulation models.

**How to**:
:doc:`Usage of Python bindings <m-simulation_model_usage>` |
:doc:`Expanding Python bindings <m-simulation_expanding_bindings>` |
:doc:`Common patterns <m-simulation_common_patterns>`

.. toctree::
    :maxdepth: 1
    :hidden:
    
    Model Usage <m-simulation_model_usage>
    Model Creation <m-simulation_expanding_bindings>
    Common Patterns <m-simulation_common_patterns>

Stubs
-----

A stub file contains a skeleton of the public interface of that Python 
module including classes, variables, functions and their types. 
They help by enabling autocompletes and type annotations. 
The package `mypy stubgen` is used to generate the stub files for the MEmilio Python Bindings 
and provide them as a separate stubs-only package.


For installing stubs you first need to install our package `memilio.simulation` 
(not in editable mode -e) and the external dependency mypy for your Python interpreter.
Then run `generate_stubs.py` to generate the stubs-only package and 
install it as memilio-stubs, e.g. from the package folder

.. code-block:: console

    python ./tools/generate_stubs.py


Testing
-------

The package provides a test suite in `memilio/simulation_test <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation/memilio/simulation_test>`_.
To run the tests, simply use the following command.

.. code-block:: console 
    
    python -m unittest
