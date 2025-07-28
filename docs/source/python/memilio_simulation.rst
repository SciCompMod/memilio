MEmilio Simulation Package
==========================

MEmilio Simulation contains Python bindings for the MEmilio C++ library. 
It enables setting up and running simulations of infectious disease dynamics from a Python interface.

The package is contained inside the folder `pycode/memilio-simulation <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation>`_.

Dependencies
------------

Required Python packages:

* scikit-build

For a successful build, the development libraries for Python need to be installed, i.e. python3.x-dev. 
Additionally, as this package builds upon the MEmilio C++ library, 
all dependencies of the main library need to be met. Read more about
the C++ dependencies at :doc:`cpp  <../getting_started>`. 

Usage
-----

For understanding the package you can follow the introductions on model usage and creation.
Additionally, multiple `examples <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation>`_ for the 
different features are provided. Lastly, this documentation provides tutorials on applications for the MEmilio Python interface,
coupling it with other libraries, showing the advantages of in language simulation models.

**Tutorials**:
:doc:`Model Usage <model_usage>` |
:doc:`Model Creation <model_creation>` |
:doc:`Common Patterns <python_bindings_patterns>`

.. toctree::
    :maxdepth: 1
    :hidden:
    
    Model Usage <model_usage>
    Model Creation <model_creation>
    Common Patterns <python_bindings_patterns>

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
