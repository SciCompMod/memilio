MEmilio Simulation Package
==========================

**Tutorials**:
:doc:`Model Usage <model_usage>` |
:doc:`Model Creation <model_creation>` |
:doc:`Common Patterns <python_bindings_patterns>`

MEmilio Simulation contains Python bindings for the MEmilio C++ library. It enables setting up and running simulations from Python code.

The package is contained inside the folder `pycode/memilio-simulation <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation>`_.


Dependencies
------------

Required python packages:

-

Stubs
-----

A stub file contains a skeleton of the public interface of that Python 
module including classes, variables, functions and their types. 
They help by enabling autocompletes and type annotations. 
mypy stubgen is used to generate the stub files for the MEmilio Python Bindings 
and provide them as a separate stubs-only package.


For installing stubs you first need to install our package memilio.simulation 
(not in editable mode -e) and the external dependency mypy for your python interpreter.
Then run generate_stubs.py to generate the stubs-only package and 
install it as memilio-stubs, e.g. from the package folder

.. code-block:: console

    python ./tools/generate_stubs.py

\setcounter{tocdepth}{1}

.. toctree::
    :maxdepth: 1
    :hidden:
    
    Model Usage <model_usage>
    Model Creation <model_creation>
    Common Patterns <python_bindings_patterns>