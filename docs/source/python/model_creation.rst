How to: Expanding Python bindings
===================================

.. toctree::
   :maxdepth: 1


This tutorial is aimed at explaining the process of adding new functionalities to the simulation package 
through writing new Python bindings. It is advantagous to have knowledge on programming in C++
as the bindings, as well as the functional source code we want to call from Python, is written in C++.

We will go over the common workflow of writing bindings for MEmilio, explaining the structure of the
bindings and providing a small step-by-step guide on adding a new model. 
For a depper look at common patterns you may encounter with the C++ MEmilio library and that are present in the bindings look into :doc:`Common patterns in Python bindings <python_bindings_patterns>`.

Expanding the simulation packages
---------------------------------

All bindings are located inside the folder `bindings <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation/memilio/simulation/bindings>`_. The simulation package is defined over modules. A main module `memilio.simulation` contains the general functionalities of the MEmilio library like the TimeSeries class, different integrators, dampings or IO functionen. 
Additionally each model is located into its own submodule, where model specific bindings are written. 

Python doesn't provide an interface for templates. Therefore, templated functions or classes need to be defined explicitly for each template argument that should be usable.
To reduce overhead, the bindings provide functions to bind templated classes and function. The binding function gets the same template arguments as the class and defines the interface for the class.
As an example, you can look at `bind_CompartmentalModel() <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation/memilio/simulation/bindings/compartments/compartmental_model.h>`_.

Pure Python additions
---------------------

The bindings can also be expanded with pure Python code by expanding the existing Python files or adding new ones underneath the `simulation <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation/memilio/simulation/>`_ folder. 
For those extensions follow the standard Python coding guidelines. 

Adding a new model
------------------

If currently a model is not available in the Python bindings or you added a new c++ model and want to bind it
then the following steps should give an overview of what needs to be done:

* Add new bindings including a model file in the `models <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation/memilio/simulation/bindings/models/>`_folder that defines the new module
* Add a Python module file similar to the other models, e.g. `osir.py <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation/memilio/simulation/osir.py>`_ (needed for the structure of the Python package) and modify getter function in `__init__.py <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation/memilio/simulation/__init__.py>`_
* Add the new module to the building process by modifying `CMakeLists.txt <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation/CMakeLists.txt>`_
* Write new tests and examples
* (Add the package to the documentation and stubs generation)

MEmilio also provides the package :doc:`memilio-generation <memilio_generation>` for automatic generation of model specific python bindings.