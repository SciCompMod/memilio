Tutorial: Expanding python bindings
===================================

.. toctree::
   :maxdepth: 1


This tutorial is aimed at explaining the process of adding new functionalities to the simulation package 
through writing new python bindings. It is advantagous to have knowledge on programming in C++
as the bindings, as well as the functional source code we want to call from python, is written in C++.

We will go over the common workflow of writing bindings for MEmilio, explaining the structure of the
bindings, providing a small step-by-step guide on adding a new model and showcasing common patterns you 
may encounter with the C++ MEmilio library.

Lets start by looking at the 

Adding a new model
------------------

If currently a model is not available in the python bindings or you added a new c++ model and want to bind it
then the following steps should give an overview of what needs to be done:

#. Add new bindings including a model file in the [models](memilio/simulation/bindings/models/) folder that defines the new module
#. Add a python module file similar to the other models, e.g. [osir.py](memilio/simulation/osir.py) (needed for the structure of the python package) and modify getter function in [__init__.py](memilio/simulation/__init__.py)
#. Add the new module to the building process by modifying [CMakeLists.txt](CMakeLists.txt)
#. Write new tests and examples
#. (Add the package to the documentation and stubs generation)

- talk about memilio generation :doc:`memilio-generation <memilio_generation>`

Expanding the simulation packages
---------------------------------


Pure python additions
---------------------

The bindings can also be expanded with pure Python code by expanding the existing python files or adding new ones underneath the [simulation](memilio/simulation/) folder. 
