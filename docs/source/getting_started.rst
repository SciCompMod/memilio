Getting Started
===============

Overview
-------------

MEmilio supports various model types including equation-based, agent-based, and hybrid graph-ODE-based models. Among the equation-based models, we provide models based on :doc:`ordinary differential equations <cpp/ode>`, :doc:`the linear chain trick, <cpp/lct>` and its :doc:`generalisation <cpp/glct>`, :doc:`integro-differential equations <cpp/ide>`, as well as :doc:`stochastic differential equations <cpp/sde>`. 

- The C++ backend powers the model and simulation components for optimal efficiency.
- Data acquisition, plotting, and machine-learned models are handled via Python.

For more details on the C++ implementation, see the sections on :doc:`model usage <cpp/model_usage>` and :doc:`model creation <cpp/model_creation>` in the documentation for the C++-interface.

We also provide several MEmilio Python packages. Via our :doc:`memilio-simulation <python/memilio_simulation>` package, you can run our C++ backend from Python; this package uses ``pybind11`` to bind our C++ model code. The :doc:`memilio-epidata <python/memilio_epidata>` package provides tools to download and structure important data such as infection or mobility data. More about the Python packages can be found in the :doc:`Python Interface Section <python/python_packages>` of this documentation.

A few things are not represented in this documentation, but part of the `github repository <https://github.com/SciCompMod/memilio>`_. In the ``data`` folder you can find some regularly used data for simulations of a pathogen's spread, currently mostly for Germany. 

Visualizations
-----------------

For visualizations we first of all recommend our :doc:`python package <python/memilio_plot>`. Apart from that we have collected some scripts that we used for visualizations in the `tools folder in our github repository <https://github.com/SciCompMod/memilio/tree/main/tools>`_.


.. note::

   This project is under active development.


Usage
__________

.. _installation:

Installation
~~~~~~~~~~~~~

To use MEmilio, first install it using ...

.. code-block:: console

   c++ init


