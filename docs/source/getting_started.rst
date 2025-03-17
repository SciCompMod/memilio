Getting Started
===============

MEmilio supports various model types including equation-based, agent-based, and hybrid graph-ODE-based models. Among the equation-based models, we provide ordinary differential equation (ODE), linear chain trick (LCT), integro-differential equation (IDE) based models.

- The C++ backend powers the model and simulation components for optimal efficiency.
- Data acquisition, plotting, and machine-learned models are handled via Python.

For more details on the C++ implementation, see the `cpp directory README <cpp/README.md>`_.

Some regularly used data for simulations of a pathogen's spread in Germany, like contact and inter-county mobility, can be found in the `data directory README <data/README.md>`_.

In `pycode`, different MEmilio Python packages are defined. Via our `memilio-simulation` package, you can run our C++ backend from Python; this package uses `pybind11` to bind our C++ model code. The `memilio-epidata` package provides tools to download and structure important data such as infection or mobility data. More about the Python packages can be found in the `Python README <pycode/README.md>`_.



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

C++ Tutorial
~~~~~~~~~~~~~~~~~~~~

TBD


For example:

>>> import memilio.epidata import progress_indicator


Other examples can be found in the :doc:`models/index` page.


