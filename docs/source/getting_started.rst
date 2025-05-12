Getting Started
===============

Overview
-------------

.. note::

   This project is under active development.


MEmilio supports various model types including equation-based, agent-based, and hybrid graph-ODE-based models. Among the equation-based models, we provide models based on :doc:`ordinary differential equations <cpp/ode>`, :doc:`the linear chain trick, <cpp/lct>` and its :doc:`generalisation <cpp/glct>`, :doc:`integro-differential equations <cpp/ide>`, as well as :doc:`stochastic differential equations <cpp/sde>`. 
The MEmilio framework is written in two languages: C++ and Python. 

- The C++ backend powers the model and simulation components for optimal efficiency.
- Data acquisition, plotting, and machine-learned models are handled via Python.

For more details on the C++ implementation, see the sections on :doc:`model usage <cpp/model_usage>` if you are interested in using or applying our models and :doc:`model creation <cpp/model_creation>` if you want to write new models inside our framework.

We also provide several MEmilio Python packages. Via our :doc:`memilio-simulation <python/memilio_simulation>` package, you can run our C++ backend from Python; this package uses ``pybind11`` to bind our C++ model code. The :doc:`memilio-epidata <python/memilio_epidata>` package provides tools to download and structure important data such as infection or mobility data. More about the Python packages can be found in the :doc:`Python Interface Section <python/python_packages>` of this documentation.

A few things are not represented in this documentation, but part of the `github repository <https://github.com/SciCompMod/memilio>`_. In the ``data`` folder you can find some regularly used data for simulations of a pathogen's spread, currently mostly for Germany. 



Usage
-----------------

.. _installation:

Installation
~~~~~~~~~~~~~

You can download the latest code version of MEmilio from our `github repository <https://github.com/SciCompMod/memilio>`_ via:

.. code-block:: console

   git clone https://github.com/SciCompMod/memilio.git

Then you can build the C++ code as described in the :doc:`C++-interface <cpp/installation>` section of this documentation. 
Installation of the Python packages is described in the :ref:`Python Interace Part <Python_Installation>` of this documentation.


Running simulations
~~~~~~~~~~~~~~~~~~~~~
You can run simulations either via the C++ interface where they are originally implemented or via the python bindings. 
For the C++ Interface you can find explanations of the models as well as guides on their usage in the :doc:`C++ model usage <cpp/model_usage>` section.
For the Python interface, you can find a short introduction in the :doc:`Python Interface <python/memilio_simulation>` section.
We also provide more examples in the ``cpp/examples`` folder of our `github repository <https://github.com/SciCompMod/memilio/tree/main/cpp/examples>`_. 
Additional explanations for our models can be found in the :doc:`Models section <models/index>` of this documentation.

Additionally we provide a python package for :doc:`surrogate models <python/memilio_surrogate>`, which can be used to create fast approximations of our models.

Loading data
~~~~~~~~~~~~~~~~~~~~~
The :doc:`memilio-epidata <python/memilio_epidata>` package provides tools to download epidemiological relevant datasets. Some 
datasets like contact matrices for Germany are also included in the ``data`` folder of the `github repository <https://github.com/SciCompMod/memilio/tree/main/data>`_ and school holidays 
(for Germany) are directly included in the `C++ code <https://github.com/SciCompMod/memilio/blob/main/cpp/memilio/geography/holiday_data.ipp>`_.  


Creating new models
~~~~~~~~~~~~~~~~~~~~~
If you want to create new models, you can do so via the C++ interface. For this, we recommend to have a look at the :doc:`C++ model creation <cpp/model_creation>` section of this documentation.


Visualizations
~~~~~~~~~~~~~~~~~~~~~

For visualizations we first of all recommend our :doc:`python package <python/memilio_plot>`. Apart from that we have collected some scripts that we used for visualizations in the `tools folder in our github repository <https://github.com/SciCompMod/memilio/tree/main/tools>`_. For the latter we don't take any responsibilities!

Further questions
~~~~~~~~~~~~~~~~~~~~~
If you have any further questions, please take a look at our :doc:`faq` and feel free to contact us via `github <https://github.com/ICB-DCM/orga/discussions/categories/q-a>`_.