Getting Started
===============

Overview
-------------

.. note:: This project is under active development.


MEmilio is an extensive framework for tasks around infectious disease modelling. It supports various :ref:`model <model-faq>` types 
including :doc:`equation-based<cpp/aggregated_models>`, :doc:`agent-based <cpp/individual_models>`, 
and :doc:`hybrid graph-ODE-based models <cpp/graph_metapop>` as well as data integration and visualizations. 
Among the equation-based models, we provide models based on :doc:`ordinary differential equations <cpp/ode>`,
:doc:`the linear chain trick, <cpp/lct>` and its :doc:`generalisation <cpp/glct>`, :doc:`integro-differential equations <cpp/ide>` 
and :doc:`stochastic differential equations <cpp/sde>`. The MEmilio framework is written in two languages: C++ and Python. 

- The C++ backend powers the model and simulation components for optimal efficiency.
- Data acquisition, plotting, and machine-learning models are handled via Python.

For more details on the C++ implementation, see the sections on :doc:`model usage <cpp/model_usage>` if you are interested 
in using or applying our models and :doc:`model creation <cpp/model_creation>` if you want to write new models inside our framework.

If you prefer using Python, you can use our :doc:`memilio-simulation <python/memilio_simulation>` package to run simulations 
in our C++ backend; this package uses ``pybind11`` to bind our C++ model code. 
The :doc:`memilio-epidata <python/memilio_epidata>` package provides tools to download and structure important data such 
as infection or mobility data. More about this and our other Python packages can be found in the :doc:`Python Interface Section <python/python_packages>` 
of this documentation.

A few things are not represented in this documentation, but part of the `github repository <https://github.com/SciCompMod/memilio>`_. 
In the `data <https://github.com/SciCompMod/memilio/tree/main/data>`_ folder you can find some regularly used data 
for simulations of a pathogen's spread, currently mostly for Germany. 


Usage
-----------------

.. _installation:

Installation
~~~~~~~~~~~~~

You can download the latest code version of MEmilio from our `github repository <https://github.com/SciCompMod/memilio>`_ 
in a Linux terminal via:

.. code-block:: console

   git clone https://github.com/SciCompMod/memilio.git

You can now build the C++ code:

.. code-block:: console

   cd memilio/cpp
   mkdir build && cd build
   cmake ..
   cmake --build .

For details on the possible compile flags, help with errors and general a more detailed instruction, see the 
:doc:`C++-interface <cpp/installation>` section of this documentation. 

For the installation of Python packages, e.g. ``memilio-epidata``, do

.. code-block:: console
   
   cd memilio/pycode
   cd memilio-epidata
   pip install .
   
For more information, we refere to the :ref:`Python Interace Part <Python_Installation>` of this documentation.


Running simulations
~~~~~~~~~~~~~~~~~~~~~
You can run simulations either via the C++ interface where they are originally implemented or via the python bindings. 
For the C++ Interface you can find explanations of the models as well as guides on their usage in the :doc:`C++ model usage <cpp/model_usage>` section.
In short, the executables for different model instatiations are build as described above and can be run via 

.. code-block:: console

   ./cpp/build/bin/<example_name>


Out of the box this works for all examples in the ``cpp/examples`` folder of our `github repository <https://github.com/SciCompMod/memilio/tree/main/cpp/examples>`_,
that do not depend on user-provided external libraries. 
Additional explanations for our models are linked at the corresponding sites of this documentation.

For the Python interface, you can find a short introduction in the :doc:`Python Interface <python/memilio_simulation>` section.

Additionally we provide a python package for :doc:`surrogate models <python/memilio_surrogate>`, which can be used to c
reate fast approximations of our models.

Loading data
~~~~~~~~~~~~~~~~~~~~~
The :doc:`memilio-epidata <python/memilio_epidata>` package provides tools to download epidemiological relevant datasets. Some 
datasets like contact matrices for Germany are also included in the ``data`` folder of the `github repository <https://github.com/SciCompMod/memilio/tree/main/data>`_ and 
school holidays (for Germany) are directly included in the `C++ code <https://github.com/SciCompMod/memilio/blob/main/cpp/memilio/geography/holiday_data.ipp>`_.  


Creating new models
~~~~~~~~~~~~~~~~~~~~~

If you want to create new models, you can do so via the C++ interface. For this, we recommend to have a look at 
the :doc:`C++ model creation <cpp/model_creation>` section of this documentation.


Visualizations
~~~~~~~~~~~~~~~~~~~~~

For visualizations we first of all recommend our :doc:`python package <python/memilio_plot>`. Apart from that we have 
collected some scripts that we used for visualizations in the `tools folder in our github repository <https://github.com/SciCompMod/memilio/tree/main/tools>`_. 
For the latter we don't take any responsibilities!

Further questions
~~~~~~~~~~~~~~~~~~~~~
If you have any further questions, please take a look at our :doc:`faq` and feel free to contact us via `github <https://github.com/ICB-DCM/orga/discussions/categories/q-a>`_.