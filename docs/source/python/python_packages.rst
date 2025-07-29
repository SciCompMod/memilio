Overview
===============

MEmilio contains a plethora of python modules containing tools to expand on the main C++.
Most of them serve their own use case:

.. grid:: 1 1 2 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: https://github.com/user-attachments/assets/b95b8803-0fdc-4d80-88cc-fcb2643c3e8f
        :text-align: center

        MEmilio Python Bindings
        ^^^

        This package provides a python interface for parts of the C++ main library,
        with the goal of exposing fast mathematical-epidemiological models to
        a bigger user base.

        +++

        .. button-ref:: memilio_simulation
            :expand:
            :color: secondary
            :click-parent:

            To the python bindings

    .. grid-item-card::
        :img-top: https://github.com/user-attachments/assets/ca098013-ec8a-4fe9-964c-f1881b7382de
        :text-align: center

        Epidata Tool
        ^^^

        This package provides tools to download and structure important 
        data such as infection or mobility data.

        +++

        .. button-ref:: memilio_epidata
            :expand:
            :color: secondary
            :click-parent:

            To the data package


.. grid:: 1 1 1 1
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: https://github.com/user-attachments/assets/064a55a0-7054-4421-a026-353f8f4cc478
        :text-align: center

        Surrogate Models
        ^^^

        This package contains machine learning based surrogate models that make predictions based on the MEmilio simulation models. 

        +++

        .. button-ref:: memilio_surrogate
            :expand:
            :color: secondary
            :click-parent:

            To the intro of surrogate models

    .. grid-item-card::
        :img-top: https://github.com/user-attachments/assets/81659df6-826c-4a34-83c0-a20c32d1d266
        :text-align: center

        Visualization
        ^^^

        Generalized visualization functions for MEmilio specific plots.

        +++

        .. button-ref:: memilio_plot
            :expand:
            :color: secondary
            :click-parent:

            To the visualization
   
    .. grid-item-card::
        :img-top: https://github.com/user-attachments/assets/0c23a9a1-5b78-477f-981e-6bb7993fc80d
        :text-align: center

        Generating Bindings
        ^^^

        Easy to use tool for helping with the creation of new bindings of C++ models.

        +++

        .. button-ref:: memilio_generation
            :expand:
            :color: secondary
            :click-parent:

            To the generation package


.. _Python_Installation:

Installation
------------

Each package provides a `setup.py` script that installs the package and its dependencies. 
The installation can be run with the following command (from the directory containing the `setup.py`)

.. code-block:: console 
    
    pip install .


For development of code use this command instead

.. code-block:: console 
    
    pip install -e .[dev]

The dependencies are denoted in the documentation of each package.