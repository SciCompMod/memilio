Python Packages
===============

MEmilio contains a plethora of python modules containing tools to expand on the main C++.
Most of them serve their own use case,

.. grid:: 1 1 2 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: ../../memilio-small.png
        :text-align: center

        MEmilio Python Bindings
        ^^^

        This module provides a python interface for parts of the C++ main library,
        with the goal of exposing fast mathematical-epidemiological models to
        a bigger user base.

        +++

        .. button-ref:: memilio_simulation
            :expand:
            :color: secondary
            :click-parent:

            To the python bindings

    .. grid-item-card::
        :img-top: ../../memilio-small.png
        :text-align: center

        Epidata Tool
        ^^^

        The memilio-epidata package provides tools to download and structure important 
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
        :text-align: center

        Surrogate Models
        ^^^

        Expanding on AI

        +++

        .. button-ref:: memilio_surrogate
            :expand:
            :color: secondary
            :click-parent:

            To the intro of surrogate models

    .. grid-item-card::
        :text-align: center

        Visualization
        ^^^

        Plot of data

        +++

        .. button-ref:: memilio_plot
            :expand:
            :color: secondary
            :click-parent:

            To the visualization
   
    .. grid-item-card::
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


For developement of code use this command instead

.. code-block:: console 
    
    pip install -e .[dev]

The dependencies are denoted in the documentation of each package.