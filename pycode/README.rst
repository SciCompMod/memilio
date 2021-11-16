MEmilio Python packages
=======================

This directory collects the python packages for the MEmilio project. Currently, it contains the following packages:

* ``memilio-epidata``: Contains scripts to download different kinds of data as RKI, John Hopkins, Population, DIVI. See `epidata README <memilio-epidata/README.rst>`_.

* ``memilio-simulation``: Contains a python/pybind11 wrapper to the MEmilio C++ library. See `simulation README <memilio-simulation/README.md>`_.

See the corresponding directory for detailed installation and usage instruction.

We recommend to use a virtual python environment to avoid dependency conflicts with other installed packages. On Linux, use the following commands to create and activate a virtual environment in directory called `epi_venv`:

.. code:: sh

    python -m venv epi_venv
    source epi_venv/bin/activate
    pip install --upgrade pip

Refer to the `Python documentation <https://docs.python.org/3/library/venv.html>`_ for more information about virtual environments.

Each package provides its own ``setup.py``, so they can be installed seperately with pip:

.. code:: sh

    cd <package directory>
    pip install .

Each package contains a test suite that can be run with

.. code:: sh

    cd <package directory>
    python -m unittest

If you are a developer, use ``pip install -e .[dev]`` instead when you are working on the code. This also installs additional dependencies used in the unit test suite or for maintenance.