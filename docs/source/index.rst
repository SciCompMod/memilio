**MEmilio** - a high performance Modular EpideMIcs simuLatIOn software
==========================================================================

.. image:: https://github.com/SciCompMod/memilio/actions/workflows/main.yml/badge.svg?branch=main
   :target: https://github.com/SciCompMod/memilio/actions/workflows/main.yml
   :alt: CI Badge

.. image:: https://codecov.io/gh/SciCompMod/memilio/branch/main/graph/badge.svg?token=DVQXIQJHBM
    :target: https://codecov.io/gh/SciCompMod/memilio


**Welcome**
MEmilio implements various models for infectious disease dynamics, ranging from simple compartmental models to complex Integro-Differential and agent-based models. Its modular design enables the combination of different models with distinct mobility patterns. Through efficient implementation and parallelization, MEmilio delivers cutting-edge and compute-intensive epidemiological models at a large scale, providing precise and high-resolution spatiotemporal infectious disease dynamics. MEmilio is continuously extended and is available open-source for community use.

If you use MEmilio, please cite our works:

- Kühn, Martin Joachim et al. (2022). *MEmilio - a High Performance Modular Epidemics Simulation Software (2022)*. Available at `GitHub <https://github.com/SciCompMod/memilio>`_ and `DLR <https://elib.dlr.de/192140/>`_.

- Koslow W, Kühn MJ, Binder S, Klitz M, Abele D, et al. (2022). *Appropriate relaxation of non-pharmaceutical interventions minimizes the risk of a resurgence in SARS-CoV-2 infections in spite of the Delta variant*. *PLOS Computational Biology* 18(5): e1010054. `https://doi.org/10.1371/journal.pcbi.1010054`

Getting Started
===============

MEmilio supports various model types including equation-based, agent-based, and hybrid graph-ODE-based models. Among the equation-based models, we provide ordinary differential equation (ODE) and integro-differential equation (IDE) based models.

- The C++ backend powers the model and simulation components for optimal efficiency.
- Data acquisition, plotting, and machine-learned models are handled via Python.

For more details on the C++ implementation, see the `cpp directory README <cpp/README.md>`_.

Regularly used data for simulations of pathogen spread in Germany, such as contact and inter-county mobility, is provided in the `data directory README <data/README.md>`_.

Python Packages
===============

- **memilio-simulation**: Interface to run the C++ backend via Python using `pybind11 <https://github.com/pybind/pybind11>`_.
- **memilio-epidata**: Tools for downloading and structuring key datasets like infection or mobility data.

More information is available in the `Python README <pycode/README.rst>`_.

Documentation
=============

Each core component is described in detail within its corresponding README file, including configuration and usage instructions for both users and developers.

Code is documented using Doxygen, and further instructions can be found in the `docs` folder.

The latest documentation for the main branch is available here: 
`https://scicompmod.github.io/memilio/documentation/index.html`

Installation, Usage, and Requirements
=====================================

Each module has its own set of requirements and usage instructions, detailed in the corresponding READMEs.

Development
===========

- `Git workflow and change process <https://github.com/SciCompMod/memilio/wiki/git-workflow>`_
- `Coding Guidelines <https://github.com/SciCompMod/memilio/wiki/coding-guidelines>`_
.. image:: https://github.com/user-attachments/assets/65af6012-106e-43c6-9e0e-c96a73aa7b1e
   :alt: Memilio_overview


**MEmilio** is a framework that ..

Check out the :doc:`usage` section for further information, including
how to :hoverxref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::
   :maxdepth: 1

   getting_started
   c++/c++
   python/python
   faq
   models/models
   pythonapi
   api/library_root

