**MEmilio** - a high performance Modular EpideMIcs simuLatIOn software
==========================================================================

.. image:: https://github.com/SciCompMod/memilio/actions/workflows/main.yml/badge.svg?branch=main
   :target: https://github.com/SciCompMod/memilio/actions/workflows/main.yml
   :alt: CI Badge

.. image:: https://codecov.io/gh/SciCompMod/memilio/branch/main/graph/badge.svg?token=DVQXIQJHBM
    :target: https://codecov.io/gh/SciCompMod/memilio


Welcome
===============

.. attention::

   This documentation is a work in progress. Some areas are already quite detailed, others are still missing completely.


MEmilio implements various models for infectious disease dynamics, ranging from simple compartmental models to complex Integro-Differential and agent-based models. Its modular design enables the combination of different models with distinct mobility patterns. Through efficient implementation and parallelization, MEmilio delivers cutting-edge and compute-intensive epidemiological models at a large scale, providing precise and high-resolution spatiotemporal infectious disease dynamics. MEmilio is continuously extended and is available open-source for community use.

.. image:: https://github.com/user-attachments/assets/65af6012-106e-43c6-9e0e-c96a73aa7b1e
   :alt: Memilio_overview

If you use MEmilio, please :doc:`cite our work<citation>`.


.. note::

   This project is under active development.

Contents
=========

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   getting_started
   citation
   references
   faq
   development

.. toctree::
   :maxdepth: 2
   :caption: C++ Interface

   cpp/overview
   cpp/interfaces
   cpp/installation
   cpp/model_usage
   cpp/model_creation
   cpp/development

.. toctree::
   :maxdepth: 2
   :caption: Python Interface

   python/python_packages
   python/memilio_simulation
   python/memilio_epidata
   python/memilio_surrogate
   python/memilio_generation
   python/memilio_plot

.. toctree::
   :maxdepth: 1
   :caption: Code API 

   pythonapi/pythonapi
   api/library_root
