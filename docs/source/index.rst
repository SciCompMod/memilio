**MEmilio** - a high performance Modular EpideMIcs simuLatIOn software
==========================================================================

.. image:: https://github.com/SciCompMod/memilio/actions/workflows/main.yml/badge.svg?branch=main
   :target: https://github.com/SciCompMod/memilio/actions/workflows/main.yml
   :alt: CI Badge

.. image:: https://codecov.io/gh/SciCompMod/memilio/branch/main/graph/badge.svg?token=DVQXIQJHBM
    :target: https://codecov.io/gh/SciCompMod/memilio


Welcome
===============

MEmilio implements various models for infectious disease dynamics, ranging from simple compartmental models to complex Integro-Differential and agent-based models. Its modular design enables the combination of different models with distinct mobility patterns. Through efficient implementation and parallelization in C++ and an easy-to-use python interface, MEmilio delivers cutting-edge and compute-intensive epidemiological models to broad range of applications and users, providing precise and high-resolution spatiotemporal infectious disease dynamics. MEmilio is continuously extended and is available open-source for community use.

.. image:: http://martinkuehn.eu/research/images/MEmilio_small.png
   :alt: Memilio Overview

If you use MEmilio, please :doc:`cite our work<citation>`.


.. note::

   This framework is under active development, and, thus, obtains constantly new features or models. If you encounter a feature not yet documented, please :ref:`contact us directly <contact>`.



.. dropdown:: :fa:`list` **Table of Contents**
   :animate: fade-in-slide-down

   .. toctree::
      :maxdepth: 1
      :caption: About

      getting_started
      citation
      references
      faq
      development
      team

   .. toctree::
      :maxdepth: 2
      :caption: C++ Interface

      cpp/overview
      cpp/installation
      cpp/model_usage
      cpp/model_creation
      cpp/development
      cpp/interfaces

   .. toctree::
      :maxdepth: 2
      :caption: Python Interface

      python/python_packages
      python/m-simulation
      python/m-epidata
      python/m-surrogate
      python/m-generation
      python/m-plot

   .. toctree::
      :maxdepth: 1
      :caption: Code API 

      pythonapi/pythonapi
      api/library_root
