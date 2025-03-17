**MEmilio** - a high performance Modular EpideMIcs simuLatIOn software
==========================================================================

.. image:: https://github.com/SciCompMod/memilio/actions/workflows/main.yml/badge.svg?branch=main
   :target: https://github.com/SciCompMod/memilio/actions/workflows/main.yml
   :alt: CI Badge

.. image:: https://codecov.io/gh/SciCompMod/memilio/branch/main/graph/badge.svg?token=DVQXIQJHBM
    :target: https://codecov.io/gh/SciCompMod/memilio


Welcome
===============

MEmilio implements various models for infectious disease dynamics, ranging from simple compartmental models to complex Integro-Differential and agent-based models. Its modular design enables the combination of different models with distinct mobility patterns. Through efficient implementation and parallelization, MEmilio delivers cutting-edge and compute-intensive epidemiological models at a large scale, providing precise and high-resolution spatiotemporal infectious disease dynamics. MEmilio is continuously extended and is available open-source for community use.

.. image:: https://github.com/user-attachments/assets/65af6012-106e-43c6-9e0e-c96a73aa7b1e
   :alt: Memilio_overview

If you use MEmilio, please cite our work

- Kühn, Martin Joachim et al. (2024). *MEmilio - a High Performance Modular Epidemics Simulation Software (2022)*. Available at `GitHub Repository <https://github.com/SciCompMod/memilio>`_ and `DLR eLib <https://elib.dlr.de/209739/>`_.

and, in particular, for

- **Ordinary differential equation-based (ODE) and Graph-ODE models**: Zunker H, Schmieding R, Kerkmann D, Schengen A, Diexer S, et al. (2024). *Novel travel time aware metapopulation models and multi-layer waning immunity for late-phase epidemic and endemic scenarios*. *PLOS Computational Biology* 20(12): e1012630. `DOI:10.1371/journal.pcbi.1012630 <https://doi.org/10.1371/journal.pcbi.1012630>`_
- **Integro-differential equation-based (IDE) models**: Wendler AC, Plötzke L, Tritzschak H, Kühn MJ. (2024). *A nonstandard numerical scheme for a novel SECIR integro differential equation-based model with nonexponentially distributed stay times*. Submitted for publication. `arXiv:2408.12228 <https://arxiv.org/abs/2408.12228>`_
- **Agent-based models (ABMs)**: Kerkmann D, Korf S, Nguyen K, Abele D, Schengen A, et al. (2024). *Agent-based modeling for realistic reproduction of human mobility and contact behavior to evaluate test and isolation strategies in epidemic infectious disease spread*. arXiv. `arXiv:2410.08050 <https://arxiv.org/abs/2410.08050>`_
- **Hybrid agent-metapopulation-based models**: Bicker J, Schmieding R, Meyer-Hermann M, Kühn MJ. (2025). *Hybrid metapopulation agent-based epidemiological models for efficient insight on the individual scale: A contribution to green computing*. *Infectious Disease Modelling* 10(2): 571-590. `DOI:10.1016/j.idm.2024.12.015 <https://doi.org/10.1016/j.idm.2024.12.015>`_
- **Graph Neural Networks**: Schmidt A, Zunker H, Heinlein A, Kühn MJ. (2024).*Towards Graph Neural Network Surrogates Leveraging Mechanistic Expert Knowledge for Pandemic Response*. arXiv. `arXiv:2411.06500 <https://arxiv.org/abs/2411.06500>`_
- **ODE-based models with Linear Chain Trick**: Plötzke L, Wendler A, Schmieding R, Kühn MJ. (2024). *Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions*. Submitted for publication. `DOI:10.48550/arXiv.2412.09140 <https://doi.org/10.48550/arXiv.2412.09140>`_
- **Behavior-based ODE models**: Zunker H, Dönges P, Lenz P, Contreras S, Kühn MJ. (2025). *Risk-mediated dynamic regulation of effective contacts de-synchronizes outbreaks in metapopulation epidemic models*. arXiv. `arXiv:2502.14428 <https://arxiv.org/abs/2502.14428>`_


Getting Started
===============

MEmilio supports various model types including equation-based, agent-based, and hybrid graph-ODE-based models. Among the equation-based models, we provide ordinary differential equation (ODE), linear chain trick (LCT), integro-differential equation (IDE) based models.

- The C++ backend powers the model and simulation components for optimal efficiency.
- Data acquisition, plotting, and machine-learned models are handled via Python.

For more details on the C++ implementation, see the `cpp directory README <cpp/README.md>`_.

Some regularly used data for simulations of a pathogen's spread in Germany, like contact and inter-county mobility, can be found in the `data directory README <data/README.md>`_.

In `pycode`, different MEmilio Python packages are defined. Via our `memilio-simulation` package, you can run our C++ backend from Python; this package uses `pybind11` to bind our C++ model code. The `memilio-epidata` package provides tools to download and structure important data such as infection or mobility data. More about the Python packages can be found in the `Python README <pycode/README.md>`_.



Development
===========

- `Git workflow and change process <https://github.com/SciCompMod/memilio/wiki/git-workflow>`_
- `Coding Guidelines <https://github.com/SciCompMod/memilio/wiki/coding-guidelines>`_

**MEmilio** is a framework that ..

Check out the :doc:`getting_started` section for further information, including
how to :hoverxref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   getting_started
   references
   faq

.. toctree::
   :maxdepth: 2
   :caption: C++ Interface

   c++/model_usage
   c++/model_creation

.. toctree::
   :maxdepth: 2
   :caption: Python Interface

   python/model_usage
   python/model_creation

.. toctree::
   :maxdepth: 1
   :caption: Code API 

   models/index
   pythonapi
   api/library_root