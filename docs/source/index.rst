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

- Kühn, Martin Joachim et al. (2024). *MEmilio - a High Performance Modular Epidemics Simulation Software (2022)*. Available at `https://github.com/SciCompMod/memilio` and `https://elib.dlr.de/209739/`.

and, in particular, for

- **Ordinary differential equation-based (ODE) and Graph-ODE models**: Zunker H, Schmieding R, Kerkmann D, Schengen A, Diexer S, et al. (2024). *Novel travel time aware metapopulation models and multi-layer waning immunity for late-phase epidemic and endemic scenarios*. *PLOS Computational Biology* 20(12): e1012630. `https://doi.org/10.1371/journal.pcbi.1012630`
- **Integro-differential equation-based (IDE) models**: Wendler AC, Plötzke L, Tritzschak H, Kühn MJ. (2024). *A nonstandard numerical scheme for a novel SECIR integro differential equation-based model with nonexponentially distributed stay times*. Submitted for publication. `https://arxiv.org/abs/2408.12228`
- **Agent-based models (ABMs)**: Kerkmann D, Korf S, Nguyen K, Abele D, Schengen A, et al. (2024). *Agent-based modeling for realistic reproduction of human mobility and contact behavior to evaluate test and isolation strategies in epidemic infectious disease spread*. arXiv. `https://arxiv.org/abs/2410.08050`
- **Hybrid agent-metapopulation-based models**: Bicker J, Schmieding R, Meyer-Hermann M, Kühn MJ. (2025). *Hybrid metapopulation agent-based epidemiological models for efficient insight on the individual scale: A contribution to green computing*. *Infectious Disease Modelling* 10(2): 571-590. `https://doi.org/10.1016/j.idm.2024.12.015`
- **Graph Neural Networks**: Schmidt A, Zunker H, Heinlein A, Kühn MJ. (2024). *Towards Graph Neural Network Surrogates Leveraging Mechanistic Expert Knowledge for Pandemic Response*. arXiv. `https://arxiv.org/abs/2411.06500`
- **ODE-based models with Linear Chain Trick**: Plötzke L, Wendler A, Schmieding R, Kühn MJ. (2024). *Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions*. Submitted for publication. `https://doi.org/10.48550/arXiv.2412.09140`
- **Behavior-based ODE models**: Zunker H, Dönges P, Lenz P, Contreras S, Kühn MJ. (2025). *Risk-mediated dynamic regulation of effective contacts de-synchronizes outbreaks in metapopulation epidemic models*. arXiv. `https://arxiv.org/abs/2502.14428`


Getting Started
===============

MEmilio supports various model types including equation-based, agent-based, and hybrid graph-ODE-based models. Among the equation-based models, we provide ordinary differential equation (ODE), linear chain trick (LCT), integro-differential equation (IDE) based models.

- The C++ backend powers the model and simulation components for optimal efficiency.
- Data acquisition, plotting, and machine-learned models are handled via Python.

For more details on the C++ implementation, see the `cpp directory README <cpp/README.md>`_.

Some regularly used data for simulations of a pathogen's spread in Germany, like contact and inter-county mobility, can be found in the `data directory README <data/README.md>`_.

In `pycode`, different MEmilio Python packages are defined. Via our `memilio-simulation` package, you can run our C++ backend from Python; this package uses `pybind11` to bind our C++ model code. The `memilio-epidata` package provides tools to download and structure important data such as infection or mobility data. More about the Python packages can be found in the `Python README <pycode/README.md>`_.

References
===========

Recently Submitted Publications
--------------------------------------

- Zunker H, Dönges P, Lenz P, Contreras S, Kühn MJ. (2025). *Risk-mediated dynamic regulation of effective contacts de-synchronizes outbreaks in metapopulation epidemic models*. arXiv. `https://arxiv.org/abs/2502.14428`
- Schmidt A, Zunker H, Heinlein A, Kühn MJ. (2024). *Towards Graph Neural Network Surrogates Leveraging Mechanistic Expert Knowledge for Pandemic Response*. arXiv. `https://arxiv.org/abs/2411.06500`
- Wendler AC, Plötzke L, Tritzschak H, Kühn MJ. (2024). *A nonstandard numerical scheme for a novel SECIR integro differential equation-based model with nonexponentially distributed stay times*. Submitted for publication. `https://arxiv.org/abs/2408.12228`
- Kerkmann D, Korf S, Nguyen K, Abele D, Schengen A, et al. (2024). *Agent-based modeling for realistic reproduction of human mobility and contact behavior to evaluate test and isolation strategies in epidemic infectious disease spread*. arXiv. `https://arxiv.org/abs/2410.08050`
- Plötzke L, Wendler A, Schmieding R, Kühn MJ. (2024). *Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions*. Submitted for publication. `https://doi.org/10.48550/arXiv.2412.09140`


Peer-Reviewed Publications
--------------------------

**2025**

- Bicker J, Schmieding R, Meyer-Hermann M, Kühn MJ. (2025). *Hybrid metapopulation agent-based epidemiological models for efficient insight on the individual scale: A contribution to green computing*. *Infectious Disease Modelling* 10(2): 571-590. `https://doi.org/10.1016/j.idm.2024.12.015`

**2024**

- Zunker H, Schmieding R, Kerkmann D, Schengen A, Diexer S, et al. (2024). *Novel travel time aware metapopulation models and multi-layer waning immunity for late-phase epidemic and endemic scenarios*. *PLOS Computational Biology* 20(12): e1012630. `https://doi.org/10.1371/journal.pcbi.1012630`

**2022**

- Kühn MJ, Abele D, Binder S, Rack K, Klitz M, et al. (2022). *Regional opening strategies with commuter testing and containment of new SARS-CoV-2 variants in Germany*. *BMC Infectious Diseases* 22(1): 333. `https://doi.org/10.1186/s12879-022-07302-9`
- Koslow W, Kühn MJ, Binder S, Klitz M, Abele D, et al. (2022). *Appropriate relaxation of non-pharmaceutical interventions minimizes the risk of a resurgence in SARS-CoV-2 infections in spite of the Delta variant*. *PLOS Computational Biology* 18(5): e1010054. `https://doi.org/10.1371/journal.pcbi.1010054`

**2021**

- Kühn MJ, Abele D, Mitra T, Koslow W, Abedi M, et al. (2021). *Assessment of effective mitigation and prediction of the spread of SARS-CoV-2 in Germany using demographic information and spatial resolution*. *Mathematical Biosciences* 108648. `https://doi.org/10.1016/j.mbs.2021.108648`


Development
===========

- `Git workflow and change process <https://github.com/SciCompMod/memilio/wiki/git-workflow>`_
- `Coding Guidelines <https://github.com/SciCompMod/memilio/wiki/coding-guidelines>`_

**MEmilio** is a framework that ..

Check out the :doc:`usage` section for further information, including
how to :hoverxref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   getting_started
   faq
   c++/c++
   python/python
   models/models
   pythonapi
   api/library_root