# MEmilio - a high performance Modular EpideMIcs simuLatIOn software #

![memilio_logo](docs/memilio-small.png)

[![CI](https://github.com/SciCompMod/memilio/actions/workflows/main.yml/badge.svg)](https://github.com/SciCompMod/memilio/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/SciCompMod/memilio/branch/main/graph/badge.svg?token=DVQXIQJHBM)](https://codecov.io/gh/SciCompMod/memilio)

MEmilio implements various models for infectious disease dynamics, from simple compartmental models through Integro-Differential equation-based models to agent- or individual-based models. Its modular design allows the combination of different models with different mobility patterns. Through efficient implementation and parallelization, MEmilio brings cutting edge and compute intensive epidemiological models to a large scale, enabling a precise and high-resolution spatiotemporal infectious disease dynamics. MEmilio will be extended continuously. It is available open-source and we encourage everyone to make use of it.

If you use MEmilio, please cite our work

- Bicker J, Kerkmann D, Korf S, Plötzke L, Schmieding R, Wendler A, Zunker H et al. (2025)  *MEmilio - a High Performance Modular Epidemics Simulation Software*. Available at `https://github.com/SciCompMod/memilio` and `https://elib.dlr.de/213614/`.

and, in particular, for

- Ordinary differential equation-based (ODE) and Graph-ODE models: Zunker H, Schmieding R, Kerkmann D, Schengen A, Diexer S, et al. (2024). *Novel travel time aware metapopulation models and multi-layer waning immunity for late-phase epidemic and endemic scenarios*. *PLOS Computational Biology* 20(12): e1012630. `https://doi.org/10.1371/journal.pcbi.1012630`
- Integro-differential equation-based (IDE) models: Wendler AC, Plötzke L, Tritzschak H, Kühn MJ. (2024). *A nonstandard numerical scheme for a novel SECIR integro differential equation-based model with nonexponentially distributed stay times*. Submitted for publication. `https://arxiv.org/abs/2408.12228`
- Agent-based models (ABMs): Kerkmann D, Korf S, Nguyen K, Abele D, Schengen A, et al. (2025). *Agent-based modeling for realistic reproduction of human mobility and contact behavior to evaluate test and isolation strategies in epidemic infectious disease spread*. *Computers in Biology and Medicine* 193: 110269. `DOI:10.1016/j.compbiomed.2025.110269 <https://doi.org/10.1016/j.compbiomed.2025.110269>`_
- Hybrid agent-metapopulation-based models: Bicker J, Schmieding R, Meyer-Hermann M, Kühn MJ. (2025). *Hybrid metapopulation agent-based epidemiological models for efficient insight on the individual scale: A contribution to green computing*. *Infectious Disease Modelling* 10(2): 571-590. `https://doi.org/10.1016/j.idm.2024.12.015`
- Graph Neural Networks: Schmidt A, Zunker H, Heinlein A, Kühn MJ. (2025). *Graph Neural Network Surrogates to leverage Mechanistic Expert Knowledge towards Reliable and Immediate Pandemic Response*. Submitted for publication. `arXiv:2411.06500 <https://arxiv.org/abs/2411.06500>`
- ODE-based models with Linear Chain Trick: Plötzke L, Wendler A, Schmieding R, Kühn MJ. (2024). *Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions*. Submitted for publication. `https://doi.org/10.48550/arXiv.2412.09140`
- Behavior-based ODE models: Zunker H, Dönges P, Lenz P, Contreras S, Kühn MJ. (2025). *Risk-mediated dynamic regulation of effective contacts de-synchronizes outbreaks in metapopulation epidemic models*. Chaos, Solitons & Fractals. `https://doi.org/10.1016/j.chaos.2025.116782`


**Getting started**


The documentation for MEmilio can be found at 
https://memilio.readthedocs.io/en/latest/index.html

**Publication simulations**

Simulations used for publications, along with their specific plotting and processing scripts, 
are available in a separate repository:
https://github.com/SciCompMod/memilio-simulations

**Development**

The coding guidelines and git workflow description can be found in the documentation at 
https://memilio.readthedocs.io/en/latest/development.html
