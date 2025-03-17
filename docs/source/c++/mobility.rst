Mobility
========

This directory contains a module to loosely couple multiple simulation instances and model the mobility between them. Each instance can represent a different geographical region. The regions are stored as nodes in a graph, with edges between them representing the mobility between the regions. Currently, only compartment models are supported as nodes, but support for other models will be extended in the future.

Simulation Process
------------------

At each time step, the simulation executes the following two phases:

1. **Advance the Simulation:**
   - Each node (region) progresses independently according to its model.

2. **Exchange of Population:**
   - People are exchanged between nodes along the edges.
   - The number of exchanged individuals is determined by coefficients.
   - The coefficient :math:`a_i` of edge :math:`e_{xy}` represents the percentage of people in compartment :math:`i` moving from node :math:`x` to node :math:`y`.
   - Coefficients may include dampings that vary over time, similar to contact matrices in compartment models.
   - During the next time step, the exchanged population will stay at their destination and influence the simulation's evolution. Afterward, individuals return to their original node.
   - The total number of people returning matches the number that left, but the compartments are adjusted to reflect changes due to the destination's epidemiological situation. For instance, some susceptible individuals may become exposed and will return in a different compartment.

Graph-Based Mobility Model
--------------------------

To utilize the novel model, assume :math:`n` different geographic units (denoted as *regions*) are given. We integrate our novel mobility model extension into a simple ODE model with the graph approach proposed in Kühn et al. (2021).

- Each region is represented by a node in the graph.
- The (multi-)edge :math:`\mathcal{E}_{ij}` between nodes :math:`\mathcal{N}_i` and :math:`\mathcal{N}_j` represents the (outgoing) mobility, with mobility-based exchange occurring twice a day (round trip).
- Each combination of sociodemographic group :math:`g_l`, where :math:`l = 1, \ldots, G`, and infection state :math:`z_l`, where :math:`l = 1, \ldots, Z`, is assigned a number of travelers. Therefore, the multi-edge consists of :math:`G \times Z` single edges.
- Return trips are mapped on the same edge, with the vector of weights on :math:`\mathcal{E}_{ij}` representing the number of outgoing travelers.

For more information, we refer to:

- Zunker H, Schmieding R, Kerkmann D, Schengen A, Diexer S, et al. (2024). *Novel travel time aware metapopulation models and multi-layer waning immunity for late-phase epidemic and endemic scenarios*. *PLOS Computational Biology* 20(12): e1012630. `https://doi.org/10.1371/journal.pcbi.1012630`
- Kühn MJ, Abele D, Mitra T, Koslow W, Abedi M, et al. (2021). *Assessment of effective mitigation and prediction of the spread of SARS-CoV-2 in Germany using demographic information and spatial resolution*. *Mathematical Biosciences* 108648. `https://doi.org/10.1016/j.mbs.2021.108648`
