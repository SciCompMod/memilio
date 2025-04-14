# MEmilio C++ mobility

This directory contains a module to loosely couple multiple simulation instances and model the mobility between them. Each instance can e.g. represent a different geographical region. The regions are stored as nodes in a graph, with edges between them representing the mobility between the regions. Currently, only compartment models are supported as nodes, but it will be extended for any model.

At each time step, the simulation executes two following phases:
1. Advance the simulation for each node independently
2. Exchange people between nodes along the edges. The number of people exchanged depends on coefficients. The coefficient `a_i` of edge `e_xy` represents the percentage of people in compartment `i` moving from node `x` to node `y`. Like the contact matrices used in compartment models, the coefficients may contain dampings that change their value over time. During the next time step, the exchanged population will stay at their destination and participate in the evolution of that simulation. Afterwards, the people return. The total number of people returning is the same as the number that left. But the number of people in each compartment is adjusted according to the epidemiological situation in the destination node, e.g. some susceptible people that went from one node to another will have been exposed, so they return in a different compartment.

See the [mobility header](metapopulation_mobility_instant.h) and the `MobilityEdge` and `SimulationNode` classes for technical details of the two phases.

Utility classes:
- Graph: Abstract class (template) that stores the simulation instances (nodes) and the connections between them (edges).
- GraphSimulation: Abstract class (template) that executes custom functions on each node and edge in each time step.

