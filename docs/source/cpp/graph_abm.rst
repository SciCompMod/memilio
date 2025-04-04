Graph-based agent-based model
================================

Introduction
-------------

This model realizes multiple instances of the mobility-based agent-based model ``abm::Model`` (see :doc:`<mobility_based_abm>` for documentation) as nodes in a directed graph. One local model represents a geographical region. The regions are connected by the graph edges. Mobility within one node and via the graph edges follows the same mobility rules which can be handed as argument to ``GraphABModel``. Therefore this graph-based agent-based model can be reduced to a single mobility-based agent-based if simulation time steps within the whole graph, i.e. the step size of each node and the step size of the edge exchange, are equal.

Simulation
-----------

The simulation runs in discrete time steps with two different step sizes. The first one is the node-internal step size with which each node advances independently its inherent model instance. The second step size is used for the regular exchange of agents between models in different nodes via the corresponding edge. The simulation therefore advances all nodes until the next exchange time step and then performs the exchange via the graph edges.
