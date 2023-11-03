This directory contains the framework used to build SIR-type compartment models.

Classes:
- CompartmentalModel: Template base class for compartmental models. Specialize the class template using a parameter set (e.g. using the [ParameterSet class](../utils/parameter_set.h)) and populations (e.g. using the [Populations class](../epidemiology/populations.h)). The population is divided into compartments (and optionally other subcategories, e.g. age groups). Derive from the class to define the change in populations.
- FlowModel: A CompartmentalModel defined by the flows between populations. Requires the additional template parameter Flows, which is a [TypeList](../utils/type_list.h) of [Flow](../utils/flow.h). Instead of defining `get_derivatives`, the function `get_flows` has to be implemented for a FlowModel. 
- Simulation: Template class that runs the simulation using a specified compartment model. Can be derived from to implement behavior that cannot be modeled inside the usual compartment flow structure.
- FlowSimulation: Template class derived from Simulation, that additionally computes the flows between populations. It requires a FlowModel.

Note that currently, only the last value of the simulation result (i.e. the TimeSeries obtained from `get_result()`) is used to advance the simulation. Thus, any changes to the population of the simulation's model will have no effect on the simulation result.


See the implemented [SEIR model](../../models/seir/README.md) for a simple example of using the classes. See the [SECIR model](../../models/secir/README.md) for an advanced example with age resolution. 
