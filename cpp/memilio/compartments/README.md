This directory contains the framework used to build SIR-type compartment models.

Classes:
- CompartmentModel: Template base class for compartment models. Specialize the class template using a parameter set (e.g. using the [ParameterSet class](../utils/parameter_set.h)) and populations (e.g. using the [Populations class](../epidemiology/populations.h)). The population is divided into compartments (and optionally other subcategories, e.g. age groups). Derive from the class to define the flows between the compartments.
- Simulation: Template class that runs the simulation using a specified compartment model. Can be derived from to implement behavior that cannot be modeled inside the usual compartment flow structure. There is an additional class called SimulationFlows, which inherits from the Simulation class itself and is also compatible with it. Its functioning closely resembles that of the Simulation class. In contrast to the Simulation class, within SimulationFlows, we focus on the flows and incorporate them into the integration process. As a result, we have the integrated flows as additional output.

It is also important to note that in our current use of the simulation, the last value of the TimeSeries in the simulation result is always used as the initial value. Thus, when the advance function is called, the model is no longer used. Consequently, any changes made in the model object, such as population changes, will not affect the results of the simulation.


See the implemented [SEIR model](../../models/seir/README.md) for a simple example of using the classes. See the [SECIR model](../../models/secir/README.md) for an advanced example with age resolution. 
