This directory contains the framework used to build SIR-type compartment models.

Classes:
- CompartmentModel: Template base class for compartment models. Specialize the class template using a parameter set (e.g. using the [ParameterSet class](../utils/parameter_set.h)) and populations (e.g. using the [Populations class](../epidemiology/populations.h)). The population is divided into compartments (and optionally other subcategories, e.g. age groups). Derive from the class to define the flows between the compartments.
- Simulation: Template class that runs the simulation using a specified compartment model.  Can be derived from to implement behavior that cannot be modeled inside the usual compartment flow structure.

See the implemented [SEIR model](../../models/seir/README.md) for a simple example of using the classes. See the [SECIR model](../../models/secir/README.md) for an advanced example with age resolution. 