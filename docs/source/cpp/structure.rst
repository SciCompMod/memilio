Structure
=========

In ``cpp/models``, create a new folder for your model. This folder will contain the following model specific files that are described below.

- ``CMakeLists.txt``: Defines your model specific library, links libraries, specifies include directories and adds compile options for your library. Note that you have to add your model specific library to the file ``cpp/CMakeLists.txt``.
- ``README.md``: Provides an overview of your model. 
- ``infection_state.h``: Defines an ``enum`` that specifies the infection states of your model.
- ``model.h``/ ``model.cpp``: Depending on the model, these files define the overall functionality of your model and, in case, of an aggregated model, the equations defining the flows between different infection states. 
- ``parameters.h``: Defines all parameters used in the model. 

Depending on your model, it may make sense to create additional files for functionality that is specific to your model. For aggregated models, the following list gives examples from already existing models but is not exhaustive.

- ``parameters_io.h``/ ``parameters_io.cpp``: Implements a scheme to initialize the model based on reported data. 
- ``simulation.h``/ ``simulation.cpp``: If your model does not derive from ``CompartmentalModel`` or ``FlowModel``, you need to provide a ``Simulation`` class that determines how your model equations are solved. 

For further information on how to create a new aggregated model, see the respective sections :doc:`above this structure <../model_creation>`. The creation of an individual-based model is much more complex and we suggest to build upon ``cpp/models/abm``. 