Structure
=========

All aggregated models follow the same overall structure. To create a new model, please proceed as follows.

In ``cpp/models``, create a new folder for your model. This folder will contain the following model specific files that are described below.

- ``CMakeLists.txt``: Defines your model specific library, links libraries, specifies include directories and adds compile options for your library. Note that you have to add your model specific library to the file ``cpp/CMakeLists.txt``.
- ``README.md``: Provides an overview of your model. 
- ``infection_state.h``: Defines an ``enum`` that specifies the infection states of your model.
- ``model.h``/ ``model.cpp``: Depending on the model, the model equations are specified or a model specific solver is implemented.
- ``parameters.h``: Defines all parameters used in the model. 

Depending on your model, it may make sense to create additional files for functionality that is specific to your model. The following list gives examples from already existing models but is not exhaustive.

- ``parameters_io.h``/ ``parameters_io.cpp``: Implements a scheme to initialize the model based on reported data. 
- ``simulation.h``/ ``simulation.cpp``: If your model does not derive from ``CompartmentalModel`` or ``FlowModel``, you need to provide a ``Simulation`` class that determines how your model equations are solved. 

For further information on how to create a new model, see the respective sections for the different types of equation-based models.