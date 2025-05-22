LCT Models
==========

Currently, a model using the Linear Chain Trick is implemented for a SECIR-type compartmental structure. Here, we want to explain how LCT models can be set up for other compartmental structures. We will do this using the example of the already implemented LCT-SECIR model. Note that this model contains different groups that can e.g. represent age groups.  

Compartments 
-------------

To set up a new LCT model, we start with defining the infection states in an enum class ``InfectionState`` in the file ``infection_state.h``. For our SECIR model this is done by

.. code-block:: cpp

    enum class InfectionState
    {
        Susceptible        = 0,
        Exposed            = 1,
        InfectedNoSymptoms = 2,
        InfectedSymptoms   = 3,
        InfectedSevere     = 4,
        InfectedCritical   = 5,
        Recovered          = 6,
        Dead               = 7,
        Count              = 8
    };

Parameters
------------

We continue by defining the necessary parameters in the file ``parameters.h``. In the case of LCT models, we require the average time spent in each compartment for the transient compartments, so in the case of the LCT-SECIR model for the compartments `Exposed`, `InfectedNoSymptoms`, `InfectedSymptoms`, `InfectedSevere` and `InfectedCritical`. These parameters can be implemented as follows, see e.g. for the average time spent in the `Exposed` compartment:

.. code-block:: cpp
    struct TimeExposed {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size)
    {
        return Type::Constant(size, 1, 1.);
    }
    static std::string name()
    {
        return "TimeExposed";
    }
    };

Here, we define a struct that contains a default value and a name for the parameter. Since the considered model can contain multiple groups, the default is a vector of dimension `size x 1` where each entry has value :math:`1` where `size` denotes the number of groups. For the remaining parameters describing the average time spent in a compartment, we can proceed analogously.  

Furthermore, we require the transition probabilities from one compartment to another where we only consider the ones that are not one by default as given by the compartmental structure. In the case of the SECIR model, we define the probabilities ``RecoveredPerInfectedNoSymptoms``, ``SeverePerInfectedSymptoms``, ``CriticalPerSevere`` and ``DeathsPerCritical``. These probabilities can be implemented as follows

.. code-block:: cpp

    struct RecoveredPerInfectedNoSymptoms {
    using Type = Eigen::VectorX<UncertainValue<ScalarType>>;
    static Type get_default(size_t size)
    {
        return Type::Constant(size, 1, 0.5);
    }
    static std::string name()
    {
        return "RecoveredPerInfectedNoSymptoms";
    }
    };

where we define a default value that is again vector of dimension `size x 1` where each entry has value `0.5` and a name for the parameter. The other transition probabilities can be defined analogously.

We also define a parameter ``ContactPatterns`` determining the contacts of the different groups by 

.. code-block:: cpp

    struct ContactPatterns {
    using Type = UncertainContactMatrix<ScalarType>;

    static Type get_default(size_t size)
    {
        mio::ContactMatrixGroup contact_matrix(1, (Eigen::Index)size);
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant((Eigen::Index)size, (Eigen::Index)size, 10.));
        return Type(contact_matrix);
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
    };

with a default contact matrix of dimension `size x size` where each entry has value :math:`10` and a name for the parameter. 

Additionally, we can determine parameters influencing the infection dynamics. In the case of the LCT-SECIR model we use the parameters ``TransmissionProbabilityOnContact``, ``RelativeTransmissionNoSymptoms``, ``RiskOfInfectionFromSymptomatic``, ``StartDay`` and ``Seasonality``. For each parameter, we need to define a default value and a name as for the above parameters. 

After having defined all parameters that are required for the model, we can define a ``ParameterSet`` containing all parameters by 

.. code-block:: cpp

    using ParametersBase =
    ParameterSet<TimeExposed, TimeInfectedNoSymptoms, TimeInfectedSymptoms, TimeInfectedSevere, TimeInfectedCritical,
                 TransmissionProbabilityOnContact, ContactPatterns, RelativeTransmissionNoSymptoms,
                 RiskOfInfectionFromSymptomatic, RecoveredPerInfectedNoSymptoms, SeverePerInfectedSymptoms,
                 CriticalPerSevere, DeathsPerCritical, StartDay, Seasonality>;

Furthermore, we define a class ``Parameters`` that inherits from this ``ParameterSet`` by 

.. code-block:: cpp

    class Parameters : public ParametersBase

This class should contain a method ``check_constraints()`` that checks if the values of the parameters are valid and a method ``deserialize()``. We will use an object of this ``Parameters`` class in the ``Model`` class (see below) so that we can use the here defined parameters within the model equations. Please check the already implemented examples for further details on the implementation.

Model equations
-----------------

Now that we have defined the compartments and parameters that we want to consider, we can define the model equations that we want to solve. 

For this, we define the class ``Model`` in the file ``model.h`` by 

.. code-block:: cpp
    template <class... LctStates>
    class Model : public CompartmentalModel<ScalarType, InfectionState, LctPopulations<ScalarType, LctStates...>, Parameters>

Note that this class has a template parameter ``LctStates`` that defines the number of subcompartments per compartment for every considered group. This class also inherits from ``CompartmentalModel``. For LCT models, the class ``CompartmentalModel`` requires the following template arguments:
    
- type of floating point type, here ``ScalarType``,
- a class ``InfectionState`` containing the compartments, see above,
- the class ``LctPopulations`` which is a class template for compartment populations of LCT models depending on the floating point type and the considered ``LctStates``
- the class ``Parameters`` containing all required parameters, see above. 

The following methods are implemented within the ``Model`` class:

- Constructor of the model
- The function ``get_derivatives()`` evaluates the right-hand-side of the ODE :math:`dydt = f(y, t)` that we want to solve.
- The function ``calculate_compartments()`` accumulates the TimeSeries containing simulation results that are divided into subcompartments to a TimeSeries that conatins the simulation results according to the infection states without subcompartments. 
- The function ``check_constraints()`` check that the model satisfies all constraints regarding parameters and populations. 

Note that you have to create a ``CMakeLists.txt`` file within your model folder where you need to create a library for your new model, link libraries, specify include directories and add compile options for your library. This new library also needs to be added to the global CMakeLists in the cpp folder. 