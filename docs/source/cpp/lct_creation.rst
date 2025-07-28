LCT model creation
==================

The mathematical model
----------------------

Before implementing a model in MEmilio, we need to do some math, in particular, define an initial value problem
given by a system of ordinary differential equations. Here, we consider a SIRD model where we 
divide the Infectious compartment into :math:`n` subcompartments. The model is defined by

.. math::  

    \begin{aligned}
        S'(t) & = -\rho\phi\ \frac{S(t)I(t)}{N_{\perp D}} \\
        I_1'(t) & = \rho\phi\ \frac{S(t)I(t)}{N_{\perp D}} - \frac{n}{T_I}I_1(t) \\
        I_j'(t) & = \frac{n}{T_I}I_{j-1}(t) - \frac{n}{T_I}I_j(t) \quad \text{for } j\in\{2,\dots,n\}\\
        R'(t) & = \frac{\mu_R}{T_I}I(t) \\
        D'(t) & = \frac{\mu_D}{T_I}I(t) \\
    \end{aligned}

with :math:`I(t) = \sum_{j=1}^n I_j(t)` and some initial values for :math:`t=0`. Here :math:`N_{\perp D} := S(t) + I(t) + R(t)`.

This type of model belongs to the class of compartmental models because the model population is represented by discrete infection
states **S**\usceptible, **I**\nfectious, **R**\ecovered, **D**\eceased, also called compartments.

Infection states
~~~~~~~~~~~~~~~~

First we define an :code:`enum class` called InfectionState in the file "infection_state.h", which contains an entry
for each infection state of the mathematical model, followed by an entry called :code:`Count`. This enumerates the 
compartments starting from 0, with Count being equal to the number of compartments. For example:

.. code-block:: cpp

    enum class InfectionState
    {
        Susceptible,
        Infectious,
        Recovered,
        Deceased,
        Count
    };

Parameters
~~~~~~~~~~

Next, we define the parameters in "parameters.h", which consist of a struct for each parameter used in the mathematical
model. This struct must define the data type, name and default value of the constant. For example, for the time a
person stays infectious, :math:`T_I`, we define a struct:

.. code-block:: cpp

    template <typename FP>
    struct TimeInfectious {
        using Type = UncertainValue<FP>;

        static Type get_default()
        {
            return Type(6.73);
        }

        static std::string name()
        {
            return "TimeInfectious";
        }
    };


We also define a parameter ``ContactPatterns`` determining the contacts of the different groups by:

.. code-block:: cpp

    struct ContactPatterns {
        using Type = UncertainContactMatrix<ScalarType>;

        static Type get_default()
        {
            mio::ContactMatrixGroup contact_matrix(1, 1);
            contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
            return Type(contact_matrix);
        }
        static std::string name()
        {
            return "ContactPatterns";
        }
    }; 

Avoid using the mathematical symbol of a constant as name for the struct. Their connection can be noted down in the
documentation of these structs.

Finally, define a type :code:`Parameters` by listing all parameter structs as template arguments of a
:code:`mio::ParameterSet`:

.. code-block:: cpp

    template <typename FP>
    using Parameters = mio::ParameterSet<TimeInfectious<FP>, RecoveryRate<FP>, LethalityRate<FP>, ContactRate<FP>,
                                         TransmissionRisk<FP>>;

For more complex models, :code:`Parameters` allows passing arguments from its constructor to the :code:`get_default`
functions. Make sure that all of these functions take the exact types as function arguments that you want to pass to
the constructor.

Population
~~~~~~~~~~

The population will be stored in a vector, with a component for each subcompartment of every infection state. We define 
it using the class ``LctPopulations``.

.. code-block:: cpp

    template <typename FP = ScalarType, class... LctStates>
    using Populations = mio::LctPopulations<FP, LctStates...>;

where ``LctStates`` contains the number of subcompartments per infection state.

Importantly, this class allows further stratifying the population vector, with the most common
example being adding age groups.

Define the model
^^^^^^^^^^^^^^^^

Now we can define the model as a **CompartmentalModel** in the file "model.h":  

.. code-block:: cpp

    template <class... LctStates>
    class Model
        : public CompartmentalModel<ScalarType, InfectionState, LctPopulations<ScalarType, LctStates...>, Parameters>
    {
    public:
        using LctStatesGroups = TypeList<LctStates...>;
        using Base = CompartmentalModel<ScalarType, InfectionState, LctPopulations<ScalarType, LctStates...>, Parameters>;
        using typename Base::ParameterSet;
        using typename Base::Populations;

        void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                             Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
        {
            const Parameters<FP>& params = this->parameters;

            const auto N = y[InfectionState::Susceptible] + y[InfectionState::Infectious] +
                           y[InfectionState::Recovered];

            dydt[InfectionState::Susceptible] = -params.template get<TransmissionRisk<FP>>() *
                                                 params.template get<ContactRate<FP>>() *
                                                 y[InfectionState::Susceptible] * y[InfectionState::Infectious] / N;
            
            . . .
        }
    };

Note that this class has a template parameter ``LctStates`` that defines the number of subcompartments per infection state. 
For LCT models, the class ``CompartmentalModel`` requires the following template arguments:
    
- type of floating point type, here ``ScalarType``,
- a class ``InfectionState`` containing the compartments, see above,
- the class ``LctPopulations`` which is a class template for compartment populations of LCT models depending on the 
floating point type and the considered ``LctStates`` and determines the type of the public member ``populations`` which contains 
the number of individuals per subcompartment and is used to pass initial conditions to the model,
- the class ``Parameters`` containing all required parameters, see above. 

The function ``get_derivatives()`` evaluates the right-hand-side of the ODE :math:`dydt = f(y, t)` that we want to solve, see above.

It is also useful to implement the following methods within the model:

- A function ``calculate_compartments()`` that accumulates the TimeSeries containing simulation results that are divided 
into subcompartments to a TimeSeries that conatins the simulation results according to the infection states without subcompartments. 
For an example, see the implementation within the LCT-SECIR model.
- A function ``check_constraints()`` that checks that the model satisfies sensible constraints regarding parameters and initial conditions. 
For an example, see the implementation within the LCT-SECIR model. 
