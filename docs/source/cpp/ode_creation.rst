ODE model creation
==================

The mathematical model
----------------------

Before implementing a model in MEmilio, we need to do a some math, in particular, define an initial value problem
given by a system of ordinary differential equations. For example we consider a SIRD model given by

.. math::  

    \begin{aligned}
        S'(t) & = -\rho\phi\ \frac{S(t)*I(t)}{N_{\perp D}} \\
        I'(t) & = \rho\phi\ \frac{S(t)*I(t)}{N_{\perp D}} - \frac{\mu_R + \mu_D}{T_I}I(t) \\
        R'(t) & = \frac{\mu_R}{T_I}I(t) \\
        D'(t) & = \frac{\mu_D}{T_I}I(t) \\
    \end{aligned}

and some initial values for :math:`t=0`. Here :math:`N_{\perp D} := S(t) + I(t) + R(t)`.

This type of model is called compartmental model, because the model population is represented by discrete infection
states **S** usceptible, **I** nfectious, **R** ecovered, **D** eceased, also called compartments.

How to define an ODE model
--------------------------

To define an ODE model in MEmilio, there are two options. You can define a CompartmentalModel or FlowModel, which
use different methods to define the right hand side of the mathematical model above. Both classes need definitions for
the infection states, population and parameters it uses. The FlowModel additionally requires a list of flows.

We start by creating a new directory for our model under "cpp/models", in this case we can call it "ode_sird". The name
must be unique and start with "ode\_", so the type of model is obvious. The rest usually contains the compartments or
other noteworthy features of the model in shortened form. All files in the following are put into this directory.

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

Next, we define the parameters in "parameters.h", which consist of a struct for each constant used in the mathematical
model. This struct must define the data type, name and default value of the constant. For example, for the time a
person stays infectious :math:`T_I` we define a struct

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

and for the contact rate :math:`\phi` a struct

.. code-block:: cpp

    template <typename FP>
    struct ContactRate {
        using Type = FP;

        static Type get_default()
        {
            return Type(10.0);
        }

        static std::string name()
        {
            return "ContactRate";
        }
    };

Avoid using the mathematical symbols of the constant as names for the struct. Their connection can be noted in the
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

The population will be stored in a vector, with a component for each infection state. We define it using the class
`mio::Population`.

.. code-block:: cpp

    template <typename FP>
    using Population = mio::Populations<FP, InfectionState>;

Importantly, this class allows further stratifying the population vector, with the most common
example being adding :code:`mio::AgeGroups`.

Define a compartmental model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we can define the model:

.. code-block:: cpp

    template <typename FP = ScalarType>
    class Model : public mio::CompartmentalModel<FP, InfectionState, Population<FP>, Parameters<FP>>
    {
    public:
        using Base = mio::CompartmentalModel<FP, InfectionState, Population<FP>, Parameters<FP>>;
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

Define a flow model
^^^^^^^^^^^^^^^^^^^

A flow model is a special case of a compartmental model, where each compartment :math:`Z_i` can be written as

.. math::

    Z_i(t) = \sum_{i \ne j} f_{Z_j \rightarrow Z_i}(t) - \sum_{i \ne j} f_{Z_i \rightarrow Z_j}(t),

where the flows :math:`f_{Z_i \rightarrow Z_j} \gt 0` are the amount of population changing from compartment
:math:`Z_i` to :math:`Z_j` at time :math:`t`. So the first sum accumulates all inflows, the second subtracts all
outflows.

The SIRD model from above can be expressed as a flow model with only three flows:

.. math::  

    \begin{aligned}
        f_{S \rightarrow I} & = \rho\phi\ \frac{S(t)*I(t)}{N_{\perp D}} \\
        f_{I \rightarrow R} & = \frac{\mu_R}{T_I}I(t) \\
        f_{I \rightarrow D} & = \frac{\mu_D}{T_I}I(t) \\
    \end{aligned}

Note that all other possible flows, like :math:`f_{I \rightarrow S}`, are constant 0.

Flows
~~~~~

To use a flow model, we need to create a list of all flows. These are used by the model to automatically assemble the
compartments. We use a :code:`mio::TypeList` with a :code:`mio::Flow` for each mathematical flow. For the SIRD model
we get:

.. code-block:: cpp

    using Flows = mio::TypeList<mio::Flow<InfectionState::Susceptible, InfectionState::Infectious>,
                                mio::Flow<InfectionState::Infectious,  InfectionState::Recovered>,
                                mio::Flow<InfectionState::Infectious,  InfectionState::Deceased>>;

Define the model
~~~~~~~~~~~~~~~~

With the flows and classes also used by the CompartmentalModel, we can define a FlowModel as such: 

.. code-block:: cpp

    template <typename FP = ScalarType>
    class Model : public mio::FlowModel<FP, InfectionState, Population<FP>, Parameters<FP>, Flows>
    {
    public:
        using Base = mio::FlowModel<FP, InfectionState, Population<FP>, Parameters<FP>, Flows>;
        using typename Base::ParameterSet;
        using typename Base::Populations;

        void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                       Eigen::Ref<Eigen::VectorX<FP>> flows) const override
        {
            const Parameters<FP>& params = this->parameters;

            const auto N = y[InfectionState::Susceptible] + y[InfectionState::Infectious] +
                           y[InfectionState::Recovered];

            flows[this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Infectious>()] =
                params.template get<TransmissionRisk<FP>>() * params.template get<ContactRate<FP>>() *
                y[InfectionState::Susceptible] * y[InfectionState::Infectious] / N;
            
            . . .
        }
    };
