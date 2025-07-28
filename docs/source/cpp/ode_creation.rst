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
states **S**\usceptible, **I**\nfectious, **R**\ecovered, **D**\eceased, also called compartments.

How to define an ODE model
--------------------------

To define an ODE model in MEmilio, there are two options. You can define a CompartmentalModel or a FlowModel, which
use different methods to define the right hand side of the mathematical model above. Both classes need definitions for
the infection states, population and parameters it uses. The FlowModel additionally requires a list of flows.

We start by creating a new directory for our model under "cpp/models", in this case we can call it "ode_sird". The name
must be unique and start with "ode\_", so the type of model is obvious. The rest usually contains the compartments or
other noteworthy features of the model in shortened form. All files in the following are put into this directory.

Infection states
~~~~~~~~~~~~~~~~

First we define an :code:`enum class` called ``InfectionState`` in the file "infection_state.h", which contains an entry
for each infection state of the mathematical model, followed by an entry called :code:`Count`. This enumerates the 
compartments starting from 0, with Count being equal to the number of compartments. In our example we have:

.. code-block:: cpp

    enum class InfectionState
    {
        Susceptible,
        Infectious,
        Recovered,
        Deceased,
        Count
    };

.. dropdown:: :fa:`gears` Expert's settings

    We use an ``enum class`` instead of a regular ``enum`` here, since it helps avoid unintended conversions, invalid
    array accesses, and "magic numbers" - you cannot simply write "2" and use it as ``InfectionState::Recovered``, since
    an ``enum class`` does not allow implicit conversions.

Parameters
~~~~~~~~~~

Next, we define the parameters in "parameters.h", which consist of a struct for each constant used in the mathematical
model. This struct must define the data type, name and default value of the constant. For example, for the time a
person stays infectious, :math:`T_I`, we define a struct

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

The template :code:`FP` and the type :code:`UncertainValue<FP>` in these examples are commonly used throughout MEmilio.
:code:`FP` is a floating point type, usually :code:`double`. An :code:`UncertainValue<FP>` holds a value of type
:code:`FP` as well as (optionally) a distribution to sample new values from, e.g. for a parameter study.

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
``mio::Populations``.

.. code-block:: cpp

    template <typename FP>
    using Population = mio::Populations<FP, InfectionState>;

Importantly, this class allows further stratifying the population vector, with the most common
example being adding :code:`mio::AgeGroups` to the template.

.. dropdown:: :fa:`gears` Expert's settings

    The type ``mio::AgeGroup`` is a typesafe ``size_t``, meaning an integer that cannot be confused with other integer
    types. So assignment, addition, etc. only works with another ``mio::AgeGroup``, not ``size_t`` or another integer
    type. This is useful for function interfaces or indexing, as it makes it (nearly) impossible to mix up, e.g., age
    groups with infection states. Check out ``mio::Index`` if you want to learn more.

    The type ``mio::Populations`` is an extension of a ``mio::CustomIndexArray``, which is a template type that manages
    a flat array. Its main purpose is to allow multidimensional indexing into this array, using typesafe indices like
    a ``mio::Index`` or a ``enum class``.

    The definition of our ``Population`` then changes to

    .. code-block:: cpp

        template <typename FP>
        using Population = mio::Populations<FP, InfectionState, AgeGroup>;

    and the access (compare with the model definition below) changes to, e.g.,

    .. code-block:: cpp

        const AgeGroup i = . . .;
        const size_t Ri = this->populations.get_flat_index({i, InfectionState::Susceptible});
        dydt[Ri] = . . . * y[Ri];

    where we use ``populations.get_flat_index`` to get the correct index in the flat state and derivative vectors.
    You may also want to change the Parameters to use age groups, check out the available ODE models as reference. 

Define a compartmental model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we can define the model:

.. code-block:: cpp

    template <typename FP>
    class Model : public mio::CompartmentalModel<FP, InfectionState, Population<FP>, Parameters<FP>>
    {
    public:
        using Base = mio::CompartmentalModel<FP, InfectionState, Population<FP>, Parameters<FP>>;
        using typename Base::ParameterSet;
        using typename Base::Populations;

        Model()
            : Base(Populations({InfectionState::Count}), ParameterSet())
        {
        }

        void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                             Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
        {
            const Parameters<FP>& params = this->parameters;

            const auto N = y[InfectionState::Susceptible] + y[InfectionState::Infectious] +
                           y[InfectionState::Recovered];

            dydt[InfectionState::Susceptible] = -params.template get<TransmissionRisk<FP>>() *
                                                params.template get<ContactRate<FP>>() *
                                                y[InfectionState::Susceptible] * pop[InfectionState::Infectious] / N;
            
            . . .
        }
    };

Here, create a new class ``Model`` that inherits from ``mio::CompartmentalModel``, which predefines some functions (like
``check_constraints``) and members (like ``populations`` and ``parameters``) for us. In the class body, we first add
a few ``using`` statements, that create shorthands for our base class and its types. Next, we define a constructor, that
creates populations and parameters with their default values, so they can be set later. Finally, we define the the right
hand side of the model equations through ``get_derivatives``. Importantly, note that the value in ``populations`` is
only used as initial value for the IVP, the current state of the model at time ``t`` is given by ``y`` instead. The
derivative of ``y`` at ``t`` is to be stored in ``dydt``.

Note that the argument ``pop`` is used once instead of ``y``. As a general rule, use ``y`` when the index matches with
the one used for ``dydt``, and ``pop`` otherwise.

Essentially, ``pop`` is the population that is interacted with, while ``y`` is the acting/changing population. Our graph
models use this to model the exchange between multiple models, e.g. by setting ``pop = y ± commuters``. Outside of graph
models both ``y`` and ``pop`` will have the same value.

Check out the page on the usage of :doc:`ODE-based models<cpp/ode>`.


Define a flow model
^^^^^^^^^^^^^^^^^^^

A flow model is a special case of a compartmental model, where the derivative of each compartment over time
:math:`Z_i'(t)` can be written as

.. math::

    Z_i'(t) = \sum_{i \ne j} f_{Z_j \rightarrow Z_i}(t) - \sum_{i \ne j} f_{Z_i \rightarrow Z_j}(t),

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

The first term in each :code:`mio::Flow` is the source compartment, the second the target. As a convention, we always
compute non-negative outflows. Hence, we only list the flow :math:`S \rightarrow I`, but not :math:`I \rightarrow S`.

Infection states, Parameters and Population
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since a Flow model is just a special case of a compartmental model, all of these are defined exactly as described above.

Define the model
~~~~~~~~~~~~~~~~

With the flows and classes also used by the CompartmentalModel, we can define a FlowModel as such: 

.. code-block:: cpp

    template <typename FP>
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
                y[InfectionState::Susceptible] * pop[InfectionState::Infectious] / N;
            
            . . .
        }
    };

This is mostly analoguous to the definition of a compartmental model, with a few important differences. First, we now
inherit from ``FlowModel``, which gets the ``Flows`` as an additional template argument. The ``Base`` alias changes
accordingly. Secondly, the function we implement is called ``get_flows`` and computes the derivative of y in terms of
its flows.

To index into the ``flows`` vector we use the function ``get_flat_flow_index``, which takes the source and target
compartments as template arguments, in that order. Indexes from further stratification (like ``mio::AgeGroup``) can be
passed as an optional function argument.
