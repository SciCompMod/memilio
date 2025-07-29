SDE model creation
==================

The mathematical model
----------------------

We assume that our stochastic models are a system of initial value problems of the following form:

.. math::

    \mathrm{d}Z_t = a(Z_t, t) \mathrm{d}t + b(Z_t, t)\mathrm{d}W_t

where :math:`Z_t \in \mathcal{R}^n` and initial values :math:`Z_0 = z_0 \in \mathcal{R}^n`. The function
:math:`a : \mathcal{R}^n \times \mathcal{R} \rightarrow \mathcal{R}^n` is called drift coefficient,
:math:`b : \mathcal{R}^n \times \mathcal{R} \rightarrow \mathcal{R}^{n \times m}` is called diffusion coefficient or
noise term for some :math:`m`, e.g. the number of flows. Both functions are deterministic, the stochasticity comes from
the (formal) derivative of the Brownian motion :math:`\mathrm{d}W_t` in :math:`\mathcal{R}^m`.


Analogous to the ODE models, :math:`Z_t` describes the sizes of compartments at time point :math:`t`, therefore the
system

.. math::

    \mathrm{d}Z_t = a(Z_t, t) \mathrm{d}t

(i.e. with :math:`b \equiv 0`) is an ODE model as described on the :doc:`ODE model creation <cpp/ode_creation>` page.

How to define an SDE model
--------------------------

In short, to define an SDE model in MEmilio, you have to implement a ``StochasticModel``, e.g. by inheriting from it.
To that end, we first need to define types that list all ``InfectionState``\s, ``Parameter``\s and initial conditions
via a ``Population``. Refer to the :doc:`ODE model creation <ode_creation>` page for more details on these types.

A valid SDE model needs to implement one function each for the deterministic part
:math:`a(Z_t, t) \mathrm{d}t` and stochastic part :math:`b(Z_t, t)\mathrm{d}W_t`. For the deterministic part, we require
a ``get_derivatives`` function, for the stochastic part a ``get_noise`` function.

Hence, you can define a ``StochasticModel`` either as a ``CompartmentalModel``

.. code:: cpp

    class Model : public mio::StochasticModel<ScalarType, InfectionState, Population<ScalarType>,
                                              Parameters<ScalarType>>
    {
    public:
        using Base = mio::StochasticModel<ScalarType, InfectionState, Population<ScalarType>, Parameters<ScalarType>>;
        using typename Base::ParameterSet;
        using typename Base::Populations;

        Model()
            : Base(Populations({InfectionState::Count}), ParameterSet())
        {
        }

        void get_derivatives(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop,
                             Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                             ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> dydt) const
        {
            . . .
        }

        void get_noise(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                       ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> noise) const
        {
            . . .
        }
    };


or alternatively as a ``FlowModel`` by additionally providing the list of ``Flows``

.. code:: cpp

    class Model : public mio::StochasticModel<ScalarType, InfectionState, Population<ScalarType>,
                                              Parameters<ScalarType>, Flows>
    {
    public:
        using Base = mio::StochasticModel<ScalarType, InfectionState, Population<ScalarType>, Parameters<ScalarType>,
                                          Flows>;
        using typename Base::ParameterSet;
        using typename Base::Populations;

        Model()
            : Base(Populations({InfectionState::Count}), ParameterSet())
        {
        }

        void get_flows(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                       ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> flows) const
        {
            . . .
        }

        void get_noise(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                       ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> noise) const
        {
            . . .
        }
    };

In both cases the computed ``noise`` vector must have the same size as the vectors ``pop`` and ``y``, i.e. the number of
compartments. For more details on how to implement the ``get_derivatives`` or ``get_flows`` methods check out the
:doc:`ODE model creation <cpp/ode_creation>` page.

.. dropdown:: :fa:`gears` Expert's knowledge

    The SDE models must work on compartments (rather than flows) due to the stochasticity being able to cause negative
    compartment values during integration, which usually make no sense in infectious disease models, so we use a
    mitigation against these negative values (see the function ``mio::map_to_nonnegative``). However, such a mitigation
    can only be applied to compartments, and, in general, propagation of changes on compartments back to flows is not
    possible.
    A ``FlowModel`` can still be used, since it defines a ``get_derivatives`` function based on the provided
    ``get_flows``.

    Note that we use ``ScalarType`` instead of an ``FP`` template. The main reason is that we are not certain that AD
    types work well with the random numbers in the model, so we recommend using ``ScalarType`` instead.

The ``StochasticModel`` base class comes with a random number generator that can be accessed via ``get_rng``, as well
as a method ``sample_standard_normal_distribution`` to draw a single random number as well as a function ``white_noise``
that returns a vector expression of independent standard normal distributed values. You can use these to implement
the ``get_noise`` function.

You may want to use a ``FlowModel`` if your noise depends on the current flow values. In that case, the noise matrix
:math:`b` may map each flow's noise contribution to its source and/or target compartment. In that case, the size of the
white noise vector :math:`m` is equal to the number of flows.

An example for a get_noise function from one of the bundled SDE models looks like this:

.. code:: cpp

    void get_noise(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                   ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> noise) const
    {
        Eigen::VectorX<ScalarType> flows(Flows::size());
        get_flows(pop, y, t, flows);
        flows = flows.array().sqrt() * Base::white_noise(Flows::size()).array();
        get_derivatives(flows, noise);
    }

Here we first compute the flows, then take the square root of each flow and multiply it by a standard normal distributed
value. The mapping from flows to compartments (that is mathematically done by a matrix multiplication) is taken care of
by the overload of ``get_derivatives``.
