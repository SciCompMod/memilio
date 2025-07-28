Stochastic metapopulation model
===============================

The stochastic metapopulation model uses a Markov process to simulate disease dynamics. Similar to the `Diffusive Agent-based Model` the Markov process is given by location and infection state changes. However, in contrast to the diffusive ABM, location changes are not given by agent-dependent diffusion processes, but by stochastic jumps between regions with the requirement that the domain is split into disjoint regions. Hence, there is no further spatial resolution within one region and locations or positions are only given by the region index. The evolution of the system state is determined by the following master equation

:math:`\partial_t p(X,Z;t) = G p(X,Z;t) + L p(X,Z;t)`.

The operator :math:`G` defines the infection state adoptions and only acts on :math:`Z`, the vector containing all subpopulations stratified by infection state. :math:`L` defines location changes, only acting on :math:`X`, the vector containing all subpopulations stratified by region. Infection state adoptions are modeled as stochastic jumps with independent Poisson processes given by adoption rate functions. Similar to the infection state dynamics, spatial transitions between regions are also modeled as stochastic jumps with independent Poisson processes given by transition rate functions. Gillespie's direct method (stochastic simulation algorithm) is used for simulation.


Infection states
----------------

The model does not have fixed infection states, but gets an enum class of infection states as template argument. Thus it can be used with any set of infection states.
Using the infection states Susceptible (S), Exposed (E), Carrier (C), Infected (I), Recovered (R) and Dead (D), this can be done as follows:

.. code-block:: cpp

    enum class InfectionState
    {
        S,
        E,
        C,
        I,
        R,
        D,
        Count

    };

    const size_t num_regions = 2;

    using Model              = mio::smm::Model<num_regions, InfectionState>;
    Model model;

Infection state transitions
---------------------------

The infection state transitions are explicitly given by the adoption rates and are therefore subject to user input. Adoption rates always depend on their source infection state. If an adoption event requires interaction of agents (e.g. disease transmission), the corresponding rate depends not only on the source infection state, but also on other infection states, the **Influence**\s. An adoption rate that only depends on the source infection state, e.g. recovery or worsening of disease symptoms, is called `first-order` adoption rate and an adoption rate that has influences is called `second-order` adoption rate. Adoption rates are region-dependent; therefore it is possible to have different rates in two regions for the same infection state transition which can be useful when having e.g. region-dependent interventions or contact behavior.

Using the infection states from above and two regions, there are five first-order adoption rates per region and one second-order adoption rate per region. In the example below, the second-order adoption rate (transition from S to E) differs between the regions:

.. code-block:: cpp

   std::vector<mio::AdoptionRate<InfectionState>> adoption_rates;

   //Set first-order adoption rates for both regions
   for (size_t r = 0; r < num_regions; ++r) {
      adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::regions::Region(r), 1.0 / 5., {}});
      adoption_rates.push_back({InfectionState::C, InfectionState::R, mio::regions::Region(r), 0.2 / 3., {}});
      adoption_rates.push_back({InfectionState::C, InfectionState::I, mio::regions::Region(r), 0.8 / 3., {}});
      adoption_rates.push_back({InfectionState::I, InfectionState::R, mio::regions::Region(r), 0.99 / 5., {}});
      adoption_rates.push_back({InfectionState::I, InfectionState::D, mio::regions::Region(r), 0.01 / 5., {}});
   }

  //Set second-order adoption rate different for the two regions
   adoption_rates.push_back({InfectionState::S, InfectionState::E, mio::regions::Region(0), 0.1, {{InfectionState::C, 1}, {InfectionState::I, 0.5}}});
   adoption_rates.push_back({InfectionState::S, InfectionState::E, mio::regions::Region(1), 0.2, {{InfectionState::C, 1}, {InfectionState::I, 0.5}}});

   //Initialize model parameter
   model.parameters.get<mio::smm::AdoptionRates<InfectionState>>()   = adoption_rates;

Sociodemographic stratification
-------------------------------

Sociodemographic stratification e.g. by age, gender or immunity can be incorporated by stratifying the set of infection states passed as template to the model.

Parameters
----------

The model has the following parameters:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Mathematical variable
     - C++ variable name
     - Description
   * - :math:`\gamma^{(k)}_{i,j}`
     - ``AdoptionRate``
     - Adoption rate in region k from infection state i to state j. Apart from the region k, the source (i) and target (j) infection state, the adoption rates get influences :math:`\tau \in \Psi` and a constant :math:`c_{i,j}`.
   * - :math:`\tau`
     - ``Influence``
     - Influence for second-order adoption rate consisting of the influencing infection state and a factor with which the population having the corresponding infection state is multiplied.
   * - :math:`\lambda^{(k,l)}_{i}`
     - ``TransitionRate``
     - Spatial transition rate for infection state i from region k to region l.

The adoption rate :math:`\gamma^{(k)}_{i,j}` at time :math:`t` is given by

:math:`\gamma^{(k)}_{i,j}(t) = c_{i,j}\frac{i^{(k)}{N}\cdot\sum_{\tau \in \Psi}\tau.factor \cdot \tau.state(t)`

and the spatial transition rate at time :math:`t` by

 :math:`\lambda^{(k,l)}_{i} = \lambda^{(k,l)}_{i}.factor\cdot i^{(k)}(t)`

with :math:`i^{(k)}` the population in region :math:`k` having infection state :math:`i`.


Initial conditions
------------------

Before running a simulation with the stochastic metapopulation model, the initial populations i.e. the number of agents per infection state for every region have to be set.
These populations have the class type **Populations** and can be set via:

.. code-block:: cpp

   double pop = 1000, numE = 0.001 * pop, numC = 0.0001 * pop, numI = 0.0001 * pop, numR = 0, numD = 0;

   //Population is distributed equally to the regions
   for (size_t r = 0; r < num_regions; ++r) {
        model.populations[{mio::regions::Region(r), InfectionState::S}] = (pop - numE - numC - numI - numR - numD) / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::E}] = numE / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::C}] = numC / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::I}] = numI / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::R}] = 0;
        model.populations[{mio::regions::Region(r), InfectionState::D}] = 0;
   }

If individuals should transition between regions, the spatial transition rates of the model have to be initialized as well.
As the spatial transition rates are dependent on infection state, region changes for specific infection states can be prevented. Below, symmetric spatial transition rates are set for every region:

.. code-block:: cpp

   std::vector<mio::smm::TransitionRate<InfectionState>> transition_rates;
   //Agents in infection state D do not transition
   for (size_t s = 0; s < static_cast<size_t>(InfectionState::D); ++s) {
      for (size_t i = 0; i < num_regions; ++i) {
         for (size_t j = 0; j < num_regions; ++j)
               if (i != j) {
            transition_rates.push_back(
               {InfectionState(s), mio::regions::Region(i), mio::regions::Region(j), 0.01});
            transition_rates.push_back(
               {InfectionState(s), mio::regions::Region(j), mio::regions::Region(i), 0.01});
         }
      }
   }

   //Initialize model parameter
   model.parameters.get<mio::smm::TransitionRates<InfectionState>>() = transition_rates;

Nonpharmaceutical interventions
--------------------------------

There are no nonpharmaceutical interventions (NPIs) explicitly implemented in the model. However, NPIs influencing the adoption or spatial transition rates can be realized by resetting the corresponding model parameters.

Simulation
----------

At the beginning of the simulation, the waiting times for all events (infection state adoptions and spatial transitions) are drawn. Then the time is advanced until the time point of the next event - which can be a spatial transition or an infection state adoption - and the event takes places. The waiting times of the other events are updated and a new waiting time for the event that just happened is drawn. The simulation saves the system state in discrete time steps.

To simulate the model from `t0` to `tmax` with given step size `dt`, a **Simulation** has to be created and advanced until `tmax`. The step size is only used to regularly save the system state during the simulation.

.. code-block:: cpp

    double t0   = 0.0;
    double dt   = 0.1;
    double tmax = 30.;

    //Pass the model, t0 and dt to the Simulation
    auto sim = mio::smm::Simulation(model, t0, dt);

    //Advance the simulation until tmax
    sim.advance(tmax);

Output
------

Subpopulations stratified by region and infection state are saved in a ``mio::TimeSeries`` object which can be accessed and printed as follows:

.. code-block:: cpp

    //Result object has size num_time_points x (num_infection_states * num_regions)
    auto result = sim.get_result();

    //Print result object to console. Infection state "Xi" with i=0,1 is the number of agents having infection state X in region i
    result.print_table({"S0", "E0", "C0", "I0", "R0", "D0", "S1", "E1", "C1", "I1", "R1", "D1"})

If one wants to interpolate the aggregated results to a ``mio::TimeSeries`` containing only full days, this can be done by

.. code-block:: cpp

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());

Examples
--------

An example of the stochastic metapopulation model with four regions can be found at: `examples/smm.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/smm.cpp>`_


Overview of the ``smm`` namespace:
-----------------------------------

.. doxygennamespace:: mio::smm
