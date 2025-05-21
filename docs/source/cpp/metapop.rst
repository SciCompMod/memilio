Metapopulation Models
=====================

Stochastic Metapopulation Model
-------------------------------

The stochastic metapopulation model uses a Markov process to simulate disease dynamics. Similar to the `Diffusive Agent-based Model` the Markov process is given by location and infection state changes. However, in contrast to the diffusive ABM, location changes are not given by agent-dependent diffusion processes, but by stochastic jumps between regions with the requirement that the domain is split into disjoint regions. Hence, there is no further spatial resolution within one region and locations or positions are only given by the region index. The evolution of the system state is determined by the following master equation

:math:`\partial_t p(X,Z;t) = G p(X,Z;t) + L p(X,Z;t)`.

The operator :math:`G` defines the infection state adoptions and only acts on :math:`Z` the vector containing all subpopulations stratified by infection state. :math:`L` defines location changes, only acting on :math:`X` the vector containing all subpopulations stratified by region. Infection state adoptions are modeled as stochastic jumps with independent Poisson processes given by adoption rate functions. Adoption rates always depend on the source infection state. If an adoption event requires interaction of agents (e.g. disease transmission), the corresponding rate depends not only on the source infection state, but also on other infection states, the `influences`. An adoption rate that only depends on the source infection state, e.g. recovery or worsening of disease symptoms, is called `first-order` adoption rate and an adoption rate that has influences is called `second-order` adoption rate. Similar to the infection state dynamics, spatial transitions between regions are also modeled as stochastic jumps with independent Poisson processes given by transition rate functions. Gillespie's direct methode (stochastic simulation algorithm) is used for simulation.

How to: Set up and run a simulation of the stochastic metapopulation model
--------------------------------------------------------------------------

The infection states used need to be defined in an enum before initializing the model and are passed - together with the number of regions - as template arguments to the model. 
Using four regions and the infection states Susceptible (S), Exposed (E), Carrier (C), Infected (I), Recovered (R) and Dead (D) this can be done as follows:

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

    const size_t num_regions = 4;

    using Model              = mio::smm::Model<num_regions, InfectionState>;
    Model model;

To initialize the model, the following parameters have to be set:

- the initial subpopulations i.e. the number of agents per infection state for every region,
- the adoption rates for infection state adoptions,
- the spatial transition rates for location (region) changes.

The example below initializes a population of size 1000 that is equally distributed to the four regions. In every region, :math:`0.1\%` of the population is Exposed, :math:`0.01\%` is Carrier, :math:`0.01\%` is Infected and the rest of the population is Susceptible. There are no Recovered or Dead agents in the beginning.

.. code-block:: cpp

   double pop = 1000, numE = 0.001 * pop, numC = 0.0001 * pop, numI = 0.0001 * pop;

   //Population is distributed equally to the four regions
   for (size_t r = 0; r < num_regions; ++r) {
        model.populations[{mio::regions::Region(r), InfectionState::S}] = (pop - numE - numC - numI - numR - numD) / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::E}] = numE / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::C}] = numC / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::I}] = numI / num_regions;
        model.populations[{mio::regions::Region(r), InfectionState::R}] = 0;
        model.populations[{mio::regions::Region(r), InfectionState::D}] = 0;
   }

As already mentioned above, infection state adoption rates can only depend on the source infection state (first-order) or they can be influenced by other infection states (second-order). An influence consists of an infection state and the corresponding multiplicative factor. Additionally, adoption rates are region-dependent. In the example below, all regions have the same adoption rates.

.. code-block:: cpp

   //Set infection state adoption rates
   std::vector<mio::AdoptionRate<InfectionState>> adoption_rates;
   for (size_t r = 0; r < num_regions; ++r) {
      //Second-order adoption rate
      adoption_rates.push_back({InfectionState::S, InfectionState::E, mio::regions::Region(r), 0.1, {{InfectionState::C, 1}, {InfectionState::I, 0.5}}});
      //First-order adoption rate
      adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::regions::Region(r), 1.0 / 5., {}});
      adoption_rates.push_back({InfectionState::C, InfectionState::R, mio::regions::Region(r), 0.2 / 3., {}});
      adoption_rates.push_back({InfectionState::C, InfectionState::I, mio::regions::Region(r), 0.8 / 3., {}});
      adoption_rates.push_back({InfectionState::I, InfectionState::R, mio::regions::Region(r), 0.99 / 5., {}});
      adoption_rates.push_back({InfectionState::I, InfectionState::D, mio::regions::Region(r), 0.01 / 5., {}});
   }

   model.parameters.get<mio::smm::AdoptionRates<InfectionState>>()   = adoption_rates;

The spatial transition rates are dependent on infection state such that location changes for specific infection states can be prevented. Below, symmetric spatial transition rates are set for every region:

.. code-block:: cpp

   //Set spatial transition rates
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

   model.parameters.get<mio::smm::TransitionRates<InfectionState>>() = transition_rates;

To simulate the model from `t0` to `tmax` with given step size `dt`, a Simulation has to be created. The step size is only used to regularly save the system state during the simulation. Advancing the simulation until `tmax` is done as follows:

.. code-block:: cpp

    double t0   = 0.0;
    double dt   = 0.1;
    double tmax = 30.;

    //Pass the model, t0 and dt to the Simulation
    auto sim = mio::smm::Simulation(model, t0, dt);

    //Advance the simulation until tmax
    sim.advance(tmax);

Subpopulations stratified by region and infection state are saved in a ``mio::TimeSeries`` object which can be accessed and printed as follows:

.. code-block:: cpp

    //Result object has size num_time_points x (num_infection_states * num_regions)
    auto result = sim.get_result();

    //Print result object to console. Infection state "Xi" with i=0,...,3 is the number of agents having infection state X in region i
    result.print_table({"S0", "E0", "C0", "I0", "R0", "D0", "S1", "E1", "C1", "I1", "R1", "D1", "S2", "E2", "C2", "I2", "R2", "D2", "S3", "E3", "C3", "I3", "R3", "D3"})

If one wants to interpolate the aggregated results to a ``mio::TimeSeries`` containing only full days, this can be done by

.. code-block:: cpp

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());
