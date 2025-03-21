Diffusive Agent-Based Model
===========================

This agent-based model uses a Markov process to simulate disease dynamics. The Markov process is given by agent movement and the evolution of disease dynamics i.e. disease transmission and progression.
The features of an agent are its position and its infection state. The evolution of the system state is determined by the following master equation

:math:`\partial_t p(X,Z;t) = G p(X,Z;t) + L p(X,Z;t)`

with :math:`X` the vector of all agents' positions and :math:`Z` the vector of all agents' infection states. The operator :math:`G` defines the infection state adoptions and only acts on :math:`Z`, while :math:`L` defines movement, i.e. location changes, only acting on :math:`X`. Infection state adoptions are modeled with independent Poisson processes given by adoption rate functions. Adoption rates always depend on the source infection state. If an adoption event requires interaction of agents (e.g. disease transmission), the corresponding rate depends not only on the source infection state, but also on other infection states, the `influences`. An adoption rate that only depends on the source infection state, e.g. recovery or worsening of disease symptoms, is called `first-order` adoption rate and an adoption rate that has influences is called `second-order` adoption rate. Movement is modeled with independent diffusion processes. A temporal Gillespie algorithm (a direct method without rejection sampling) is used for simulation. Therefore, :math:`G` and :math:`L` are not implemented explicitly, instead their effects are sampled via the `move` and `adoption_rate` functions, respectively.

The Model class needs an Implementation class as template argument which provides the domain agents move and interact in. A quadwell potential given in the class ``QuadWellModel`` is implemented, but any other suitable potential can be used as implementation. 

Simulation
-----------

The simulation runs in discrete time steps. In every step the model is advanced until the next infection state adoption event. Then the corresponding agent's infection state is adopted and a new waiting time until the next adoption event is drawn. If the waiting time until the next adoption event is bigger than the remaining time in the time step, we advance the model until the end of the time step.

For a detailed description and application of the model, see:

- Bicker J, Schmieding R, et al. (2025) Hybrid metapopulation agent-based epidemiological models for efficient insight on the individual scale: A contribution to green computing. Infectious Disease Modelling, Volume 10, Issue 2. https://doi.org/10.1016/j.idm.2024.12.015

How to: Set up and run a simulation of the diffusive ABM
---------------------------------------------------------

The infection states used need to be defined in an enum before initializing the model and are passed as template argument to the model. 
Using the infection states Susceptible (S), Exposed (E), Carrier (C), Infected (I), Recovered (R) and Dead (D) this can be done as follows:

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

    using Model = mio::dabm::Model<QuadWellModel<InfectionState>>;

To initialze the model using the quad well potential, the following inputs need to be passed to the model constructor:

- a vector containing the agents,
- a vector containing the adoption rates for infection state adoptions,
- the interaction radius of agents,
- the noise term for the diffusion process.

Agents have two attributes: A position on the domain and an infection state. The example below initializes 100 agents with an agent's position sampled uniformly from :math:`\left[-2,2\right]\times\left[-2,2\right]` and its infection state sampled from a discrete distribution with probabilities given by :math:`98\%` (S), :math:`1\%` (E), :math:`0.5\%` (C), :math:`0.5\%` (I), :math:`0\%` (R), :math:`0\%` (D).  

.. code-block:: cpp

    std::vector<Model::Agent> agents(100);

    //Random variables for initialization of agents' position and infection state
    auto& pos_rng = mio::UniformDistribution<double>::get_instance();
    auto& sta_rng = mio::DiscreteDistribution<size_t>::get_instance();

    //Infection state distribution
    std::vector<double> pop_dist{0.98, 0.01, 0.005, 0.005, 0., 0.};

    for (auto& a : agents) {
        //Agents are uniformly distributed in [-2,2]x[-2,2]
        a.position = Eigen::Vector2d{pos_rng(mio::thread_local_rng(), -2., 2.), pos_rng(mio::thread_local_rng(), -2., 2.)};
        a.status = static_cast<InfectionState>(sta_rng(mio::thread_local_rng(), pop_dist));
    }

As already mentioned above, infection state adoption rates can only depend on the source infection state (first-order) or they can be influenced by other infection states (second-order). An influence consists of an infection state and the corresponding multiplicative factor. Additionally, adoption rates are region-dependent with the regions corresponding to the potential wells for the implemented quad well potential. That means the quad well potential has four regions given by :math:`\Omega_0=\left(-\infty,0\right)\times\left(0,\infty\right), \Omega_1=\left(0,\infty\right)\times\left(0,\infty\right), \Omega_2=\left(-\infty,0\right)\times\left(-\infty,0\right)` and :math:`\Omega_3=\left(0,\infty\right)\times\left(-\infty,0\right)`.

In the example below, regions 0-2 have the same adoption rates while region 3 has a two times higher transmission rate (adoption from S to E):  

.. code-block:: cpp

    std::vector<mio::AdoptionRate<InfectionState>> adoption_rates;

    //First-order adoption rates
    for (size_t region = 0; region < 4; ++region) {
        adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::regions::Region(region), 1.0 / 5., {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::R, mio::regions::Region(region), 0.2 / 3., {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::I, mio::regions::Region(region), 0.8 / 3., {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::R, mio::regions::Region(region), 0.99 / 5., {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::D, mio::regions::Region(region), 0.01 / 5., {}});
    }

    //Second-order adoption rate for regions 0-2
    for (size_t region = 0; region < 3; ++region) {        
        adoption_rates.push_back({InfectionState::S, InfectionState::E, mio::regions::Region(region), 0.1, {{InfectionState::C, 1}, {InfectionState::I, 0.5}}});
    }

    //Second-order adoption rate for region 3
    adoption_rates.push_back({InfectionState::S, InfectionState::E, mio::regions::Region(3), 0.2, {{InfectionState::C, 1}, {InfectionState::I, 0.5}}});

Choosing an interaction radius of 0.5 and a noise term of 0.4, the model is initialized by

.. code-block:: cpp

    double interaction_radius = 0.5;
    double noise = 0.4;

    Model model(agents, adoption_rates, interaction_radius, noise);

There is the option to pass a vector of non-moving infection states to the model. Defining that agents in infection state D do not move, the model initialization would be

.. code-block:: cpp

    Model model(agents, adoption_rates, interaction_radius, noise, {InfectionState::D});

To simulate the model from `t0` to `tmax` with given step size `dt`, a Simulation has to be created and advanced until `tmax`, which is done as follows:

.. code-block:: cpp

    double t0   = 0.0;
    double dt   = 0.1;
    double tmax = 30.;

    //Pass the model, t0 and dt to the Simulation
    auto sim = mio::dabm::Simulation(model, t0, dt);

    //Advance the simulation until tmax
    sim.advance(tmax);

The result vector containing all agents can be accessed via 

.. code-block:: cpp

    sim.get_model().populations

Additionally, the agents are automatically aggregated by region and infection state in a ``mio::TimeSeries`` object which can be accessed and printed as follows:

.. code-block:: cpp

    //Result object has size num_time_points x (num_infection_states * num_regions)
    auto result = sim.get_result();

    //Print result object to console. Infection state "Xi" with i=0,...,3 is the number of agents having infection state X in region i
    result.print_table({"S0", "E0", "C0", "I0", "R0", "D0", "S1", "E1", "C1", "I1", "R1", "D1", "S2", "E2", "C2", "I2", "R2", "D2", "S3", "E3", "C3", "I3", "R3", "D3"})

If one wants to interpolate the aggregated results to a ``mio::TimeSeries`` containing only full days, this can be done by

.. code-block:: cpp

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());
