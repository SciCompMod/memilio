Graph-based metapopulation model
================================

Introduction
-------------

Todo Julia

In order to extend the local models with spatial resolution, we provide a module to loosely couple multiple simulation instances and model the mobility between them. Each instance can represent a different geographical region. The regions are stored as nodes in a multi-edge graph as proposed in Kühn et al. (2021), with edges between them representing the mobility between the regions. 
We assume that :math:`n` different geographic units (denoted as *regions*) are given. The mobility extension can contain a arbitrary number of regions. The number of edges is limited by the number of regions squared. 
A pair of nodes is connected via multi edges, where each edge represents the mobility of a specific group of individuals. The number of edges between a pair of nodes is given as the product of the number of age groups and number of compartments in the (local) model.


Currently, only ODE models are supported as nodes; support for other models will be added soon.

For further details, please refer to:

- Zunker H, Schmieding R, Kerkmann D, Schengen A, Diexer S, et al. (2024). *Novel travel time aware metapopulation models and multi-layer waning immunity for late-phase epidemic and endemic scenarios*. *PLOS Computational Biology* 20(12): e1012630. `<https://doi.org/10.1371/journal.pcbi.1012630>`_
- Kühn MJ, Abele D, Mitra T, Koslow W, Abedi M, et al. (2021). *Assessment of effective mitigation and prediction of the spread of SARS-CoV-2 in Germany using demographic information and spatial resolution*. *Mathematical Biosciences* 108648. `<https://doi.org/10.1016/j.mbs.2021.108648>`_

Simulation
-----------
Todo Julia general

Instant mobility - Henrik

Detailed Mobility - Henrik

Stochastic Mobility - Julia
The mobility can contain a different scale of complexity.  One option is the "instant" mobility scheme, where we implement two daily exchanges of commuters that occur instantaneously, at the beginning of the day and after half a day. With this approach,
commuters can change their infection state during the half-day stay. Theres also another more "detailed" mobility scheme, where travel times and stay durations are considered; see Zunker et al. (2024) for more details.

Basicallly, the graph simulation consists of two phases:

1. **Advance the Simulation:**
   
   - Each node progresses independently according to its model.

2. **Exchange of Population:**
   
   - Individuals are exchanged between nodes along the graph edges.
   - The number of exchanged individuals is determined by coefficients.
   - Number of daily commuters may include time-dependent damping factors, similar to the dampings used in contact matrices.
   - While located in another node, the individuals can change their infection state. Therefore, we need to estimate the infection state of the commuters at the time of the return.

How to: Set up a graph and run a graph simulation
-------------------------------------------------

Henrik How to for normal graph

Julia adds for stochastic

The graph simulation couples multiple simulation instances representing different geographic regions via a graph-based approach. In this example, a compartment model based on ODEs (ODE-SECIR) is used for each region. The simulation proceeds by advancing each region independently and then exchanging portions of the population between regions along the graph edges.

The following steps detail how to configure and execute a graph simulation:

1. **Initialize the local compartment model:**

   First, set up the compartment model by initializing the populations and parameters. In this example, we use a single age group and start with 10,000 susceptible individuals. Additionally, the necessary epidemiological parameters (e.g. time periods, transmission probabilities) are set.

   .. code-block:: cpp

           const size_t num_groups = 1;
           mio::osecir::Model model(num_groups);
           model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = 10000;
           model.parameters.set<mio::osecir::StartDay>(0);
           model.parameters.set<mio::osecir::Seasonality<ScalarType>>(0.2);
       
           model.parameters.get<mio::osecir::TimeExposed<ScalarType>>()            = 3.2;
           model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>() = 2.0;
           model.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()   = 5.8;
           model.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()     = 9.5;
           model.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()   = 7.1;
       
           model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()  = 0.1;
           model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()    = 0.7;
           model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()    = 0.09;
           model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()    = 0.25;
           model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<ScalarType>>() = 0.45;
           model.parameters.get<mio::osecir::TestAndTraceCapacity<ScalarType>>()              = 35;
           model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()         = 0.2;
           model.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()                 = 0.25;
           model.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()                 = 0.3;
       
           mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
           contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));

2. **Create simulation groups and adjust contact patterns:**

   To represent different geographic regions, clone the base model into separate model groups. In this example, two model groups are created. The first group is modified by applying a contact damping to simulate contact restrictions.

   .. code-block:: cpp

           // Create two mostly identical groups
           auto model_group1 = model;
           auto model_group2 = model;
       
           // Apply contact restrictions to model_group1
           mio::ContactMatrixGroup& contact_matrix_m1 =
               model_group1.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
           contact_matrix_m1[0].add_damping(0.7, mio::SimulationTime(15.));
       
           // Initialize infection in group 1
           model_group1.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = 9990;
           model_group1.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]     = 100;

3. **Define compartments to save from edges:**

   To extract the mobility results, define the compartments to save from the edges. In this example, the compartments for infected individuals with and without symptoms are saved for each region.

   .. code-block:: cpp

           // Get indices of INS and ISy compartments.
           std::vector<std::vector<size_t>> indices_save_edges(2);
           for (auto& vec : indices_save_edges) {
               vec.reserve(2 * num_groups);
           }
           for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(num_groups); ++i) {
               indices_save_edges[0].emplace_back(
                   model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedNoSymptoms}));
               indices_save_edges[0].emplace_back(
                   model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}));
               indices_save_edges[1].emplace_back(
                   model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedSymptoms}));
               indices_save_edges[1].emplace_back(
                   model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}));
           }

4. **Construct the Mobility Graph:**

   Build a graph where each node represents a simulation and each edge represents mobility between a pair of nodes. Mobility coefficients (here, 0.1 for all compartments) determine the fraction of the population exchanged between nodes.

   .. code-block:: cpp
    
           const auto t0   = 0.;
           mio::Graph<mio::SimulationNode<mio::osecir::Simulation<>>, mio::MobilityEdge<ScalarType>> g;
           g.add_node(1001, model_group1, t0);
           g.add_node(1002, model_group2, t0);
           g.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)mio::osecir::InfectionState::Count, 0.1), indices_save_edges);
           g.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)mio::osecir::InfectionState::Count, 0.1), indices_save_edges);

5. **Initialize and Advance the Mobility Simulation:**

   With the graph constructed, initialize the simulation with the starting time and time step. Then, advance the simulation until the final time :math:`t_{max}`.

   .. code-block:: cpp
            
           const auto tmax = 30.;
           const auto dt   = 0.5; // time step or Mobility (daily mobility occurs every second step)
           auto sim = mio::make_mobility_sim(t0, dt, std::move(g));
           sim.advance(tmax);

6. **Access and Display Mobility Results:**

   After the simulation, the mobility results can be extracted from a specific edge. In this example, the results for the edge from node 1 to node 0 are printed.

   .. code-block:: cpp

           auto& edge_1_0 = sim.get_graph().edges()[1];
           auto& results  = edge_1_0.property.get_mobility_results();
           results.print_table({"Commuter INS", "Commuter ISy", "Commuter Total"});
       
           return 0;
       }
