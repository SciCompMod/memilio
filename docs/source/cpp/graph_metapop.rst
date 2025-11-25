.. include:: ../literature.rst

Graph-based metapopulation model
================================

Introduction
-------------

This model realizes a metapopulation model using a directed graph with instances of an ODE-based model as graph nodes. 
Any ODE-based model implemented in MEmilio can be used in the graph nodes, see :doc:`ode`, representing one geographical region. 
Currently, only ODE-based models are supported as nodes, but support for other models will be added soon. 
Mobility between graph nodes is modelled via graph edges. There are two different types of mobility edges implemented 
in MEmilio: ``mio::MobilityEdge`` and ``mio::MobilityEdgeStochastic``. The first one implements a deterministic instant 
mobility approach where subpopulations are exchanged between graph nodes instantaneously in fixed time intervals. 
The latter implements mobility as stochastic jumps between graph nodes using a temporal Gillespie algorithm for simulation.

An overview of nonstandard but often used data types can be found under :doc:`data_types`.


Simulation
-----------

During the simulation, graph nodes are advanced independently from each other according to their local solver until the next time point where mobility via the edges takes place. The concrete exchange of individuals or subpopulations via the graph edges depends on the chosen edge type.

**Instant mobility approach**

The instant mobility approach was introduced by Kühn et al. (2021) and mainly represents daily instant commuting activities between different regions. Using realistic commuting data, such as the absolute number of commuters between regions, we can model realistic exchanges of individuals in the metapopulation context.

In this approach, the commuters are exchanged twice daily in an instantaneous manner, allowing the integration of age-specific mobility restrictions. The implementation also allows for the restriction of mobility activities for certain infection states, such as symptomatic infected individuals. An important feature implemented in Kühn et al. (2022) is the testing of commuters, which enables the detection and potential isolation of asymptomatic or symptomatic infected individuals before they can infect people in other regions.

The instant mobility approach assumes instant exchange between regions without considering travel times. The length of stay is fixed to half a day.
Since all regions exchange commuters simultaneously, this scheme offers great properties for parallelization and has fewer constraints on the integration time step compared to more detailed mobility approaches.

For further details, please refer to:

- |Assessment_of_effective_mitigation|
- |Regional_opening_strategies_with|

**Detailed mobility approach**

The detailed mobility scheme was introduced in Zunker et al. (2024) and further develops the instant approach. The detailed mobility model addresses some of the constraints of the instant approach by adding realistic travel and stay times.

The core idea of this detailed approach is that each region has a second local model, which can be parameterized independently from the existing model, to represent mobility within the region. This allows for more detailed modeling of population movement.

To have a realistic representation of the exchange process, we calculate the centroids of each region, which is represented by a polygon. When modeling travel between two non-adjacent regions (nodes :math:`k` and :math:`l`), the model creates a path connecting their centroids, which may pass through other regions. The predefined travel time between regions is then distributed across all regions along this path.
Commuters move along these paths and pass through various mobility models during their trip. This allows for potential interactions with other commuters during transit. The detailed mobility approach enables modeling of scenarios that cannot be addressed with the instant approach, such as the impact of different mask type usage in public transport. However, this increased realism comes at the cost of greater computational effort and requires more detailed input data.

For further details, please refer to:

- |Novel_travel_time_aware_metapopulation_models|

**Stochastic mobility approach**

In the stochastic mobility approach, transitions of individuals between regions, i.e. graph nodes, are modeled as stochastic jumps. The frequency of individuals transitioning from region :math:`i` to region :math:`j` is determined by a rate :math:`\lambda_{ij}` which can be interpreted as the (expected) fraction of the population in region :math:`i` that commutes to region :math:`j` within one day. In contrast to the instant and detailed mobility approach, the time points and the number of transitions between two regions are stochastic.


Behavior-based ODE models
-------------------------

The graph-based simulation can be combined with the feedback mechanism described in the :doc:`ODE models documentation <ode>`.
This allows for modeling behavioral changes in response to the epidemiological situation not only on a local level but also considering regional and global information. The ``FeedbackGraphSimulation`` extends the local feedback by incorporating a weighted blend of local, regional, and global ICU occupancy to calculate the perceived risk. This blended risk then influences the contact patterns in each node, allowing for a more nuanced and realistic simulation of public health responses across a metapopulation.

The extension and inclusion of regional and global information is based on the following paper:

- |Risk-mediated_dynamic_regulation|


How to: Set up a graph and run a graph simulation
-------------------------------------------------

The graph simulation couples multiple simulation instances representing different geographic regions via a graph-based approach. In this example, a compartmental model based on ODEs (ODE-SECIR) is used for each region. The simulation proceeds by advancing each region independently and then exchanging portions of the population between regions along the graph edges.

The following steps detail how to configure and execute a graph simulation:

1. **Initialize the local compartmental model:**

   First, set up the compartmental model by initializing the parameters that are equal in every region. In this example, we use a single age group and set the necessary epidemiological parameters (e.g. time periods, transmission probabilities).

   .. code-block:: cpp

           const size_t num_groups = 1;
           mio::osecir::Model model(num_groups);
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

   To represent different geographic regions, clone the base model into separate model groups. In this example, two model groups are created. The first group is modified by applying a contact damping to simulate contact restrictions. Additionally, set the populations in all models.

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

           model_group2.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = 10000;

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


   For the stochastic mobility, ``mio::MobilityEdgeStochastic`` has to be used as edge type for the graph. The rates or mobility coefficients can be set as follows:

    .. code-block:: cpp

        mio::Graph<mio::SimulationNode<mio::Simulation<double, mio::osecir::Model<double>>>, mio::MobilityEdgeStochastic> graph;
        graph.add_node(1001, model_group1, t0);
        graph.add_node(1002, model_group2, t0);

        auto transition_rates = mio::MobilityCoefficients(model.populations.numel());

        for (auto compartment = mio::Index<mio::osecir::InfectionState>(0); compartment < model.populations.size<mio::osecir::InfectionState>(); compartment++) {
            auto coeff_idx = model.populations.get_flat_index({mio::AgeGroup(0), compartment});
            transition_rates.get_baseline()[coeff_idx] = 0.01;
        }

        graph.add_edge(0, 1, std::move(transition_rates));
        graph.add_edge(1, 0, std::move(transition_rates));

.. dropdown:: :fa:`gears` Working with large graphs

    When working with very large graphs, i.e. starting from a few thousand edges, it will be faster to not use the standard ``add_edge`` function.
    For this case, we provide a ``GraphBuilder``. There you can add all edges without checking for uniqueness and sorting, thus improving the speed.
    The edges will be sorted when the graph is generated:

    .. code-block:: cpp

        mio::GraphBuilder<mio::SimulationNode<mio::Simulation<double, mio::osecir::Model<double>>>, mio::MobilityEdgeStochastic> builder;
        builder.add_node(1001, model_group1, t0);
        builder.add_node(1002, model_group2, t0);
        builder.add_edge(0, 1, std::move(transition_rates));
        builder.add_edge(1, 0, std::move(transition_rates));
        auto graph = builder.build();


    Usually, there should be no duplicate edges in your input. If this is not certain, the ``GraphBuilder`` can also remove duplicates.
    Here, duplicate means that the start and end node are the same. The parameters in the edge will not be compared. 
    When duplicates are found, only the **last** inserted edge is kept, all previous edges are discarded: 

    .. code-block:: cpp

        mio::GraphBuilder<Int, Int> builder;
        builder.add_node(1001, 100);
        builder.add_node(1002, 100);
        builder.add_edge(0, 1, 100);
        builder.add_edge(1, 0, 100);
        builder.add_edge(0, 1, 200);
        auto graph = builder.build(true);
        // graph contains the edges (0, 1, 100) and (1, 0, 100)


5. **Initialize and Advance the Mobility Simulation:**

   With the graph constructed, initialize the simulation with the starting time and time step. Then, advance the simulation until the final time :math:`t_{max}`.

   .. code-block:: cpp
            
           const auto tmax = 30.;
           const auto dt   = 0.5; // time step for Mobility (daily mobility occurs every second step)
           auto sim = mio::make_mobility_sim(t0, dt, std::move(g));
           sim.advance(tmax);

6. **Access and Display Mobility Results:**

   After the simulation, the mobility results can be extracted from a specific edge. In this example, the results for the edge from node 1 to node 0 are printed.

   .. code-block:: cpp

           auto& edge_1_0 = sim.get_graph().edges()[1];
           auto& results  = edge_1_0.property.get_mobility_results();
           results.print_table({"Commuter INS", "Commuter ISy", "Commuter Total"});

