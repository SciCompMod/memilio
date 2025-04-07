Graph-based agent-based model
================================

Introduction
-------------

This model realizes multiple instances of the mobility-based agent-based model ``abm::Model`` (see :doc:`<mobility_based_abm>` for documentation) as nodes in a directed graph. One local model represents a geographical region. The regions are connected by the graph edges. Mobility within one node and via the graph edges follows the same mobility rules that can be handed as argument to ``mio::GraphABModel``. Therefore this graph-based agent-based model can be reduced to a single mobility-based agent-based if simulation time steps within the whole graph, i.e. the step size of each node and the step size of the edge exchange, are equal.

Simulation
-----------

The simulation runs in discrete time steps with two different step sizes. The first one is the node-internal step size with which each node advances independently its inherent model instance. The second step size is used for the regular exchange of agents between models in different nodes via the corresponding edge. The simulation therefore advances all nodes until the next exchange time step and then performs the exchange via the graph edges.

How to: Set up and run a simulation of the graph-based agent-based model
------------------------------------------------------------------------

First, the models (corresponding to the graph nodes) have to be created and initialized with the number of age groups and a unique node id:

.. code-block:: cpp

    //This example considers two age groups representing children and adults
    size_t num_age_groups         = 2;
    const auto age_group_children = mio::AgeGroup(0);
    const auto age_group_adults   = mio::AgeGroup(1);

    //Initialize the models with id 0 and 1
    auto model1 = mio::GraphABModel(num_age_groups, 0);
    auto model2 = mio::GraphABModel(num_age_groups, 1);

For all models/nodes, the infection parameters have to be set. If this is not done, the default values are used. For information on what parameters the mobility-based agent-based model ``abm::Model`` has and how to set them, see :doc:`<mobility_based_abm>`. Below, only the age groups that go to school and work are set:

.. code-block:: cpp

    //Age group 0 goes to school and age group 1 goes to work
    model1.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_children] = true;
    model1.parameters.get<mio::abm::AgeGroupGotoWork>()[age_group_adults]     = true;
    model2.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_children] = true;
    model2.parameters.get<mio::abm::AgeGroupGotoWork>()[age_group_adults]     = true;

In the example below, only one home location is added for each model. Model 1 has one adult and one child assigned to that home location and model 2 has two adults assigned to the home location:

.. code-block:: cpp

    //Create home location and persons for node 1
    auto home_id1    = model1.add_location(mio::abm::LocationType::Home);
    auto& home1      = model1.get_location(home_id1);
    auto child1_id   = model1.add_person(home_id1, mio::AgeGroup(0));
    auto adult1_id   = model1.add_person(home_id1, mio::AgeGroup(1));
    auto& child1     = model1.get_person(child1_id);
    auto& adult1     = model1.get_person(adult1_id);
    //Assign home1
    child1.set_assigned_location(home1.get_type(), home1.get_id(), home1.get_model_id());
    adult1.set_assigned_location(home1.get_type(), home1.get_id(), home1.get_model_id());

    //Create home location and persons for node 2
    auto home_id2    = model2.add_location(mio::abm::LocationType::Home);
    auto& home2      = model2.get_location(home_id2);
    auto adult2_id   = model2.add_person(home_id2, mio::AgeGroup(1));
    auto adult3_id   = model2.add_person(home_id2, mio::AgeGroup(1));
    auto& adult2     = model2.get_person(adult2_id);
    auto& adult3     = model2.get_person(adult3_id);
    //Assign home2
    adult2.set_assigned_location(home2.get_type(), home2.get_id(), home2.get_model_id());
    adult3.set_assigned_location(home2.get_type(), home2.get_id(), home2.get_model_id());

Next, for all models in the graph, locations have to be added. This can be done as follows:

.. code-block:: cpp

    //Add an event and a shop to both models
    auto event1 = model1.add_location(mio::abm::LocationType::SocialEvent);
    auto event2 = model2.add_location(mio::abm::LocationType::SocialEvent);
    auto shop1  = model1.add_location(mio::abm::LocationType::BasicsShop);
    auto shop2  = model2.add_location(mio::abm::LocationType::BasicsShop);
    //Add a school, a hospital and an ICU only to model 1
    auto school   = model1.add_location(mio::abm::LocationType::School);
    auto hospital = model1.add_location(mio::abm::LocationType::Hospital);
    auto icu      = model1.add_location(mio::abm::LocationType::ICU);
    //Add a work place only to model2
    auto work = model2.add_location(mio::abm::LocationType::Work);

Assigning infection states and locations to persons in all models can be done via

.. code-block:: cpp

    //Simulation start date
    auto start_date = mio::abm::TimePoint(0);

    //Add infection to persons in home1
    auto rng_child1 = mio::abm::PersonalRandomNumberGenerator(child1);
    child1.add_new_infection(mio::abm::Infection(rng_child1, mio::abm::VirusVariant::Wildtype, child1.get_age(),
                                                         model1.parameters, start_date, mio::abm::InfectionState::InfectedNoSymptoms));
    auto rng_adult1 = mio::abm::PersonalRandomNumberGenerator(adult1);
    adult1.add_new_infection(mio::abm::Infection(rng_adult1, mio::abm::VirusVariant::Wildtype, adult1.get_age(),
                                                         model1.parameters, start_date, mio::abm::InfectionState::Exposed));

    //Assign Event, Shop, Hospital and ICU to all persons, school only to the child and work to the adults
    //Event
    child1.set_assigned_location(mio::abm::LocationType::SocialEvent, event1, model1.get_id());
    adult1.set_assigned_location(mio::abm::LocationType::SocialEvent, event1, model1.get_id());
    adult2.set_assigned_location(mio::abm::LocationType::SocialEvent, event2, model2.get_id());
    adult3.set_assigned_location(mio::abm::LocationType::SocialEvent, event2, model2.get_id());
    //Shop
    child1.set_assigned_location(mio::abm::LocationType::BasicsShop, shop1, model1.get_id());
    adult1.set_assigned_location(mio::abm::LocationType::BasicsShop, shop1, model1.get_id());
    adult2.set_assigned_location(mio::abm::LocationType::BasicsShop, shop2, model2.get_id());
    adult3.set_assigned_location(mio::abm::LocationType::BasicsShop, shop2, model2.get_id());
    //Hospital
    child1.set_assigned_location(mio::abm::LocationType::Hospital, hospital, model1.get_id());
    adult1.set_assigned_location(mio::abm::LocationType::Hospital, hospital, model1.get_id());
    adult2.set_assigned_location(mio::abm::LocationType::Hospital, hospital, model1.get_id());
    adult3.set_assigned_location(mio::abm::LocationType::Hospital, hospital, model1.get_id());
    //ICU
    child1.set_assigned_location(mio::abm::LocationType::ICU, icu, model1.get_id());
    adult1.set_assigned_location(mio::abm::LocationType::ICU, icu, model1.get_id());
    adult2.set_assigned_location(mio::abm::LocationType::ICU, icu, model1.get_id());
    adult3.set_assigned_location(mio::abm::LocationType::ICU, icu, model1.get_id());
    //School
    child1.set_assigned_location(mio::abm::LocationType::School, school, model1.get_id());
    //Work
    adult1.set_assigned_location(mio::abm::LocationType::Work, work, model2.get_id());
    adult2.set_assigned_location(mio::abm::LocationType::Work, work, model2.get_id());
    adult3.set_assigned_location(mio::abm::LocationType::Work, work, model2.get_id());

For initializing the graph nodes and edges a ``mio::Graph`` is created which gets ``mio::ABMSimulationNode`` and ``mio::ABMMobilityEdge`` as templates. Additionally, every node needs a ``mio::History`` object to log its results during the simulation. See :doc:`<io>` for information on how to use ``mio::History``. Below, ``mio::abm::LogInfectionState`` is used as logger.

.. code-block:: cpp

    //Define history type
    using HistoryType = mio::History<mio::DataWriterToMemory, mio::abm::LogInfectionState>;
    //Create graph and add nodes and edges
    mio::Graph<mio::ABMSimulationNode<HistoryType>, mio::ABMMobilityEdge<HistoryType>> graph;
    graph.add_node(model1.get_id(), HistoryType{}, start_date, std::move(model1));
    graph.add_node(model2.get_id(), HistoryType{}, start_date, std::move(model2));
    graph.add_edge(model1.get_id(), model2.get_id());
    graph.add_edge(model2.get_id(), model1.get_id());

To simulate the model from `start_date` to `end_date` with given graph step size `exchange_time_span`, a GraphSimulation has to be created. The step size is used to regularly exchange agents via the graph edges. Advancing the simulation until `end_date` is done as follows:

.. code-block:: cpp

    //Simulation end date
    auto end_date   = start_date + mio::abm::days(30);

    //Agents are exchanged via the graph edges every 12 hours
    auto exchange_time_span = mio::abm::hours(12);
    //Create GraphSimulation and advance until end_date
    auto sim                = mio::make_abm_graph_sim<HistoryType>(start_date, exchange_time_span, std::move(graph));
    sim.advance(end_date);
