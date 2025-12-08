#pragma once

#include "boost/filesystem.hpp"

#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "models/ode_secir/model.h"
#include "memilio/io/result_io.h"
#include "memilio/io/epi_data.h"

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <any>
#include <typeindex>
#include <mutex>

template <class FP>
mio::IOResult<Eigen::MatrixX<FP>> read_mobility_plain_template(const std::string& mobility_data_file)
{
    Eigen::MatrixX<ScalarType> matrix = mio::read_mobility_plain<ScalarType>(mobility_data_file).value();
    Eigen::MatrixX<FP> matrix_FP      = matrix.template cast<FP>();
    return mio::success(matrix_FP);
}

class OptimizationModel
{
public:
    OptimizationModel(const boost::filesystem::path& data_directory,

                      const boost::filesystem::path& optimal_control_settings_directory,

                      size_t simulation_days, size_t num_age_groups, double npis_duration, double npis_interval,
                      double npis_base_value);

    boost::filesystem::path data_directory() const;
    boost::filesystem::path optimal_control_settings_directory() const;
    size_t simulation_days() const;
    size_t num_age_groups() const;
    double npis_duration() const;
    double npis_interval() const;
    double npis_base_value() const;

    // Retrieves (or generates) a graph model of the simulation.
    // Uses template-based caching since reading in models can be time-consuming.
    template <typename FP>
    mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>> get_graph_model() const;

private:
    boost::filesystem::path m_data_directory;
    boost::filesystem::path m_optimal_control_settings_directory;
    size_t m_simulation_days;
    size_t m_num_age_groups;
    double m_npis_duration;
    double m_npis_interval;
    double m_npis_base_value;

    // Creates an artificial (synthetic) graph model with preset parameters and population data.
    template <typename FP>
    mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>> create_graph_model() const;

    struct Cache {
        std::unordered_map<std::type_index, std::any> data;
        std::mutex mutex;
    };
    std::shared_ptr<Cache> m_cache;
};

// Retrieves a graph model from the cache or builds a new one if necessary.
template <typename FP>
auto OptimizationModel::get_graph_model() const
    -> mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>>
{
    const auto key = std::type_index(typeid(FP));
    {
        // Check if cached
        std::scoped_lock lock(m_cache->mutex);
        auto it = m_cache->data.find(key);
        if (it != m_cache->data.end()) {
            return std::any_cast<
                mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>>>(it->second);
        }
    }
    auto graph = create_graph_model<FP>();
    {
        // Cache the result
        std::scoped_lock lock(m_cache->mutex);
        m_cache->data.emplace(key, graph);
    }
    return graph;
}

template <typename FP>
mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>>
OptimizationModel::create_graph_model() const
{
    /* --------------------------- */
    /* Define the model parameters */
    /* --------------------------- */
    constexpr size_t num_age_groups = 1;
    mio::osecir::Parameters<FP> model_parameters(num_age_groups);

    model_parameters.template get<mio::osecir::StartDay<FP>>()    = 0.0;
    model_parameters.template get<mio::osecir::Seasonality<FP>>() = 0.2;
    model_parameters.template get<mio::osecir::ICUCapacity<FP>>() = std::numeric_limits<FP>::max();

    auto set_all_groups = [&](auto Tag, FP value) {
        auto& current_parameter = model_parameters.template get<decltype(Tag)>();
        for (mio::AgeGroup i = 0; i < model_parameters.get_num_groups(); i++) {
            current_parameter[i] = value;
        }
    };

    // Read in custom model parameters
    std::string samples_data_path = mio::path_join(m_optimal_control_settings_directory.string(), "samples.json");
    mio::IOResult<Json::Value> result_parameter_list = mio::read_json(samples_data_path);
    if (!result_parameter_list) {
        std::cerr << "Error: mio::read_json failed. IOStatus = " << static_cast<int>(result_parameter_list.error())
                  << std::endl;
        throw std::runtime_error("Failed to read_json(samples_data_path) in create_graph_model()");
    }
    Json::Value parameter_list = result_parameter_list.value();

    // — times (all groups same value)
    set_all_groups(mio::osecir::TimeExposed<FP>{}, parameter_list["t_E"][0].asDouble());
    set_all_groups(mio::osecir::TimeInfectedNoSymptoms<FP>{}, 5.2 - parameter_list["t_E"][0].asDouble());
    set_all_groups(mio::osecir::TimeInfectedSymptoms<FP>{}, parameter_list["t_ISy"][0].asDouble());
    set_all_groups(mio::osecir::TimeInfectedSevere<FP>{}, parameter_list["t_ISev"][0].asDouble());
    set_all_groups(mio::osecir::TimeInfectedCritical<FP>{}, parameter_list["t_Cr"][0].asDouble());
    // — probabilities (all groups same value)
    set_all_groups(mio::osecir::TransmissionProbabilityOnContact<FP>{},
                   parameter_list["transmission_prob"][0].asDouble());
    set_all_groups(mio::osecir::RelativeTransmissionNoSymptoms<FP>{}, 1.0);
    set_all_groups(mio::osecir::RiskOfInfectionFromSymptomatic<FP>{}, 1.0); // Default value
    set_all_groups(mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>{}, 0.0); // Default value
    set_all_groups(mio::osecir::RecoveredPerInfectedNoSymptoms<FP>{}, parameter_list["mu_CR"][0].asDouble());
    set_all_groups(mio::osecir::SeverePerInfectedSymptoms<FP>{}, parameter_list["mu_IH"][0].asDouble());
    set_all_groups(mio::osecir::CriticalPerSevere<FP>{}, parameter_list["mu_HU"][0].asDouble());
    set_all_groups(mio::osecir::DeathsPerCritical<FP>{}, parameter_list["mu_UD"][0].asDouble());

    // --------------------------- //
    // Define the contact matrices //
    // Contact locations for the edges
    enum class ContactLocation
    {
        Work  = 0,
        Count = 1
    };
    constexpr size_t contact_locations_size = static_cast<size_t>(ContactLocation::Count);
    Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum;
    baseline.setConstant(7.95);
    minimum.setZero();
    mio::ContactMatrixGroup<FP> contact_matrices(contact_locations_size, num_age_groups);
    contact_matrices[0].get_baseline() = baseline;
    contact_matrices[0].get_minimum()  = minimum;
    model_parameters.template get<mio::osecir::ContactPatterns<FP>>() =
        mio::UncertainContactMatrix<FP>(std::move(contact_matrices));

    // --------------------------- //
    // DynamicNPIsInfectedSymptoms //
    mio::DynamicNPIs<FP> dynamic_npis;
    dynamic_npis.set_duration(mio::SimulationTime<FP>{m_npis_duration});
    dynamic_npis.set_interval(mio::SimulationTime<FP>{m_npis_interval});
    dynamic_npis.set_base_value(m_npis_base_value);
    model_parameters.template get<mio::osecir::DynamicNPIsInfectedSymptoms<FP>>() = dynamic_npis;

    model_parameters.template get<mio::osecir::DynamicNPIsImplementationDelay<FP>>() = 0.0; // Default value
    // Or is t&t std::numeric_limits<FP>::max() ???
    // model_parameters.template get<mio::osecir::TestAndTraceCapacity<FP>>()        = 0.0; // Set implicity via set_nodes
    model_parameters.template get<mio::osecir::TestAndTraceCapacity<FP>>() =
        std::numeric_limits<FP>::max(); // Set implicity via set_nodes

    model_parameters.template get<mio::osecir::TestAndTraceCapacityMaxRisk<FP>>() = 5.0; // Default value

    model_parameters.check_constraints();
    model_parameters.apply_constraints();

    /* ------------------------- */
    /* Create vector of node ids */
    /* ------------------------- */
    bool is_node_for_county                         = true;
    bool rki_age_groups                             = true;
    mio::IOResult<std::vector<int>> result_node_ids = mio::get_node_ids(
        mio::path_join(m_data_directory.string(), "Germany", "pydata", "county_current_population.json"),
        is_node_for_county, rki_age_groups);

    if (!result_node_ids) {
        std::cerr << "Error: mio::get_node_ids failed. IOStatus = " << static_cast<int>(result_node_ids.error())
                  << std::endl;
        throw std::runtime_error("Failed to get_node_ids() in create_graph_model()");
    }
    std::vector<int> node_ids = result_node_ids.value();
    size_t num_graph_nodes    = node_ids.size();

    /* --------------------------- */
    /* Read in the population data */
    mio::IOResult<std::vector<mio::SimulationResult<ScalarType>>> result_graph_population_data =
        mio::read_result<ScalarType>(mio::path_join(m_optimal_control_settings_directory.string(), "results_run0.h5"));
    if (!result_graph_population_data) {
        std::cerr << "Error: mio::set_nodes failed. IOStatus = "
                  << static_cast<int>(result_graph_population_data.error()) << std::endl;
        throw std::runtime_error("Failed to read_result(population_data_path) in create_graph_model()");
    }
    std::vector<mio::SimulationResult<ScalarType>> graph_population_data = result_graph_population_data.value();

    /* -------------------- */
    /* Set graph population */
    std::vector<mio::Populations<FP, mio::AgeGroup, mio::osecir::InfectionState>> graph_population(
        num_graph_nodes, mio::Populations<FP, mio::AgeGroup, mio::osecir::InfectionState>(
                             {model_parameters.get_num_groups(), mio::osecir::InfectionState::Count}));

    size_t num_infection_states = static_cast<size_t>(mio::osecir::InfectionState::Count);
    for (size_t node_idx = 0; node_idx < num_graph_nodes; ++node_idx) {
        const auto& sim_result                                 = graph_population_data[node_idx];
        const mio::TimeSeries<ScalarType>& groups              = sim_result.get_groups();
        Eigen::Ref<const Eigen::VectorX<ScalarType>> last_step = groups.get_last_value();
        for (size_t group_idx = 0; group_idx < num_age_groups; ++group_idx) {
            for (size_t state_idx = 0; state_idx < num_infection_states; ++state_idx) {
                mio::AgeGroup age(group_idx);
                mio::osecir::InfectionState state        = static_cast<mio::osecir::InfectionState>(state_idx);
                size_t idx                               = group_idx * num_infection_states + state_idx;
                graph_population[node_idx][{age, state}] = static_cast<FP>(last_step[idx]);
            }
        }
    }

    /* ------------------------- */
    /* Construct parameter graph */
    mio::Graph<mio::osecir::Model<FP>, mio::MobilityParameters<FP>> parameter_graph;
    for (size_t node_idx = 0; node_idx < num_graph_nodes; ++node_idx) {
        mio::osecir::Model<FP> node_model(graph_population[node_idx], model_parameters);
        parameter_graph.add_node(node_ids[node_idx], node_model);
    }

    /* -------------------- */
    /* Sets the graph edges */
    /* -------------------- */

    // Compartments that commute
    std::initializer_list<mio::osecir::InfectionState> mobile_compartments = {
        mio::osecir::InfectionState::Susceptible, mio::osecir::InfectionState::Exposed,
        mio::osecir::InfectionState::InfectedNoSymptoms, mio::osecir::InfectionState::InfectedSymptoms,
        mio::osecir::InfectionState::Recovered};
    // File that contains the commuting matrix of size num_graph_nodes x num_graph_nodes
    std::string mobility_data_file =
        mio::path_join(m_data_directory.string(), "Germany", "mobility", "commuter_mobility_2022.txt");
    // Vector with a commuting weight for every AgeGroup.
    std::vector<FP> commuting_weights = {1.0};
    // Vector of vectors with indices of the compartments that should be saved on the edges.
    std::vector<std::vector<size_t>> indices_of_saved_edges = {};

    // states statt county
    mio::IOResult<void> result_set_edges =
        mio::set_edges<FP, ContactLocation, mio::osecir::Model<FP>, mio::MobilityParameters<FP>,
                       mio::MobilityCoefficientGroup<FP>, mio::osecir::InfectionState,
                       decltype(read_mobility_plain_template<FP>)>(
            mobility_data_file, parameter_graph, mobile_compartments, contact_locations_size,
            read_mobility_plain_template<FP>, commuting_weights, indices_of_saved_edges);

    if (!result_set_edges) {
        std::cerr << "Error: mio::set_edges failed. IOStatus = " << static_cast<int>(result_set_edges.error())
                  << std::endl;
        throw std::runtime_error("Failed to set edges in create_graph_model()");
    }

    mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>> graph_simulation;

    // Add nodes (regional SECIR models) to the simulation graph
    for (auto& node : parameter_graph.nodes()) {
        if (14628 == node.id)
            graph_simulation.add_node(node.id, node.property);
    }
    // Add edges (mobility links) to the simulation graph
    for (auto& edge : parameter_graph.edges()) {
        /* We don't add any edges, because otherwise the simulation gets stuck in an infinite loop. */
        // graph_simulation.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
    }

    return graph_simulation;
}
