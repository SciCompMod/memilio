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
mio::IOResult<Eigen::MatrixX<FP>> read_mobility_plain_FP(const std::string& mobility_data_file)
{
    Eigen::MatrixX<ScalarType> matrix = mio::read_mobility_plain(mobility_data_file).value();
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
    constexpr size_t num_age_groups = 1;
    mio::osecir::Model<FP> model(num_age_groups);
    auto& params = model.parameters;

    params.template get<mio::osecir::StartDay<FP>>()    = 334.0;
    params.template get<mio::osecir::Seasonality<FP>>() = 0.2;
    /* ICU Capactiy default: std::numeric_limits<FP>::max() */
    // params.template get<mio::osecir::ICUCapacity<FP>>() = 100;

    auto set_all_groups = [&](auto Tag, FP value) {
        auto& age_group_params = params.template get<decltype(Tag)>();
        for (size_t age_group = 0; age_group < num_age_groups; age_group++) {
            age_group_params[mio::AgeGroup(age_group)] = value;
        }
    };

    // -------------------------------- //
    // Read sample parameters from JSON //
    // -------------------------------- //
    std::string samples_data_path = mio::path_join(m_optimal_control_settings_directory.string(), "samples.json");
    mio::IOResult<Json::Value> result_parameter_list = mio::read_json(samples_data_path);
    if (!result_parameter_list) {
        std::cerr << "Error: mio::read_json failed. IOStatus = " << static_cast<int>(result_parameter_list.error())
                  << std::endl;
        throw std::runtime_error("Failed to read_json(samples_data_path) in create_graph_model()");
    }
    Json::Value parameter_list = result_parameter_list.value();

    // ---------------- //
    // Transition Times //
    set_all_groups(mio::osecir::TimeExposed<FP>{}, parameter_list["t_E"][0].asDouble());
    set_all_groups(mio::osecir::TimeInfectedNoSymptoms<FP>{}, 5.2 - parameter_list["t_E"][0].asDouble());
    set_all_groups(mio::osecir::TimeInfectedSymptoms<FP>{}, parameter_list["t_ISy"][0].asDouble());
    set_all_groups(mio::osecir::TimeInfectedSevere<FP>{}, parameter_list["t_ISev"][0].asDouble());
    set_all_groups(mio::osecir::TimeInfectedCritical<FP>{}, parameter_list["t_Cr"][0].asDouble());
    // ---------------------- //
    // Transition Probability //
    /* TransmissionProbabilityOnContact default: 1.0 */
    set_all_groups(mio::osecir::TransmissionProbabilityOnContact<FP>{},
                   parameter_list["transmission_prob"][0].asDouble());
    /* RelativeTransmissionNoSymptoms default: 1.0 */
    set_all_groups(mio::osecir::RelativeTransmissionNoSymptoms<FP>{}, 1.0);
    /* RiskOfInfectionFromSymptomatic default: 1.0 */
    // set_all_groups(mio::osecir::RiskOfInfectionFromSymptomatic<FP>{}, 0.2);
    /* MaxRiskOfInfectionFromSymptomatic default: 0.0 */
    // set_all_groups(mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>{}, 0.45);
    /* RecoveredPerInfectedNoSymptoms default: 0.0 */
    set_all_groups(mio::osecir::RecoveredPerInfectedNoSymptoms<FP>{}, parameter_list["mu_CR"][0].asDouble());
    /* RecoveredPerInfectedNoSymptoms default: 0.0 */
    set_all_groups(mio::osecir::SeverePerInfectedSymptoms<FP>{}, parameter_list["mu_IH"][0].asDouble());
    /* RecoveredPerInfectedNoSymptoms default: 0.0 */
    set_all_groups(mio::osecir::CriticalPerSevere<FP>{}, parameter_list["mu_HU"][0].asDouble());
    /* RecoveredPerInfectedNoSymptoms default: 0.0 */
    set_all_groups(mio::osecir::DeathsPerCritical<FP>{}, parameter_list["mu_UD"][0].asDouble());

    // --------------------------- //
    // Define the contact matrices //
    constexpr size_t num_contact_locations = 1;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline;
    Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum;
    baseline.setConstant(7.95);
    minimum.setZero();
    mio::ContactMatrixGroup<FP> contact_matrices(num_contact_locations, num_age_groups);
    contact_matrices[0].get_baseline() = baseline;
    contact_matrices[0].get_minimum()  = minimum;
    params.template get<mio::osecir::ContactPatterns<FP>>() =
        mio::UncertainContactMatrix<FP>(std::move(contact_matrices));

    // DynamicNPIsInfectedSymptoms
    mio::DynamicNPIs<FP> dynamic_npis;
    dynamic_npis.set_duration(mio::SimulationTime<FP>{m_npis_duration});
    dynamic_npis.set_interval(mio::SimulationTime<FP>{m_npis_interval});
    dynamic_npis.set_base_value(m_npis_base_value);
    params.template get<mio::osecir::DynamicNPIsInfectedSymptoms<FP>>() = dynamic_npis;

    // Percentage of infected commuters that are not detected.
    params.get_commuter_nondetection() = 0.0;
    // Time in simulation before which no infected commuters are detected.
    params.get_start_commuter_detection() = 0.0;
    // Time in simulation after which no infected commuters are detected.
    params.get_end_commuter_detection() = 0.0;
    // Time in simulation after which no dynamic NPIs are applied.
    params.get_end_dynamic_npis() = std::numeric_limits<FP>::max();

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

    size_t num_graph_nodes = result_node_ids.value().size();

    /* ---------------------------- */
    /* Create SECIR parameter graph */
    /* ---------------------------- */
    mio::Graph<mio::osecir::Model<FP>, mio::MobilityParameters<FP>> parameter_graph(result_node_ids.value(), model);

    std::string population_data_path = mio::path_join(m_optimal_control_settings_directory.string(), "results_run0.h5");
    mio::IOResult<std::vector<mio::SimulationResult>> result_read_population_data =
        mio::read_result(population_data_path);
    if (!result_read_population_data) {
        std::cerr << "Error: mio::set_nodes failed. IOStatus = "
                  << static_cast<int>(result_read_population_data.error()) << std::endl;
        throw std::runtime_error("Failed to read_result(population_data_path) in create_graph_model()");
    }
    std::cout << "Population data loaded for " << num_graph_nodes << " nodes." << std::endl;

    // Set population data
    size_t num_infection_states = static_cast<size_t>(mio::osecir::InfectionState::Count);
    for (size_t node_idx = 0; node_idx < num_graph_nodes; node_idx++) {
        auto& node                                             = parameter_graph.nodes()[node_idx];
        mio::SimulationResult simulation_result                = result_read_population_data.value()[node_idx];
        const mio::TimeSeries<ScalarType>& group_data          = simulation_result.get_groups();
        Eigen::Ref<const Eigen::VectorX<ScalarType>> last_step = group_data.get_last_value();

        for (size_t state_idx = 0; state_idx < num_infection_states; state_idx++) {
            mio::osecir::InfectionState infection_state = static_cast<mio::osecir::InfectionState>(state_idx);
            node.property.populations[{mio::AgeGroup(0), infection_state}] = static_cast<FP>(last_step[state_idx]);
        }
    }

    // Set Dampings from day 45
    for (size_t node_idx = 0; node_idx < parameter_graph.nodes().size(); ++node_idx) {
        auto& node                                  = parameter_graph.nodes()[node_idx];
        auto& node_params                           = node.property.parameters;
        int state                                   = static_cast<int>(mio::regions::get_state_id(node.id)) - 1;
        double damping_value                        = parameter_list["damping_values"][state][2].asDouble();
        mio::ContactMatrixGroup<FP>& contact_matrix = node_params.template get<mio::osecir::ContactPatterns<FP>>();
        contact_matrix[0].add_damping(Eigen::MatrixX<FP>::Constant(1, 1, damping_value), mio::SimulationTime<FP>(0.0));
    }

    for (size_t node_idx = 0; node_idx < parameter_graph.nodes().size(); ++node_idx) {
        auto& node        = parameter_graph.nodes()[node_idx];
        auto& node_params = node.property.parameters;
        /* DynamicNPIsImplementationDelay default: 0.0 */
        node_params.template get<mio::osecir::DynamicNPIsImplementationDelay<FP>>() = 0.0;
        /* TestAndTraceCapacity default: std::numeric_limits<FP>::max() */
        FP test_and_trace_capacity_factor = 0.0 / 100000.0;
        node_params.template get<mio::osecir::TestAndTraceCapacity<FP>>() =
            node.property.populations.get_total() * test_and_trace_capacity_factor;
        /* TestAndTraceCapacityMaxRisk default: 5.0 */
        node_params.template get<mio::osecir::TestAndTraceCapacityMaxRisk<FP>>() = 5.0;
    }

    /* -------------------- */
    /* Sets the graph edges */
    /* -------------------- */
    enum class ContactLocation
    {
        Work  = 1,
        Count = 2,
    };

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

    size_t contact_locations_size = static_cast<size_t>(ContactLocation::Count);

    // states statt county
    mio::IOResult<void> result_set_edges =
        mio::set_edges<FP, ContactLocation, mio::osecir::Model<FP>, mio::MobilityParameters<FP>,
                       mio::MobilityCoefficientGroup<FP>, mio::osecir::InfectionState,
                       decltype(read_mobility_plain_FP<FP>)>(mobility_data_file, parameter_graph, mobile_compartments,
                                                             contact_locations_size, read_mobility_plain_FP<FP>,
                                                             commuting_weights, indices_of_saved_edges);

    if (!result_set_edges) {
        std::cerr << "Error: mio::set_edges failed. IOStatus = " << static_cast<int>(result_set_edges.error())
                  << std::endl;
        throw std::runtime_error("Failed to set edges in create_graph_model()");
    }

    // ---------------------- //
    // Define the graph model //
    mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>> simulation_graph;

    // Add nodes (regional SECIR models) to the simulation graph
    for (auto& node : parameter_graph.nodes()) {
        simulation_graph.add_node(node.id, node.property);
    }
    // Add edges (mobility links) to the simulation graph
    for (auto& edge : parameter_graph.edges()) {
        simulation_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
    }

    return simulation_graph;
}
