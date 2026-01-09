#pragma once

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstddef>
#include <functional>

#include "control_parameters.h"

#include "models/ode_secir/model.h"

#include "../optimization_settings/secir_optimization.h"
#include "../helpers/make_time_grid.h"

enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Count
};
enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Holidays,
    Count
};
enum class ContactLocation
{
    Home,
    School,
    Work,
    Other,
    Count
};

template <typename FP>
mio::DampingSampling<FP> make_school_closure_damping(mio::SimulationTime<FP> time, FP value, size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::Main));
    const mio::DampingType type(static_cast<size_t>(Intervention::SchoolClosure));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::School)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, 1.0);

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}

template <typename FP>
mio::DampingSampling<FP> make_home_office_damping(mio::SimulationTime<FP> time, FP value, size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::Main));
    const mio::DampingType type(static_cast<size_t>(Intervention::HomeOffice));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::Work)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, 1.0);

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}

template <typename FP>
mio::DampingSampling<FP> make_physical_distancing_school_damping(mio::SimulationTime<FP> time, FP value,
                                                                 size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::PhysicalDistanceAndMasks));
    const mio::DampingType type(static_cast<size_t>(Intervention::PhysicalDistanceAndMasks));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::School)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, 1.0);

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}

template <typename FP>
mio::DampingSampling<FP> make_physical_distancing_work_damping(mio::SimulationTime<FP> time, FP value,
                                                               size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::PhysicalDistanceAndMasks));
    const mio::DampingType type(static_cast<size_t>(Intervention::PhysicalDistanceAndMasks));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::Work)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, 1.0);

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}

template <typename FP>
mio::DampingSampling<FP> make_physical_distancing_other_damping(mio::SimulationTime<FP> time, FP value,
                                                                size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::PhysicalDistanceAndMasks));
    const mio::DampingType type(static_cast<size_t>(Intervention::PhysicalDistanceAndMasks));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::Other)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, 1.0);

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}

template <typename FP>
void set_dynamic_NPIs(
    const SecirOptimization& settings,
    mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>>& graph_model,
    const std::vector<std::pair<FP, FP>>& dynamic_NPI_values)
{
    size_t num_graph_nodes = graph_model.nodes().size();

    for (size_t node_index = 0; node_index < num_graph_nodes; node_index++) {
        mio::SimulationNode<FP, mio::osecir::Simulation<FP>>& node = graph_model.nodes()[node_index].property;
        mio::osecir::Simulation<FP>& node_simulation               = node.get_simulation();
        mio::osecir::Model<FP>& node_model                         = node_simulation.get_model();

        mio::SimulationTime<FP> start_time(0.0);
        size_t num_age_groups = static_cast<size_t>(node_model.parameters.get_num_groups());

        auto effectiveness = [&](const std::string& name) {
            size_t control_index = static_cast<size_t>(string_to_control(name));
            return settings.dynamic_NPI_parameters()[control_index].effectiveness();
        };
        auto strength = [&](const std::string& name) {
            size_t control_index = static_cast<size_t>(string_to_control(name));
            return dynamic_NPI_values[control_index].second;
        };
        auto threshold = [&](const std::string& name) {
            size_t control_index = static_cast<size_t>(string_to_control(name));
            return dynamic_NPI_values[control_index].first;
        };

        auto damping_school_closure = make_school_closure_damping<FP>(
            start_time, effectiveness("SchoolClosure") * strength("SchoolClosure"), num_age_groups);

        auto damping_home_office = make_home_office_damping<FP>(
            start_time, effectiveness("HomeOffice") * strength("HomeOffice"), num_age_groups);

        auto damping_physical_distancing_school = make_physical_distancing_school_damping<FP>(
            start_time, effectiveness("PhysicalDistancingSchool") * strength("PhysicalDistancingSchool"),
            num_age_groups);

        auto damping_physical_distancing_work = make_physical_distancing_work_damping<FP>(
            start_time, effectiveness("PhysicalDistancingWork") * strength("PhysicalDistancingWork"), num_age_groups);

        auto damping_physical_distancing_other = make_physical_distancing_other_damping<FP>(
            start_time, effectiveness("PhysicalDistancingOther") * strength("PhysicalDistancingOther"), num_age_groups);

        // We asume that there are no other dynamic_npis stored before!
        mio::DynamicNPIs<FP>& dynamic_npis =
            node_model.parameters.template get<mio::osecir::DynamicNPIsInfectedSymptoms<FP>>();

        // Collect all dampings with their thresholds
        std::vector<std::pair<FP, mio::DampingSampling<FP>>> dampings_with_thresholds = {
            {threshold("SchoolClosure"), damping_school_closure},
            {threshold("HomeOffice"), damping_home_office},
            {threshold("PhysicalDistancingSchool"), damping_physical_distancing_school},
            {threshold("PhysicalDistancingWork"), damping_physical_distancing_work},
            {threshold("PhysicalDistancingOther"), damping_physical_distancing_other}};

        // Sort dampings by threshold ascending
        std::sort(dampings_with_thresholds.begin(), dampings_with_thresholds.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

        // Map thresholds to cumulative dampings
        std::map<FP, std::vector<mio::DampingSampling<FP>>> threshold_to_dampings;

        // Accumulate dampings
        for (const auto& [thresh, damping] : dampings_with_thresholds) {
            // Add this damping to all thresholds >= current damping
            for (auto& [key, vec] : threshold_to_dampings) {
                if (key >= thresh) {
                    vec.push_back(damping);
                }
            }
            // Also ensure this threshold exists in map
            threshold_to_dampings[thresh].push_back(damping);
        }

        // Now thresholds in ascending order
        std::vector<FP> thresholds;
        for (const auto& [thresh, _] : threshold_to_dampings) {
            thresholds.push_back(thresh);
        }
        std::sort(thresholds.begin(), thresholds.end());

        // Set thresholds in dynamic_npis
        std::vector<mio::DampingSampling<FP>> cumulative_dampings;
        for (FP thresh : thresholds) {
            // Start with previous cumulative dampings
            cumulative_dampings = {};
            for (const auto& [d_thresh, dampings] : dampings_with_thresholds) {
                if (d_thresh <= thresh) {
                    cumulative_dampings.push_back(dampings);
                }
            }
            dynamic_npis.set_threshold(thresh, cumulative_dampings);
        }
    }
}

// template <typename FP>
// void set_commuter_testing(
//     const SecirOptimization& settings,
//     mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>>& graph_model,
//     const std::vector<FP>& commuter_testing_values)
// {
//     size_t num_graph_nodes = graph_model.nodes().size();

//     for (size_t node_index = 0; node_index < num_graph_nodes; node_index++) {
//         mio::SimulationNode<FP, mio::osecir::Simulation<FP>>& node = graph_model.nodes()[node_index].property;
//         mio::osecir::Simulation<FP>& node_simulation               = node.get_simulation();
//         mio::osecir::Model<FP>& node_model                         = node_simulation.get_model();

//         node_model.parameters.get_commuter_nondetection() = commuter_testing_values[node_index];
//     }
// }

template <typename FP>
void set_commuter_testing(
    const SecirOptimization& settings,
    mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>>& graph_model,
    const std::vector<FP>& commuter_testing_values)
{
    size_t num_graph_nodes = graph_model.nodes().size();

    for (size_t node_index = 0; node_index < num_graph_nodes; ++node_index) {
        auto& node            = graph_model.nodes()[node_index].property;
        auto& node_simulation = node.get_simulation();
        auto& node_model      = node_simulation.get_model();

        // Compute Bundesland index from node_model.id
        size_t state_index = static_cast<size_t>(std::floor(graph_model.nodes()[node_index].id / 1000.0)) - 1;

        // Safety check: ensure state_index is within bounds [0, 15]
        if (state_index >= commuter_testing_values.size()) {
            throw std::out_of_range("Invalid state index derived from node_model.id: " + std::to_string(state_index));
        }

        // Apply state-level commuter testing value to this node
        node_model.parameters.get_commuter_nondetection() = commuter_testing_values[state_index];
    }
}
