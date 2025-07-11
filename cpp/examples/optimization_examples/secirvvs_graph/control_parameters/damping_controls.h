#pragma once

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstddef>
#include <functional>

#include "control_parameters.h"

#include "models/ode_secirvvs/model.h"

#include "../optimization_settings/secirvvs_optimization.h"
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
void set_control_dampings(const SecirvvsOptimization& settings,
                          mio::Graph<mio::SimulationNode<FP, mio::Simulation<FP, mio::osecirvvs::Model<FP>>>,
                                     mio::MobilityEdge<FP>>& graph_model,
                          const std::vector<std::vector<FP>>& control_parameters)
{
    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    size_t num_graph_nodes = graph_model.nodes().size();

    for (size_t node_index = 0; node_index < num_graph_nodes; node_index++) {
        mio::SimulationNode<FP, mio::Simulation<FP, mio::osecirvvs::Model<FP>>>& node =
            graph_model.nodes()[node_index].property;
        mio::Simulation<FP, mio::osecirvvs::Model<FP>>& node_simulation = node.get_simulation();
        mio::osecirvvs::Model<FP>& node_model                           = node_simulation.get_model();

        mio::UncertainContactMatrix<FP>& contacts =
            node_model.parameters.template get<mio::osecirvvs::ContactPatterns<FP>>();
        std::vector<mio::DampingSampling<FP>>& contact_dampings = contacts.get_dampings();

        for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {
            mio::SimulationTime<FP> time(time_steps[control_interval * settings.pc_resolution()]);

            auto param_at = [&](const std::string& name) {
                size_t control_index = static_cast<size_t>(string_to_control(name));
                return control_parameters[node_index]
                                         [control_index + control_interval * settings.num_control_parameters()];
            };

            auto effectiveness = [&](const std::string& name) {
                size_t control_index = static_cast<size_t>(string_to_control(name));
                return settings.control_parameters()[control_index].effectiveness();
            };

            size_t num_age_groups = static_cast<size_t>(node_model.parameters.get_num_groups());

            auto damping_school_closure = make_school_closure_damping<FP>(
                time, effectiveness("SchoolClosure") * param_at("SchoolClosure"), num_age_groups);

            auto damping_home_office = make_home_office_damping<FP>(
                time, effectiveness("HomeOffice") * param_at("HomeOffice"), num_age_groups);

            auto damping_physical_distancing_school = make_physical_distancing_school_damping<FP>(
                time, effectiveness("PhysicalDistancingSchool") * param_at("PhysicalDistancingSchool"), num_age_groups);

            auto damping_physical_distancing_work = make_physical_distancing_work_damping<FP>(
                time, effectiveness("PhysicalDistancingWork") * param_at("PhysicalDistancingWork"), num_age_groups);

            auto damping_physical_distancing_other = make_physical_distancing_other_damping<FP>(
                time, effectiveness("PhysicalDistancingOther") * param_at("PhysicalDistancingOther"), num_age_groups);

            contact_dampings.push_back(damping_school_closure);
            contact_dampings.push_back(damping_home_office);
            contact_dampings.push_back(damping_physical_distancing_school);
            contact_dampings.push_back(damping_physical_distancing_work);
            contact_dampings.push_back(damping_physical_distancing_other);
        }
        contacts.make_matrix();
    }
}
