/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Henrik Zunker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef METAPOPULATION_MOBILITY_DETAILED_H
#define METAPOPULATION_MOBILITY_DETAILED_H

#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "memilio/math/eigen_util.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/math/euler.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/dynamic_npis.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/date.h"
#include "memilio/mobility/graph.h"
#include "memilio/io/mobility_io.h"

#include "boost/filesystem.hpp"

#include <cassert>

namespace mio
{

template <class NodePropertyT>
struct NodeDetailed : Node<NodePropertyT> {
    template <class... Args>
    NodeDetailed(int node_id, Args&&... args)
        : Node<NodePropertyT>(node_id, std::forward<Args>(args)...)
        , stay_duration(0.5)
        , mobility(std::forward<Args>(args)...)
    {
    }

    template <typename Model>
    NodeDetailed(int node_id, double duration, Model property_arg, Model mobility_arg, double m_t0,
                 double m_dt_integration)
        : Node<NodePropertyT>(node_id, property_arg, m_t0, m_dt_integration)
        , stay_duration(duration)
        , mobility(mobility_arg, m_t0, m_dt_integration)
    {
    }

    template <class... Args>
    NodeDetailed(int node_id, double duration, Args&&... args)
        : Node<NodePropertyT>(node_id, std::forward<Args>(args)...)
        , stay_duration(duration)
        , mobility(std::forward<Args>(args)...)
    {
    }

    NodeDetailed(int node_id, double duration, NodePropertyT property_arg, NodePropertyT mobility_pt_arg)
        : Node<NodePropertyT>(node_id, property_arg)
        , stay_duration(duration)
        , mobility(mobility_pt_arg)
    {
    }

    double stay_duration;
    NodePropertyT mobility;
};

template <class EdgePropertyT>
struct EdgeDetailed : Edge<EdgePropertyT> {
    template <class... Args>
    EdgeDetailed(size_t start, size_t end, Args&&... args)
        : Edge<EdgePropertyT>(start, end, std::forward<Args>(args)...)
        , traveltime(0.)
        , path{static_cast<int>(start), static_cast<int>(end)}
    {
    }

    template <class... Args>
    EdgeDetailed(size_t start, size_t end, double t_travel, Args&&... args)
        : Edge<EdgePropertyT>(start, end, std::forward<Args>(args)...)
        , traveltime(t_travel)
        , path{static_cast<int>(start), static_cast<int>(end)}
    {
    }

    template <class... Args>
    EdgeDetailed(size_t start, size_t end, double t_travel, std::vector<int> path_mobility, Args&&... args)
        : Edge<EdgePropertyT>(start, end, std::forward<Args>(args)...)
        , traveltime(t_travel)
        , path(path_mobility)
    {
    }
    double traveltime;
    std::vector<int> path;
};

template <class NodePropertyT, class EdgePropertyT>
class GraphDetailed : public Graph<NodePropertyT, EdgePropertyT>
{
public:
    using Graph<NodePropertyT, EdgePropertyT>::Graph;

    template <class... Args>
    NodeDetailed<NodePropertyT>& add_node(int id, double duration_stay, Args&&... args)
    {
        this->m_nodes.emplace_back(id, duration_stay, std::forward<Args>(args)...);
        return this->m_nodes.back();
    }

    template <class ModelType>
    NodeDetailed<NodePropertyT>& add_node(int id, double duration_stay, ModelType& model1, ModelType& model2)
    {
        this->m_nodes.emplace_back(id, duration_stay, model1, model2);
        return this->m_nodes.back();
    }

    template <class ModelType>
    NodeDetailed<NodePropertyT>& add_node(int id, double duration_stay, ModelType& model1, ModelType& model2,
                                          double m_t0, double m_dt_integration)
    {
        this->m_nodes.emplace_back(id, duration_stay, model1, model2, m_t0, m_dt_integration);
        return this->m_nodes.back();
    }

    template <class... Args>
    EdgeDetailed<EdgePropertyT>& add_edge(size_t start_node_idx, size_t end_node_idx, double traveltime, Args&&... args)
    {
        assert(m_nodes.size() > start_node_idx && m_nodes.size() > end_node_idx);
        return *insert_sorted_replace(
            this->m_edges, Edge<EdgePropertyT>(start_node_idx, end_node_idx, traveltime, std::forward<Args>(args)...),
            [](auto&& e1, auto&& e2) {
                return e1.start_node_idx == e2.start_node_idx ? e1.end_node_idx < e2.end_node_idx
                                                              : e1.start_node_idx < e2.start_node_idx;
            });
    }

    template <class... Args>
    EdgeDetailed<EdgePropertyT>& add_edge(size_t start_node_idx, size_t end_node_idx, double traveltime,
                                          std::vector<int> path, Args&&... args)
    {
        assert(m_nodes.size() > start_node_idx && m_nodes.size() > end_node_idx);
        return *insert_sorted_replace(
            this->m_edges,
            Edge<EdgePropertyT>(start_node_idx, end_node_idx, traveltime, path, std::forward<Args>(args)...),
            [](auto&& e1, auto&& e2) {
                return e1.start_node_idx == e2.start_node_idx ? e1.end_node_idx < e2.end_node_idx
                                                              : e1.start_node_idx < e2.start_node_idx;
            });
    }
};

/**
 * @brief Sets the graph nodes for counties or districts.
 * Reads the node ids which could refer to districts or counties and the epidemiological
 * data from json files and creates one node for each id. Every node contains a model.
 * @param[in] params Model Parameters that are used for every node.
 * @param[in] start_date Start date for which the data should be read.
 * @param[in] end_data End date for which the data should be read.
 * @param[in] data_dir Directory that contains the data files.
 * @param[in] population_data_path Path to json file containing the population data.
 * @param[in] is_node_for_county Specifies whether the node ids should be county ids (true) or district ids (false).
 * @param[in, out] params_graph Graph whose nodes are set by the function.
 * @param[in] read_func Function that reads input data for german counties and sets Model compartments.
 * @param[in] node_func Function that returns the county ids.
 * @param[in] scaling_factor_inf Factor of confirmed cases to account for undetected cases in each county.
 * @param[in] scaling_factor_icu Factor of ICU cases to account for underreporting.
 * @param[in] tnt_capacity_factor Factor for test and trace capacity.
 * @param[in] num_days Number of days to be simulated; required to load data for vaccinations during the simulation.
 * @param[in] export_time_series If true, reads data for each day of simulation and writes it in the same directory as the input files.
 */
template <class TestAndTrace, class ContactPattern, class Model, class MigrationParams, class Parameters,
          class ReadFunction, class NodeIdFunction>
IOResult<void>
set_nodes(const Parameters& params, Date start_date, Date end_date, const fs::path& data_dir,
          const std::string& population_data_path, bool is_node_for_county, Graph<Model, MigrationParams>& params_graph,
          ReadFunction&& read_func, NodeIdFunction&& node_func, const std::vector<double>& scaling_factor_inf,
          double scaling_factor_icu, double tnt_capacity_factor, int num_days = 0, bool export_time_series = false)
{
    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);

    auto read_duration = mio::read_duration_stay(durations_path);
    if (!read_duration) {
        std::cout << read_duration.error().formatted_message() << '\n';
    }
    auto duration_stay = read_duration.value();

    std::string path = "/localdata1/test/memilio/data/pydata/Germany/county_population.json";
    auto read_ids    = node_func(path, true);
    auto node_ids    = read_ids.value();

    std::vector<Model> nodes(node_ids.size(), Model(int(size_t(params.get_num_groups()))));

    for (auto& node : nodes) {
        node.parameters = params;
    }
    auto read_node = read_input_data_county(nodes, start_date, node_ids, scaling_factor_inf, scaling_factor_icu,
                                            "/localdata1/test/memilio/data", num_days);

    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {

        auto tnt_capacity = nodes[node_idx].populations.get_total() * tnt_capacity_factor;

        //local parameters
        auto& tnt_value = nodes[node_idx].parameters.template get<TestAndTrace>();
        tnt_value       = mio::UncertainValue(0.5 * (1.2 * tnt_capacity + 0.8 * tnt_capacity));
        tnt_value.set_distribution(mio::ParameterDistributionUniform(0.8 * tnt_capacity, 1.2 * tnt_capacity));

        //holiday periods
        auto id              = int(mio::regions::CountyId(node_ids[node_idx]));
        auto holiday_periods = mio::regions::get_holidays(mio::regions::get_state_id(id), start_date, end_date);
        auto& contacts       = nodes[node_idx].parameters.template get<mio::osecirvvs::ContactPatterns>();
        contacts.get_school_holidays() =
            std::vector<std::pair<mio::SimulationTime, mio::SimulationTime>>(holiday_periods.size());
        std::transform(
            holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(), [=](auto& period) {
                return std::make_pair(mio::SimulationTime(mio::get_offset_in_days(period.first, start_date)),
                                      mio::SimulationTime(mio::get_offset_in_days(period.second, start_date)));
            });

        //uncertainty in populations
        for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = mio::Index<typename mio::osecirvvs::Model::Compartments>(0);
                 j < mio::osecirvvs::Model::Compartments::Count; ++j) {
                auto& compartment_value = nodes[node_idx].populations[{i, j}];
                compartment_value =
                    mio::UncertainValue(0.5 * (1.1 * double(compartment_value) + 0.9 * double(compartment_value)));
                compartment_value.set_distribution(mio::ParameterDistributionUniform(0.9 * double(compartment_value),
                                                                                     1.1 * double(compartment_value)));
            }
        }

        // Add mobility node
        auto mobility = nodes[node_idx];
        mobility.populations.set_total(0);

        // reduce transmission on contact due to mask obligation in mobility node
        // first age group not able to  (properly) wear masks
        if (masks) {

            const double fact_surgical_mask = 0.1;
            const double fact_ffp2          = 0.001;

            double factor_mask[] = {1, fact_surgical_mask, fact_ffp2, fact_ffp2, fact_ffp2, fact_ffp2};
            if (!ffp2) {
                std::cout << "surgical masks"
                          << "\n";
                // multipliziere alle außer den ersten EIntrag mit 40
                for (size_t j = 1; j < 6; j++) {
                    factor_mask[j] = factor_mask[j] * fact_surgical_mask / fact_ffp2;
                }
            }

            double fac_variant                           = 1.4; //1.94; //https://doi.org/10.7554/eLife.78933
            double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                            0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

            double transmissionProbabilityOnContactMax[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                                            0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};
            for (int i = 0; i < 6; i++) {
                transmissionProbabilityOnContactMin[i] = transmissionProbabilityOnContactMin[i] * factor_mask[i];
                transmissionProbabilityOnContactMax[i] = transmissionProbabilityOnContactMax[i] * factor_mask[i];
            }
            array_assign_uniform_distribution(
                mobility.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>(),
                transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
        }

        for (size_t t_idx = 0; t_idx < num_days; ++t_idx) {
            auto t = mio::SimulationDay((size_t)t_idx);
            for (auto j = mio::AgeGroup(0); j < params.get_num_groups(); j++) {
                mobility.parameters.template get<mio::osecirvvs::DailyFirstVaccination>()[{j, t}] = 0;
                mobility.parameters.template get<mio::osecirvvs::DailyFullVaccination>()[{j, t}]  = 0;
            }
        }

        auto& params_mobility = mobility.parameters;

        params_graph.add_node(node_ids[node_idx], duration_stay((Eigen::Index)node_idx), nodes[node_idx], mobility);
    }
    return success();
}

class MigrationEdgeDetailed : public MigrationEdge
{
public:
    template <class Sim>
    void apply_migration(double t, double dt, SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to, int mode);
};

/**
 * @brief number of migrated people when they return according to the model.
 * E.g. during the time in the other node, some people who left as susceptible will return exposed.
 * Implemented for general compartmentmodel simulations, overload for your custom model if necessary
 * so that it can be found with argument-dependent lookup, i.e. in the same namespace as the model.
 * @param migrated number of people that migrated as input, number of people that return as output
 * @param sim Simulation that is used for the migration
 * @param integrator Integrator that is used for the estimation. Has to be a one-step scheme.
 * @param total total population in the node that the people migrated to.
 * @param t time of migration
 * @param dt timestep
 */
template <class Sim, class = std::enable_if_t<is_compartment_model_simulation<Sim>::value>>
void update_status_migrated(Eigen::Ref<TimeSeries<double>::Vector> migrated, const Sim& sim, IntegratorCore& integrator,
                            Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt)
{
    auto y0 = migrated.eval();
    auto y1 = migrated;
    integrator.step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            sim.get_model().get_derivatives(total, y, t_, dydt);
        },
        y0, t, dt, y1);
}

template <typename FP>
using Vector = Eigen::Matrix<FP, Eigen::Dynamic, 1>;

/*
    * move migrated people from one node to another.
    * @param migrated number of people that migrated 
    * @param results_from current state of node that people migrated from
    * @param results_to current state of node that people migrated to
*/
template <typename FP>
void move_migrated(Eigen::Ref<Vector<FP>> migrated, Eigen::Ref<Vector<FP>> results_from,
                   Eigen::Ref<Vector<FP>> results_to)
{
    // The performance of the low order integrator often fails for small compartments i.e. compartments are sometimes negative.
    // We handel this by checking for values below 0 in m_migrated and add them to the highest value
    Eigen::Index max_index_migrated;
    for (size_t i = 0; i < migrated.size(); ++i) {
        migrated.maxCoeff(&max_index_migrated);
        if (migrated(i) < 0) {
            migrated(max_index_migrated) -= migrated(i);
            migrated(i) = 0;
        }
        if (results_from(i) < migrated(i)) {
            migrated(max_index_migrated) -= (migrated(i) - results_from(i));
            migrated(i) = results_from(i);
        }
        results_from(i) -= migrated(i);
        results_to(i) += migrated(i);
    }
}

template <class Sim>
void MigrationEdgeDetailed::apply_migration(double t, double dt, SimulationNode<Sim>& node_from,
                                            SimulationNode<Sim>& node_to, int mode)
{

    if (mode == 0) {
        //normal daily migration
        m_migrated.add_time_point(
            t, (node_from.get_last_state().array() * m_parameters.get_coefficients().get_matrix_at(t).array() *
                get_migration_factors(node_from, t, node_from.get_last_state()).array())
                   .matrix());
        m_return_times.add_time_point(t + dt);
        move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
                      node_to.get_result().get_last_value());
    }
    // change county of migrated
    else if (mode == 1) {
        // update status of migrated before moving to next county
        IntegratorCore& integrator_node = node_from.get_simulation().get_integrator();
        update_status_migrated(m_migrated.get_last_value(), node_from.get_simulation(), integrator_node,
                               node_from.get_result().get_last_value(), t, dt);
        move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
                      node_to.get_result().get_last_value());
    }
    // option for last time point to remove time points
    else if (mode == 2) {
        Eigen::Index idx                = m_return_times.get_num_time_points() - 1;
        IntegratorCore& integrator_node = node_from.get_simulation().get_integrator();
        update_status_migrated(m_migrated[idx], node_from.get_simulation(), integrator_node,
                               node_from.get_result().get_last_value(), t, dt);

        move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
                      node_to.get_result().get_last_value());

        for (Eigen::Index i = m_return_times.get_num_time_points() - 1; i >= 0; --i) {
            if (m_return_times.get_time(i) <= t) {
                m_migrated.remove_time_point(i);
                m_return_times.remove_time_point(i);
            }
        }
    }
    // just update status of migrated
    else if (mode == 3) {
        Eigen::Index idx                = m_return_times.get_num_time_points() - 1;
        IntegratorCore& integrator_node = node_from.get_simulation().get_integrator();
        update_status_migrated(m_migrated[idx], node_from.get_simulation(), integrator_node,
                               node_from.get_result().get_last_value(), t, dt);

        move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
                      node_from.get_result().get_last_value());
    }
    else {
        std::cout << "Invalid input mode. Should be 0 or 1."
                  << "\n";
    }
}

/**
 * edge functor for migration simulation.
 * @see MigrationEdge::apply_migration
 */
template <class Sim>
void apply_migration(double t, double dt, MigrationEdge& migrationEdge, SimulationNode<Sim>& node_from,
                     SimulationNode<Sim>& node_to, int mode)
{
    migrationEdge.apply_migration(t, dt, node_from, node_to, mode);
}
} // namespace mio

#endif //METAPOPULATION_MOBILITY_DETAILED_H
