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

#include "boost/filesystem.hpp"

#include <cassert>

namespace mio
{

//  Erbe von MigrationEdge aus metapopulation_mobility_instant.h und überschreibe die apply_migration Funktion indem diese zusätzlich das argument mode bekommt
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
