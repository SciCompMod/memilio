/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef METAPOPULATION_MOBILITY_INSTANT_H
#define METAPOPULATION_MOBILITY_INSTANT_H

#include "memilio/config.h"
#include "memilio/mobility/graph_simulation.h"
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
#include <vector>

namespace mio
{

/**
 * represents the simulation in one node of the graph.
 */
template <class Sim>
class SimulationNode
{
public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Sim, Args...>::value, void>>
    SimulationNode(Args&&... args)
        : m_simulation(std::forward<Args>(args)...)
        , m_last_state(m_simulation.get_result().get_last_value())
        , m_t0(m_simulation.get_result().get_last_time())
    {
    }

    // when this is called with property(property_arg, m_t0, m_dt_integration). Write the correct constructor
    SimulationNode(const Sim& property_arg, double t0, double dt_integration)
        : m_simulation(property_arg, t0, dt_integration)
        , m_last_state(m_simulation.get_result().get_last_value())
        , m_t0(m_simulation.get_result().get_last_time())
    {
    }

    // Copy-Konstruktor
    SimulationNode(const SimulationNode& other)
        : m_simulation(other.m_simulation)
        , m_last_state(other.m_last_state)
        , m_t0(other.m_t0)
    {
    }

    /**
     * get the result of the simulation in this node.
     * @{
     */
    decltype(auto) get_result() const
    {
        return m_simulation.get_result();
    }
    decltype(auto) get_result()
    {
        return m_simulation.get_result();
    }
    /**@}*/

    /**
     * get the the simulation in this node.
     * @{
     */
    Sim& get_simulation()
    {
        return m_simulation;
    }
    const Sim& get_simulation() const
    {
        return m_simulation;
    }
    /**@}*/

    Eigen::Ref<const Eigen::VectorXd> get_last_state() const
    {
        return m_last_state;
    }

    double get_t0() const
    {
        return m_t0;
    }

    void evolve(double t, double dt)
    {
        m_simulation.advance(t + dt);
        m_last_state = m_simulation.get_result().get_last_value();
    }

private:
    Sim m_simulation;
    Eigen::VectorXd m_last_state;
    double m_t0;
};

/**
 * time dependent migration coefficients.
 */
using MigrationCoefficients = DampingMatrixExpression<VectorDampings>;

/**
 * sum of time dependent migration coefficients.
 * differentiate between sources of migration.
 */
using MigrationCoefficientGroup = DampingMatrixExpressionGroup<MigrationCoefficients>;

/**
 * parameters that influence migration.
 */
class MigrationParameters
{
public:
    /**
     * constructor from migration coefficients.
     * @param coeffs migration coefficients
     */
    MigrationParameters(const MigrationCoefficientGroup& coeffs)
        : m_coefficients(coeffs)
    {
    }

    /**
     * constructor from migration coefficients.
     * @param coeffs migration coefficients
     */
    MigrationParameters(const Eigen::VectorXd& coeffs)
        : m_coefficients({MigrationCoefficients(coeffs)})
    {
    }

    /** 
     * equality comparison operators
     */
    //@{
    bool operator==(const MigrationParameters& other) const
    {
        return m_coefficients == other.m_coefficients;
    }
    bool operator!=(const MigrationParameters& other) const
    {
        return m_coefficients != other.m_coefficients;
    }
    //@}

    /**
     * Get/Setthe migration coefficients.
     * The coefficients represent the (time-dependent) percentage of people migrating 
     * from one node to another by age and infection compartment. 
     * @{
     */
    /**
     * @return the migration coefficients.
     */
    const MigrationCoefficientGroup& get_coefficients() const
    {
        return m_coefficients;
    }
    MigrationCoefficientGroup& get_coefficients()
    {
        return m_coefficients;
    }
    /**
     * @param coeffs the migration coefficients.
     */
    void set_coefficients(const MigrationCoefficientGroup& coeffs)
    {
        m_coefficients = coeffs;
    }
    /** @} */

    /**
     * Get/Set dynamic NPIs that are implemented when relative infections exceed thresholds.
     * This feature is optional. The simulation model needs to overload the get_infected_relative function.
     * @{
     */
    /**
     * @return dynamic NPIs for relative infections.
     */
    const DynamicNPIs& get_dynamic_npis_infected() const
    {
        return m_dynamic_npis;
    }
    DynamicNPIs& get_dynamic_npis_infected()
    {
        return m_dynamic_npis;
    }
    /**
     * @param v dynamic NPIs for relative infections.
     */
    void set_dynamic_npis_infected(const DynamicNPIs& v)
    {
        m_dynamic_npis = v;
    }
    /** @} */

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("MigrationParameters");
        obj.add_element("Coefficients", m_coefficients);
        obj.add_element("DynamicNPIs", m_dynamic_npis);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<MigrationParameters> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("MigrationParameters");
        auto c   = obj.expect_element("Coefficients", Tag<MigrationCoefficientGroup>{});
        auto d   = obj.expect_element("DynamicNPIs", Tag<DynamicNPIs>{});
        return apply(
            io,
            [](auto&& c_, auto&& d_) {
                MigrationParameters params(c_);
                params.set_dynamic_npis_infected(d_);
                return params;
            },
            c, d);
    }

private:
    MigrationCoefficientGroup m_coefficients; //one per group and compartment
    DynamicNPIs m_dynamic_npis;
};

/** 
 * represents the migration between two nodes.
 */
class MigrationEdge
{
public:
    /**
     * create edge with coefficients.
     * @param coeffs % of people in each group and compartment that migrate in each time step.
     */
    MigrationEdge(const MigrationParameters& params)
        : m_parameters(params)
        , m_migrated(params.get_coefficients().get_shape().rows())
        , m_return_times(0)
    {
    }

    /**
     * create edge with coefficients.
     * @param coeffs % of people in each group and compartment that migrate in each time step.
     */
    MigrationEdge(const Eigen::VectorXd& coeffs)
        : m_parameters(coeffs)
        , m_migrated(coeffs.rows())
        , m_return_times(0)
    {
    }

    /**
     * get the migration parameters.
     */
    const MigrationParameters& get_parameters() const
    {
        return m_parameters;
    }

    /**
     * compute migration from node_from to node_to.
     * migration is based on coefficients.
     * migrants are added to the current state of node_to, subtracted from node_from.
     * on return, migrants (adjusted for infections) are subtracted from node_to, added to node_from.
     * @param t current time
     * @param dt last time step (fixed to 0.5 for migration model)
     * @param node_from node that people migrated from, return to
     * @param node_to node that people migrated to, return from
     */
    template <class Sim>
    void apply_migration(double t, double dt, SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to, int mode);

private:
    MigrationParameters m_parameters;
    TimeSeries<double> m_migrated;
    TimeSeries<double> m_return_times;
    double m_t_last_dynamic_npi_check               = -std::numeric_limits<double>::infinity();
    std::pair<double, SimulationTime> m_dynamic_npi = {-std::numeric_limits<double>::max(), SimulationTime(0)};
};

/**
 * adjust number of migrated people when they return according to the model.
 * E.g. during the time in the other node, some people who left as susceptible will return exposed.
 * Implemented for general compartmentmodel simulations, overload for your custom model if necessary
 * so that it can be found with argument-dependent lookup, i.e. in the same namespace as the model.
 * @param[inout] migrated number of people that migrated as input, number of people that return as output
 * @param params parameters of model in the node that the people migrated to.
 * @param total total population in the node that the people migrated to.
 * @param t time of migration
 * @param dt time between migration and return
 */
template <class Sim, class = std::enable_if_t<is_compartment_model_simulation<Sim>::value>>
void update_status_migrated(Eigen::Ref<TimeSeries<double>::Vector> migrated, Sim& sim, IntegratorCore& integrator,
                            Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt)
{
    // 1) get flows in mobility comp
    // 2) assure that this does not change the local model before return. (könnte den integrator zur diff nutzen)
    auto& model = sim.get_model();
    auto y0     = migrated.eval();
    auto y1     = migrated;
    EulerIntegratorCore().step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            sim.get_model().get_derivatives(total, y, t_, dydt);
        },
        y0, t, dt, y1);

    // if integrator type == EulerIntegratorCore
    if (dynamic_cast<EulerIntegratorCore*>(&integrator) != nullptr) {
        auto flows_step                  = model.get_flow_values();
        std::vector<int> indx_infections = {0,   33,  17,  53,  86,  70,  106, 139, 123,
                                            159, 192, 176, 212, 245, 229, 265, 298, 282};
        auto sum_infections =
            std::accumulate(indx_infections.begin(), indx_infections.end(), 0.0, [&flows_step](double sum, int i) {
                return sum + flows_step(i);
            });
        sim.get_flows().get_last_value()(0) += sum_infections;

        std::vector<int> indx_symp = {20, 36, 56, 73, 89, 109, 126, 142, 162, 179, 195, 215, 232, 248, 268, 285, 301};
        auto sum_symp = std::accumulate(indx_symp.begin(), indx_symp.end(), 0.0, [&flows_step](double sum, int i) {
            return sum + flows_step(i);
        });
        sim.get_flows().get_last_value()(1) += sum_symp;
    }
}

template <typename FP>
using Vector = Eigen::Matrix<FP, Eigen::Dynamic, 1>;

/*
    * move migrated people from one node to another.
    * @param migrated number of people that migrated 
    * @param results_from current state of node that people migrated from
    * @param results_to current state of node that people migrated to
*/
// template <typename FP>
// void move_migrated(Eigen::Ref<Vector<FP>> migrated, Eigen::Ref<Vector<FP>> results_from,
//                    Eigen::Ref<Vector<FP>> results_to)
// {
//     double sum_migrated = migrated.sum();
//     for (size_t i = 0; i < migrated.size(); ++i) {
//         // Ensure that migrated(i) is not greater than results_from(i)
//         if (migrated(i) > results_from(i)) {
//             // Distribute the excess to other elements in migrated
//             FP excess   = migrated(i) - results_from(i);
//             migrated(i) = results_from(i);
//             for (size_t j = 0; j < migrated.size(); ++j) {
//                 if (j != i && migrated(j) + excess <= results_from(j)) {
//                     migrated(j) += excess;
//                     break;
//                 }
//             }
//         }
//     }
//     double sum_migrated_new = migrated.sum();
//     if (std::abs(sum_migrated - sum_migrated_new) > 1e-10) {
//         std::cout << "Sum of migration is not the same after migration. Difference: " << sum_migrated - sum_migrated_new
//                   << "\n";
//     }

//     for (size_t i = 0; i < migrated.size(); ++i) {
//         results_from(i) -= migrated(i);
//         results_to(i) += migrated(i);
//     }
// }

inline std::vector<double> calc_sum_groupss(Eigen::Ref<Eigen::Matrix<double, -1, 1>> y)
{
    const int num_age_groups = 6;
    std::vector<double> sum_groups(num_age_groups);
    auto const num_comparts = y.size() / num_age_groups;
    for (int age_group = 0; age_group < num_age_groups; ++age_group) {
        sum_groups[age_group] = y.segment(age_group * num_comparts, num_comparts).sum();
    }
    return sum_groups;
}

template <typename FP>
void move_migrated(Eigen::Ref<Vector<FP>> migrated, Eigen::Ref<Vector<FP>> results_from,
                   Eigen::Ref<Vector<FP>> results_to)
{
    const auto group        = 6;
    const auto num_comparts = results_to.size() / group;

    // double sum_migrated = migrated.sum();
    // double sum_from     = results_from.sum();
    // double sum_to       = results_to.sum();

    // auto groups_from  = calc_sum_groupss(results_from);
    // auto groups_to    = calc_sum_groupss(results_to);
    // auto groups_moved = calc_sum_groupss(migrated);

    // std::vector<double> migrated_as_vector(migrated.data(), migrated.data() + migrated.size());
    // std::vector<double> results_from_as_vector(results_from.data(), results_from.data() + results_from.size());
    // std::vector<double> results_to_as_vector(results_to.data(), results_to.data() + results_to.size());

    for (Eigen::Index j = 0; j < migrated.size(); ++j) {
        if (migrated(j) < -1e-10) {
            auto compart        = j % num_comparts;
            auto curr_age_group = int(j / num_comparts);
            auto indx_begin     = curr_age_group * num_comparts;
            auto indx_end       = (curr_age_group + 1) * num_comparts;
            // calculate max index in indx boundaries
            Eigen::Index max_index = indx_begin;
            for (Eigen::Index i = indx_begin; i < indx_end; ++i) {
                if (migrated(i) > migrated(max_index)) {
                    max_index = i;
                }
            }
            // we assume that the solution from migrated is bettter because there is contact with other nodes
            migrated(max_index) = migrated(max_index) + migrated(j);
            migrated(j)         = migrated(j) - migrated(j);
        }
    }

    // // check if a entry of migrated_as_vector is negative
    // for (size_t i = 0; i < migrated_as_vector.size(); ++i) {
    //     if (migrated_as_vector[i] < 0) {
    //         std::cout << "Negative entry in migrated_as_vector at index " << i << "\n";
    //     }
    // }

    Eigen::VectorXd remaining_after_return = (results_from - migrated).eval();
    // auto remaining_after_return_as_vector  = std::vector<double>(
    //     remaining_after_return.data(), remaining_after_return.data() + remaining_after_return.size());
    for (Eigen::Index j = 0; j < results_to.size(); ++j) {
        if (remaining_after_return(j) < -1e-10) {
            auto compart        = j % num_comparts;
            auto curr_age_group = int(j / num_comparts);
            auto indx_begin     = curr_age_group * num_comparts;
            auto indx_end       = (curr_age_group + 1) * num_comparts;
            // calculate max index in indx boundaries
            Eigen::Index max_index = indx_begin;
            for (Eigen::Index i = indx_begin; i < indx_end; ++i) {
                if (remaining_after_return(i) > remaining_after_return(max_index)) {
                    max_index = i;
                }
            }

            // its possible that the max value in the boundaries is not enough to fill the negative value.
            // Therefore we have to find multiple max values
            while (remaining_after_return(max_index) + remaining_after_return(j) < -1e-10) {

                // calculate sum between indx_begin and indx_end
                double result_from_sum_group      = 0;
                double result_migration_sum_group = 0;
                for (Eigen::Index i = indx_begin; i < indx_end; ++i) {
                    result_from_sum_group += results_from(i);
                    result_migration_sum_group += migrated(i);
                }
                auto diff = result_from_sum_group - result_migration_sum_group;
                if (diff < -1e-10) {
                    std::cout << "Sum of results_from is smaller than sum of migrated. Diff is "
                              << result_from_sum_group - result_migration_sum_group << "\n";
                    // std::cout << "diff overall: " << sum_from - sum_migrated << "\n";
                }

                results_from(j)         = results_from(j) + remaining_after_return(max_index);
                results_from(max_index) = results_from(max_index) - remaining_after_return(max_index);
                remaining_after_return  = (results_from - migrated).eval();

                max_index = indx_begin;
                for (Eigen::Index i = indx_begin; i < indx_end; ++i) {
                    if (remaining_after_return(i) > remaining_after_return(max_index)) {
                        max_index = i;
                    }
                }
                if (max_index == indx_begin && remaining_after_return(max_index) == 0) {
                    std::cout << "Fixing negative Value in migration not possible."
                              << "\n";
                }
            }

            // we assume that the solution from migrated is bettter because there is contact with other nodes
            results_from(j)         = results_from(j) - remaining_after_return(j);
            results_from(max_index) = results_from(max_index) + remaining_after_return(j);
            remaining_after_return  = (results_from - migrated).eval();
        }
    }

    // Subtrahiere migrated von results_from
    results_from -= migrated;

    // Addiere migrated zu results_to
    results_to += migrated;

    // again as vectors
    // std::vector<double> migrated_as_vector_new(migrated.data(), migrated.data() + migrated.size());
    // std::vector<double> results_from_as_vector_new(results_from.data(), results_from.data() + results_from.size());
    // std::vector<double> results_to_as_vector_new(results_to.data(), results_to.data() + results_to.size());

    // // check if a entry of migrated_as_vector_new is negative
    // for (size_t i = 0; i < migrated_as_vector_new.size(); ++i) {
    //     if (migrated_as_vector_new[i] < -1e-8) {
    //         std::cout << "Negative entry in migrated_as_vector_new at index " << i << " with value "
    //                   << migrated_as_vector_new[i] << "\n";
    //     }
    // }
    // // same for results_from_as_vector_new and results_to_as_vector_new
    // for (size_t i = 0; i < results_from_as_vector_new.size(); ++i) {
    //     if (results_from_as_vector_new[i] < -1e-8) {
    //         std::cout << "Negative entry in results_from_as_vector_new at index " << i << "with value "
    //                   << results_from_as_vector_new[i] << "\n";
    //     }
    // }
    // for (size_t i = 0; i < results_to_as_vector_new.size(); ++i) {
    //     if (results_to_as_vector_new[i] < -1e-8) {
    //         std::cout << "Negative entry in results_to_as_vector_new at index " << i << "with value "
    //                   << results_to_as_vector_new[i] << "\n";
    //     }
    // }

    // // teste if sum_from und sum_to gleich sind
    // if (std::abs(sum_from - sum_to) > 1e-7) {
    //     double sum_migrated_new = migrated.sum();
    //     double sum_from_new     = results_from.sum();
    //     double sum_to_new       = results_to.sum();
    //     if (std::abs(sum_migrated - sum_migrated_new) > 1e-7) {
    //         std::cout << "Sum of migration is not the same after migration. Difference: "
    //                   << sum_migrated - sum_migrated_new << "\n";
    //     }
    //     if (std::abs(sum_from - sum_migrated - sum_from_new) > 1e-7) {
    //         std::cout << "Sum of from is not the same after migration. Difference: " << sum_from - sum_from_new << "\n";
    //     }
    //     if (std::abs(sum_to + sum_migrated - sum_to_new) > 1e-7) {
    //         std::cout << "Sum of to is not the same after migration. Difference: " << sum_to - sum_to_new << "\n";
    //     }
    // }

    // // check if the sum of the groups is the same
    // auto groups_from_new  = calc_sum_groupss(results_from);
    // auto groups_to_new    = calc_sum_groupss(results_to);
    // auto groups_moved_new = calc_sum_groupss(migrated);

    // for (size_t i = 0; i < groups_from_new.size(); ++i) {
    //     if (std::abs(groups_from[i] + groups_to[i] - groups_from_new[i] - groups_to_new[i]) > 1e-7) {
    //         std::cout << "Sum of group " << i << " in from is not the same after migration. Difference: "
    //                   << groups_from[i] + groups_to[i] - groups_from_new[i] - groups_to_new[i] << "\n";
    //     }
    //     if (std::abs(groups_moved[i] - groups_moved_new[i]) > 1e-7) {
    //         std::cout << "Sum of group " << i << " in moved is not the same after migration. Difference: "
    //                   << groups_moved[i] - groups_moved_new[i] << "\n";
    //     }
    // }

} // namespace mio

/**
 * detect a get_infections_relative function for the Model type.
 */
template <class Sim>
using get_infections_relative_expr_t = decltype(get_infections_relative(
    std::declval<const Sim&>(), std::declval<double>(), std::declval<const Eigen::Ref<const Eigen::VectorXd>&>()));

/**
 * get the percantage of infected people of the total population in the node
 * If dynamic NPIs are enabled, there needs to be an overload of get_infections_relative(model, y)
 * for the Model type that can be found with argument-dependent lookup. Ideally define get_infections_relative 
 * in the same namespace as the Model type.
 * @param node a node of a migration graph.
 * @param y the current value of the simulation.
 * @param t the current simulation time
 */
template <class Sim,
          std::enable_if_t<!is_expression_valid<get_infections_relative_expr_t, Sim>::value, void*> = nullptr>
double get_infections_relative(const SimulationNode<Sim>& /*node*/, double /*t*/,
                               const Eigen::Ref<const Eigen::VectorXd>& /*y*/)
{
    assert(false && "Overload get_infections_relative for your own model/simulation if you want to use dynamic NPIs.");
    return 0;
}
template <class Sim, std::enable_if_t<is_expression_valid<get_infections_relative_expr_t, Sim>::value, void*> = nullptr>
double get_infections_relative(const SimulationNode<Sim>& node, double t, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    return get_infections_relative(node.get_simulation(), t, y);
}

/**
 * detect a get_migration_factors function for the Model type.
 */
template <class Sim>
using get_migration_factors_expr_t = decltype(get_migration_factors(
    std::declval<const Sim&>(), std::declval<double>(), std::declval<const Eigen::Ref<const Eigen::VectorXd>&>()));

/**
 * Get an additional migration factor.
 * The absolute migration for each compartment is computed by c_i * y_i * f_i, wher c_i is the coefficient set in 
 * MigrationParameters, y_i is the current compartment population, f_i is the factor returned by this function.
 * This factor is optional, default 1.0. If you need to adjust migration in that way, overload get_migration_factors(model, t, y) 
 * for your Model type so that can be found with argument-dependent lookup.
 * @param node a node of a migration graph.
 * @param y the current value of the simulation.
 * @param t the current simulation time
 * @return a vector expression, same size as y, with the factor for each compartment.
 */
template <class Sim, std::enable_if_t<!is_expression_valid<get_migration_factors_expr_t, Sim>::value, void*> = nullptr>
auto get_migration_factors(const SimulationNode<Sim>& /*node*/, double /*t*/,
                           const Eigen::Ref<const Eigen::VectorXd>& y)
{
    return Eigen::VectorXd::Ones(y.rows());
}
template <class Sim, std::enable_if_t<is_expression_valid<get_migration_factors_expr_t, Sim>::value, void*> = nullptr>
auto get_migration_factors(const SimulationNode<Sim>& node, double t, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    return get_migration_factors(node.get_simulation(), t, y);
}

/**
 * detect a get_migration_factors function for the Model type.
 */
template <class Sim>
using test_commuters_expr_t = decltype(
    test_commuters(std::declval<Sim&>(), std::declval<Eigen::Ref<const Eigen::VectorXd>&>(), std::declval<double>()));

/**
 * Test persons when migrating from their source node.
 * May transfer persons between compartments, e.g., if an infection was detected.
 * This feature is optional, default implementation does nothing.
 * In order to support this feature for your model, implement a test_commuters overload 
 * that can be found with argument-dependent lookup.
 * @param node a node of a migration graph.
 * @param migrated mutable reference to vector of persons per compartment that migrate.
 * @param t the current simulation time.
 */
template <class Sim, std::enable_if_t<!is_expression_valid<test_commuters_expr_t, Sim>::value, void*> = nullptr>
void test_commuters(SimulationNode<Sim>& /*node*/, Eigen::Ref<Eigen::VectorXd> /*migrated*/, double /*time*/)
{
}
template <class Sim, std::enable_if_t<is_expression_valid<test_commuters_expr_t, Sim>::value, void*> = nullptr>
void test_commuters(SimulationNode<Sim>& node, Eigen::Ref<Eigen::VectorXd> migrated, double time)
{
    return test_commuters(node.get_simulation(), migrated, time);
}

template <class Sim>
void MigrationEdge::apply_migration(double t, double dt, SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to,
                                    int mode)
{
    //check dynamic npis
    if (m_t_last_dynamic_npi_check == -std::numeric_limits<double>::infinity()) {
        m_t_last_dynamic_npi_check = node_from.get_t0();
    }

    auto& dyn_npis = m_parameters.get_dynamic_npis_infected();
    if (dyn_npis.get_thresholds().size() > 0 &&
        floating_point_greater_equal(t, m_t_last_dynamic_npi_check + dyn_npis.get_interval().get())) {
        auto inf_rel = get_infections_relative(node_from, t, node_from.get_last_state()) * dyn_npis.get_base_value();
        auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
        if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
            (exceeded_threshold->first > m_dynamic_npi.first ||
             t > double(m_dynamic_npi.second))) { //old NPI was weaker or is expired
            auto t_end    = SimulationTime(t + double(dyn_npis.get_duration()));
            m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
            implement_dynamic_npis(
                m_parameters.get_coefficients(), exceeded_threshold->second, SimulationTime(t), t_end, [this](auto& g) {
                    return make_migration_damping_vector(m_parameters.get_coefficients().get_shape(), g);
                });
        }
        m_t_last_dynamic_npi_check = t;
    }

    // //returns
    // for (Eigen::Index i = m_return_times.get_num_time_points() - 1; i >= 0; --i) {
    //     if (mode == 1) {
    //         auto v0 = find_value_reverse(node_to.get_result(), m_migrated.get_time(i), 1e-10, 1e-10);
    //         assert(v0 != node_to.get_result().rend() && "unexpected error.");
    //         IntegratorCore& integrator_node = node_to.get_simulation().get_integrator();
    //         update_status_migrated(m_migrated[i], node_to.get_simulation(), integrator_node, *v0,
    //                                m_migrated.get_time(i), dt);

    //         //the lower-order return calculation may in rare cases produce negative compartments,
    //         //especially at the beginning of the simulation.
    //         //fix by subtracting the supernumerous returns from the biggest compartment of the age group.
    //         Eigen::VectorXd remaining_after_return = (node_to.get_result().get_last_value() - m_migrated[i]).eval();
    //         for (Eigen::Index j = 0; j < node_to.get_result().get_last_value().size(); ++j) {
    //             if (remaining_after_return(j) < 0) {
    //                 auto num_comparts = (Eigen::Index)Sim::Model::Compartments::Count;
    //                 auto group        = Eigen::Index(j / num_comparts);
    //                 auto compart      = j % num_comparts;
    //                 log(remaining_after_return(j) < -1e-3 ? LogLevel::warn : LogLevel::info,
    //                     "Underflow during migration returns at time {}, compartment {}, age group {}: {}", t, compart,
    //                     group, remaining_after_return(j));
    //                 Eigen::Index max_index;
    //                 slice(remaining_after_return, {group * num_comparts, num_comparts}).maxCoeff(&max_index);
    //                 log_info("Transferring to compartment {}", max_index);
    //                 max_index += group * num_comparts;
    //                 m_migrated[i](max_index) -= remaining_after_return(j);
    //                 m_migrated[i](j) += remaining_after_return(j);
    //             }
    //         }
    //         node_from.get_result().get_last_value() += m_migrated[i];
    //         node_to.get_result().get_last_value() -= m_migrated[i];
    //         m_migrated.remove_time_point(i);
    //         m_return_times.remove_time_point(i);
    //     }
    // }

    if (mode == 0) {
        //normal daily migration
        m_migrated.add_time_point(
            t, (node_from.get_last_state().array() * m_parameters.get_coefficients().get_matrix_at(t).array() *
                get_migration_factors(node_from, t, node_from.get_last_state()).array())
                   .matrix());
        m_return_times.add_time_point(t + dt);
        // printe twas, wenn die summe von m_migration kleiner als 1e-2 ist
        if (m_migrated.get_last_value().sum() < 1e-2) {
            std::cout << "Sum of migration is smaller than 1e-2. Sum: " << m_migrated.get_last_value().sum() << "\n";
        }
        move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
                      node_to.get_result().get_last_value());
    }
    // change county of migrated
    else if (mode == 1) {
        // update status of migrated before moving to next county
        Eigen::Index idx = m_return_times.get_num_time_points() - 1;

        auto node_from_as_vector = std::vector<double>(node_from.get_result().get_last_value().data(),
                                                       node_from.get_result().get_last_value().data() +
                                                           node_from.get_result().get_last_value().size());
        auto node_to_as_vector   = std::vector<double>(node_to.get_result().get_last_value().data(),
                                                     node_to.get_result().get_last_value().data() +
                                                         node_to.get_result().get_last_value().size());

        IntegratorCore& integrator_node = node_from.get_simulation().get_integrator();
        // when the sum of both is equal, we can assume that only the migrated group is currently in the node.
        // Therefore, we can use the last value calculated by the node integration.

        if (std::abs(m_migrated[idx].sum() - node_from.get_result().get_last_value().sum()) > 1e-9) {
            update_status_migrated(m_migrated[idx], node_from.get_simulation(), integrator_node,
                                   node_from.get_result().get_last_value(), t, dt);
        }
        else {
            for (size_t i = 0; i < m_migrated[idx].size(); ++i) {
                m_migrated[idx](i) = node_from.get_result().get_last_value()(i);
            }
        }

        // teste, ob die abweichung im bereich 1e-4 bis 1e-12 liegt
        double diff_from_migrated = std::abs(m_migrated[idx].sum() - node_from.get_result().get_last_value().sum());
        if (diff_from_migrated < 1e-2 && diff_from_migrated > 1e-9) {
            if (m_migrated[idx].sum() - node_from.get_result().get_last_value().sum() > 0)
                std::cout << "Difference between migrated and node_from is in the range 1e-4 to 1e-10. Difference: "
                          << diff_from_migrated << "\n";
        }
        if (m_migrated[idx].sum() - node_from.get_result().get_last_value().sum() > 0) {
            auto dif = m_migrated[idx].sum() - node_from.get_result().get_last_value().sum();
            std::cout << "Difference between migrated and node_from bad. Difference: " << dif << "\n";
        }

        move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
                      node_to.get_result().get_last_value());
    }
    // option for last time point to remove time points
    else if (mode == 2) {
        // Eigen::Index idx = m_return_times.get_num_time_points() - 1;

        // if (std::abs(m_migrated[idx].sum() - node_from.get_result().get_last_value().sum()) > 1e-9) {
        //     update_status_migrated(m_migrated[idx], node_from.get_simulation(),
        //                            node_from.get_simulation().get_integrator(), node_from.get_result().get_last_value(),
        //                            t, dt);
        // }
        // else {
        //     for (size_t i = 0; i < m_migrated[idx].size(); ++i) {
        //         m_migrated[idx](i) = node_from.get_result().get_last_value()(i);
        //     }
        // }

        // double diff_from_migrated = std::abs(m_migrated[idx].sum() - node_from.get_result().get_last_value().sum());
        // if (diff_from_migrated < 1e-2 && diff_from_migrated > 1e-9) {
        //     if (m_migrated[idx].sum() - node_from.get_result().get_last_value().sum() > 0)
        //         std::cout << "Difference between migrated and node_from is in the range 1e-4 to 1e-9. Difference: "
        //                   << diff_from_migrated << "\n";
        // }
        // if (m_migrated[idx].sum() - node_from.get_result().get_last_value().sum() > 0) {
        //     auto dif = m_migrated[idx].sum() - node_from.get_result().get_last_value().sum();<s
        //     std::cout << "Difference between migrated and node_from bad. Difference: " << dif << "\n";
        // }

        // auto node_from_as_vector = std::vector<double>(node_from.get_result().get_last_value().data(),
        //                                                node_from.get_result().get_last_value().data() +
        //                                                    node_from.get_result().get_last_value().size());

        // move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
        //               node_to.get_result().get_last_value());

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

        if (std::abs(m_migrated[idx].sum() - node_from.get_result().get_last_value().sum()) > 1e-12) {
            update_status_migrated(m_migrated[idx], node_from.get_simulation(), integrator_node,
                                   node_from.get_result().get_last_value(), t, dt);
        }
        else {
            for (size_t i = 0; i < m_migrated[idx].size(); ++i) {
                // if theres only one group located in the node, we take the result from the high order integration.
                m_migrated[idx](i) = node_from.get_result().get_last_value()(i);
            }
        }
        move_migrated(m_migrated.get_last_value(), node_from.get_result().get_last_value(),
                      node_from.get_result().get_last_value());
    }

    else {
        std::cout << "Invalid input mode. Should be 0, 1 or 2."
                  << "\n";
    }
}

/**
 * edge functor for migration simulation.
 * @see SimulationNode::evolve
 */
template <class Sim>
void evolve_model(double t, double dt, SimulationNode<Sim>& node)
{
    node.evolve(t, dt);
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
/**
 * create a migration simulation.
 * After every second time step, for each edge a portion of the population corresponding to the coefficients of the edge
 * moves from one node to the other. In the next timestep, the migrated population return to their "home" node. 
 * Returns are adjusted based on the development in the target node. 
 * @param t0 start time of the simulation
 * @param dt time step between migrations
 * @param graph set up for migration simulation
 * @{
 */
template <class Sim>
GraphSimulation<Graph<SimulationNode<Sim>, MigrationEdge>>
make_migration_sim(double t0, double dt, const Graph<SimulationNode<Sim>, MigrationEdge>& graph)
{
    return make_graph_sim(t0, dt, graph, &evolve_model<Sim>, &apply_migration<Sim>);
}

template <class Sim>
GraphSimulation<Graph<SimulationNode<Sim>, MigrationEdge>>
make_migration_sim(double t0, double dt, Graph<SimulationNode<Sim>, MigrationEdge>&& graph)
{
    return make_graph_sim(t0, dt, std::move(graph), &evolve_model<Sim>, &apply_migration<Sim>);
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_INSTANT_H
