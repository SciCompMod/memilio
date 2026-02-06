/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Martin J. Kuehn, Daniel Abele
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
#ifndef MIO_EPI_DYNAMIC_LOCKDOWN_H
#define MIO_EPI_DYNAMIC_LOCKDOWN_H

#include "memilio/epidemiology/damping_sampling.h"
#include "memilio/utils/stl_util.h"

namespace mio
{

/**
 * represents non-pharmaceutical interventions (NPI) that are activated during the simulation if
 * some value (e.g. infections) exceeds specified thresholds.
 */
template <typename FP>
class DynamicNPIs
{
public:
    /**
     * set a threshold and the NPIs that should be enacted.
     * @param threshold the threshold that may be exceeded.
     * @param dampings the NPIs
     */
    void set_threshold(FP threshold, const std::vector<DampingSampling<FP>>& dampings)
    {
        insert_sorted_replace(m_thresholds, std::make_pair(threshold, dampings), [](auto& t1, auto& t2) {
            return t1.first > t2.first;
        });
    }

    /**
     * find the highest threshold that is smaller than the value.
     * @param value value to compare against the thresholds.
     * @return iterator (see get_thresholds()) to the exceeded threshold if found.
     *         get_thresholds().end() otherwise
     */
    auto get_max_exceeded_threshold(FP value)
    {
        //thresholds are sorted by value descending, so upper_bound returns the first threshold that is smaller using binary search
        auto iter_max_exceeded_threshold =
            std::upper_bound(m_thresholds.begin(), m_thresholds.end(), value, [](auto& val, auto& t2) {
                return val > t2.first;
            });
        return iter_max_exceeded_threshold;
    }

    /**
     * range of pairs of threshold values and NPIs.
     * thresholds are sorted by value in descending order.
     * @returns the range of set thresholds.
     * @{
     */
    auto get_thresholds() const
    {
        return make_range(m_thresholds.begin(), m_thresholds.end());
    }
    auto get_thresholds()
    {
        return make_range(m_thresholds.begin(), m_thresholds.end());
    }
    /** @} */

    /**
     * get the duration of the NPIs.
     * @return the duration of the NPIs.
     */
    SimulationTime<FP> get_duration() const
    {
        return m_duration;
    }

    /**
     * set the duration of the NPIs.
     * @param v the duration of the NPIs.
     */
    void set_duration(SimulationTime<FP> v)
    {
        m_duration = v;
    }

    /**
     * Get/Set the interval at which the NPIs are checked.
     * @{
     */
    /**
     * @return the interval at which the NPIs are checked.
     */
    SimulationTime<FP> get_interval() const
    {
        return m_interval;
    }
    /**
     * @param interval The interval at which the NPIs are checked.
     */
    void set_interval(SimulationTime<FP> interval)
    {
        m_interval = interval;
    }
    /**@}*/

    /**
     * Get/Set the base value of the thresholds.
     * The base value determines the unit of the threshold values.
     * E.g. If the base value is X, the thresholds should be interpreted as cases per X people.
     * @{
     */
    /**
     * @return The base value of the thresholds.
     */
    FP get_base_value() const
    {
        return m_base;
    }
    /**
     * @return The base value of the thresholds.
     */
    void set_base_value(FP v)
    {
        m_base = v;
    }
    /**@}*/

    /**
     * draw a random sample from the damping distributions
     */
    void draw_sample()
    {
        for (auto&& t : m_thresholds) {
            for (auto&& d : t.second) {
                d.draw_sample();
            }
        }
    }

    /**
     * serialize this.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("DynamicNPIs");
        obj.add_list("Thresholds", get_thresholds().begin(), get_thresholds().end());
        obj.add_element("Duration", get_duration());
        obj.add_element("Interval", get_interval());
        obj.add_element("BaseValue", get_base_value());
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<DynamicNPIs> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("DynamicNPIs");
        auto t   = obj.expect_list("Thresholds", Tag<std::pair<FP, std::vector<DampingSampling<FP>>>>{});
        auto d   = obj.expect_element("Duration", Tag<SimulationTime<FP>>{});
        auto i   = obj.expect_element("Interval", Tag<SimulationTime<FP>>{});
        auto b   = obj.expect_element("BaseValue", Tag<FP>{});
        return apply(
            io,
            [](auto&& t_, auto&& d_, auto&& i_, auto&& b_) {
                auto npis = DynamicNPIs();
                npis.set_duration(d_);
                npis.set_interval(i_);
                npis.set_base_value(b_);
                for (auto&& e : t_) {
                    npis.set_threshold(e.first, e.second);
                }
                return npis;
            },
            t, d, i, b);
    }

private:
    std::vector<std::pair<FP, std::vector<DampingSampling<FP>>>> m_thresholds;
    SimulationTime<FP> m_duration{14.0};
    SimulationTime<FP> m_interval{3.0};
    FP m_base{1.0};
};

/**
 * Get a list of indices of specified dampings.
 * Returns the indices of dampings that match the given type and level and that become active in the specified
 * time span (excluding the particular interval boundaries, begin and end).
 * Utility for implementation of dynamic NPIs.
 * @param damping_expr some matrix expression that contains dampings, e.g. a ContactMatrix.
 * @param lvl damping level to match
 * @param type damping type to match
 * @param begin beginning of the time span that contains the dampings
 * @param end end of the time span that contains the dampings.
 * @return list of indices in range damping_expr.get_dampings()
 */
template <typename FP, class DampingExpr>
std::vector<size_t> get_damping_indices(const DampingExpr& damping_expr, DampingLevel lvl, DampingType type,
                                        SimulationTime<FP> begin, SimulationTime<FP> end)
{
    std::vector<size_t> indices;
    for (size_t i = 0; i < damping_expr.get_dampings().size(); ++i) {
        const auto d = damping_expr.get_dampings()[i];
        if (d.get_level() == lvl && d.get_type() == type && d.get_time() > begin && d.get_time() < end) {
            indices.push_back(i);
        }
    }
    return indices;
}

/**
 * Get the value of the damping that matches the given type and level and that is active at the specified time.
 * If no damping is found, returns a zero matrix of the correct shape.
 * Utility for implementation of dynamic NPIs.
 * @param damping_expr some matrix expression that contains dampings, e.g. a ContactMatrix.
 * @param lvl damping level to match
 * @param type damping type to match
 * @param time time where the damping is active
 * @return matrix of damping coefficients if active damping is found.
 *         zero matrix otherwise.
 */
template <typename FP, class DampingExpr>
Eigen::Ref<const typename DampingExpr::Matrix> get_active_damping(const DampingExpr& damping_expr, DampingLevel lvl,
                                                                  DampingType type, SimulationTime<FP> t)
{
    auto ub =
        std::find_if(damping_expr.get_dampings().rbegin(), damping_expr.get_dampings().rend(), [lvl, type, t](auto& d) {
            return d.get_level() == lvl && d.get_type() == type && d.get_time() <= t;
        });
    if (ub != damping_expr.get_dampings().rend()) {
        return ub->get_coeffs();
    }
    return DampingExpr::Matrix::Zero(damping_expr.get_shape().rows(), damping_expr.get_shape().cols());
}

/**
 * implement dynamic NPIs for a time span.
 * Adds or removes dampings to ensure that the active dampings during the specified
 * time span is at least as big as the specified dynamic dampings.
 * If another damping of the same type and level is active at the beginning of the time span or
 * becomes active during the time span, the coefficient wise maximum of the new damping and the existing damping
 * is used.
 * At the end of the time span, another set of dampings may be added that restores the dampings on each level and type as they
 * would have been without the dynamic npis that have just been implemented.
 * Examples:
 * a) no damping exists yet, dynamic npi of value `d`:
 *     one damping is added at the beginning of the time span that has the value `d`,
 *     another damping is added at the end of the time span that has a value zero.
 * b) damping of value `a` is active before the beginning of the time span, dynamic npi of value `d` is added:
 *     one damping is added at the beginning of the time span that has the value `max(a, d)`,
 *     another damping is added at the end of the time span that has the value a
 * b) damping of value `a` becomes active at a time `t_a` between the beginning of the time span and the end, dynamic npi of value `d` is added:
 *     one damping is added at the beginning of the time span that has the value `d`,
 *     the value of the damping at time `t_a` is set to `max(d, a)`,
 *     another damping is added at the end of the time span that has the value a
 * @param damping_expr_group a group of matrix expressions that contains dampings, e.g. a ContactMatrixGroup.
 * @param dynamic_npis the NPIs to be implemented
 * @param begin beginning of the time span that the NPIs will be active for.
 * @param end end of the time span that the NPIs will be active for.
 * @param make_matrix function to make a matrix of the same shape as the damping expression, see e.g. make_contact_damping_matrix
 */
template <typename FP, class DampingExprGroup, class MakeMatrix>
void implement_dynamic_npis(DampingExprGroup& damping_expr_group, const std::vector<DampingSampling<FP>>& npis,
                            SimulationTime<FP> begin, SimulationTime<FP> end, MakeMatrix&& make_matrix)
{
    for (auto& npi : npis) {
        for (auto& mat_idx : npi.get_matrix_indices()) {
            auto type          = npi.get_type();
            auto level         = npi.get_level();
            auto& damping_expr = damping_expr_group[mat_idx];

            auto active     = get_active_damping<FP>(damping_expr, level, type, begin);
            auto active_end = get_active_damping<FP>(damping_expr, level, type, end)
                                  .eval(); //copy because it may be removed or changed
            auto value = make_matrix(npi.get_value().value() * npi.get_group_weights());

            auto npi_implemented = false;

            //add begin of npi if not already bigger
            if ((active.array() < value.array()).any()) {
                damping_expr.add_damping(max(value, active), level, type, begin);
                npi_implemented = true;
            }

            //replace dampings during the new npi
            auto damping_indices = get_damping_indices<FP>(damping_expr, level, type, begin, end);
            for (auto& i : damping_indices) {
                auto& d = damping_expr.get_dampings()[i];
                damping_expr.add_damping(max(d.get_coeffs(), value), level, type, d.get_time());
                npi_implemented = true;
            }

            //add end of npi to restore active dampings if any change was made
            if (npi_implemented) {
                damping_expr.add_damping(active_end, level, type, end);
            }
        }
    }

    //remove duplicates that accumulated because of dampings that become active during the time span
    //a damping is obsolete if the most recent damping of the same type and level has the same value
    for (auto& damping_expr : damping_expr_group) {
        //go from the back so indices aren't invalidated when dampings are removed
        //use indices to loop instead of reverse iterators because removing invalidates the current iterator
        for (auto i = int(0); i < int(damping_expr.get_dampings().size()) - 1; ++i) {
            auto it = damping_expr.get_dampings().rbegin() + i;

            //look for previous damping of the same type/level
            auto it_prev = std::find_if(it + 1, damping_expr.get_dampings().rend(), [&di = *it](auto& dj) {
                return di.get_level() == dj.get_level() && di.get_type() == dj.get_type();
            });

            //remove if match is found and has same value
            if (it_prev != damping_expr.get_dampings().rend() && it->get_coeffs() == it_prev->get_coeffs()) {
                damping_expr.remove_damping(damping_expr.get_dampings().size() - 1 - i);
            }
        }
    }
}

} // namespace mio

#endif // MIO_EPI_DYNAMIC_LOCKDOWN_H
