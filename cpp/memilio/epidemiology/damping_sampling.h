/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef EPI_SECIR_DAMPING_SAMPLING_H
#define EPI_SECIR_DAMPING_SAMPLING_H

#include "memilio/epidemiology/damping.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/uncertain_value.h"
#include <memory>

namespace mio
{

/**
 * randomly sample dampings for e.g. contact matrices.
 * All coefficients of the damping matrix depend on a single random value.
 * The damping applies to one or more of the matrices of a DampingExpressionGroup (e.g. ContactMatrixGroup).
 * The damping value is weighted by group (e.g. age) to be able to e.g. construct dampings that only
 * apply to specific groups.
 */
template <typename FP = double>
class DampingSampling
{
public:
    /**
     * Creates a DampingSampling.  
     * @param value random value that all matrix coefficients depend on.
     * @param level level of the damping.
     * @param type type of the damping.
     * @param time time of the damping.
     * @param matrices list of matrix indices that the damping applies to.
     * @param groups weights of age groups.
     */
    template <class V>
    DampingSampling(const UncertainValue<FP>& value, DampingLevel level, DampingType type, SimulationTime time,
                    const std::vector<size_t> matrices, const Eigen::MatrixBase<V>& groups)
        : m_value(value)
        , m_level(level)
        , m_type(type)
        , m_time(time)
        , m_matrices(matrices)
        , m_groups(groups)
    {
    }

    /**
     * Get the random value.
     * @return the random value.
     * @{
     */
    const UncertainValue<FP>& get_value() const
    {
        return m_value;
    }
    UncertainValue<FP>& get_value()
    {
        return m_value;
    }
    /**@}*/

    /**
     * Set the random value.
     * @param v random value.
     */
    void set_value(const UncertainValue<FP>& v)
    {
        m_value = v;
    }

    /**
     * Get the damping level.
     * @return damping level.
     */
    DampingLevel get_level() const
    {
        return m_level;
    }

    /**
     * Set the damping level.
     * @param l the damping level
     */
    void set_level(DampingLevel l)
    {
        m_level = l;
    }

    /**
     * Get the damping type.
     * @return damping type.
     */
    DampingType get_type() const
    {
        return m_type;
    }

    /**
     * Set the damping type.
     * @param t the damping type.
     */
    void set_type(DampingType t)
    {
        m_type = t;
    }

    /**
     * Get the time the damping becomes active.
     * @return the damping time.
     */
    SimulationTime get_time() const
    {
        return m_time;
    }

    /**
     * Set the time the damping becomes active.
     * @param t the damping time.
     */
    void set_time(SimulationTime t)
    {
        m_time = t;
    }

    /**
     * Get a list of matrix indices that the damping applies to.
     * The indices correspond to the indices of matrix in a DampingExpressionGroup (e.g. ContactMatrixGroup).
     * @return list of matrix indices.
     */
    const std::vector<size_t>& get_matrix_indices() const
    {
        return m_matrices;
    }

    /**
     * Set a list of matrix indices that the damping applies to.
     * @return list of matrix indices.
     */
    void set_matrix_indices(const std::vector<size_t>& v)
    {
        m_matrices = v;
    }

    /**
     * Get the group weights.
     * The groups correspond to e.g. age groups in the SECIR model.
     * @return weights of groups.
     */
    const Eigen::VectorXd& get_group_weights() const
    {
        return m_groups;
    }
    /**
     * Set the group weights.
     * @param v a vector expression of group weights.
     * @tparam V Eigen3 vector expression type.
     */
    template <class V>
    void set_group_weights(const Eigen::MatrixBase<V>& v)
    {
        m_groups = v;
    }

    /**
     * draw a value from the distribution.
     */
    void draw_sample()
    {
        m_value.draw_sample();
    }

    /**
     * equality comparison operators. 
     * @{
     */
    bool operator==(const DampingSampling& other) const
    {
        return m_value == other.m_value && m_level == other.m_level && m_type == other.m_type &&
               m_time == other.m_time && m_matrices == other.m_matrices && m_groups == other.m_groups;
    }
    bool operator!=(const DampingSampling& other) const
    {
        return !(*this == other);
    }
    /**@}*/

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("DampingSampling");
        obj.add_element("Time", get_time());
        obj.add_element("Type", get_type());
        obj.add_element("Level", get_level());
        obj.add_element("Value", get_value());
        obj.add_list("MatrixIndices", get_matrix_indices().begin(), get_matrix_indices().end());
        obj.add_element("GroupWeights", get_group_weights());
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<DampingSampling> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("DampingSampling");
        auto ti  = obj.expect_element("Time", Tag<SimulationTime>{});
        auto ty  = obj.expect_element("Type", Tag<DampingType>{});
        auto l   = obj.expect_element("Level", Tag<DampingLevel>{});
        auto v   = obj.expect_element("Value", Tag<UncertainValue<FP>>{});
        auto m   = obj.expect_list("MatrixIndices", Tag<size_t>{});
        auto g   = obj.expect_element("GroupWeights", Tag<Eigen::VectorXd>{});
        return apply(
            io,
            [](auto&& ti_, auto&& ty_, auto&& l_, auto&& v_, auto&& m_, auto&& g_) {
                return DampingSampling(v_, l_, ty_, ti_, m_, g_);
            },
            ti, ty, l, v, m, g);
    }

private:
    UncertainValue<FP> m_value;
    DampingLevel m_level;
    DampingType m_type;
    SimulationTime m_time;
    std::vector<size_t> m_matrices;
    Eigen::VectorXd m_groups;
};

/**
 * add sampled dampings to a damping expression.
 * does not draw new random value, just adds dampings.
 * @param damping_expression e.g. contact matrix group.
 * @param dampings sampled dampings.
 * @param make_mask functor that creates a matrix from damping value weighted by group.
 */
template <class DampingExpression, class DampingSamplings, class F>
void apply_dampings(DampingExpression& damping_expression, const DampingSamplings& dampings, F make_matrix)
{
    damping_expression.set_automatic_cache_update(false);
    for (auto& d : dampings) {
        for (auto& i : d.get_matrix_indices()) {
            auto m = make_matrix(double(d.get_value()) * d.get_group_weights());
            damping_expression[i].add_damping(m, d.get_level(), d.get_type(), d.get_time());
        }
    }
    damping_expression.set_automatic_cache_update(true);
}

/**
 * Make a contact damping matrix from dampings by group.
 * Maps a vector of dampings by group onto a contact damping matrix according to the formula
 * d_ij = 1 - sqrt((1 - g_i) * (1 - g_j))
 * where d_ij is a coefficient of the matrix
 * and g_i,g_j are coefficients of the group vector.
 * For diagonal elements (i.e. contacts of group with itself): d_ii = g_i; 
 * the damping of the corresponding group is applied directly. 
 * For off diagonal elements (i.e. contacts of group with other group): d_ij between g_i and g_j; 
 * the dampings of both groups are combined and applied equally.
 * @param groups damping value weighted by group.
 * @return square matrix expression of damping coefficients.
 */
template <class V>
auto make_contact_damping_matrix(V&& groups)
{
    auto ones_v = std::decay_t<V>::PlainObject::Constant(groups.size(), 1.0);
    auto prod   = (ones_v - groups) * (ones_v.transpose() - groups.transpose());
    auto ones_m = decltype(prod)::PlainObject::Constant(groups.size(), groups.size(), 1.0);
    return ones_m - sqrt(prod.array()).matrix();
}

/**
 * Make mobility coefficient damping vector from dampings by group.
 * Maps the vector of dampings by group onto a mobility coefficient damping vector
 * [g_0, g_0, ..., g_1, g_1, ..., g_2, ...].
 * @param shape shape (i.e. size) of the mobility coefficient vector.
 * @param groups damping value weighted by group.
 * @return vector expression of mobility coefficient damping.
 */
template <class V>
auto make_mobility_damping_vector(ColumnVectorShape shape, V&& groups)
{
    return Eigen::VectorXd::NullaryExpr(shape.size(), [shape, groups = std::forward<V>(groups)](Eigen::Index i) {
        auto num_groups       = groups.size();
        auto num_compartments = size_t(shape.size()) / num_groups;
        return groups[size_t(i) / num_compartments];
    });
}

} // namespace mio

#endif //EPI_SECIR_DAMPING_SAMPLING_H
