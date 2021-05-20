#ifndef EPI_SECIR_DAMPING_SAMPLING_H
#define EPI_SECIR_DAMPING_SAMPLING_H

#include "epidemiology/secir/damping.h"
#include "epidemiology/utils/uncertain_value.h"
#include <memory>

namespace epi
{

/**
 * randomly sample dampings for e.g. contact matrices.
 * All coefficients of the damping matrix depend on a single random value.
 * The damping applies to one or more of the matrices of a DampingExpressionGroup (e.g. ContactMatrixGroup).
 * The damping value is weighted by group (e.g. age) to be able to e.g. construct dampings that only
 * apply to specific groups.
 */
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
    DampingSampling(const UncertainValue& value, DampingLevel level, DampingType type, SimulationTime time,
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
    const UncertainValue& get_value() const
    {
        return m_value;
    }
    UncertainValue& get_value()
    {
        return m_value;
    }
    /**@}*/

    /**
     * Set the random value.
     * @param v random value.
     */
    void set_value(const UncertainValue& v)
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
     * make the damping matrix.
     * @param make_mask function that makes a matrix from group weights.
     * @return the full damping matrix in the shape of the mask
     */
    template <class F>
    auto make_matrix(F&& make_mask) const
    {
        return double(m_value) * make_mask(m_groups);
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

private:
    UncertainValue m_value;
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
 * @param make_mask functor that creates a weight matrix from group weights.
 */
template <class DampingExpression, class DampingSamplings, class F>
void apply_dampings(DampingExpression& damping_expression, const DampingSamplings& dampings, F make_mask)
{
    for (auto& d : dampings) {
        for (auto& i : d.get_matrix_indices()) {
            auto m = d.make_matrix(make_mask);
            damping_expression[i].add_damping(m, d.get_level(), d.get_type(), d.get_time());
        }
    }
}

/**
 * make a mask for contact matrices.
 * @param groups group weights.
 * @return square matrix expression.
 */
template <class V>
auto make_contact_damping_sampling_mask(V&& groups)
{
    auto ones_v = std::decay_t<V>::PlainObject::Constant(groups.size(), 1.0);
    auto prod   = (ones_v - groups) * (ones_v.transpose() - groups.transpose());
    auto ones_m = decltype(prod)::PlainObject::Constant(groups.size(), groups.size(), 1.0);
    return ones_m - prod;
}

/**
 * make a mask for migration coefficients.
 * @param shape shape (i.e. size) of the vector.
 * @param groups group weights.
 * @return vector expression.
 */
template <class V>
auto make_migration_damping_sampling_mask(ColumnVectorShape shape, V&& groups)
{
    return Eigen::VectorXd::NullaryExpr(shape.size(), [shape, groups = std::forward<V>(groups)](Eigen::Index i) {
        auto num_groups       = groups.size();
        auto num_compartments = size_t(shape.size()) / num_groups;
        return groups[size_t(i) / num_compartments];
    });
}

} // namespace epi

#endif //EPI_SECIR_DAMPING_SAMPLING_H