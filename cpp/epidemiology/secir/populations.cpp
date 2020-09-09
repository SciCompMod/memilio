#include "epidemiology/secir/populations.h"
#include "epidemiology/utils/tensor_helpers.h"

#include <numeric>

namespace epi
{

Populations::Populations(std::vector<size_t> const& category_sizes)
    : m_category_sizes(category_sizes)
{
    size_t prod = std::accumulate(category_sizes.begin(), category_sizes.end(), 1, std::multiplies<size_t>());
    m_y         = std::vector<UncertainValue>(prod, 0.0);
}

size_t Populations::get_num_compartments() const
{
    return m_y.size();
}

std::vector<size_t> const& Populations::get_category_sizes() const
{
    return m_category_sizes;
}

Eigen::VectorXd Populations::get_compartments() const
{
    Eigen::VectorXd m_y_eigen(m_y.size());
    for (auto i = 0; i < m_y.size(); i++) {
        m_y_eigen[i] = m_y[i];
    }

    return m_y_eigen;
}

UncertainValue const& Populations::get(std::vector<size_t> const& indices) const
{
#ifndef NDEBUG
    assert(indices.size() == m_category_sizes.size());
    size_t i = 0;
    for (auto idx : indices) {
        assert(idx < m_category_sizes[i++]);
    }
#endif
    return m_y[flatten_index(indices, m_category_sizes)];
}

UncertainValue& Populations::get(std::vector<size_t> const& indices)
{
#ifndef NDEBUG
    assert(indices.size() == m_category_sizes.size());
    size_t i = 0;
    for (auto idx : indices) {
        assert(idx < m_category_sizes[i++]);
    }
#endif
    return m_y[flatten_index(indices, m_category_sizes)];
}

double Populations::get_group_total(size_t category_idx, size_t group_idx) const
{
    //TODO maybe implement an iterator/visitor pattern rather than calculating indices?

    double sum   = 0;
    auto indices = get_slice_indices(category_idx, group_idx, m_category_sizes);
    for (auto i : indices) {
        sum += m_y[i];
    }
    return sum;
}

void Populations::set_group_total(size_t category_idx, size_t group_idx, double value)
{
    //TODO slice indices are calcualated twice...
    double current_population = get_group_total(category_idx, group_idx);
    auto indices              = get_slice_indices(category_idx, group_idx, m_category_sizes);

    if (fabs(current_population) < 1e-12) {
        for (auto i : indices) {
            m_y[i] = value / indices.size();
        }
    }
    else {
        for (auto i : indices) {
            m_y[i] *= value / current_population;
        }
    }
}

void Populations::set_difference_from_group_total(std::vector<size_t> const& indices, size_t category_idx,
                                                  size_t group_idx, double total_group_population)
{
    // is the given index part of the group?
    assert(indices[category_idx] == group_idx);

    double current_population = get_group_total(category_idx, group_idx);
    size_t idx                = get_flat_index(indices);
    current_population -= m_y[idx];

    assert(current_population <= total_group_population);

    m_y[idx] = total_group_population - current_population;
}

double Populations::get_total() const
{
    double sum = 0;
    for (auto i = 0; i < m_y.size(); i++) {
        sum += m_y[i];
    }
    return sum;
}

void Populations::set(std::vector<size_t> const& indices, UncertainValue const& value)
{
#ifndef NDEBUG
    assert(indices.size() == m_category_sizes.size());
    size_t i = 0;
    for (auto idx : indices) {
        assert(idx < m_category_sizes[i++]);
    }
#endif
    m_y[get_flat_index(indices)] = value;
}

void Populations::set(std::vector<size_t> const& indices, double value)
{
#ifndef NDEBUG
    assert(indices.size() == m_category_sizes.size());
    size_t i = 0;
    for (auto idx : indices) {
        assert(idx < m_category_sizes[i++]);
    }
#endif
    m_y[get_flat_index(indices)] = value;
}

void Populations::set(std::vector<size_t> const& indices, ParameterDistribution const& dist)
{
#ifndef NDEBUG
    assert(indices.size() == m_category_sizes.size());
    size_t i = 0;
    for (auto idx : indices) {
        assert(idx < m_category_sizes[i++]);
    }
#endif
    m_y[get_flat_index(indices)].set_distribution(dist);
}

void Populations::set_total(double value)
{
    double current_population = get_total();
    if (fabs(current_population) < 1e-12) {
        double ysize = m_y.size();
        for (auto i = 0; i < m_y.size(); i++) {
            m_y[i] = value / ysize;
        }
    }
    else {
        for (auto i = 0; i < m_y.size(); i++) {
            m_y[i] *= value / current_population;
        }
    }
}

void Populations::set_difference_from_total(std::vector<size_t> const& indices, double total_population)
{
    double current_population = get_total();
    size_t idx                = get_flat_index(indices);
    current_population -= m_y[idx];

    assert(current_population <= total_population);

    m_y[idx] = total_population - current_population;
}

size_t Populations::get_flat_index(std::vector<size_t> const& indices) const
{
    return flatten_index(indices, m_category_sizes);
}

void Populations::check_constraints()
{
    for (auto i = 0; i < m_y.size(); i++) {
        if (m_y[i] < 0) {
            log_warning("Constraint check: Compartment size {:d} changed from {:.4f} to {:d}", i, m_y[i], 0);
            m_y[i] = 0;
        }
    }
}

} // namespace epi
