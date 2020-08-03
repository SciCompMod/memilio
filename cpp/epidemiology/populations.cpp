#include "populations.h"
#include "tensor_helpers.h"

#include <numeric>

namespace epi
{

Populations::Populations(std::vector<size_t> const& category_sizes)
    : m_category_sizes(category_sizes)
{
    size_t prod = std::accumulate(category_sizes.begin(), category_sizes.end(), 1, std::multiplies<size_t>());
    m_y         = Eigen::VectorXd::Zero(prod);
}

size_t Populations::get_num_compartments() const
{
    return m_y.size();
}

std::vector<size_t> const& Populations::get_category_sizes() const
{
    return m_category_sizes;
}

Eigen::VectorXd const& Populations::get_compartments() const
{
    return m_y;
}

double Populations::get(std::vector<size_t> const& indices) const
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

double Populations::get_group_population(size_t category_idx, size_t group_idx) const
{
    //TODO maybe implement an iterator/visitor pattern rather than calculating indices?

    double sum   = 0;
    auto indices = get_slice_indices(category_idx, group_idx, m_category_sizes);
    for (auto i : indices) {
        sum += m_y[i];
    }
    return sum;
}

void Populations::set_group_population(size_t category_idx, size_t group_idx, double value)
{
    //TODO slice indices are calcualated twice...
    double current_population = get_group_population(category_idx, group_idx);
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

double Populations::get_total() const
{
    return m_y.sum();
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

void Populations::set_total(double value)
{
    double current_population = get_total();
    if (fabs(current_population) < 1e-12) {
        m_y.fill(value / m_y.size());
    }
    else {
        m_y *= value / current_population;
    }
}

void Populations::set_remaining(std::vector<size_t> const& indices, size_t total_population)
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

} // namespace epi
