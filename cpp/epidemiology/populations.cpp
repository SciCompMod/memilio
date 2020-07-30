#include "populations.h".h "
#include "common.h"

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

double Populations::get(std::initializer_list<size_t> const& indices) const
{
    return m_y[flatten_index(indices, m_category_sizes)];
}

double Populations::get_group_population(size_t category_idx, size_t group_idx) const
{
    return 0;
}

double Populations::get_total_population() const
{
    return m_y.sum();
}

void Populations::set(std::initializer_list<size_t> const& indices, double value)
{
    m_y[flatten_index(indices, m_category_sizes)] = value;
}

void Populations::set_total_populaton(double value)
{
    double current_population = get_total_population();
    if (fabs(current_population) < 1e-12) {
        m_y.fill(value / m_y.size());
    }
    else {
        m_y *= value / current_population;
    }
}

} // namespace epi
