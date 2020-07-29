#include "compartmental_model.h"
#include "common.h"

namespace epi
{

size_t CompartmentalModel::get_num_compartments() const
{
    return m_y.size();
}

std::vector<size_t> const& CompartmentalModel::get_category_sizes() const
{
    return m_category_sizes;
}

Eigen::VectorXd const& CompartmentalModel::get_compartments() const
{
    return m_y;
}

CompartmentalModel::CompartmentalModel(std::vector<size_t> const& category_sizes)
    : m_category_sizes(category_sizes)
{
    size_t prod = std::accumulate(category_sizes.begin(), category_sizes.end(), 1, std::multiplies<size_t>());
    m_y         = Eigen::VectorXd::Zero(prod);
}

double CompartmentalModel::get(std::vector<size_t> const& indices) const
{
    return m_y[flatten_index(indices, m_category_sizes)];
}

double CompartmentalModel::get_total_population() const
{
    return m_y.sum();
}

void CompartmentalModel::set(std::initializer_list<size_t> const& indices, double value)
{
    m_y[flatten_index(indices, m_category_sizes)] = value;
}

void CompartmentalModel::set_total_populaton(double value)
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
