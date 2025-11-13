#include "control_parameters.h"

DampingControlParameter::DampingControlParameter(const std::string& name, std::pair<double, double> allowed_range,
                                                 double effectiveness_value, double cost_value, mio::DampingLevel level,
                                                 mio::DampingType type, std::vector<size_t> locations,
                                                 const Eigen::VectorXd& group_weights)
    : m_name(name)
    , m_range(std::move(allowed_range))
    , m_effectiveness(effectiveness_value)
    , m_cost(cost_value)
    , m_level(std::move(level))
    , m_type(std::move(type))
    , m_locations(std::move(locations))
    , m_group_weights(group_weights)
{
    if (m_range.first > m_range.second) {
        throw std::invalid_argument("Invalid range: min must not exceed max.");
    }
    if (m_group_weights.size() == 0) {
        throw std::invalid_argument("Group weights vector must not be empty.");
    }
}

std::string DampingControlParameter::name() const
{
    return m_name;
}

std::pair<double, double> DampingControlParameter::range() const
{
    return m_range;
}

double DampingControlParameter::min() const
{
    return m_range.first;
}

double DampingControlParameter::max() const
{
    return m_range.second;
}

double DampingControlParameter::effectiveness() const
{
    return m_effectiveness;
}

double DampingControlParameter::cost() const
{
    return m_cost;
}

const mio::DampingLevel& DampingControlParameter::level() const
{
    return m_level;
}

const mio::DampingType& DampingControlParameter::type() const
{
    return m_type;
}

const std::vector<size_t>& DampingControlParameter::locations() const
{
    return m_locations;
}

const Eigen::VectorXd& DampingControlParameter::group_weights() const
{
    return m_group_weights;
}

void DampingControlParameter::set_name(const std::string& new_name)
{
    m_name = new_name;
}

void DampingControlParameter::set_range(const std::pair<double, double>& new_range)
{
    if (new_range.first > new_range.second) {
        throw std::invalid_argument("Invalid range: min must not exceed max.");
    }
    m_range = new_range;
}

void DampingControlParameter::set_effectiveness(double new_effectiveness)
{
    m_effectiveness = new_effectiveness;
}

void DampingControlParameter::set_cost(double new_cost)
{
    m_cost = new_cost;
}

void DampingControlParameter::set_level(const mio::DampingLevel& new_level)
{
    m_level = new_level;
}

void DampingControlParameter::set_type(const mio::DampingType& new_type)
{
    m_type = new_type;
}

void DampingControlParameter::set_locations(const std::vector<size_t>& new_locations)
{
    m_locations = new_locations;
}

void DampingControlParameter::set_group_weights(const Eigen::VectorXd& new_group_weights)
{
    if (new_group_weights.size() == 0) {
        throw std::invalid_argument("Group weights vector must not be empty.");
    }
    m_group_weights = new_group_weights;
}
