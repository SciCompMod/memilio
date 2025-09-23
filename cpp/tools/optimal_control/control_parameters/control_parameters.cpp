#include "tools/optimal_control/control_parameters/control_parameters.h"

template <class V>
ControlParameter::ControlParameter(std::string name, std::pair<double, double> allowed_range,
                                   double effectiveness_value, double cost_value, mio::DampingLevel level, mio::DampingType type, 
                                   const std::vector<size_t> locations, const Eigen::MatrixBase<V>& group_weights)
    : m_name(std::move(name))
    , m_range(allowed_range)
    , m_cost(cost_value)
    , m_damping(mio::UncertainValue<double>(effectiveness_value), level, type, mio::SimulationTime<double>(0), locations, group_weights)
{
}

ControlParameter::ControlParameter(std::string name, std::pair<double, double> allowed_range,
                                   double effectiveness_value, double cost_value, mio::DampingLevel level, mio::DampingType type, 
                                   const std::vector<size_t> locations, size_t num_age_groups)
    : m_name(std::move(name))
    , m_range(allowed_range)
    , m_cost(cost_value)
    , m_damping(mio::UncertainValue<double>(effectiveness_value), level, type, mio::SimulationTime<double>(0), locations, Eigen::VectorX<double>::Constant(num_age_groups, 1.0))
{
}

const std::string& ControlParameter::name() const
{
    return m_name;
}

std::pair<double, double> ControlParameter::range() const
{
    return m_range;
}

double ControlParameter::min() const
{
    return m_range.first;
}

double ControlParameter::max() const
{
    return m_range.second;
}

double ControlParameter::effectiveness() const
{
    return m_damping.get_value().value();
}

double ControlParameter::cost() const
{
    return m_cost;
}

void ControlParameter::set_name(const std::string& new_name)
{
    m_name = new_name;
}

void ControlParameter::set_range(const std::pair<double, double>& new_range)
{
    m_range = new_range;
}

void ControlParameter::set_effectiveness(double new_effectiveness)
{
    m_damping.get_value() = new_effectiveness;
}

void ControlParameter::set_cost(double new_cost)
{
    m_cost = new_cost;
}