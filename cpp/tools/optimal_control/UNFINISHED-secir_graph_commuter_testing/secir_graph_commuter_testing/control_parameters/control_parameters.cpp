#include "control_parameters.h"

ControlParameter::ControlParameter(std::string name, std::pair<double, double> allowed_range,
                                   double effectiveness_value, double cost_value)
    : m_name(std::move(name))
    , m_range(allowed_range)
    , m_effectiveness(effectiveness_value)
    , m_cost(cost_value)
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
    return m_effectiveness;
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
    m_effectiveness = new_effectiveness;
}

void ControlParameter::set_cost(double new_cost)
{
    m_cost = new_cost;
}
