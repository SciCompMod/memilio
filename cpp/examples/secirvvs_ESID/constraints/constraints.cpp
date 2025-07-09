#include "constraints.h"

Constraint::Constraint(std::string name, std::pair<double, double> allowed_range)
    : m_name(std::move(name))
    , m_range(allowed_range)
{
}

const std::string& Constraint::name() const
{
    return m_name;
}

std::pair<double, double> Constraint::range() const
{
    return m_range;
}

double Constraint::min() const
{
    return m_range.first;
}

double Constraint::max() const
{
    return m_range.second;
}

void Constraint::set_name(const std::string& new_name)
{
    m_name = new_name;
}

void Constraint::set_range(const std::pair<double, double>& new_range)
{
    m_range = new_range;
}
