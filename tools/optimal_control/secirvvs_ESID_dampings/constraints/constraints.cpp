#include "constraints.h"

Constraint::Constraint(std::string name, std::pair<double, double> allowed_range, std::optional<size_t> node_id)
    : m_name(std::move(name))
    , m_range(allowed_range)
    , m_node_id(node_id)
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

std::optional<size_t> Constraint::node_id() const
{
    return m_node_id;
}

void Constraint::set_name(const std::string& new_name)
{
    m_name = new_name;
}

void Constraint::set_range(const std::pair<double, double>& new_range)
{
    m_range = new_range;
}

void Constraint::set_node_id(std::optional<size_t> new_node_id)
{
    m_node_id = new_node_id;
}
