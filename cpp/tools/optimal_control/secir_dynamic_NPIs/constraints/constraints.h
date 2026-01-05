#pragma once

#include <string>
#include <utility>
#include <optional>

class Constraint
{
public:
    Constraint(std::string name, std::pair<double, double> allowed_range, std::optional<size_t> node_id);

    const std::string& name() const;
    std::pair<double, double> range() const;
    double min() const;
    double max() const;
    std::optional<size_t> node_id() const;

    void set_name(const std::string& new_name);
    void set_range(const std::pair<double, double>& new_range);
    void set_node_id(std::optional<size_t> new_node_id);

private:
    std::string m_name;
    std::pair<double, double> m_range;
    std::optional<size_t> m_node_id;
};
