#pragma once

#include <string>
#include <utility>

class Constraint
{
public:
    Constraint(std::string name, std::pair<double, double> allowed_range);

    const std::string& name() const;
    std::pair<double, double> range() const;
    double min() const;
    double max() const;

    void set_name(const std::string& new_name);
    void set_range(const std::pair<double, double>& new_range);

private:
    std::string m_name;
    std::pair<double, double> m_range;
};
