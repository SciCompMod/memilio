#pragma once

#include <string>
#include <map> 
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

class ControlParameter
{
public:
    ControlParameter(std::string name, std::pair<double, double> allowed_range, double effectiveness_value,
                     double cost_value);

    const std::string& name() const;
    std::pair<double, double> range() const;
    double min() const;
    double max() const;
    double effectiveness() const;
    double cost() const;

    void set_name(const std::string& new_name);
    void set_range(const std::pair<double, double>& new_range);
    void set_effectiveness(double new_effectiveness);
    void set_cost(double new_cost);

private:
    std::string m_name;
    std::pair<double, double> m_range;
    double m_effectiveness;
    double m_cost;
};

enum class ControlType
{
    SchoolClosure,
    HomeOffice,
    PhysicalDistancingSchool,
    PhysicalDistancingWork,
    PhysicalDistancingOther,
    Count
};

inline ControlType string_to_control(const std::string& control_name)
{
    static const std::map<std::string, ControlType> control_pairs = {
        {"SchoolClosure", ControlType::SchoolClosure},
        {"HomeOffice", ControlType::HomeOffice},
        {"PhysicalDistancingSchool", ControlType::PhysicalDistancingSchool},
        {"PhysicalDistancingWork", ControlType::PhysicalDistancingWork},
        {"PhysicalDistancingOther", ControlType::PhysicalDistancingOther},
    };

    try
    {
        return control_pairs.at(control_name);
    } catch(const std::out_of_range& ex)
    {
        throw std::runtime_error("Invalid control name: " + control_name);
    }
}
