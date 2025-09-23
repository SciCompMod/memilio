#pragma once

#include <string>
#include <map> 
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "memilio/epidemiology/damping_sampling.h"

class ControlParameter
{
public:
    template <class V>
    ControlParameter(std::string name, std::pair<double, double> allowed_range, double effectiveness_value,
                     double cost_value, mio::DampingLevel level, mio::DampingType type, const std::vector<size_t> matrices, const Eigen::MatrixBase<V>& groups);

    ControlParameter(std::string name, std::pair<double, double> allowed_range, double effectiveness_value,
                     double cost_value, mio::DampingLevel level, mio::DampingType type, const std::vector<size_t> matrices, size_t num_age_groups);
    
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
    double m_cost;
    mio::DampingSampling<double> m_damping;
};
