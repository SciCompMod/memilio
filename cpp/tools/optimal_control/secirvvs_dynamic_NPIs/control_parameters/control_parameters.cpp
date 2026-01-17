#include "control_parameters.h"

DynamicNPIControlParameter::DynamicNPIControlParameter(const std::string& name, double threshold,
                                                       const std::vector<std::string>& damping_names,
                                                       const std::vector<mio::DampingSampling<double>>& dampings)
    : m_name(name)
    , m_threshold(threshold)
    , m_damping_names(damping_names)
    , m_dampings(dampings)
{
    if (m_damping_names.size() != m_dampings.size()) {
        throw std::invalid_argument("Number of damping names must equal number of dampings.");
    }
}

const std::string& DynamicNPIControlParameter::name() const
{
    return m_name;
}

double DynamicNPIControlParameter::threshold() const
{
    return m_threshold;
}

const std::vector<mio::DampingSampling<double>>& DynamicNPIControlParameter::dampings() const
{
    return m_dampings;
}

const std::vector<std::string>& DynamicNPIControlParameter::damping_names() const
{
    return m_damping_names;
}
