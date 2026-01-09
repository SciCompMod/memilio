#pragma once

#include "memilio/epidemiology/damping_sampling.h"
#include "memilio/math/eigen.h"

#include <string>
#include <utility>
#include <vector>

class DynamicNPIControlParameter
{
public:
    DynamicNPIControlParameter(const std::string& name, double threshold,
                               const std::vector<std::string>& dampings_names,
                               const std::vector<mio::DampingSampling<double>>& dampings);

    const std::string& name() const;
    double threshold() const;
    const std::vector<mio::DampingSampling<double>>& dampings() const;
    const std::vector<std::string>& damping_names() const;

    template <class FP>
    std::vector<mio::DampingSampling<FP>> get_dampings() const;

private:
    std::string m_name;
    double m_threshold;
    std::vector<std::string> m_damping_names;
    std::vector<mio::DampingSampling<double>> m_dampings;
};

template <class FP>
std::vector<mio::DampingSampling<FP>> DynamicNPIControlParameter::get_dampings() const
{
    std::vector<mio::DampingSampling<FP>> result;
    result.reserve(m_dampings.size());

    for (const auto& d : m_dampings) {
        result.emplace_back(mio::UncertainValue<FP>(static_cast<FP>(d.get_value())), mio::DampingLevel(d.get_level()),
                            mio::DampingType(d.get_type()), mio::SimulationTime<FP>(d.get_time().get()),
                            d.get_matrix_indices(), d.get_group_weights().template cast<FP>());
    }

    return result;
}
