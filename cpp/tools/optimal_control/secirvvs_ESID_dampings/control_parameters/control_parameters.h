#pragma once

#include "memilio/math/eigen.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/damping_sampling.h"

#include <utility>
#include <vector>
#include <stdexcept>

class DampingControlParameter
{
public:
    DampingControlParameter(const std::string& name, std::pair<double, double> allowed_range,
                            double effectiveness_value, double cost_value, mio::DampingLevel level,
                            mio::DampingType type, std::vector<size_t> locations, const Eigen::VectorXd& group_weights);

    std::string name() const;
    std::pair<double, double> range() const;
    double min() const;
    double max() const;
    double effectiveness() const;
    double cost() const;
    const mio::DampingLevel& level() const;
    const mio::DampingType& type() const;
    const std::vector<size_t>& locations() const;
    const Eigen::VectorX<double>& group_weights() const;

    void set_name(const std::string& new_name);
    void set_range(const std::pair<double, double>& new_range);
    void set_effectiveness(double new_effectiveness);
    void set_cost(double new_cost);
    void set_level(const mio::DampingLevel& new_level);
    void set_type(const mio::DampingType& new_type);
    void set_locations(const std::vector<size_t>& new_locations);
    void set_group_weights(const Eigen::VectorX<double>& new_group_weights);

    template <class FP>
    mio::DampingSampling<FP> make_damping_sampling(const mio::UncertainValue<FP>& strength,
                                                   mio::SimulationTime<FP> time) const;

private:
    std::string m_name;
    std::pair<double, double> m_range;
    double m_effectiveness;
    double m_cost;
    mio::DampingLevel m_level;
    mio::DampingType m_type;
    std::vector<size_t> m_locations;
    Eigen::VectorX<double> m_group_weights;
};

template <class FP>
mio::DampingSampling<FP> DampingControlParameter::make_damping_sampling(const mio::UncertainValue<FP>& strength,
                                                                        mio::SimulationTime<FP> time) const
{
    return mio::DampingSampling<FP>(m_effectiveness * strength, m_level, m_type, time, m_locations,
                                    m_group_weights.cast<FP>());
}
