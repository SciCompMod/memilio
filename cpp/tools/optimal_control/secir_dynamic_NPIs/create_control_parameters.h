#pragma once

#include "boost/filesystem.hpp"

#include "memilio/epidemiology/damping.h"

#include "control_parameters/control_parameters.h"

#include <vector>
#include <string>
#include <set>

enum class Intervention
{
    PhysicalDistanceAndMasks = 0,
    Count                    = 1
};

enum class InterventionLevel
{
    Main  = 0,
    Count = 1
};

enum class ContactLocation
{
    Work  = 0,
    Count = 1
};

std::vector<DynamicNPIControlParameter> create_control_parameters(size_t num_age_groups)
{
    std::vector<DynamicNPIControlParameter> control_parameters;

    std::vector<mio::DampingSampling<double>> dampings_A;
    // physical_distancing_work
    dampings_A.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::PhysicalDistanceAndMasks)), mio::SimulationTime<double>(0.0),
        std::vector<size_t>{static_cast<size_t>(ContactLocation::Work)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));

    control_parameters.emplace_back("DynamicNPI A", 50.0, std::vector<std::string>{"physical_distancing_work"},
                                    dampings_A);

    std::vector<mio::DampingSampling<double>> dampings_B;
    // physical_distancing_work
    dampings_B.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::PhysicalDistanceAndMasks)), mio::SimulationTime<double>(0.0),
        std::vector<size_t>{static_cast<size_t>(ContactLocation::Work)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));

    control_parameters.emplace_back("DynamicNPI B", 250.0, std::vector<std::string>{"physical_distancing_work"},
                                    dampings_B);

    return control_parameters;
}
