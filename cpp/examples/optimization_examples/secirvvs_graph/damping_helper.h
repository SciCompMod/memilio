#pragma once

#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/damping.h"

#include <vector>
#include <functional>

enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Count
};
enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Holidays,
    Count
};
enum class ContactLocation
{
    Home,
    School,
    Work,
    Other,
    Count
};

template <typename FP>
mio::DampingSampling<FP> set_school_closure(mio::SimulationTime<FP> time, FP value, size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::Main));
    const mio::DampingType type(static_cast<size_t>(Intervention::SchoolClosure));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::School)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, FP(1.0));

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}

template <typename FP>
mio::DampingSampling<FP> set_home_office(mio::SimulationTime<FP> time, FP value, size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::Main));
    const mio::DampingType type(static_cast<size_t>(Intervention::HomeOffice));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::Work)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, FP(1.0));

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}

template <typename FP>
mio::DampingSampling<FP> set_physical_distancing_school(mio::SimulationTime<FP> time, FP value, size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::PhysicalDistanceAndMasks));
    const mio::DampingType type(static_cast<size_t>(Intervention::PhysicalDistanceAndMasks));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::School)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, FP(1.0));

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}

template <typename FP>
mio::DampingSampling<FP> set_physical_distancing_work(mio::SimulationTime<FP> time, FP value, size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::PhysicalDistanceAndMasks));
    const mio::DampingType type(static_cast<size_t>(Intervention::PhysicalDistanceAndMasks));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::Work)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, FP(1.0));

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}

template <typename FP>
mio::DampingSampling<FP> set_physical_distancing_other(mio::SimulationTime<FP> time, FP value, size_t num_age_groups)
{
    const mio::UncertainValue<FP> damping_value(value);
    const mio::DampingLevel level(static_cast<size_t>(InterventionLevel::PhysicalDistanceAndMasks));
    const mio::DampingType type(static_cast<size_t>(Intervention::PhysicalDistanceAndMasks));
    const std::vector<size_t> locations    = {static_cast<size_t>(ContactLocation::Other)};
    const Eigen::VectorX<FP> group_weights = Eigen::VectorX<FP>::Constant(num_age_groups, FP(1.0));

    return mio::DampingSampling<FP>(damping_value, level, type, time, locations, group_weights);
}