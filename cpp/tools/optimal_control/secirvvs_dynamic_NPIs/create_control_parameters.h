#pragma once

#include "boost/filesystem.hpp"

#include "memilio/epidemiology/damping.h"

#include "control_parameters/control_parameters.h"

#include <vector>
#include <string>
#include <set>

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

std::vector<DynamicNPIControlParameter> create_control_parameters(size_t num_age_groups)
{
    std::vector<DynamicNPIControlParameter> control_parameters;

    std::vector<mio::DampingSampling<double>> dampings_A;

    // contacts_at_home
    dampings_A.emplace_back(mio::UncertainValue<double>(1.0),
                            mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
                            mio::DampingType(static_cast<size_t>(Intervention::Home)), mio::SimulationTime<double>(0.0),
                            std::vector<size_t>{}, Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // school closure
    dampings_A.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::SchoolClosure)), mio::SimulationTime<double>(0.0),
        std::vector<size_t>{static_cast<size_t>(ContactLocation::School)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // social_events
    dampings_A.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::GatheringBanFacilitiesClosure)),
        mio::SimulationTime<double>(0.0), std::vector<size_t>{static_cast<size_t>(ContactLocation::Other)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // social_events_work
    dampings_A.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::GatheringBanFacilitiesClosure)),
        mio::SimulationTime<double>(0.0), std::vector<size_t>{static_cast<size_t>(ContactLocation::Work)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // physical_distancing_work
    dampings_A.emplace_back(mio::UncertainValue<double>(1.0),
                            mio::DampingLevel(static_cast<size_t>(InterventionLevel::PhysicalDistanceAndMasks)),
                            mio::DampingType(static_cast<size_t>(Intervention::GatheringBanFacilitiesClosure)),
                            mio::SimulationTime<double>(0.0),
                            std::vector<size_t>{static_cast<size_t>(ContactLocation::Work)},
                            Eigen::VectorX<double>::Constant(num_age_groups, 1.0));

    control_parameters.emplace_back("DynamicNPI A", 50.0,
                                    std::vector<std::string>{"contacts_at_home", "school closure", "social_events",
                                                             "social_events_work", "physical_distancing_work"},
                                    dampings_A);

    // // ------------ //
    // // Dampings 200 //
    // // ------------ //
    // std::vector<mio::DampingSampling<double>> dampings_B;

    // // contacts_at_home
    // dampings_B.emplace_back(mio::UncertainValue<double>(1.0),
    //                         mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
    //                         mio::DampingType(static_cast<size_t>(Intervention::Home)), mio::SimulationTime<double>(0.0),
    //                         std::vector<size_t>{static_cast<size_t>(ContactLocation::Home)},
    //                         Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // // school closures
    // dampings_B.emplace_back(
    //     mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
    //     mio::DampingType(static_cast<size_t>(Intervention::SchoolClosure)), mio::SimulationTime<double>(0.0),
    //     std::vector<size_t>{static_cast<size_t>(ContactLocation::School)},
    //     Eigen::VectorX<double>::Constant(num_age_groups, 1.0));

    // control_parameters.emplace_back("DynamicNPI B", 200.0,
    //                                 std::vector<std::string>{"contacts_at_home", "school closure"}, dampings_B);

    std::vector<mio::DampingSampling<double>> dampings_B;

    // contacts_at_home
    dampings_B.emplace_back(mio::UncertainValue<double>(1.0),
                            mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
                            mio::DampingType(static_cast<size_t>(Intervention::Home)), mio::SimulationTime<double>(0.0),
                            std::vector<size_t>{}, Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // school closure
    dampings_B.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::SchoolClosure)), mio::SimulationTime<double>(0.0),
        std::vector<size_t>{static_cast<size_t>(ContactLocation::School)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // social_events
    dampings_B.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::GatheringBanFacilitiesClosure)),
        mio::SimulationTime<double>(0.0), std::vector<size_t>{static_cast<size_t>(ContactLocation::Other)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // social_events_work
    dampings_B.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::GatheringBanFacilitiesClosure)),
        mio::SimulationTime<double>(0.0), std::vector<size_t>{static_cast<size_t>(ContactLocation::Work)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // physical_distancing_work
    dampings_B.emplace_back(mio::UncertainValue<double>(1.0),
                            mio::DampingLevel(static_cast<size_t>(InterventionLevel::PhysicalDistanceAndMasks)),
                            mio::DampingType(static_cast<size_t>(Intervention::GatheringBanFacilitiesClosure)),
                            mio::SimulationTime<double>(0.0),
                            std::vector<size_t>{static_cast<size_t>(ContactLocation::Work)},
                            Eigen::VectorX<double>::Constant(num_age_groups, 1.0));

    control_parameters.emplace_back("DynamicNPI B", 200.0,
                                    std::vector<std::string>{"contacts_at_home", "school closure", "social_events",
                                                             "social_events_work", "physical_distancing_work"},
                                    dampings_B);

    std::vector<mio::DampingSampling<double>> dampings_C;

    // contacts_at_home
    dampings_C.emplace_back(mio::UncertainValue<double>(1.0),
                            mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
                            mio::DampingType(static_cast<size_t>(Intervention::Home)), mio::SimulationTime<double>(0.0),
                            std::vector<size_t>{}, Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // school closure
    dampings_C.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::SchoolClosure)), mio::SimulationTime<double>(0.0),
        std::vector<size_t>{static_cast<size_t>(ContactLocation::School)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // social_events
    dampings_C.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::GatheringBanFacilitiesClosure)),
        mio::SimulationTime<double>(0.0), std::vector<size_t>{static_cast<size_t>(ContactLocation::Other)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));
    // social_events_work
    dampings_C.emplace_back(
        mio::UncertainValue<double>(1.0), mio::DampingLevel(static_cast<size_t>(InterventionLevel::Main)),
        mio::DampingType(static_cast<size_t>(Intervention::GatheringBanFacilitiesClosure)),
        mio::SimulationTime<double>(0.0), std::vector<size_t>{static_cast<size_t>(ContactLocation::Work)},
        Eigen::VectorX<double>::Constant(num_age_groups, 1.0));

    // physical_distancing_work
    dampings_C.emplace_back(mio::UncertainValue<double>(1.0),
                            mio::DampingLevel(static_cast<size_t>(InterventionLevel::PhysicalDistanceAndMasks)),
                            mio::DampingType(static_cast<size_t>(Intervention::GatheringBanFacilitiesClosure)),
                            mio::SimulationTime<double>(0.0),
                            std::vector<size_t>{static_cast<size_t>(ContactLocation::Work)},
                            Eigen::VectorX<double>::Constant(num_age_groups, 1.0));

    control_parameters.emplace_back("DynamicNPI C", 500.0,
                                    std::vector<std::string>{"contacts_at_home", "school closure", "social_events",
                                                             "social_events_work", "physical_distancing_work"},
                                    dampings_C);

    return control_parameters;
}
