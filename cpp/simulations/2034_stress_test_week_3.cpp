/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: 
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameters.h"
#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/parameter_space.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/stl_util.h"
#include "memilio/utils/uncertain_value.h"
#include "boost/filesystem.hpp"
#include <boost/fusion/functional/invocation/invoke.hpp>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace fs = boost::filesystem;

/**
 * indices of contact matrix corresponding to locations where contacts occur.
 */
enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

/**
 * different types of NPI, used as DampingType.
 */
enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
};

/**
 * different level of NPI, used as DampingLevel.
 */
enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Holidays,
};

/**
 * Set a value and distribution of an UncertainValue.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution.
 * @param p uncertain value to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void assign_uniform_distribution(mio::UncertainValue& p, double min, double max)
{
    p = mio::UncertainValue(0.5 * (max + min));
    p.set_distribution(mio::ParameterDistributionUniform(min, max));
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min[i] and max[i] as a value and UNIFORM(min[i], max[i]) as a distribution for
 * each element i of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution for each element of array.
 * @param max minimum of distribution for each element of array.
 */
template <size_t N>
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N])
{
    assert(N == array.numel());
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(N); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)]);
    }
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution to every element of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array, double min,
                                       double max)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

/**
 * Set epidemiological parameters of Covid19.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
//TODO: set parameters according to information from Andre
mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters& params)
{
    //times
    const double timeExposedMin            = 2.;
    const double timeExposedMax            = 4.;
    const double timeInfectedNoSymptomsMin = 1.;
    const double timeInfectedNoSymptomsMax = 1.;

    const double timeInfectedSymptomsMin[] = {7.0, 7.0, 7.0, 7.0, 7.0, 7.0};
    const double timeInfectedSymptomsMax[] = {7.0, 7.0, 7.0, 7.0, 7.0, 7.0};
    //52-86% are discharged from the hospital within 2-12 days (median 5) -> including time in ICU
    //8-34% of these are in ICU where they spend 1-28 days (we set median here to 10)
    //Calculation for median of time infected severe: 5 - 0.2 (percentage ICU) * 10 (median days ICU) = 3
    const double timeInfectedSevereMin[]   = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    const double timeInfectedSevereMax[]   = {6.0, 6.0, 6.0, 6.0, 6.0, 6.0};
    const double timeInfectedCriticalMin[] = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
    const double timeInfectedCriticalMax[] = {12.0, 12.0, 12.0, 12.0, 12.0, 12.0};

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeExposed>(), timeExposedMin, timeExposedMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedNoSymptoms>(), timeInfectedNoSymptomsMin,
                                      timeInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSymptoms>(), timeInfectedSymptomsMin,
                                      timeInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSevere>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedCritical>(), timeInfectedCriticalMin,
                                      timeInfectedCriticalMax);

    //probabilities
    double fac_variant                                 = 1.0;
    const auto transmission_prob                       = 3.5 / 80.85552683191771;
    const double transmissionProbabilityOnContactMin[] = {
        transmission_prob * fac_variant, transmission_prob * fac_variant, transmission_prob * fac_variant,
        transmission_prob * fac_variant, transmission_prob * fac_variant, transmission_prob * fac_variant};

    const double transmissionProbabilityOnContactMax[] = {
        transmission_prob * fac_variant, transmission_prob * fac_variant, transmission_prob * fac_variant,
        transmission_prob * fac_variant, transmission_prob * fac_variant, transmission_prob * fac_variant};
    const double relativeTransmissionNoSymptomsMin = 1.0;
    const double relativeTransmissionNoSymptomsMax = 1.0;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    const double riskOfInfectionFromSymptomaticMin    = 0.6; //Cd
    const double riskOfInfectionFromSymptomaticMax    = 0.8; //Cd
    const double maxRiskOfInfectionFromSymptomaticMin = 0.6;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.8;
    const double recoveredPerInfectedNoSymptomsMin[]  = {0., 0., 0., 0., 0., 0.};
    const double recoveredPerInfectedNoSymptomsMax[]  = {0., 0., 0., 0., 0., 0.};
    const double severePerInfectedSymptomsMin[]       = {0.25, 0.1, 0.1, 0.1, 0.1, 0.3};
    const double severePerInfectedSymptomsMax[]       = {0.25, 0.1, 0.1, 0.1, 0.1, 0.3};
    const double criticalPerSevereMin[]               = {0.18, 0.18, 0.18, 0.18, 0.18, 0.18};
    const double criticalPerSevereMax[]               = {0.24, 0.24, 0.24, 0.24, 0.24, 0.24};
    const double deathsPerCriticalMin[]               = {0., 0., 0., 0., 0., 0.};
    const double deathsPerCriticalMax[]               = {0., 0., 0., 0., 0., 0.};

    const double reducExposedPartialImmunityMin                     = 1.0;
    const double reducExposedPartialImmunityMax                     = 1.0;
    const double reducExposedImprovedImmunityMin                    = 1.0;
    const double reducExposedImprovedImmunityMax                    = 1.0;
    const double reducInfectedSymptomsPartialImmunityMin            = 1.0;
    const double reducInfectedSymptomsPartialImmunityMax            = 1.0;
    const double reducInfectedSymptomsImprovedImmunityMin           = 1.0;
    const double reducInfectedSymptomsImprovedImmunityMax           = 1.0;
    const double reducInfectedSevereCriticalDeadPartialImmunityMin  = 1.0;
    const double reducInfectedSevereCriticalDeadPartialImmunityMax  = 1.0;
    const double reducInfectedSevereCriticalDeadImprovedImmunityMin = 1.0;
    const double reducInfectedSevereCriticalDeadImprovedImmunityMax = 1.0;

    const double reducTimeInfectedMild = 1.0;

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TransmissionProbabilityOnContact>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RelativeTransmissionNoSymptoms>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SeverePerInfectedSymptoms>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::CriticalPerSevere>(), criticalPerSevereMin,
                                      criticalPerSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::DeathsPerCritical>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax);

    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedPartialImmunity>(),
                                      reducExposedPartialImmunityMin, reducExposedPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedImprovedImmunity>(),
                                      reducExposedImprovedImmunityMin, reducExposedImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>(),
                                      reducInfectedSymptomsPartialImmunityMin, reducInfectedSymptomsPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>(),
                                      reducInfectedSymptomsImprovedImmunityMin,
                                      reducInfectedSymptomsImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>(),
                                      reducInfectedSevereCriticalDeadPartialImmunityMin,
                                      reducInfectedSevereCriticalDeadPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>(),
                                      reducInfectedSevereCriticalDeadImprovedImmunityMin,
                                      reducInfectedSevereCriticalDeadImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducTimeInfectedMild>(), reducTimeInfectedMild,
                                      reducTimeInfectedMild);

    //sasonality
    const double seasonality_min = 0.2; //0.1;
    const double seasonality_max = 0.2; //0.3;

    assign_uniform_distribution(params.get<mio::osecirvvs::Seasonality>(), seasonality_min, seasonality_max);

    // no waning of immunity
    params.get<mio::osecirvvs::TimeWaningPartialImmunity>()  = std::numeric_limits<double>::max();
    params.get<mio::osecirvvs::TimeWaningImprovedImmunity>() = std::numeric_limits<double>::max();
    params.get<mio::osecirvvs::TimeTemporaryImmunityPI>()    = std::numeric_limits<double>::max();
    params.get<mio::osecirvvs::TimeTemporaryImmunityII>()    = std::numeric_limits<double>::max();

    // set vaccinations to zero
    params.get<mio::osecirvvs::DailyPartialVaccination>().resize(mio::SimulationDay(1000));
    params.get<mio::osecirvvs::DailyFullVaccination>().resize(mio::SimulationDay(1000));
    params.get<mio::osecirvvs::DailyBoosterVaccination>().resize(mio::SimulationDay(1000));
    for (size_t i = 0; i < 1000; ++i) {
        for (auto age_group = mio::AgeGroup(0); age_group < mio::AgeGroup(6); ++age_group) {
            params.get<mio::osecirvvs::DailyPartialVaccination>()[{(mio::AgeGroup)age_group, mio::SimulationDay(i)}] =
                0.;
            params.get<mio::osecirvvs::DailyFullVaccination>()[{(mio::AgeGroup)age_group, mio::SimulationDay(i)}] = 0.;
            params.get<mio::osecirvvs::DailyBoosterVaccination>()[{(mio::AgeGroup)age_group, mio::SimulationDay(i)}] =
                0.;
        }
    }

    return mio::success();
}

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

/**
 * Set contact matrices.
 * Reads contact matrices from files in the data directory.
 * @param data_dir data directory.
 * @param params Object that the contact matrices will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirvvs::Parameters& params)
{
    //TODO: io error handling
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);
    }
    params.get<mio::osecirvvs::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrices);

    return mio::success();
}

/**
 * @brief Sets the NPIs for a model based on the intervention parameter.
 *
 * @param params A reference to the parameters of the simulation. This object is modified by the function.
 * @param intervention An integer that determines the type of npi to apply.
 *                     If intervention is 0, a fixed damping is applied with a 20% reduction.
 *                     If intervention is 1, a dynamic damping is applied with a 60% reduction.
 *                     If intervention is any other value, no intervention is applied.
 */
mio::IOResult<void> set_interventions(mio::osecirvvs::Parameters& params, const int intervention)
{
    auto group_weights_all    = Eigen::VectorXd::Constant(size_t(params.get_num_groups()), 1.0);
    auto& contacts            = params.get<mio::osecirvvs::ContactPatterns>();
    auto& contact_dampings    = contacts.get_dampings();
    auto& dynamic_npis        = params.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms>();
    auto dynamic_npi_dampings = std::vector<mio::DampingSampling>();
    auto start_day            = mio::SimulationTime(14); //Cd

    auto physical_distancing_school = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                    mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                    {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto physical_distancing_work = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                    mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                    {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto physical_distancing_other = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                    mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                    {size_t(ContactLocation::Other)}, group_weights_all);
    };

    switch (intervention) {
    case 0:
        mio::log_info("Fixed Damping (20% reduction).");
        contact_dampings.push_back(physical_distancing_school(start_day, 0.2, 0.2));
        contact_dampings.push_back(physical_distancing_work(start_day, 0.2, 0.2));
        contact_dampings.push_back(physical_distancing_other(start_day, 0.2, 0.2));
        break;
    case 1:
        mio::log_info("Dynamic Damping (60% reduction).");
        dynamic_npi_dampings.push_back(physical_distancing_school(start_day, 0.4, 0.4));
        dynamic_npi_dampings.push_back(physical_distancing_work(start_day, 0.4, 0.4));
        dynamic_npi_dampings.push_back(physical_distancing_other(start_day, 0.4, 0.4));

        dynamic_npis.set_interval(mio::SimulationTime(1.0)); // how often we check if we need to activate the NPI
        dynamic_npis.set_duration(mio::SimulationTime(14.0)); // duration of the NPI
        dynamic_npis.set_base_value(100'000);
        dynamic_npis.set_threshold(1., dynamic_npi_dampings); // activation when incidence is above 1.0
        break;
    case 2:
        mio::log_info("Fixed Damping (40% reduction).");
        contact_dampings.push_back(physical_distancing_school(start_day, 0.4, 0.4));
        contact_dampings.push_back(physical_distancing_work(start_day, 0.4, 0.4));
        contact_dampings.push_back(physical_distancing_other(start_day, 0.4, 0.4));
        break;
    case 3:
        mio::log_info("Fixed Damping (60% reduction).");
        contact_dampings.push_back(physical_distancing_school(start_day, 0.6, 0.6));
        contact_dampings.push_back(physical_distancing_work(start_day, 0.6, 0.6));
        contact_dampings.push_back(physical_distancing_other(start_day, 0.6, 0.6));
        break;
    default:
        mio::log_info("No intervention is applied.");
    }

    return mio::success();
}

//TODO: start initialization?
mio::IOResult<void> set_population(std::vector<mio::osecirvvs::Model>& nodes, const std::string& path_pop,
                                   const std::vector<int>& node_ids)
{
    BOOST_OUTCOME_TRY(num_population, mio::osecirvvs::details::read_population_data(path_pop, node_ids));
    const auto num_groups = nodes[0].parameters.get_num_groups();

    assert(node_ids.size() == nodes.size());

    for (size_t i = 0; i < node_ids.size(); ++i) {
        auto& population = nodes[i].populations;
        for (auto age_group = mio::AgeGroup(0); age_group < num_groups; ++age_group) {
            population[{age_group, mio::osecirvvs::InfectionState::SusceptibleNaive}] =
                num_population[i][(size_t)age_group];
        }
        // assert that size from population and num_population is equal. So no other compartments are > 0.
        assert(population.get_total() == std::accumulate(num_population[i].begin(), num_population[i].end(), 0.0));
    }

    return mio::success();
}

std::vector<int> node_ids_global;

mio::IOResult<
    std::pair<mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>, std::vector<mio::osecirvvs::Model>>>
get_graph(mio::Date start_date, const fs::path& data_dir, const int intervention)
{
    // global parameters
    const int num_age_groups = 6;
    mio::osecirvvs::Parameters params(num_age_groups);
    params.get<mio::osecirvvs::StartDay>() = mio::get_day_in_year(start_date);
    params.get_end_dynamic_npis()          = 32;
    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));
    BOOST_OUTCOME_TRY(set_interventions(params, intervention));

    auto population_data_path =
        mio::path_join((data_dir / "pydata" / "Germany").string(), "county_current_population.json");
    BOOST_OUTCOME_TRY(node_ids, mio::get_node_ids(population_data_path, true));

    node_ids_global = node_ids;

    std::vector<mio::osecirvvs::Model> nodes(node_ids.size(),
                                             mio::osecirvvs::Model(int(size_t(params.get_num_groups()))));
    for (auto& node : nodes) {
        node.parameters = params;
    }

    BOOST_OUTCOME_TRY(set_population(nodes, population_data_path, node_ids));

    mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> params_graph;
    for (size_t i = 0; i < node_ids.size(); ++i) {
        params_graph.add_node(node_ids[i], nodes[i]);
    }

    auto migrating_compartments     = {mio::osecirvvs::InfectionState::SusceptibleNaive,
                                       mio::osecirvvs::InfectionState::ExposedNaive,
                                       mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive,
                                       mio::osecirvvs::InfectionState::InfectedSymptomsNaive,
                                       mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity,
                                       mio::osecirvvs::InfectionState::SusceptiblePartialImmunity,
                                       mio::osecirvvs::InfectionState::ExposedPartialImmunity,
                                       mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity,
                                       mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity,
                                       mio::osecirvvs::InfectionState::ExposedImprovedImmunity,
                                       mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity,
                                       mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity,
                                       mio::osecirvvs::InfectionState::TemporaryImmunPartialImmunity,
                                       mio::osecirvvs::InfectionState::TemporaryImmunImprovedImmunity};
    const auto& read_function_edges = mio::read_mobility_plain;
    const auto& set_edge_function =
        mio::set_edges<ContactLocation, mio::osecirvvs::Model, mio::MigrationParameters, mio::MigrationCoefficientGroup,
                       mio::osecirvvs::InfectionState, decltype(read_function_edges)>;
    BOOST_OUTCOME_TRY(set_edge_function(data_dir, params_graph, migrating_compartments, contact_locations.size(),
                                        read_function_edges, std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.}));

    return mio::success(std::pair(params_graph, nodes));
}

/**
 * Different modes for running the parameter study.
 */
enum class RunMode
{
    Load,
    Save,
};

const std::vector<int> region_ids;

static const std::map<int, std::string> intervention_mapping = {
    {-1, "No_intervention"}, {0, "20p_reduc"}, {1, "Dynamic_NPI"}, {2, "40p_reduc"}, {3, "60p_reduc"}};

/**
 * Run the parameter study.
 * Load a previously stored graph or create a new one from data.
 * The graph is the input for the parameter study.
 * A newly created graph is saved and can be reused.
 * @param mode Mode for running the parameter study.
 * @param data_dir data directory. Not used if mode is RunMode::Load.
 * @param save_dir directory where the graph is loaded from if mode is RunMOde::Load or save to if mode is RunMode::Save.
 * @param result_dir directory where all results of the parameter study will be stored.
 * @param save_single_runs [Default: true] Defines if single run results are written to the disk.
 * @returns any io error that occurs during reading or writing of files.
 */
mio::IOResult<void> run(const fs::path& data_dir, std::string result_dir)
{
    const auto start_date                   = mio::Date(2034, 5, 1);
    constexpr auto num_days_sim             = 44.0;
    constexpr auto num_runs_per_scenario    = 5;
    constexpr auto num_runs_per_param_study = 50;
    std::vector<int> interventions{-1, 0, 2, 3};

    for (int intervention : interventions) {

        auto result_dir_run_copy = result_dir + "/" + intervention_mapping.at(intervention);

        //create or load graph
        BOOST_OUTCOME_TRY(graph_pair, get_graph(start_date, data_dir, intervention));

        auto& params_graph = graph_pair.first;
        auto& nodes        = graph_pair.second;

        auto nodes_copy = nodes;

        std::vector<int> county_ids(params_graph.nodes().size());
        std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
            return n.id;
        });

        double initially_infected = 10.0 / 100'000.0;
        auto rng                  = mio::RandomNumberGenerator();
        rng.seed({3236549026, 3706391501, 886432438, 190527773, 3180356503, 1314521368});
        if (mio::mpi::is_root()) {
            printf("Seeds regions: ");
            for (auto s : rng.get_seeds()) {
                printf("%u, ", s);
            }
            printf("\n");
        }
        auto result_dir_run = result_dir_run_copy;

        // 20% of all counties have infected individuals
        size_t num_regions = static_cast<size_t>(0.2 * params_graph.nodes().size());

        for (size_t run_per_scenario = 0; run_per_scenario < num_runs_per_scenario; run_per_scenario++) {

            size_t node_indx = 0;
            for (auto& node : params_graph.nodes()) {
                node.property.populations = nodes_copy[node_indx].populations;
                node_indx++;
            }

            auto result_dir_scenario = result_dir_run + "/run_" + std::to_string(run_per_scenario);

            std::vector<int> hotspot_ids;
            //draw regions
            for (size_t region = 0; region < num_regions; ++region) {
                int region_id = node_ids_global[mio::UniformIntDistribution<size_t>::get_instance()(
                    rng, size_t(0), node_ids_global.size())];
                bool new_id   = std::find(hotspot_ids.begin(), hotspot_ids.end(), region_id) == hotspot_ids.end();
                while (!new_id) {
                    region_id = node_ids_global[mio::UniformIntDistribution<size_t>::get_instance()(
                        rng, size_t(0), node_ids_global.size())];
                    new_id    = std::find(hotspot_ids.begin(), hotspot_ids.end(), region_id) == hotspot_ids.end();
                }

                hotspot_ids.push_back(region_id);

                //find correct node and set number of infected
                node_indx = 0;
                for (auto& node : params_graph.nodes()) {
                    node_indx++;
                    if (node.id == region_id) {
                        auto percentage_symptomatic =
                            node.property.parameters.get<mio::osecirvvs::TimeInfectedSymptoms>()[{mio::AgeGroup(3)}] /
                            (node.property.parameters.get<mio::osecirvvs::TimeInfectedSymptoms>()[{mio::AgeGroup(3)}] +
                             node.property.parameters
                                 .get<mio::osecirvvs::TimeInfectedNoSymptoms>()[{mio::AgeGroup(3)}] +
                             node.property.parameters.get<mio::osecirvvs::TimeExposed>()[{mio::AgeGroup(3)}]);
                        auto percentage_exposed =
                            node.property.parameters.get<mio::osecirvvs::TimeExposed>()[{mio::AgeGroup(3)}] /
                            (node.property.parameters.get<mio::osecirvvs::TimeInfectedSymptoms>()[{mio::AgeGroup(3)}] +
                             node.property.parameters
                                 .get<mio::osecirvvs::TimeInfectedNoSymptoms>()[{mio::AgeGroup(3)}] +
                             node.property.parameters.get<mio::osecirvvs::TimeExposed>()[{mio::AgeGroup(3)}]);
                        auto infected = initially_infected * node.property.populations.get_total();
                        node.property
                            .populations[{mio::AgeGroup(3), mio::osecirvvs::InfectionState::InfectedSymptomsNaive}] =

                            percentage_symptomatic * infected;
                        node.property
                            .populations[{mio::AgeGroup(3), mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}] =
                            (1 - percentage_exposed - percentage_symptomatic) * infected;

                        node.property.populations[{mio::AgeGroup(3), mio::osecirvvs::InfectionState::ExposedNaive}] =
                            percentage_exposed * infected;

                        node.property
                            .populations[{mio::AgeGroup(3), mio::osecirvvs::InfectionState::SusceptibleNaive}] -=
                            infected;
                    }
                }
            }

            //run parameter study
            auto parameter_study =
                mio::ParameterStudy<mio::osecirvvs::Simulation<mio::FlowSimulation<mio::osecirvvs::Model>>>{
                    params_graph, 0.0, num_days_sim, 0.5, num_runs_per_param_study};

            parameter_study.get_rng().seed({1744668429, 3100904884, 949309539, 3730340632, 1029148146,
                                            3502301618}); //set seeds, e.g., for debugging
            if (mio::mpi::is_root()) {
                printf("Seeds parameter study: ");
                for (auto s : parameter_study.get_rng().get_seeds()) {
                    printf("%u, ", s);
                }
                printf("\n");
            }

            auto save_single_run_result = mio::IOResult<void>(mio::success());
            auto ensemble               = parameter_study.run(
                [&](auto&& graph) {
                    return draw_sample(graph);
                },
                [&](auto results_graph, auto&& run_idx) {
                    auto interpolated_result = mio::interpolate_simulation_result(results_graph);
                    auto params              = std::vector<mio::osecirvvs::Model>();
                    params.reserve(results_graph.nodes().size());
                    std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                                                 std::back_inserter(params), [](auto&& node) {
                                       return node.property.get_simulation().get_model();
                                   });

                    auto flows = std::vector<mio::TimeSeries<ScalarType>>{};
                    flows.reserve(results_graph.nodes().size());
                    std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                                                 std::back_inserter(flows), [](auto&& node) {
                                       auto& flow_node         = node.property.get_simulation().get_flows();
                                       auto interpolated_flows = mio::interpolate_simulation_result(flow_node);
                                       return interpolated_flows;
                                   });

                    std::cout << "run " << run_idx << " complete." << std::endl;
                    return std::make_tuple(std::move(interpolated_result), std::move(params), std::move(flows));
                });

            if (mio::mpi::is_root()) {
                boost::filesystem::path dir(result_dir_scenario);
                bool created = boost::filesystem::create_directories(dir);

                if (created) {
                    mio::log_info("Directory '{:s}' was created.", dir.string());
                }
                printf("Saving results to \"%s\".\n", result_dir_scenario.c_str());
            }

            if (ensemble.size() > 0) {
                auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
                ensemble_results.reserve(ensemble.size());
                auto ensemble_params = std::vector<std::vector<mio::osecirvvs::Model>>{};
                ensemble_params.reserve(ensemble.size());
                auto ensemble_flows = std::vector<std::vector<mio::TimeSeries<double>>>{};
                ensemble_flows.reserve(ensemble.size());
                for (auto&& run : ensemble) {
                    ensemble_results.emplace_back(std::move(std::get<0>(run)));
                    ensemble_params.emplace_back(std::move(std::get<1>(run)));
                    ensemble_flows.emplace_back(std::move(std::get<2>(run)));
                }

                BOOST_OUTCOME_TRY(save_single_run_result);
                BOOST_OUTCOME_TRY(
                    save_results(ensemble_results, ensemble_params, county_ids, result_dir_scenario, false));

                auto result_dir_run_flows = result_dir_scenario + "/flows";
                if (mio::mpi::is_root()) {
                    boost::filesystem::path dir(result_dir_run_flows);
                    bool created = boost::filesystem::create_directories(dir);

                    if (created) {
                        mio::log_info("Directory '{:s}' was created.", dir.string());
                    }
                    printf("Saving Flow results to \"%s\".\n", result_dir_run_flows.c_str());
                }
                BOOST_OUTCOME_TRY(
                    save_results(ensemble_flows, ensemble_params, county_ids, result_dir_run_flows, false));
            }
        }
    }

    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::err);
    mio::mpi::init();

    std::string save_dir;
    std::string data_dir;
    std::string result_dir;
    if (argc == 4) {
        data_dir   = argv[1];
        save_dir   = argv[2];
        result_dir = argv[3];
    }

    else if (argc == 1) {
        data_dir   = "/localdata1/test/memilio/data";
        save_dir   = "/localdata1/test/memilio/test";
        result_dir = "/localdata1/test/memilio/test";
    }
    else {
        if (mio::mpi::is_root()) {
            printf("Given arguments are not valid");
        }
        mio::mpi::finalize();
        return 0;
    }

    if (mio::mpi::is_root()) {
        boost::filesystem::path dir(result_dir);
        bool created = boost::filesystem::create_directories(dir);

        if (created) {
            mio::log_info("Directory '{:s}' was created.", dir.string());
        }
        printf("Saving results to \"%s\".\n", result_dir.c_str());
    }

    auto result = run(data_dir, result_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        mio::mpi::finalize();
        return -1;
    }
    mio::mpi::finalize();
    return 0;
}
