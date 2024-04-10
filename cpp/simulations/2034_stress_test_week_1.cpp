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
#include "boost/filesystem.hpp"
#include <boost/fusion/functional/invocation/invoke.hpp>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <string>

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
    const double timeExposedMin            = 3.; // 2.;
    const double timeExposedMax            = 3.; //4.;
    const double timeInfectedNoSymptomsMin = 1.;
    const double timeInfectedNoSymptomsMax = 1.;

    const double timeInfectedSymptomsMin[] = {7.0, 7.0, 7.0, 7.0, 7.0, 7.0};
    const double timeInfectedSymptomsMax[] = {7.0, 7.0, 7.0, 7.0, 7.0, 7.0};
    const double timeInfectedSevereMin[]   = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
    const double timeInfectedSevereMax[]   = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
    const double timeInfectedCriticalMin[] = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    const double timeInfectedCriticalMax[] = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0};

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
    const auto transmission_prob                       = 3 / 80.85552683191771;
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
    const double riskOfInfectionFromSymptomaticMin    = 1.0;
    const double riskOfInfectionFromSymptomaticMax    = 1.0;
    const double maxRiskOfInfectionFromSymptomaticMin = 1.0;
    const double maxRiskOfInfectionFromSymptomaticMax = 1.0;
    const double recoveredPerInfectedNoSymptomsMin[]  = {0., 0., 0., 0., 0., 0.};
    const double recoveredPerInfectedNoSymptomsMax[]  = {0., 0., 0., 0., 0., 0.};
    const double severePerInfectedSymptomsMin[]       = {0.25, 0.1, 0.1, 0.1, 0.1, 0.3};
    const double severePerInfectedSymptomsMax[]       = {0.25, 0.1, 0.1, 0.1, 0.1, 0.3};
    const double criticalPerSevereMin[]               = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    const double criticalPerSevereMax[]               = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    const double deathsPerCriticalMin[]               = {0.5, 0.5, 0.5, 0.5, 0.5, 1.};
    const double deathsPerCriticalMax[]               = {0.5, 0.5, 0.5, 0.5, 0.5, 1.};

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
    auto start_day            = mio::SimulationTime(0);

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
        mio::log_info("Fixed Damping (Masks - 20% reduction).");
        contact_dampings.push_back(physical_distancing_school(start_day, 0.2, 0.2));
        contact_dampings.push_back(physical_distancing_work(start_day, 0.2, 0.2));
        contact_dampings.push_back(physical_distancing_other(start_day, 0.2, 0.2));
        break;
    case 1:
        mio::log_info("Dynamic Damping (40% reduction).");
        dynamic_npi_dampings.push_back(physical_distancing_school(start_day, 0.4, 0.4));
        dynamic_npi_dampings.push_back(physical_distancing_work(start_day, 0.4, 0.4));
        dynamic_npi_dampings.push_back(physical_distancing_other(start_day, 0.4, 0.4));

        dynamic_npis.set_interval(mio::SimulationTime(1.0)); // how often we check if we need to activate the NPI
        dynamic_npis.set_duration(mio::SimulationTime(14.0)); // duration of the NPI
        dynamic_npis.set_base_value(100'000);
        dynamic_npis.set_threshold(1., dynamic_npi_dampings); // activation when incidence is above 1.0
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

    // print indices flow for transmission
    // mio::osecirvvs::Model modele(num_age_groups);
    // for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(6); i++) {
    //     auto indx_naive   = modele.get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptibleNaive,
    //                                                  mio::osecirvvs::InfectionState::ExposedNaive>({i});
    //     auto indx_partial = modele.get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptiblePartialImmunity,
    //                                                    mio::osecirvvs::InfectionState::ExposedPartialImmunity>({i});
    //     auto indx_ii      = modele.get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity,
    //                                               mio::osecirvvs::InfectionState::ExposedImprovedImmunity>({i});

    //     std::cout << indx_naive << "," << indx_partial << "," << indx_ii << std::endl;
    // }

    // std::exit(0);

    auto population_data_path =
        mio::path_join((data_dir / "pydata" / "Germany").string(), "county_current_population.json");
    BOOST_OUTCOME_TRY(node_ids, mio::get_node_ids(population_data_path, true));

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

/* 1: "Metropole, Regiopole und Großstädte"
   2: "Mittelstädte, städtischer Raum und kleinstädtischer, dörflicher Raum einer Ländlichen Region"
*/
static const std::map<std::string, std::vector<int>> region_mapping = {
    {"Metropolis",
     {2000, 4011, 5111,  5112,  5113,  5315,  5913,  6412,  8111,  8222,  9162,  9564, 11000, 14612, 14713,
      1002, 1003, 3101,  3102,  3103,  3403,  3404,  4012,  5114,  5116,  5117,  5119, 5120,  5122,  5124,
      5314, 5316, 5334,  5512,  5513,  5515,  5711,  5911,  5914,  5915,  5916,  6411, 6413,  6414,  6611,
      7111, 7211, 7312,  7314,  7315,  8121,  8212,  8221,  8231,  8311,  8421,  9161, 9362,  9562,  9563,
      9663, 9761, 10041, 12052, 12054, 13003, 13004, 14511, 15002, 15003, 16051, 16053}},
    {"Rural area",
     {1055,  3153,  3251,  3252,  3257,  3351,  3358,  3452,  3453,  3454,  3455,  3457,  3460,  5366,  5374,
      5554,  5570,  5770,  5958,  5966,  5974,  6437,  6531,  6532,  6533,  6534,  6631,  6632,  6634,  7316,
      7320,  8117,  8126,  8127,  8128,  8135,  8136,  8235,  8237,  8325,  8327,  8426,  8435,  8436,  9171,
      9173,  9180,  9181,  9182,  9187,  9190,  9271,  9373,  9376,  9473,  9475,  9478,  9479,  9671,  9676,
      9677,  9774,  9776,  9780,  10042, 10046, 12066, 14521, 14523, 14625, 14626, 15001, 15081, 15084, 15087,
      15090, 16062, 16064, 16065, 16066, 16070, 16072, 16073, 16076, 16077, 15082, 1051,  1054,  1058,  1059,
      1061,  3255,  3256,  3354,  3357,  3360,  3462,  5762,  6535,  6635,  6636,  7131,  7133,  7134,  7135,
      7140,  7141,  7143,  7231,  7232,  7233,  7331,  7333,  7336,  7337,  7340,  8225,  8337,  8417,  8437,
      9183,  9189,  9272,  9274,  9275,  9276,  9277,  9278,  9279,  9371,  9372,  9374,  9377,  9471,  9472,
      9476,  9477,  9571,  9575,  9577,  9672,  9673,  9674,  9678,  9773,  9777,  9778,  9779,  12062, 12068,
      12070, 12073, 13071, 13073, 13076, 16061, 16063, 16069, 16075, 13075, 15085, 15091}},
};

static const std::map<int, std::string> intervention_mapping = {
    {-1, "No_intervention"}, {0, "Fixed_Damping"}, {1, "Dynamic_NPI"}};

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
    const auto start_date                   = mio::Date(2034, 3, 1);
    constexpr auto num_days_sim             = 30.0;
    constexpr auto num_runs_per_scenario    = 5;
    constexpr auto num_runs_per_param_study = 1;
    constexpr int intervention              = 0;

    result_dir += "/" + intervention_mapping.at(intervention);

    //create or load graph
    BOOST_OUTCOME_TRY(graph_pair, get_graph(start_date, data_dir, intervention));

    auto& params_graph = graph_pair.first;
    auto& nodes        = graph_pair.second;

    auto nodes_copy = nodes;

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    std::vector<double> num_infected_range = {1.0, 10.0};

    for (const auto& region : region_mapping) {
        auto rng = mio::RandomNumberGenerator();
        rng.seed({3236549026, 3706391501, 886432438, 190527773, 3180356503, 1314521368});
        if (mio::mpi::is_root()) {
            printf("Seeds regions: ");
            for (auto s : rng.get_seeds()) {
                printf("%u, ", s);
            }
            printf("\n");
        }
        for (size_t run_per_scenario = 0; run_per_scenario < num_runs_per_scenario; run_per_scenario++) {
            int region_id =
                region
                    .second[mio::UniformIntDistribution<size_t>::get_instance()(rng, size_t(0), region.second.size())];
            for (double num_infected : num_infected_range) {

                std::string result_dir_run = result_dir + "/" + region.first;
                result_dir_run += "/" + std::to_string((int)num_infected) + "_Infected";

                result_dir_run += "/" + std::to_string(region_id);

                auto nodes_tmp = nodes_copy;

                //find correct node and set number of infected
                size_t node_indx = 0;
                for (auto& node : params_graph.nodes()) {
                    // reset population
                    node.property.populations = nodes_tmp[node_indx].populations;
                    node_indx++;
                    if (node.id == region_id) {
                        node.property
                            .populations[{mio::AgeGroup(3), mio::osecirvvs::InfectionState::InfectedSymptomsNaive}] =
                            num_infected;
                        node.property
                            .populations[{mio::AgeGroup(3), mio::osecirvvs::InfectionState::SusceptibleNaive}] -=
                            num_infected;
                    }
                }

                //run parameter study
                auto parameter_study =
                    mio::ParameterStudy<mio::osecirvvs::Simulation<mio::FlowSimulation<mio::osecirvvs::Model>>>{
                        params_graph, 0.0, num_days_sim, 0.5, num_runs_per_param_study};

                // parameter_study.get_rng().seed(
                //    {1744668429, 3100904884, 949309539, 3730340632, 1029148146, 3502301618}); //set seeds, e.g., for debugging
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

                        // auto& model_node_berlin = results_graph.nodes()[324].property.get_simulation().get_model();
                        // auto results_berlin     = interpolated_result[324];
                        // for (auto t_indx = 0; t_indx < results_berlin.get_num_time_points(); t_indx++) {
                        //     double timm_pi = 0.0;
                        //     double timm_ii = 0.0;
                        //     for (mio::AgeGroup i = 0; i < mio::AgeGroup(6); i++) {
                        //         timm_pi += results_berlin.get_value(t_indx)[model_node_berlin.populations.get_flat_index(
                        //             {i, mio::osecirvvs::InfectionState::TemporaryImmunPartialImmunity})];
                        //         timm_ii += results_berlin.get_value(t_indx)[model_node_berlin.populations.get_flat_index(
                        //             {i, mio::osecirvvs::InfectionState::TemporaryImmunImprovedImmunity})];
                        //     }
                        //     printf("t=%i, timm_pi=%f, timm_ii=%f\n", int(results_berlin.get_time(t_indx)), timm_pi, timm_ii);
                        // }

                        std::cout << "run " << run_idx << " complete." << std::endl;
                        return std::make_tuple(std::move(interpolated_result), std::move(params), std::move(flows));
                    });

                if (mio::mpi::is_root()) {
                    boost::filesystem::path dir(result_dir_run);
                    bool created = boost::filesystem::create_directories(dir);

                    if (created) {
                        mio::log_info("Directory '{:s}' was created.", dir.string());
                    }
                    printf("Saving results to \"%s\".\n", result_dir_run.c_str());
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
                        save_results(ensemble_results, ensemble_params, county_ids, result_dir_run, false));

                    auto result_dir_run_flows = result_dir_run + "/flows";
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
    }

    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::warn);
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
