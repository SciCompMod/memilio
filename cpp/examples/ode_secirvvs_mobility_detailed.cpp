/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Henrik Zunker
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
#include "memilio/config.h"
#include "memilio/data/analyze_result.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/logging.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/parameters.h"
#include "ode_secirvvs/parameters_io.h"
#include "memilio/mobility/metapopulation_mobility_detailed.h"
#include "memilio/compartments/simulation.h"
#include "memilio/io/result_io.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"

#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/mobility_io.h"
#include "memilio/utils/date.h"
#include "ode_secirvvs/parameter_space.h"
#include "memilio/utils/stl_util.h"
#include "boost/filesystem.hpp"
#include <cstddef>
#include <cstdio>
#include <iomanip>

#include <iostream>
#include <string>
#include <vector>

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
mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters& params)
{
    //times
    const double incubationTime          = 3.1; // doi.org/10.3201/eid2806.220158
    const double serialIntervalMin       = 2.38; // doi.org/10.1016/j.lanepe.2022.100446
    const double serialIntervalMax       = 2.38; // doi.org/10.1016/j.lanepe.2022.100446
    const double timeInfectedSymptomsMin = 6.58; //https://doi.org/10.1016/S0140-6736(22)00327-0
    const double timeInfectedSymptomsMax = 7.16; //https://doi.org/10.1016/S0140-6736(22)00327-0
    const double timeInfectedSevereMin[] = {1.8, 1.8, 1.8, 2.5, 3.5, 4.91}; // doi.org/10.1186/s12879-022-07971-6
    const double timeInfectedSevereMax[] = {2.3, 2.3, 2.3, 3.67, 5, 7.01}; // doi.org/10.1186/s12879-022-07971-6

    const ScalarType fact_icu              = 1.;
    const double timeInfectedCriticalMin[] = {fact_icu * 4.95,  fact_icu * 4.95, fact_icu * 4.86,
                                              fact_icu * 14.14, fact_icu * 14.4, fact_icu * 10.};
    const double timeInfectedCriticalMax[] = {fact_icu * 8.95,  fact_icu * 8.95, fact_icu * 8.86,
                                              fact_icu * 20.58, fact_icu * 19.8, fact_icu * 13.2};

    array_assign_uniform_distribution(params.get<mio::osecirvvs::IncubationTime>(), incubationTime, incubationTime);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SerialInterval>(), serialIntervalMin,
                                      serialIntervalMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSymptoms>(), timeInfectedSymptomsMin,
                                      timeInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSevere>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedCritical>(), timeInfectedCriticalMin,
                                      timeInfectedCriticalMax);

    double fac_variant                                 = 1.94; //https://doi.org/10.7554/eLife.78933
    const double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                          0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

    const double transmissionProbabilityOnContactMax[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                                          0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};

    const double relativeTransmissionNoSymptomsMin    = 0.5;
    const double relativeTransmissionNoSymptomsMax    = 0.5;
    const double riskOfInfectionFromSymptomaticMin    = 0.0; // beta (depends on incidence and test and trace capacity)
    const double riskOfInfectionFromSymptomaticMax    = 0.2;
    const double maxRiskOfInfectionFromSymptomaticMin = 0.4;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.5;

    const double recoveredPerInfectedNoSymptomsMin[] = {0.2, 0.25,  0.2,
                                                        0.2, 0.175, 0.1}; // doi.org/10.1101/2022.05.05.22274697
    const double recoveredPerInfectedNoSymptomsMax[] = {0.4, 0.45, 0.35,
                                                        0.3, 0.25, 0.15}; // doi.org/10.1101/2022.05.05.22274697

    // doi:10.1136/bmjgh-2023-0123

    // Factors from https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00462-7/fulltext
    const double severePerInfectedSymptomsMin[] = {1 * 0.006,   0.8 * 0.006, 0.4 * 0.015,
                                                   0.3 * 0.049, 0.25 * 0.15, 0.35 * 0.2};
    const double severePerInfectedSymptomsMax[] = {1 * 0.009,   0.8 * 0.009, 0.4 * 0.023, 0.3 * 0.074,
                                                   0.25 * 0.18, 0.35 * 0.25}; // quelle 2021 paper + factors

    // delta paper
    // risk of icu admission https://doi.org/10.1177/14034948221108548
    const double fac_icu                = 0.52;
    const double criticalPerSevereMin[] = {0.05 * fac_icu, 0.05 * fac_icu, 0.05 * fac_icu,
                                           0.10 * fac_icu, 0.25 * fac_icu, 0.35 * fac_icu};
    const double criticalPerSevereMax[] = {0.10 * fac_icu, 0.10 * fac_icu, 0.10 * fac_icu,
                                           0.20 * fac_icu, 0.35 * fac_icu, 0.45 * fac_icu};

    // 61% less risk doi:10.1136/bmjgh-2023-0123
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10347449/pdf/bmjgh-2023-012328.pdf
    const double fac_dead               = 0.39;
    const double deathsPerCriticalMin[] = {fac_dead * 0.00, fac_dead * 0.00, fac_dead * 0.10,
                                           fac_dead * 0.10, fac_dead * 0.30, fac_dead * 0.5}; // 2021 paper
    const double deathsPerCriticalMax[] = {fac_dead * 0.10, fac_dead * 0.10, fac_dead * 0.18,
                                           fac_dead * 0.18, fac_dead * 0.50, fac_dead * 0.7};

    const double reducExposedPartialImmunityMin = 0.569; // doi.org/10.1136/bmj-2022-071502
    const double reducExposedPartialImmunityMax = 0.637; // doi.org/10.1136/bmj-2022-071502
    const double reducExposedImprovedImmunityMin =
        0.34 * reducExposedPartialImmunityMin; // https://jamanetwork.com/journals/jama/fullarticle/2788487 0.19346
    const double reducExposedImprovedImmunityMax =
        0.34 * reducExposedPartialImmunityMax; // https://jamanetwork.com/journals/jama/fullarticle/2788487 0.21658

    const double reducInfectedSymptomsPartialImmunityMin  = 0.746; // doi.org/10.1056/NEJMoa2119451
    const double reducInfectedSymptomsPartialImmunityMax  = 0.961; // doi.org/10.1056/NEJMoa2119451
    const double reducInfectedSymptomsImprovedImmunityMin = 0.295; // doi.org/10.1056/NEJMoa2119451
    const double reducInfectedSymptomsImprovedImmunityMax = 0.344; // doi.org/10.1056/NEJMoa2119451

    const double reducInfectedSevereCriticalDeadPartialImmunityMin =
        0.52; // www.assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf
    const double reducInfectedSevereCriticalDeadPartialImmunityMax =
        0.82; // www.assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf
    const double reducInfectedSevereCriticalDeadImprovedImmunityMin = 0.1; // doi.org/10.1136/bmj-2022-071502
    const double reducInfectedSevereCriticalDeadImprovedImmunityMax = 0.19; // doi.org/10.1136/bmj-2022-071502
    const double temp_reducTimeInfectedMild                         = 0.5; // doi.org/10.1101/2021.09.24.21263978

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
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducTimeInfectedMild>(), temp_reducTimeInfectedMild,
                                      temp_reducTimeInfectedMild);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecirvvs::Seasonality>(), seasonality_min, seasonality_max);
    return mio::success();
}

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
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirvvs::Parameters& params,
                                         ScalarType scale_contacts = 1.0, ScalarType share_staying = 1.0)
{
    auto contact_transport_status = mio::read_mobility_plain(data_dir.string() + "//contacts//contacts_transport.txt");
    auto contact_matrix_transport = contact_transport_status.value();
    auto contact_matrices         = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        if (contact_location.second == "other") {
            contact_matrices[size_t(contact_location.first)].get_baseline() =
                (1 - share_staying) * abs(baseline - contact_matrix_transport) / scale_contacts +
                share_staying * abs(baseline - contact_matrix_transport);
            contact_matrices[size_t(contact_location.first)].get_minimum() = Eigen::MatrixXd::Zero(6, 6);
        }
        else {
            contact_matrices[size_t(contact_location.first)].get_baseline() =
                (1 - share_staying) * baseline / scale_contacts + share_staying * baseline;
            contact_matrices[size_t(contact_location.first)].get_minimum() = Eigen::MatrixXd::Zero(6, 6);
        }
    }
    params.get<mio::osecirvvs::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrices);

    return mio::success();
}

mio::IOResult<void> set_contact_matrices_transport(const fs::path& data_dir, mio::osecirvvs::Parameters& params,
                                                   ScalarType scale_contacts = 1.0)
{
    auto contact_transport_status = mio::read_mobility_plain(data_dir.string() + "//contacts//contacts_transport.txt");
    auto contact_matrix_transport = contact_transport_status.value();
    auto contact_matrices         = mio::ContactMatrixGroup(1, size_t(params.get_num_groups()));
    // ScalarType const polymod_share_contacts_transport = 1 / 0.2770885028949545;

    contact_matrices[0].get_baseline() = contact_matrix_transport / scale_contacts;
    contact_matrices[0].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);

    params.get<mio::osecirvvs::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrices);

    return mio::success();
}

// reset population in graph
void init_pop_cologne_szenario(mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>& graph)
{
    std::vector<std::vector<double>> immunity = {{0.04, 0.61, 0.35}, {0.04, 0.61, 0.35},   {0.075, 0.62, 0.305},
                                                 {0.08, 0.62, 0.3},  {0.035, 0.58, 0.385}, {0.01, 0.41, 0.58}};

    for (auto& node : graph.nodes()) {
        for (auto age = mio::AgeGroup(0); age < mio::AgeGroup(6); age++) {
            auto pop_age = 0.0;
            for (auto inf_state = mio::Index<mio::osecirvvs::InfectionState>(0);
                 inf_state < mio::osecirvvs::InfectionState::Count; ++inf_state) {
                pop_age += node.property.populations[{age, inf_state}];
                node.property.populations[{age, inf_state}] = 0.0;
            }
            size_t immunity_index = static_cast<size_t>(age);
            node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptibleNaive}] =
                pop_age * immunity[immunity_index][0];
            node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}] =
                pop_age * immunity[immunity_index][1];
            node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}] =
                pop_age * immunity[immunity_index][2];
        }
        if (node.id == 5315) {
            // infect p% of population
            ScalarType p = 0.05;
            for (auto age = mio::AgeGroup(0); age < graph.nodes()[0].property.parameters.get_num_groups(); age++) {
                node.property.populations[{mio::AgeGroup(age), mio::osecirvvs::InfectionState::InfectedSymptomsNaive}] =
                    node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptibleNaive}] * p;
                node.property.populations[{age, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}] =
                    node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}] * p;
                node.property.populations[{mio::AgeGroup(age),
                                           mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}] =
                    node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}] * p;
                node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptibleNaive}] *= 0.99;
                node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}] *= 0.99;
                node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}] *= 0.99;
            }
        }
    }
}

/**
 * Create the input graph for the parameter study.
 * Reads files from the data directory.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @returns created graph or any io errors that happen during reading of the files.
 */
mio::IOResult<mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>>
get_graph(const int num_days, bool masks, bool ffp2, bool szenario_cologne, bool edges)
{
    std::string data_dir         = "/localdata1/test/memilio//data";
    std::string traveltimes_path = "/localdata1/code/memilio/test/travel_times_pathes.txt";
    std::string durations_path   = "/localdata1/code/memilio/test/activity_duration_work.txt";
    auto mobility_data_commuter =
        mio::read_mobility_plain(("/localdata1/code/memilio/test/commuter_migration_with_locals.txt"));
    auto mob_data   = mobility_data_commuter.value();
    auto start_date = mio::Date(2022, 8, 1);
    auto end_date   = mio::Date(2022, 11, 1);

    // global parameters
    const int num_age_groups = 6;
    mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> params_graph;

    // Nodes
    mio::osecirvvs::Parameters params(num_age_groups);

    params.get<mio::osecirvvs::StartDay>() = mio::get_day_in_year(start_date);
    auto params_status                     = set_covid_parameters(params);
    auto contacts_status                   = set_contact_matrices(data_dir, params);

    params.get<mio::osecirvvs::StartDay>() = mio::get_day_in_year(start_date);

    // Set nodes
    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;
    auto tnt_capacity_factor     = 1.43 / 100000.;

    auto read_duration = mio::read_duration_stay(durations_path);
    if (!read_duration) {
        std::cout << read_duration.error().formatted_message() << '\n';
    }
    auto duration_stay = read_duration.value();

    std::string path = "/localdata1/test/memilio/data/pydata/Germany/county_population.json";
    auto read_ids    = mio::get_node_ids(path, true);
    auto node_ids    = read_ids.value();

    std::vector<mio::osecirvvs::Model> nodes(node_ids.size(),
                                             mio::osecirvvs::Model(int(size_t(params.get_num_groups()))));
    std::vector<double> scaling_factor_inf(6, 1.0);

    for (auto& node : nodes) {
        node.parameters = params;
    }
    auto read_node = read_input_data_county(nodes, start_date, node_ids, scaling_factor_inf, scaling_factor_icu,
                                            "/localdata1/test/memilio/data", num_days);

    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {

        auto tnt_capacity = nodes[node_idx].populations.get_total() * tnt_capacity_factor;

        //local parameters
        auto& tnt_value = nodes[node_idx].parameters.template get<mio::osecirvvs::TestAndTraceCapacity>();
        tnt_value       = mio::UncertainValue(0.5 * (1.2 * tnt_capacity + 0.8 * tnt_capacity));
        tnt_value.set_distribution(mio::ParameterDistributionUniform(0.8 * tnt_capacity, 1.2 * tnt_capacity));

        //holiday periods
        auto id              = int(mio::regions::CountyId(node_ids[node_idx]));
        auto holiday_periods = mio::regions::get_holidays(mio::regions::get_state_id(id), start_date, end_date);
        auto& contacts       = nodes[node_idx].parameters.template get<mio::osecirvvs::ContactPatterns>();
        contacts.get_school_holidays() =
            std::vector<std::pair<mio::SimulationTime, mio::SimulationTime>>(holiday_periods.size());
        std::transform(
            holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(), [=](auto& period) {
                return std::make_pair(mio::SimulationTime(mio::get_offset_in_days(period.first, start_date)),
                                      mio::SimulationTime(mio::get_offset_in_days(period.second, start_date)));
            });

        //uncertainty in populations
        for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = mio::Index<typename mio::osecirvvs::Model::Compartments>(0);
                 j < mio::osecirvvs::Model::Compartments::Count; ++j) {
                auto& compartment_value = nodes[node_idx].populations[{i, j}];
                compartment_value =
                    mio::UncertainValue(0.5 * (1.1 * double(compartment_value) + 0.9 * double(compartment_value)));
                compartment_value.set_distribution(mio::ParameterDistributionUniform(0.9 * double(compartment_value),
                                                                                     1.1 * double(compartment_value)));
            }
        }

        // Add mobility node
        auto mobility = nodes[node_idx];
        mobility.populations.set_total(0);

        // reduce transmission on contact due to mask obligation in mobility node
        // first age group not able to  (properly) wear masks
        if (masks) {

            const double fact_surgical_mask = 0.1;
            const double fact_ffp2          = 0.001;

            double factor_mask[] = {1, fact_surgical_mask, fact_ffp2, fact_ffp2, fact_ffp2, fact_ffp2};
            if (!ffp2) {
                std::cout << "surgical masks"
                          << "\n";
                // multipliziere alle außer den ersten EIntrag mit 40
                for (size_t j = 1; j < 6; j++) {
                    factor_mask[j] = factor_mask[j] * fact_surgical_mask / fact_ffp2;
                }
            }

            double fac_variant                           = 1.4; //1.94; //https://doi.org/10.7554/eLife.78933
            double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                            0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

            double transmissionProbabilityOnContactMax[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                                            0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};
            for (int i = 0; i < 6; i++) {
                transmissionProbabilityOnContactMin[i] = transmissionProbabilityOnContactMin[i] * factor_mask[i];
                transmissionProbabilityOnContactMax[i] = transmissionProbabilityOnContactMax[i] * factor_mask[i];
            }
            array_assign_uniform_distribution(
                mobility.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>(),
                transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
        }

        for (size_t t_idx = 0; t_idx < num_days; ++t_idx) {
            auto t = mio::SimulationDay((size_t)t_idx);
            for (auto j = mio::AgeGroup(0); j < params.get_num_groups(); j++) {
                mobility.parameters.template get<mio::osecirvvs::DailyFirstVaccination>()[{j, t}] = 0;
                mobility.parameters.template get<mio::osecirvvs::DailyFullVaccination>()[{j, t}]  = 0;
            }
        }

        auto& params_mobility = mobility.parameters;

        params_graph.add_node(node_ids[node_idx], duration_stay((Eigen::Index)node_idx), nodes[node_idx], mobility);
    }

    if (szenario_cologne) {
        init_pop_cologne_szenario(params_graph);
    }

    // Edges
    if (edges) {
        ScalarType theshold_edges   = 4e-5;
        auto migrating_compartments = {
            mio::osecirvvs::InfectionState::SusceptibleNaive,
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
        };

        // mobility between nodes
        auto read_travel_times   = mio::read_mobility_plain(traveltimes_path);
        auto travel_times        = read_travel_times.value();
        auto read_paths_mobility = mio::read_path_mobility(
            "/localdata1/code/memilio/test/wegketten_ohne_komma.txt"); // wegketten_ohne_komma.txt
        auto path_mobility = read_paths_mobility.value();

        for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
            for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
                // if (county_idx_i == county_idx_j)
                //     continue;
                auto& populations = nodes[county_idx_i].populations;
                // mobility coefficients have the same number of components as the contact matrices.
                // so that the same NPIs/dampings can be used for both (e.g. more home office => fewer commuters)
                auto mobility_coeffs = mio::MigrationCoefficientGroup(contact_locations.size(), populations.numel());

                //commuters
                auto working_population = 0.0;
                auto min_commuter_age   = mio::AgeGroup(2);
                auto max_commuter_age   = mio::AgeGroup(4); //this group is partially retired, only partially commutes
                for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                    working_population += populations.get_group_total(age) * (age == max_commuter_age ? 0.33 : 1.0);
                }
                auto commuter_coeff_ij = mob_data(county_idx_i, county_idx_j) /
                                         working_population; //data is absolute numbers, we need relative
                for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                    for (auto compartment : migrating_compartments) {
                        auto coeff_index = populations.get_flat_index({age, compartment});
                        mobility_coeffs[size_t(ContactLocation::Work)].get_baseline()[coeff_index] =
                            commuter_coeff_ij * (age == max_commuter_age ? 0.33 : 1.0);
                    }
                }

                auto path = path_mobility[county_idx_i][county_idx_j];
                if (path[0] != county_idx_i || path[path.size() - 1] != county_idx_j)
                    std::cout << "falsch"
                              << "\n";
                if (commuter_coeff_ij > theshold_edges) {
                    params_graph.add_edge(county_idx_i, county_idx_j, travel_times(county_idx_i, county_idx_j),
                                          path_mobility[county_idx_i][county_idx_j], std::move(mobility_coeffs));
                }
            }
        }

        // average travel time
        const ScalarType avg_traveltime = std::accumulate(params_graph.edges().begin(), params_graph.edges().end(), 0.0,
                                                          [](double sum, const auto& e) {
                                                              return sum + e.traveltime;
                                                          }) /
                                          params_graph.edges().size();

        // average share of commuters for all counties relative to the total population
        const ScalarType total_population = std::accumulate(params_graph.nodes().begin(), params_graph.nodes().end(),
                                                            0.0, [](double sum, const auto& n) {
                                                                return sum + n.property.populations.get_total();
                                                            });
        const ScalarType num_commuters = std::accumulate(params_graph.edges().begin(), params_graph.edges().end(), 0.0,
                                                         [mob_data](double sum, const auto& e) {
                                                             auto start_node = e.start_node_idx;
                                                             auto end_node   = e.end_node_idx;
                                                             return sum + mob_data(start_node, end_node);
                                                         });
        const ScalarType avg_commuter_share = num_commuters / total_population;

        if (mio::mpi::is_root()) {
            std::cout << "avg commuter share = " << avg_commuter_share << "\n";
            std::cout << "avg travel time = " << avg_traveltime << "\n";
        }

        // scale all contact matrices
        for (auto& node : params_graph.nodes()) {
            // mobility node
            auto contacts_status_mobility =
                set_contact_matrices_transport(data_dir, node.mobility.parameters, avg_traveltime * avg_commuter_share);

            // local node
            auto contacts_status =
                set_contact_matrices(data_dir, node.property.parameters, 1 - avg_traveltime, 1 - avg_commuter_share);
        }
        mio::ContactMatrixGroup const& contact_matrix3 =
            params_graph.nodes()[0].mobility.parameters.get<mio::osecirvvs::ContactPatterns>();
        std::cout << "Mobility node: contact matrix at t = 0\n";
        for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(6); i++) {
            for (auto j = mio::AgeGroup(0); j < mio::AgeGroup(6); j++) {
                std::cout << contact_matrix3.get_matrix_at(0)(static_cast<Eigen::Index>((size_t)i),
                                                              static_cast<Eigen::Index>((size_t)j))
                          << " ";
            }
            std::cout << "\n";
        }
    }

    if (mio::mpi::is_root()) {
        std::cout << "Anzahl kanten = " << params_graph.edges().size() << "\n";
    }
    return params_graph;
}
mio::IOResult<void> run()
{
    // mio::set_log_level(mio::LogLevel::critical);
    const auto num_days = 70;
    const int num_runs  = 150;
    const bool masks    = true;
    const bool ffp2     = true;
    const bool edges    = true;

    // wenn masks false und ffp2 true, dann error ausgeben
    if (!masks && ffp2) {
        mio::log_error("ffp2 only possible with masks");
    }

    bool szenario_cologne = false;

    // auto params_graph = get_graph(num_days, masks);
    BOOST_OUTCOME_TRY(created, get_graph(num_days, masks, ffp2, szenario_cologne, edges));
    auto params_graph = created;

    std::string res_dir = "/localdata1/code/memilio/results_paper/mask_" + std::to_string(masks);

    // if (ffp2)
    //     res_dir = res_dir + "_ffp2";
    // if (szenario_cologne)
    //     res_dir = res_dir + "_cologne";
    // if (!edges)
    //     res_dir = res_dir + "_no_edges";
    if (!edges && (ffp2 || masks))
        mio::log_error("no edges only possible without masks and ffp2");

    res_dir += std::string(ffp2 ? "_ffp2" : "") + std::string(szenario_cologne ? "_cologne" : "") +
               std::string(!edges ? "_no_edges" : "");

    // check if boths dir exist, otherwise create them
    if (!fs::exists(res_dir)) {
        fs::create_directory(res_dir);
    }
    auto write_graph_status = write_graph(params_graph, "/localdata1/code/memilio/save_graph");

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    // parameter study
    auto parameter_study =
        mio::ParameterStudy<mio::FlowSimulation<mio::osecirvvs::Model>>(params_graph, 0.0, num_days, 0.01, num_runs);
    if (mio::mpi::is_root()) {
        printf("Seeds: ");
        for (auto s : parameter_study.get_rng().get_seeds()) {
            printf("%u, ", s);
        }
        printf("\n");
    }
    auto save_single_run_result = mio::IOResult<void>(mio::success());
    auto ensemble               = parameter_study.run(
        [&](auto&& graph) {
            return draw_sample(graph, false);
        },
        [&](auto results_graph, auto&& run_idx) {
            auto interpolated_result = mio::interpolate_simulation_result(results_graph);

            auto params = std::vector<mio::osecirvvs::Model>();
            // params.reserve(results_graph.nodes().size());
            // std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
            //                [](auto&& node) {
            //                    return node.property.get_simulation().get_model();
            //                });
            // const auto size_vec = interpolated_result[0].get_num_time_points();
            // std::vector<ScalarType> num_transmissions_in_mobility_nodes(size_vec, 0.0);
            // std::vector<std::vector<ScalarType>> num_transmission_per_county(results_graph.nodes().size(),
            //                                                                  std::vector<ScalarType>(size_vec, 0.0));
            // std::vector<std::vector<ScalarType>> num_symptomatic_per_county(results_graph.nodes().size(),
            //                                                                 std::vector<ScalarType>(size_vec, 0.0));

            // // vaccinations
            // // std::vector<std::vector<ScalarType>> vacc_county_N(results_graph.nodes().size());
            // // std::vector<std::vector<ScalarType>> vacc_county_PI(results_graph.nodes().size());
            // // std::vector<std::vector<ScalarType>> vacc_county_II(results_graph.nodes().size());

            // std::vector<ScalarType> waning_partial_immunity(size_vec, 0.0);
            // std::vector<ScalarType> waning_improved_immunity(size_vec, 0.0);
            // auto node_idx = 0;
            // for (auto& node : results_graph.nodes()) {
            //     auto flows             = mio::interpolate_simulation_result(node.mobility.get_simulation().get_flows());
            //     auto flows_local_model = mio::interpolate_simulation_result(node.property.get_simulation().get_flows());
            //     for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(6); ++i) {
            //         auto& model = node.mobility.get_simulation().get_model();
            //         auto flow_index_N =
            //             model.template get_flat_flow_index<mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive,
            //                                                mio::osecirvvs::InfectionState::InfectedSymptomsNaive>({i});
            //         auto flow_index_N_confirmed = model.template get_flat_flow_index<
            //             mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed,
            //             mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed>({i});
            //         auto flow_index_PI = model.template get_flat_flow_index<
            //             mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity,
            //             mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity>({i});
            //         auto flow_index_PI_confirmed = model.template get_flat_flow_index<
            //             mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,
            //             mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed>({i});
            //         auto flow_index_II = model.template get_flat_flow_index<
            //             mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity,
            //             mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity>({i});
            //         auto flow_index_II_confirmed = model.template get_flat_flow_index<
            //             mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed,
            //             mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed>({i});

            //         // Flows from Susceptible to Exposed
            //         auto flow_index_SN =
            //             model.template get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptibleNaive,
            //                                                mio::osecirvvs::InfectionState::ExposedNaive>({i});
            //         auto flow_index_SII =
            //             model.template get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity,
            //                                                mio::osecirvvs::InfectionState::ExposedImprovedImmunity>(
            //                 {i});
            //         auto flow_index_SPI =
            //             model.template get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptiblePartialImmunity,
            //                                                mio::osecirvvs::InfectionState::ExposedPartialImmunity>({i});

            //         auto flow_index_waning_pi =
            //             model.template get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptiblePartialImmunity,
            //                                                mio::osecirvvs::InfectionState::SusceptibleNaive>({i});
            //         auto flow_index_waning_ii =
            //             model.template get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity,
            //                                                mio::osecirvvs::InfectionState::SusceptiblePartialImmunity>(
            //                 {i});
            //         for (Eigen::Index time_idx = 0; time_idx < flows.get_num_time_points(); ++time_idx) {
            //             const auto& flow_values = flows.get_value(time_idx);
            //             auto time               = flows.get_time(time_idx);

            //             auto indx_vector = static_cast<size_t>(time_idx);
            //             // num_transmissions_in_mobility_nodes[indx_vector] +=
            //             //     flow_values(flow_index_N) + flow_values(flow_index_PI) + flow_values(flow_index_II) +
            //             //     flow_values(flow_index_N_confirmed) + flow_values(flow_index_PI_confirmed) +
            //             //     flow_values(flow_index_II_confirmed);
            //             num_transmissions_in_mobility_nodes[indx_vector] +=
            //                 flow_values(flow_index_SN) + flow_values(flow_index_SPI) + flow_values(flow_index_SII);
            //             waning_partial_immunity[indx_vector] += flow_values(flow_index_waning_pi);
            //             waning_improved_immunity[indx_vector] += flow_values(flow_index_waning_ii);
            //             num_symptomatic_per_county[node_idx][indx_vector] +=
            //                 flow_values(flow_index_N) + flow_values(flow_index_PI) + flow_values(flow_index_II) +
            //                 flow_values(flow_index_N_confirmed) + flow_values(flow_index_PI_confirmed) +
            //                 flow_values(flow_index_II_confirmed);
            //             num_transmission_per_county[node_idx][indx_vector] +=
            //                 flow_values(flow_index_SN) + flow_values(flow_index_SPI) + flow_values(flow_index_SII);
            //         }
            //         for (Eigen::Index time_idx = 0; time_idx < flows_local_model.get_num_time_points(); ++time_idx) {
            //             const auto& flow_values_local_model = flows_local_model.get_value(time_idx);
            //             auto time                           = flows_local_model.get_time(time_idx);
            //             auto indx_vector                    = static_cast<size_t>(time_idx);
            //             waning_partial_immunity[indx_vector] += flow_values_local_model(flow_index_waning_pi);
            //             waning_improved_immunity[indx_vector] += flow_values_local_model(flow_index_waning_ii);
            //             num_transmission_per_county[node_idx][indx_vector] += flow_values_local_model(flow_index_SN) +
            //                                                                   flow_values_local_model(flow_index_SPI) +
            //                                                                   flow_values_local_model(flow_index_SII);
            //             num_symptomatic_per_county[node_idx][indx_vector] +=
            //                 flow_values_local_model(flow_index_N) + flow_values_local_model(flow_index_PI) +
            //                 flow_values_local_model(flow_index_II) + flow_values_local_model(flow_index_N_confirmed) +
            //                 flow_values_local_model(flow_index_PI_confirmed) +
            //                 flow_values_local_model(flow_index_II_confirmed);
            //         }
            //     }

            //     // // vaccinations
            //     // vacc_county_N[node_idx].resize(flows_local_model.get_num_time_points(), 0.0);
            //     // vacc_county_PI[node_idx].resize(flows_local_model.get_num_time_points(), 0.0);
            //     // vacc_county_II[node_idx].resize(flows_local_model.get_num_time_points(), 0.0);
            //     // for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(6); ++i) {
            //     //     auto& model = node.mobility.get_simulation().get_model();
            //     //     auto vacc_N =
            //     //         model.template get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptibleNaive,
            //     //                                       mio::osecirvvs::InfectionState::TemporaryImmunity1>({i});
            //     //     auto flow_index_PI =
            //     //         model.template get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptiblePartialImmunity,
            //     //                                       mio::osecirvvs::InfectionState::TemporaryImmunity2>({i});
            //     //     auto flow_index_II =
            //     //         model.template get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity,
            //     //                                       mio::osecirvvs::InfectionState::TemporaryImmunity2>({i});

            //     //     for (Eigen::Index time_idx = 0; time_idx < flows_local_model.get_num_time_points(); ++time_idx) {
            //     //         const auto& flow_values_local_model = flows_local_model.get_value(time_idx);
            //     //         auto time                           = flows_local_model.get_time(time_idx);
            //     //         auto indx_vector                    = static_cast<size_t>(time_idx);
            //     //         vacc_county_N[node_idx][indx_vector] += flow_values_local_model(vacc_N);
            //     //         vacc_county_PI[node_idx][indx_vector] += flow_values_local_model(flow_index_PI);
            //     //         vacc_county_II[node_idx][indx_vector] += flow_values_local_model(flow_index_II);
            //     //     }
            //     // }
            //     node_idx++;
            // }
            // std::string results_path = "/localdata1/code/memilio/results_paper/mask_" + std::to_string(masks);
            // if (ffp2)
            //     results_path = results_path + "_ffp2";

            // if (szenario_cologne)
            //     results_path = results_path + "_cologne";

            // if (!edges)
            //     results_path = results_path + "_no_edges";
            // const std::string results_path_mb = results_path + "/flows_mb";
            // if (!fs::exists(results_path_mb)) {
            //     fs::create_directory(results_path_mb);
            // }
            // const std::string results_filename = results_path_mb + "/results_run_" + std::to_string(run_idx) + ".txt";
            // std::ofstream outfile(results_filename);

            // outfile << "TimeIndex Num_Transmissions_in_Mobility_Nodes Waning_Partial_Immunity "
            //            "Waning_Improved_Immunity\n";
            // size_t vec_size = num_transmissions_in_mobility_nodes.size();
            // for (size_t i = 0; i < vec_size; ++i) {
            //     outfile << i << " " << num_transmissions_in_mobility_nodes[i] << " " << waning_partial_immunity[i]
            //             << " " << waning_improved_immunity[i] << "\n";
            // }
            // outfile.close();

            // // erzeuge ordner flows_county
            // const std::string results_path_flows = results_path + "/flows_county";
            // if (!fs::exists(results_path_flows)) {
            //     fs::create_directory(results_path_flows);
            // }
            // const std::string results_filename_county =
            //     results_path_flows + "/transmission_county_run_" + std::to_string(run_idx) + ".txt";
            // std::ofstream outfile_county(results_filename_county);
            // for (size_t i = 0; i < results_graph.nodes().size(); ++i) {
            //     outfile_county << i << " ";
            //     for (size_t j = 0; j < vec_size; ++j) {
            //         outfile_county << num_transmission_per_county[i][j] << " ";
            //     }
            //     outfile_county << "\n";
            // }
            // outfile_county.close();

            // const std::string results_filename_county_symptomatic =
            //     results_path_flows + "/symptomatic_county_run_" + std::to_string(run_idx) + ".txt";
            // std::ofstream outfile_county_symptomatic(results_filename_county_symptomatic);
            // for (size_t i = 0; i < results_graph.nodes().size(); ++i) {
            //     outfile_county_symptomatic << i << " ";
            //     for (size_t j = 0; j < vec_size; ++j) {
            //         outfile_county_symptomatic << num_symptomatic_per_county[i][j] << " ";
            //     }
            //     outfile_county_symptomatic << "\n";
            // }
            // outfile_county_symptomatic.close();

            // // // speichere vacc_county_N, vacc_county_PI, vacc_county_II in jeweils eigene Datei
            // // const std::string results_path_vacc = results_path + "/vacc_county";
            // // // check if dir exist, otherwise create them
            // // if (!fs::exists(results_path_vacc)) {
            // //     fs::create_directory(results_path_vacc);
            // // }
            // // const std::string results_filename_vacc_N =
            // //     results_path_vacc + "/vacc_county_N_run_" + std::to_string(run_idx) + ".txt";
            // // const std::string results_filename_vacc_PI =
            // //     results_path_vacc + "/vacc_county_PI_run_" + std::to_string(run_idx) + ".txt";
            // // const std::string results_filename_vacc_II =
            // //     results_path_vacc + "/vacc_county_II_run_" + std::to_string(run_idx) + ".txt";

            // // std::ofstream outfile_vacc_N(results_filename_vacc_N);
            // // std::ofstream outfile_vacc_PI(results_filename_vacc_PI);
            // // std::ofstream outfile_vacc_II(results_filename_vacc_II);

            // // for (size_t i = 0; i < results_graph.nodes().size(); ++i) {
            // //     outfile_vacc_N << i << " ";
            // //     outfile_vacc_PI << i << " ";
            // //     outfile_vacc_II << i << " ";

            // //     // Abfrage der Größe des aktuellen inneren Vektors:
            // //     size_t current_vec_size_N  = vacc_county_N[i].size();
            // //     size_t current_vec_size_PI = vacc_county_PI[i].size();
            // //     size_t current_vec_size_II = vacc_county_II[i].size();

            // //     for (size_t j = 0; j < current_vec_size_N; ++j) {
            // //         outfile_vacc_N << vacc_county_N[i][j] << " ";
            // //     }
            // //     for (size_t j = 0; j < current_vec_size_PI; ++j) {
            // //         outfile_vacc_PI << vacc_county_PI[i][j] << " ";
            // //     }
            // //     for (size_t j = 0; j < current_vec_size_II; ++j) {
            // //         outfile_vacc_II << vacc_county_II[i][j] << " ";
            // //     }

            // //     outfile_vacc_N << "\n";
            // //     outfile_vacc_PI << "\n";
            // //     outfile_vacc_II << "\n";
            // // }
            // // outfile_vacc_N.close();
            // // outfile_vacc_PI.close();
            // // outfile_vacc_II.close();

            // std::cout << "Run " << run_idx << " complete." << std::endl;
            return std::make_pair(interpolated_result, params);
        });

    if (ensemble.size() > 0) {
        auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_results.reserve(ensemble.size());
        auto ensemble_params = std::vector<std::vector<mio::osecirvvs::Model>>{};
        ensemble_params.reserve(ensemble.size());
        for (auto&& run : ensemble) {
            ensemble_results.emplace_back(std::move(run.first));
            ensemble_params.emplace_back(std::move(run.second));
        }
        BOOST_OUTCOME_TRY(save_results(ensemble_results, ensemble_params, county_ids, res_dir, false));
    }

    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::warn);
    mio::mpi::init();

    auto result = run();
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        mio::mpi::finalize();
        return -1;
    }
    mio::mpi::finalize();
    return 0;
}