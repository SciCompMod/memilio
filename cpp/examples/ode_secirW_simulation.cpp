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
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/parameters.h"
#include "ode_secirvvs/parameters_io.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"
#include "memilio/io/result_io.h"
#include "memilio/compartments/parameter_studies.h"

#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/mobility_io.h"
#include "memilio/utils/date.h"
#include "ode_secirvvs/parameter_space.h"
#include "memilio/utils/stl_util.h"
#include "boost/filesystem.hpp"
#include "memilio/math/eigen.h"
#include <cstdio>
#include <iomanip>

#include <iostream>
#include <string>

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
    const double timeInfectedSymptomsMin = 6.37; // doi.org/10.1136/bmj.o922
    const double timeInfectedSymptomsMax = 7.37; // doi.org/10.1136/bmj.o922
    // const double t_inf_hosp_min[] = {7, 7, 7, 7, 7, 7}; // doi.org/10.1016/S0140-6736(20)30183-5
    // const double t_inf_hosp_max[] = {7, 7, 7, 7, 7, 7}; // doi.org/10.1016/S0140-6736(20)30183-5
    const double timeInfectedSevereMin[] = {1.8, 1.8, 1.8, 2.5, 3.5, 4.91}; // doi.org/10.1101/2022.03.16.22271361
    const double timeInfectedSevereMax[] = {2.3, 2.3, 2.3, 3.67, 5, 7.01}; // doi.org/10.1101/2022.03.16.22271361
    // const double t_hosp_icu_min[] = {0.67, 0.67, 0.67, 1.54, 1.7, 1.83}; // doi.org/10.1101/2022.03.16.22271361
    // const double t_hosp_icu_max[] = {0.97, 0.97, 0.97, 2.08, 2.2, 2.42}; // doi.org/10.1101/2022.03.16.22271361
    const double timeInfectedCriticalMin[] = {5, 5, 5, 14, 14, 10}; // improvable doi.org/10.1371/journal.pcbi.1010054
    const double timeInfectedCriticalMax[] = {9, 9, 9, 21, 21, 15}; // improvable doi.org/10.1371/journal.pcbi.1010054
    // const double t_icu_dead_min[]          = {4, 4, 4, 15, 15, 10}; // improvable doi.org/10.1371/journal.pcbi.1010054
    // const double t_icu_dead_max[]          = {8, 8, 8, 18, 18, 12}; // improvable doi.org/10.1371/journal.pcbi.1010054

    array_assign_uniform_distribution(params.get<mio::osecirvvs::IncubationTime>(), incubationTime, incubationTime);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SerialInterval>(), serialIntervalMin,
                                      serialIntervalMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSymptoms>(), timeInfectedSymptomsMin,
                                      timeInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSevere>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedCritical>(), timeInfectedCriticalMin,
                                      timeInfectedCriticalMax);

    const double transmissionProbabilityOnContactMin[] = {0.14, 0.15, 0.12,
                                                          0.1,  0.08, 0.16}; // doi.org/10.21203/rs.3.rs-1729679/v1
    const double transmissionProbabilityOnContactMax[] = {0.15, 0.18, 0.16,
                                                          0.13, 0.16, 0.18}; // doi.org/10.21203/rs.3.rs-1729679/v1

    const double relativeTransmissionNoSymptomsMin[]  = {0.6, 0.55, 0.65,
                                                        0.7, 0.75, 0.85}; // doi.org/10.1101/2022.05.05.22274697
    const double relativeTransmissionNoSymptomsMax[]  = {0.8, 0.75,  0.8,
                                                        0.8, 0.825, 0.9}; // doi.org/10.1101/2022.05.05.22274697
    const double riskOfInfectionFromSymptomaticMin    = 0.0; // beta (depends on incidence and test and trace capacity)
    const double riskOfInfectionFromSymptomaticMax    = 0.2;
    const double maxRiskOfInfectionFromSymptomaticMin = 0.4;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.5;

    const double recoveredPerInfectedNoSymptomsMin[] = {0.2, 0.25,  0.2,
                                                        0.2, 0.175, 0.1}; // doi.org/10.1101/2022.05.05.22274697
    const double recoveredPerInfectedNoSymptomsMax[] = {0.4, 0.45, 0.35,
                                                        0.3, 0.25, 0.15}; // doi.org/10.1101/2022.05.05.22274697
    const double severePerInfectedSymptomsMin[]      = {0.054, 0.005, 0.01, 0.013,
                                                   0.076, 0.251}; // doi.org/10.2807/1560-7917.ES.2022.27.22.2200396
    const double severePerInfectedSymptomsMax[]      = {0.054, 0.005, 0.01, 0.013,
                                                   0.076, 0.251}; // doi.org/10.2807/1560-7917.ES.2022.27.22.2200396

    const double criticalPerSevereMin[] = {
        0.0511, 0.0686, 0.0491, 0.114,
        0.1495, 0.0674}; // www.sozialministerium.at/dam/jcr:f472e977-e1bf-415f-95e1-35a1b53e608d/Factsheet_Coronavirus_Hospitalisierungen.pdf
    const double criticalPerSevereMax[] = {
        0.0511, 0.0686, 0.0491, 0.114,
        0.1495, 0.0674}; // www.sozialministerium.at/dam/jcr:f472e977-e1bf-415f-95e1-35a1b53e608d/Factsheet_Coronavirus_Hospitalisierungen.pdf
    const double deathsPerCriticalMin[] = {0.00, 0.00, 0.10,
                                           0.10, 0.30, 0.5}; // improvable doi.org/10.1371/journal.pcbi.1010054
    const double deathsPerCriticalMax[] = {0.10, 0.10, 0.18,
                                           0.18, 0.50, 0.7}; // improvable doi.org/10.1371/journal.pcbi.1010054

    const double reducExposedPartialImmunityMin           = 0.569; // doi.org/10.1136/bmj-2022-071502
    const double reducExposedPartialImmunityMax           = 0.637; // doi.org/10.1136/bmj-2022-071502
    const double reducExposedImprovedImmunityMin          = 0.46; // doi.org/10.1136/bmj-2022-071502
    const double reducExposedImprovedImmunityMax          = 0.57; // doi.org/10.1136/bmj-2022-071502
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
    // else {
    //     double vacc_add{0.25};
    //     const double reduc_partial_exp_min  = 0.569 + vacc_add; // doi.org/10.1136/bmj-2022-071502
    //     const double reduc_partial_exp_max  = 0.637 + vacc_add; // doi.org/10.1136/bmj-2022-071502
    //     const double reduc_improved_exp_min = 0.46 + vacc_add; // doi.org/10.1136/bmj-2022-071502
    //     const double reduc_improved_exp_max = 0.57 + vacc_add; // doi.org/10.1136/bmj-2022-071502
    //     const double reduc_partial_inf_min  = 0.746 + vacc_add; // doi.org/10.1056/NEJMoa2119451
    //     const double reduc_partial_inf_max  = 1.; // doi.org/10.1056/NEJMoa2119451
    //     const double reduc_improved_inf_min = 0.295 + vacc_add; // doi.org/10.1056/NEJMoa2119451
    //     const double reduc_improved_inf_max = 0.344 + vacc_add; // doi.org/10.1056/NEJMoa2119451
    //     const double reduc_vacc_hosp_min =
    //         0.52 +
    //         vacc_add; // www.assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf
    //     const double reduc_vacc_hosp_max =
    //         1.; // www.assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf
    //     const double reduc_improved_hosp_min  = 0.1 + vacc_add; // doi.org/10.1136/bmj-2022-071502
    //     const double reduc_improved_hosp_max  = 0.19 + vacc_add; // doi.org/10.1136/bmj-2022-071502
    //     const double temp_reduc_mild_rec_time = 0.5 + vacc_add; // doi.org/10.1101/2021.09.24.21263978
    // }

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

    // const double ImmunityInterval1Min = 30;
    // const double ImmunityInterval1Max = 30;

    // array_assign_uniform_distribution(params.get<mio::osecirvvs::ImmunityInterval1>(), ImmunityInterval1Min,
    //                                   ImmunityInterval1Max);

    // const double ImmunityInterval2Min = 30;
    // const double ImmunityInterval2Max = 30;

    // array_assign_uniform_distribution(params.get<mio::osecirvvs::ImmunityInterval2>(), ImmunityInterval2Min,
    //                                   ImmunityInterval2Max);

    // params.get<mio::osecirvvs::ImmunityInterval1>() = 30.0;
    // params.get<mio::osecirvvs::ImmunityInterval2>() = 60.0;

    // params.get<mio::osecirvvs::WaningPartialImmunity>()  = 30.0;
    // params.get<mio::osecirvvs::WaningImprovedImmunity>() = 60.0;

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
 * Create the input graph for the parameter study.
 * Reads files from the data directory.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @returns created graph or any io errors that happen during reading of the files.
 */
mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> get_graph(const int num_days)
{
    std::string data_dir         = "/localdata1/code/memilio//data";
    std::string traveltimes_path = "/localdata1/test/memilio/traveltimes_100kmh.txt";
    std::string durations_path   = "/localdata1/code/memilio/activity_duration_work.txt";
    auto start_date              = mio::Date(2021, 6, 1);
    auto end_date                = mio::Date(2022, 1, 1);

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
        // for (size_t node_idx = 0; node_idx < 5; ++node_idx) {

        auto tnt_capacity = nodes[node_idx].populations.get_total() * tnt_capacity_factor;

        //local parameters
        auto& tnt_value = nodes[node_idx].parameters.template get<mio::osecirvvs::TestAndTraceCapacity>();
        tnt_value       = mio::UncertainValue(0.5 * (1.2 * tnt_capacity + 0.8 * tnt_capacity));
        tnt_value.set_distribution(mio::ParameterDistributionUniform(0.8 * tnt_capacity, 1.2 * tnt_capacity));

        auto id = int(mio::regions::CountyId(node_ids[node_idx]));

        //holiday periods
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

        params_graph.add_node(node_ids[node_idx], duration_stay((Eigen::Index)node_idx), nodes[node_idx]);
        // params_graph.add_node(node_ids[node_idx], nodes[node_idx]);
    }

    // Edges
    auto migrating_compartments = {mio::osecirvvs::InfectionState::SusceptibleNaive,
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
                                   mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity};

    // mobility between nodes
    auto mobility_data_commuter =
        mio::read_mobility_plain(("/localdata1/test/memilio/data/mobility/commuter_migration_scaled.txt"));
    auto mob_data = mobility_data_commuter.value();
    auto mobility_data_twitter =
        mio::read_mobility_plain(("/localdata1/test/memilio/data/mobility/twitter_scaled_1252.txt"));
    auto read_travel_times = mio::read_mobility_plain(traveltimes_path);
    auto travel_times      = read_travel_times.value();

    auto read_paths_mobility =
        mio::read_path_mobility("/localdata1/test/memilio/wegketten_ohne_komma.txt"); // wegketten_ohne_komma.txt
    auto path_mobility = read_paths_mobility.value();

    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            if (county_idx_i == county_idx_j)
                continue;
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
            auto commuter_coeff_ij =
                mob_data(county_idx_i, county_idx_j) / working_population; //data is absolute numbers, we need relative
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
            if (commuter_coeff_ij > 4e-5) {
                params_graph.add_edge(county_idx_i, county_idx_j, travel_times(county_idx_i, county_idx_j),
                                      path_mobility[county_idx_i][county_idx_j], std::move(mobility_coeffs));
            }
        }
    }

    int num_edges = 0;
    for (auto& e : params_graph.edges()) {
        num_edges += 1;
    }
    std::cout << "Anzahl kanten = " << num_edges << "\n";

    return params_graph;
}

int main(int argc, char** argv)
{
    // TODO: No vacc currently
    mio::set_log_level(mio::LogLevel::critical);
    const auto t0       = 0.;
    const auto num_days = 15.0;
    const auto dt       = 0.5;
    const int num_runs  = 2;

    auto params_graph = get_graph(num_days);

    // auto write_graph_status = write_graph(params_graph, "/localdata1/code/memilio/save_graph");

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });
    // parameter study

    auto parameter_study =
        mio::ParameterStudy<mio::osecirvvs::Simulation<>>(params_graph, 0.0, num_days, 0.5, num_runs);
    auto save_single_run_result = mio::IOResult<void>(mio::success());
    auto ensemble               = parameter_study.run(
        [&](auto&& graph) {
            return draw_sample(graph, false);
        },
        [&](auto results_graph, auto&& run_idx) {
            auto interpolated_result = mio::interpolate_simulation_result(results_graph);
            auto params              = std::vector<mio::osecirvvs::Model>();
            params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
                           [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            save_single_run_result = save_result_with_params(interpolated_result, params, county_ids,
                                                             "/localdata1/test/memilio/test", run_idx);

            std::cout << "run " << run_idx << " complete." << std::endl;
            return std::make_pair(interpolated_result, params);
        });
    auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
    ensemble_results.reserve(ensemble.size());
    auto ensemble_params = std::vector<std::vector<mio::osecirvvs::Model>>{};
    ensemble_params.reserve(ensemble.size());
    for (auto&& run : ensemble) {
        ensemble_results.emplace_back(std::move(run.first));
        ensemble_params.emplace_back(std::move(run.second));
    }
    auto save_results_status =
        save_results(ensemble_results, ensemble_params, county_ids, "/localdata1/test/memilio/test2", false);

    // std::cout << "Landkreis 0" << std::endl;
    // std::cout << gr.nodes()[0].property.get_result().get_last_value() << std::endl;
    // std::cout << "Landkreis 1" << std::endl;
    // std::cout << gr.nodes()[1].property.get_result().get_last_value() << std::endl;

    // auto region0_result = gr.nodes()[0].property.get_result();
    // auto region1_result = gr.nodes()[1].property.get_result();

    // auto num_points = static_cast<size_t>(region0_result.get_num_time_points());
    // for (size_t i = 0; i < num_points; i++) {
    //     std::cout << region0_result.get_time(i);
    //     std::cout << region0_result.get_value(i);
    // }

    // std::vector<mio::TimeSeries<double>> results_from_sim = {region0_result, region1_result};
    // std::vector<int> ids                                  = {0, 1};

    // auto save_result_status = mio::save_result(results_from_sim, ids, (int)(size_t)1, "test_result_1.h5");

    return 0;
}
