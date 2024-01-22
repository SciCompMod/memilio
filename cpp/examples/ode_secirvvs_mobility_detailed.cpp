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
mio::IOResult<mio::GraphDetailed<mio::osecirvvs::Model, mio::MigrationParameters>>
get_graph(const std::string data_dir, const int num_days, bool masks, bool ffp2, bool szenario_cologne, bool edges)
{
    auto mobility_data_commuter =
        mio::read_mobility_plain(("/localdata1/code/memilio/test/commuter_migration_with_locals.txt"));
    auto mob_data   = mobility_data_commuter.value();
    auto start_date = mio::Date(2021, 8, 1);
    auto end_date   = mio::Date(2021, 11, 1);

    // global parameters
    const int num_age_groups = 6;
    mio::GraphDetailed<mio::osecirvvs::Model, mio::MigrationParameters> params_graph;

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

    const auto& read_function_nodes = mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model>;
    const auto& read_function_edges = mio::read_mobility_plain;
    const auto& node_id_function    = mio::get_node_ids;

    const auto& set_node_function =
        mio::set_nodes_detailed<mio::osecirvvs::TestAndTraceCapacity, mio::osecirvvs::ContactPatterns,
                                mio::osecirvvs::Model, mio::MigrationParameters, mio::osecirvvs::Parameters,
                                decltype(read_function_nodes), decltype(node_id_function)>;
    BOOST_OUTCOME_TRY(set_node_function(
        params, start_date, end_date, data_dir, data_dir + "//pydata//Germany//county_current_population.json",
        mio::path_join(data_dir, "mobility", "activity_duration_work.txt"), true, params_graph, read_function_nodes,
        node_id_function, scaling_factor_infected, scaling_factor_icu, tnt_capacity_factor,
        mio::get_offset_in_days(end_date, start_date), false));

    for (auto& node : params_graph.nodes()) {
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

            double fac_variant                           = 1.94; //https://doi.org/10.7554/eLife.78933
            double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                            0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

            double transmissionProbabilityOnContactMax[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                                            0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};
            for (int i = 0; i < 6; i++) {
                transmissionProbabilityOnContactMin[i] = transmissionProbabilityOnContactMin[i] * factor_mask[i];
                transmissionProbabilityOnContactMax[i] = transmissionProbabilityOnContactMax[i] * factor_mask[i];
            }
            array_assign_uniform_distribution(
                node.mobility.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>(),
                transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
        }

        for (size_t t_idx = 0; t_idx < num_days; ++t_idx) {
            auto t = mio::SimulationDay((size_t)t_idx);
            for (auto j = mio::AgeGroup(0); j < params.get_num_groups(); j++) {
                node.mobility.parameters.template get<mio::osecirvvs::DailyFirstVaccination>()[{j, t}] = 0;
                node.mobility.parameters.template get<mio::osecirvvs::DailyFullVaccination>()[{j, t}]  = 0;
            }
        }
    }

    if (szenario_cologne) {
        init_pop_cologne_szenario(params_graph);
    }

    // Edges
    if (edges) {
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

        const auto& read_function_edges = mio::read_mobility_plain;
        const auto& set_edge_function =
            mio::set_edges_detailed<ContactLocation, mio::osecirvvs::Model, mio::MigrationParameters,
                                    mio::MigrationCoefficientGroup, mio::osecirvvs::InfectionState>;
        BOOST_OUTCOME_TRY(set_edge_function(mio::path_join(data_dir, "mobility", "travel_times_pathes.txt"),
                                            data_dir + "//mobility//commuter_migration_with_locals.txt",
                                            mio::path_join(data_dir, "mobility", "wegketten_ohne_komma.txt"),
                                            params_graph, migrating_compartments, contact_locations.size(),
                                            std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.}, false));

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
    std::string data_dir = "/localdata1/test/memilio/data";

    bool szenario_cologne = false;

    BOOST_OUTCOME_TRY(created, get_graph(data_dir, num_days, masks, ffp2, szenario_cologne, edges));
    auto params_graph = created;

    std::string res_dir = "/localdata1/code/memilio/results_paper/mask_" + std::to_string(masks);

    if (ffp2)
        res_dir = res_dir + "_ffp2";
    if (szenario_cologne)
        res_dir = res_dir + "_cologne";
    if (!edges)
        res_dir = res_dir + "_no_edges";
    if (!edges && (ffp2 || masks))
        mio::log_error("no edges only possible without masks and ffp2");

    res_dir += std::string(ffp2 ? "_ffp2" : "") + std::string(szenario_cologne ? "_cologne" : "") +
               std::string(!edges ? "_no_edges" : "");

    // check if boths dir exist, otherwise create them
    if (!fs::exists(res_dir)) {
        fs::create_directory(res_dir);
    }

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    // parameter study
    auto parameter_study = mio::ParameterStudyDetailed<mio::FlowSimulation<mio::osecirvvs::Model>>(
        params_graph, 0.0, num_days, 0.01, num_runs);
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
            return draw_sample(graph);
        },
        [&](auto results_graph, auto&& run_idx) {
            auto interpolated_result = mio::interpolate_simulation_result(results_graph);
            auto params              = std::vector<mio::osecirvvs::Model>();

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

int main(int, char**)
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