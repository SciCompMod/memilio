/*
* Copyright (C) 2020-2024 MEmilio
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

#include "memilio/compartments/parameter_studies.h"
#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/parameter_space.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/metapopulation_mobility_detailed.h"
#include "memilio/utils/stl_util.h"
#include "boost/filesystem.hpp"
#include <cstdio>
#include <iomanip>

namespace fs = boost::filesystem;

/**
 * Set a value and distribution of an UncertainValue.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution.
 * @param p uncertain value to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void assign_uniform_distribution(mio::UncertainValue<double>& p, double min, double max)
{
    p = mio::UncertainValue<double>(0.5 * (max + min));
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
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
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
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       double min, double max)
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
mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters<double>& params)
{
    //times
    // doi.org/10.1016/j.lanepe.2022.100446 , doi.org/10.3201/eid2806.220158
    const double timeExposedMin            = 1.66;
    const double timeExposedMax            = 1.66;
    const double timeInfectedNoSymptomsMin = 1.44;
    const double timeInfectedNoSymptomsMax = 1.44;

    const double timeInfectedSymptomsMin = 6.58; //https://doi.org/10.1016/S0140-6736(22)00327-0
    const double timeInfectedSymptomsMax = 7.16; //https://doi.org/10.1016/S0140-6736(22)00327-0
    const double timeInfectedSevereMin[] = {1.8, 1.8, 1.8, 2.5, 3.5, 4.91}; // doi.org/10.1186/s12879-022-07971-6
    const double timeInfectedSevereMax[] = {2.3, 2.3, 2.3, 3.67, 5, 7.01}; // doi.org/10.1186/s12879-022-07971-6

    const double timeInfectedCriticalMin[] = {9.29,   9.29,  9.29,
                                              10.842, 11.15, 11.07}; // https://doi.org/10.1186/s12879-022-07971-6
    const double timeInfectedCriticalMax[] = {10.57, 10.57, 10.57,
                                              12.86, 13.23, 13.25}; // https://doi.org/10.1186/s12879-022-07971-6

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeExposed<double>>(), timeExposedMin,
                                      timeExposedMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>(),
                                      timeInfectedNoSymptomsMin, timeInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSymptoms<double>>(),
                                      timeInfectedSymptomsMin, timeInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSevere<double>>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedCritical<double>>(),
                                      timeInfectedCriticalMin, timeInfectedCriticalMax);

    //probabilities
    double fac_variant                                 = 1.5; //https://doi.org/10.7554/eLife.78933
    const double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                          0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

    const double transmissionProbabilityOnContactMax[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                                          0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};

    const double relativeTransmissionNoSymptomsMin = 0.5;

    //{0.6, 0.55, 0.65,0.7, 0.75, 0.85}; // DOI: 10.1097/INF.0000000000003791
    const double relativeTransmissionNoSymptomsMax = 0.5;
    // {0.8, 0.75,  0.8,0.8, 0.825, 0.9}; // DOI: 10.1097/INF.0000000000003791
    const double riskOfInfectionFromSymptomaticMin    = 0.0; // beta (depends on incidence and test and trace capacity)
    const double riskOfInfectionFromSymptomaticMax    = 0.2;
    const double maxRiskOfInfectionFromSymptomaticMin = 0.4;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.5;

    // DOI: 10.1097/INF.0000000000003791 geht hier auch. aber ähnliche werte
    const double recoveredPerInfectedNoSymptomsMin[] = {0.2, 0.25,  0.2,
                                                        0.2, 0.175, 0.1}; // doi.org/10.1101/2022.05.05.22274697
    const double recoveredPerInfectedNoSymptomsMax[] = {0.4, 0.45, 0.35,
                                                        0.3, 0.25, 0.15}; // doi.org/10.1101/2022.05.05.22274697

    // 56% weniger risiko ins krankenhaus doi:10.1136/bmjgh-2023-0123
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10347449/pdf/bmjgh-2023-012328.pdf
    // alternativ: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9321237/pdf/vaccines-10-01049.pdf

    // Faktoren aus https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00462-7/fulltext
    const double severePerInfectedSymptomsMin[] = {1 * 0.006,   0.8 * 0.006, 0.4 * 0.015,
                                                   0.3 * 0.049, 0.25 * 0.15, 0.35 * 0.2}; // 2021 paper
    const double severePerInfectedSymptomsMax[] = {1 * 0.009,   0.8 * 0.009, 0.4 * 0.023, 0.3 * 0.074,
                                                   0.25 * 0.18, 0.35 * 0.25}; // quelle 2021 paper + factors

    // const double criticalPerSevereMin[] = {
    //     0.0511, 0.0686, 0.0491, 0.114,
    //     0.1495, 0.0674}; // www.sozialministerium.at/dam/jcr:f472e977-e1bf-415f-95e1-35a1b53e608d/Factsheet_Coronavirus_Hospitalisierungen.pdf
    // const double criticalPerSevereMax[] = {
    //     0.0511, 0.0686, 0.0491, 0.114,
    //     0.1495, 0.0674}; // www.sozialministerium.at/dam/jcr:f472e977-e1bf-415f-95e1-35a1b53e608d/Factsheet_Coronavirus_Hospitalisierungen.pdf

    // delta paper
    // risk of icu admission https://doi.org/10.1177/14034948221108548
    const double fac_icu                = 0.52;
    const double criticalPerSevereMin[] = {0.05 * fac_icu, 0.05 * fac_icu, 0.05 * fac_icu,
                                           0.10 * fac_icu, 0.25 * fac_icu, 0.35 * fac_icu};
    const double criticalPerSevereMax[] = {0.10 * fac_icu, 0.10 * fac_icu, 0.10 * fac_icu,
                                           0.20 * fac_icu, 0.35 * fac_icu, 0.45 * fac_icu};

    // 61% weniger risiko zu sterben doi:10.1136/bmjgh-2023-0123
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10347449/pdf/bmjgh-2023-012328.pdf
    const double fac_dead               = 0.39;
    const double deathsPerCriticalMin[] = {fac_dead * 0.00, fac_dead * 0.00, fac_dead * 0.10,
                                           fac_dead * 0.10, fac_dead * 0.30, fac_dead * 0.5}; // 2021 paper
    const double deathsPerCriticalMax[] = {fac_dead * 0.10, fac_dead * 0.10, fac_dead * 0.18,
                                           fac_dead * 0.18, fac_dead * 0.50, fac_dead * 0.7};

    // alternative https://doi.org/10.1136/bmj-2022-070695
    // const double fac_dead_u59           = 0.14;
    // const double fac_dead_p59           = 0.44;
    // const double deathsPerCriticalMin[] = {fac_dead_u59 * 0.00, fac_dead_u59 * 0.00, fac_dead_u59 * 0.10,
    //                                        fac_dead_u59 * 0.10, fac_dead_p59 * 0.30, fac_dead_p59 * 0.5}; // 2021 paper
    // const double deathsPerCriticalMax[] = {fac_dead_u59 * 0.10, fac_dead_u59 * 0.10, fac_dead_u59 * 0.18,
    //                                        fac_dead_u59 * 0.18, fac_dead_p59 * 0.50, fac_dead_p59 * 0.7};

    const double reducExposedPartialImmunityMin = 1.0; //0.569; // doi.org/10.1136/bmj-2022-071502
    const double reducExposedPartialImmunityMax = 1.0; // 0.637; // doi.org/10.1136/bmj-2022-071502
    // const double reducExposedImprovedImmunityMin = 0.36; // https://doi.org/10.1038/s41591-023-02219-5
    // const double reducExposedImprovedImmunityMax = 0.66; // https://doi.org/10.1038/s41591-023-02219-5

    const double reducExposedImprovedImmunityMin = 1.0;
    //0.34 * reducExposedPartialImmunityMin; // https://jamanetwork.com/journals/jama/fullarticle/2788487 0.19346
    const double reducExposedImprovedImmunityMax = 1.0;
    //0.34 * reducExposedPartialImmunityMax; // https://jamanetwork.com/journals/jama/fullarticle/2788487 0.21658

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
    const double reducTimeInfectedMild                              = 0.5; // doi.org/10.1101/2021.09.24.21263978

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<double>>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::CriticalPerSevere<double>>(), criticalPerSevereMin,
                                      criticalPerSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::DeathsPerCritical<double>>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax);

    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>(),
                                      reducExposedPartialImmunityMin, reducExposedPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>(),
                                      reducExposedImprovedImmunityMin, reducExposedImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>(),
                                      reducInfectedSymptomsPartialImmunityMin, reducInfectedSymptomsPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>(),
                                      reducInfectedSymptomsImprovedImmunityMin,
                                      reducInfectedSymptomsImprovedImmunityMax);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(),
        reducInfectedSevereCriticalDeadPartialImmunityMin, reducInfectedSevereCriticalDeadPartialImmunityMax);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(),
        reducInfectedSevereCriticalDeadImprovedImmunityMin, reducInfectedSevereCriticalDeadImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducTimeInfectedMild<double>>(),
                                      reducTimeInfectedMild, reducTimeInfectedMild);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecirvvs::Seasonality<double>>(), seasonality_min, seasonality_max);

    // Delta specific parameter
    params.get<mio::osecirvvs::StartDayNewVariant>() = mio::get_day_in_year(mio::Date(2021, 6, 6));

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
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirvvs::Parameters<double>& params,
                                         ScalarType scale_contacts = 1.0, ScalarType share_staying = 1.0)
{
    auto contact_transport_status = mio::read_mobility_plain(data_dir.string() + "//contacts//contacts_transport.txt");
    auto contact_matrix_transport = contact_transport_status.value();
    auto contact_matrices         = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
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
    params.get<mio::osecirvvs::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

mio::IOResult<void> set_contact_matrices_transport(const fs::path& data_dir, mio::osecirvvs::Parameters<double>& params,
                                                   ScalarType scale_contacts = 1.0)
{
    auto contact_transport_status = mio::read_mobility_plain(data_dir.string() + "//contacts//contacts_transport.txt");
    auto contact_matrix_transport = contact_transport_status.value();
    auto contact_matrices         = mio::ContactMatrixGroup(1, size_t(params.get_num_groups()));
    // ScalarType const polymod_share_contacts_transport = 1 / 0.2770885028949545;

    contact_matrices[0].get_baseline() = contact_matrix_transport / scale_contacts;
    contact_matrices[0].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);

    params.get<mio::osecirvvs::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

// reset population in graph
void init_pop_cologne_szenario(mio::Graph<mio::osecirvvs::Model<double>, mio::MigrationParameters<double>>& graph,
                               const int id_cologne)
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
        if (node.id == id_cologne) {
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
                node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptibleNaive}] *= (1 - p);
                node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}] *= (1 - p);
                node.property.populations[{age, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}] *=
                    (1 - p);
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
mio::IOResult<mio::ExtendedGraph<mio::osecirvvs::Model<double>>>
get_graph(const int num_days, const std::string& data_dir, bool masks, bool ffp2, bool szenario_cologne, bool edges)
{
    mio::unused(num_days, masks, ffp2, szenario_cologne, edges);
    std::string travel_times_dir = mio::path_join(data_dir, "mobility", "travel_times_pathes.txt");
    std::string durations_dir    = mio::path_join(data_dir, "mobility", "activity_duration_work.txt");
    auto start_date              = mio::Date(2022, 8, 1);
    auto end_date                = mio::Date(2022, 11, 1);

    // global parameters
    const int num_age_groups = 6;
    mio::osecirvvs::Parameters params(num_age_groups);
    params.get<mio::osecirvvs::StartDay>() = mio::get_day_in_year(start_date);
    auto params_status                     = set_covid_parameters(params);
    auto contacts_status                   = set_contact_matrices(data_dir, params);
    params.get<mio::osecirvvs::StartDay>() = mio::get_day_in_year(start_date);

    // create graph
    mio::ExtendedGraph<mio::osecirvvs::Model<double>> params_graph;

    // set nodes
    auto scaling_factor_infected    = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu         = 1.0;
    auto tnt_capacity_factor        = 1.43 / 100000.;
    const auto& read_function_nodes = mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model<double>>;
    const auto& node_id_function    = mio::get_node_ids;

    auto set_nodes_status = mio::set_nodes<mio::osecirvvs::TestAndTraceCapacity<double>,
                                           mio::osecirvvs::ContactPatterns<double>, mio::osecirvvs::Model<double>,
                                           mio::MigrationParameters<double>, mio::osecirvvs::Parameters<double>>(
        params, start_date, end_date, data_dir,
        mio::path_join(data_dir, "pydata", "Germany", "county_current_population.json"), durations_dir, true,
        params_graph, read_function_nodes, node_id_function, scaling_factor_infected, scaling_factor_icu,
        tnt_capacity_factor, num_days, false);

    if (!set_nodes_status) {
        return set_nodes_status.error();
    }

    // set edges
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
    auto set_edges_status =
        mio::set_edges<ContactLocation, mio::osecirvvs::Model<double>, mio::MigrationParameters<double>,
                       mio::MigrationCoefficientGroup, mio::osecirvvs::InfectionState>(
            travel_times_dir, mio::path_join(data_dir, "mobility", "commuter_migration_with_locals.txt"),
            mio::path_join(data_dir, "mobility", "wegketten_ohne_komma.txt"), params_graph, migrating_compartments,
            contact_locations.size(), std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.});

    if (!set_edges_status) {
        return set_edges_status.error();
    }

    return params_graph;
}
mio::IOResult<void> run(const std::string data_dir)
{
    // mio::set_log_level(mio::LogLevel::critical);
    const auto num_days = 90;
    // const int num_runs  = 12;
    const bool masks = true;
    const bool ffp2  = true;
    const bool edges = true;

    // wenn masks false und ffp2 true, dann error ausgeben
    if (!masks && ffp2) {
        mio::log_error("ffp2 only possible with masks");
    }

    bool szenario_cologne = false;

    // auto params_graph = get_graph(num_days, masks);
    BOOST_OUTCOME_TRY(auto&& created, get_graph(num_days, data_dir, masks, ffp2, szenario_cologne, edges));
    auto params_graph = created;

    // std::string res_dir = "/localdata1/code/memilio/results_paper/t90_f145_mask_" + std::to_string(masks);

    // res_dir += std::string(ffp2 ? "_ffp2" : "") + std::string(szenario_cologne ? "_cologne" : "") +
    //            std::string(!edges ? "_no_edges" : "");

    // if (mio::mpi::is_root())
    //     std::cout << "res_dir = " << res_dir << "\n";

    // // check if boths dir exist, otherwise create them
    // if (!fs::exists(res_dir)) {
    //     fs::create_directory(res_dir);
    // }
    // auto write_graph_status = mio::write_graph(params_graph, "/localdata1/code/memilio/save_graph");

    // std::vector<int> county_ids(params_graph.nodes().size());
    // std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
    //     return n.id;
    // });

    // // parameter study
    // auto parameter_study = mio::ParameterStudy<double, mio::FlowSimulation<mio::osecirvvs::Model<double>>>(
    //     params_graph, 0.0, num_days, 0.01, num_runs);
    // if (mio::mpi::is_root()) {
    //     printf("Seeds: ");
    //     for (auto s : parameter_study.get_rng().get_seeds()) {
    //         printf("%u, ", s);
    //     }
    //     printf("\n");
    // }
    // auto save_single_run_result = mio::IOResult<void>(mio::success());
    // auto ensemble               = parameter_study.run(
    //     [&](auto&& graph) {
    //         return draw_sample(graph, false);
    //     },
    //     [&](auto results_graph, auto&& run_idx) {
    //         auto interpolated_result = mio::interpolate_simulation_result(results_graph);

    //         auto params = std::vector<mio::osecirvvs::Model>();
    //         params.reserve(results_graph.nodes().size());
    //         std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
    //                        [](auto&& node) {
    //                            return node.property.get_simulation().get_model();
    //                        });

    //         auto tp_mob_node = results_graph.nodes()[0].node_pt.get_simulation().get_flows().get_num_time_points();

    //         // transmission mobility models
    //         const std::string results_path_mb = res_dir + "/flows_mb";
    //         if (!fs::exists(results_path_mb)) {
    //             fs::create_directory(results_path_mb);
    //         }
    //         const std::string results_filename =
    //             results_path_mb + "/transmission_mobility_run_" + std::to_string(run_idx) + ".txt";
    //         std::ofstream outfile(results_filename);

    //         outfile << "Node_from Node_to Infections_mobility_node \n";
    //         for (size_t i = 0; i < results_graph.edges().size(); ++i) {
    //             auto& edge = results_graph.edges()[i];
    //             outfile << edge.start_node_idx << " " << edge.end_node_idx << " ";
    //             const auto size_vec = edge.infecions_mobility.size();
    //             for (size_t j = 0; j < size_vec; ++j) {
    //                 outfile << edge.infecions_mobility[j] << " ";
    //             }
    //             outfile << "\n";
    //         }
    //         outfile.close();

    //         // transmission commuter in local models
    //         const std::string results_filename_local_transmission =
    //             results_path_mb + "/transmission_local_run_" + std::to_string(run_idx) + ".txt";
    //         std::ofstream outfile_transmission_local(results_filename_local_transmission);

    //         outfile_transmission_local << "Node_from Node_to Infections_mobility_node \n";
    //         for (size_t i = 0; i < results_graph.edges().size(); ++i) {
    //             auto& edge = results_graph.edges()[i];
    //             outfile_transmission_local << edge.start_node_idx << " " << edge.end_node_idx << " ";
    //             const auto size_vec = edge.infecions_local_model.size();
    //             for (size_t j = 0; j < size_vec; ++j) {
    //                 outfile_transmission_local << edge.infecions_local_model[j] << " ";
    //             }
    //             outfile_transmission_local << "\n";
    //         }
    //         outfile_transmission_local.close();

    //         // symptomatic infections commuter in mobility models
    //         const std::string results_filename_mobility_symp =
    //             results_path_mb + "/symptomps_mobility_run_" + std::to_string(run_idx) + ".txt";
    //         std::ofstream outfile_symp_mobility(results_filename_mobility_symp);

    //         outfile_symp_mobility << "Node_from Node_to Infections_mobility_node \n";
    //         for (size_t i = 0; i < results_graph.edges().size(); ++i) {
    //             auto& edge = results_graph.edges()[i];
    //             outfile_symp_mobility << edge.start_node_idx << " " << edge.end_node_idx << " ";
    //             const auto size_vec = edge.infecions_symptomatic_mobility.size();
    //             for (size_t j = 0; j < size_vec; ++j) {
    //                 outfile_symp_mobility << edge.infecions_symptomatic_mobility[j] << " ";
    //             }
    //             outfile_symp_mobility << "\n";
    //         }
    //         outfile_symp_mobility.close();

    //         // symptomatic infections commuter in local models
    //         const std::string results_filename_local_symp =
    //             results_path_mb + "/symptomps_local_run_" + std::to_string(run_idx) + ".txt";
    //         std::ofstream outfile_symp_local(results_filename_local_symp);

    //         outfile_symp_local << "Node_from Node_to Infections_mobility_node \n";
    //         for (size_t i = 0; i < results_graph.edges().size(); ++i) {
    //             auto& edge = results_graph.edges()[i];
    //             outfile_symp_local << edge.start_node_idx << " " << edge.end_node_idx << " ";
    //             const auto size_vec = edge.infecions_symptomatic_local_model.size();
    //             for (size_t j = 0; j < size_vec; ++j) {
    //                 outfile_symp_local << edge.infecions_symptomatic_local_model[j] << " ";
    //             }
    //             outfile_symp_local << "\n";
    //         }
    //         outfile_symp_local.close();

    //         // symptomatic infections total in local models
    //         const std::string results_filename_total_symp =
    //             results_path_mb + "/total_symptomps_local_run_" + std::to_string(run_idx) + ".txt";
    //         std::ofstream outfile_symp_local_total(results_filename_total_symp);
    //         std::vector<int> indx_symp = {20,  36,  56,  73,  89,  109, 126, 142, 162,
    //                                       179, 195, 215, 232, 248, 268, 285, 301};
    //         outfile_symp_local_total << "Node Infections_mobility_node \n";
    //         for (size_t i = 0; i < results_graph.nodes().size(); ++i) {
    //             auto& node = results_graph.nodes()[i];
    //             auto flows = mio::interpolate_simulation_result(node.property.get_simulation().get_flows());
    //             outfile_symp_local_total << i << " ";
    //             for (auto idx = 0; idx < flows.get_num_time_points(); ++idx) {
    //                 auto flows_step = flows.get_value(idx);
    //                 auto sum_infections =
    //                     std::accumulate(indx_symp.begin(), indx_symp.end(), 0.0, [&flows_step](double sum, int i) {
    //                         return sum + flows_step(i);
    //                     });
    //                 outfile_symp_local_total << sum_infections << " ";
    //             }
    //             outfile_symp_local_total << "\n";
    //         }
    //         outfile_symp_local_total.close();

    //         // transmission infections total in local models
    //         const std::string results_filename_total_transmission =
    //             results_path_mb + "/total_transmission_local_run_" + std::to_string(run_idx) + ".txt";
    //         std::ofstream outfile_tra_local_total(results_filename_total_transmission);
    //         std::vector<int> indx_infections = {0,   33,  17,  53,  86,  70,  106, 139, 123,
    //                                             159, 192, 176, 212, 245, 229, 265, 298, 282};
    //         outfile_tra_local_total << "Node Infections_mobility_node \n";
    //         for (size_t i = 0; i < results_graph.nodes().size(); ++i) {
    //             auto& node = results_graph.nodes()[i];
    //             auto flows = mio::interpolate_simulation_result(node.property.get_simulation().get_flows());
    //             outfile_tra_local_total << i << " ";
    //             for (auto idx = 0; idx < flows.get_num_time_points(); ++idx) {
    //                 auto flows_step     = flows.get_value(idx);
    //                 auto sum_infections = std::accumulate(indx_infections.begin(), indx_infections.end(), 0.0,
    //                                                       [&flows_step](double sum, int i) {
    //                                                           return sum + flows_step(i);
    //                                                       });
    //                 outfile_tra_local_total << sum_infections << " ";
    //             }
    //             outfile_tra_local_total << "\n";
    //         }
    //         outfile_tra_local_total.close();

    //         std::cout << "Run " << run_idx << " complete." << std::endl;

    //         return std::make_pair(interpolated_result, params);
    //     });

    // if (ensemble.size() > 0) {
    //     auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
    //     ensemble_results.reserve(ensemble.size());
    //     auto ensemble_params = std::vector<std::vector<mio::osecirvvs::Model<double>>>{};
    //     ensemble_params.reserve(ensemble.size());
    //     for (auto&& run : ensemble) {
    //         ensemble_results.emplace_back(std::move(run.first));
    //         ensemble_params.emplace_back(std::move(run.second));
    //     }
    //     BOOST_OUTCOME_TRY(mio::save_results(ensemble_results, ensemble_params, county_ids, res_dir, false));
    // }

    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    mio::mpi::init();

    const std::string data_dir = "/localdata1/code/memilio/data";

    auto result = run(data_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        mio::mpi::finalize();
        return -1;
    }
    mio::mpi::finalize();
    return 0;
}