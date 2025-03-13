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
#include "ode_secirts/parameters.h"
#include "ode_secirts/parameters_io.h"
#include "ode_secirts/parameter_space.h"
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
mio::IOResult<void> set_covid_parameters(mio::osecirts::Parameters<double>& params)
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

    array_assign_uniform_distribution(params.get<mio::osecirts::TimeExposed<double>>(), timeExposedMin, timeExposedMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::TimeInfectedNoSymptoms<double>>(),
                                      timeInfectedNoSymptomsMin, timeInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::TimeInfectedSymptoms<double>>(),
                                      timeInfectedSymptomsMin, timeInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::TimeInfectedSevere<double>>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::TimeInfectedCritical<double>>(),
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

    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10347449/pdf/bmjgh-2023-012328.pdf
    const double fac_dead               = 0.39;
    const double deathsPerCriticalMin[] = {fac_dead * 0.00, fac_dead * 0.00, fac_dead * 0.10,
                                           fac_dead * 0.10, fac_dead * 0.30, fac_dead * 0.5}; // 2021 paper
    const double deathsPerCriticalMax[] = {fac_dead * 0.10, fac_dead * 0.10, fac_dead * 0.18,
                                           fac_dead * 0.18, fac_dead * 0.50, fac_dead * 0.7};

    const double reducExposedPartialImmunityMin = 1.0;
    const double reducExposedPartialImmunityMax = 1.0;

    const double reducExposedImprovedImmunityMin = 1.0;
    const double reducExposedImprovedImmunityMax = 1.0;

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

    array_assign_uniform_distribution(params.get<mio::osecirts::TransmissionProbabilityOnContact<double>>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::RelativeTransmissionNoSymptoms<double>>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<double>>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::SeverePerInfectedSymptoms<double>>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::CriticalPerSevere<double>>(), criticalPerSevereMin,
                                      criticalPerSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::DeathsPerCritical<double>>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax);

    array_assign_uniform_distribution(params.get<mio::osecirts::ReducExposedPartialImmunity<double>>(),
                                      reducExposedPartialImmunityMin, reducExposedPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::ReducExposedImprovedImmunity<double>>(),
                                      reducExposedImprovedImmunityMin, reducExposedImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>(),
                                      reducInfectedSymptomsPartialImmunityMin, reducInfectedSymptomsPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>(),
                                      reducInfectedSymptomsImprovedImmunityMin,
                                      reducInfectedSymptomsImprovedImmunityMax);
    array_assign_uniform_distribution(
        params.get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(),
        reducInfectedSevereCriticalDeadPartialImmunityMin, reducInfectedSevereCriticalDeadPartialImmunityMax);
    array_assign_uniform_distribution(
        params.get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(),
        reducInfectedSevereCriticalDeadImprovedImmunityMin, reducInfectedSevereCriticalDeadImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirts::ReducTimeInfectedMild<double>>(), reducTimeInfectedMild,
                                      reducTimeInfectedMild);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecirts::Seasonality<double>>(), seasonality_min, seasonality_max);

    // Delta specific parameter
    params.get<mio::osecirts::StartDayNewVariant>() = mio::get_day_in_year(mio::Date(2022, 6, 6));

    const double ImmunityInterval1Min = 60; ///https://doi.org/10.1016/S1473-3099(22)00801-5
    const double ImmunityInterval1Max = 60;

    array_assign_uniform_distribution(params.get<mio::osecirts::TimeTemporaryImmunityPI<double>>(),
                                      ImmunityInterval1Min, ImmunityInterval1Max);

    const double ImmunityInterval2Min = 60; // https://doi.org/10.1016/S1473-3099(22)00801-5
    const double ImmunityInterval2Max = 60;

    array_assign_uniform_distribution(params.get<mio::osecirts::TimeTemporaryImmunityII<double>>(),
                                      ImmunityInterval2Min, ImmunityInterval2Max);

    // https://doi.org/10.1016/S1473-3099(22)00801-5
    params.get<mio::osecirts::TimeWaningPartialImmunity<double>>()  = 365.0;
    params.get<mio::osecirts::TimeWaningImprovedImmunity<double>>() = 365.0;

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
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirts::Parameters<double>& params,
                                         ScalarType avg_tt = 0.0, ScalarType share_staying = 1.0)
{
    auto contact_transport_status = mio::read_mobility_plain(data_dir.string() + "//contacts//contacts_transport.txt");
    auto contact_matrix_transport = contact_transport_status.value();
    auto contact_matrices         = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        // the other contact matrix containts also the transport contact matrix. Therefore, we need to subtract it.
        if (contact_location.second == "other") {
            baseline = abs(baseline - contact_matrix_transport);
        }
        // Scale the number of contacts
        auto scaled_baseline = (1 - share_staying) * baseline / (1 - avg_tt) + share_staying * baseline;
        // Because we only model the acitvity from commuters, some age group should just have the unscaled contacts
        Eigen::MatrixXd baseline_adjusted = Eigen::MatrixXd::Zero(scaled_baseline.rows(), scaled_baseline.cols());
        for (auto i = 0; i < 6; i++) {
            for (auto j = 0; j < 6; j++) {
                if ((i >= 2 && i <= 5) || (j >= 2 && j <= 5)) {
                    baseline_adjusted(i, j) = scaled_baseline(i, j);
                }
                else {
                    baseline_adjusted(i, j) = baseline(i, j);
                }
            }
        }
        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline_adjusted;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);
    }
    params.get<mio::osecirts::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

mio::IOResult<void> set_contact_matrices_transport(const fs::path& data_dir, mio::osecirts::Parameters<double>& params)
{
    auto contact_transport_status = mio::read_mobility_plain(data_dir.string() + "//contacts//contacts_transport.txt");
    auto contact_matrix_transport = contact_transport_status.value();
    auto contact_matrices         = mio::ContactMatrixGroup(1, size_t(params.get_num_groups()));
    // ScalarType const polymod_share_contacts_transport = 1 / 0.2770885028949545;

    contact_matrices[0].get_baseline() = contact_matrix_transport;
    contact_matrices[0].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);

    params.get<mio::osecirts::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

mio::IOResult<void> scale_contacts_local(mio::ExtendedGraph<mio::osecirts::Model<double>>& params_graph,
                                         const fs::path& data_dir, const std::string mobility_data_dir,
                                         const std::vector<ScalarType> commuting_weights)
{
    BOOST_OUTCOME_TRY(auto&& mobility_data_commuter, mio::read_mobility_plain(mobility_data_dir));
    // average travel time
    const ScalarType avg_traveltime = std::accumulate(params_graph.edges().begin(), params_graph.edges().end(), 0.0,
                                                      [](double sum, const auto& e) {
                                                          return sum + e.property.travel_time;
                                                      }) /
                                      params_graph.edges().size();

    // average share of commuters for all counties relative to the total population
    const ScalarType total_population = std::accumulate(
        params_graph.nodes().begin(), params_graph.nodes().end(), 0.0, [commuting_weights](double sum, const auto& n) {
            return sum + n.property.base_sim.populations.get_group_total(mio::AgeGroup(0)) * commuting_weights[0] +
                   n.property.base_sim.populations.get_group_total(mio::AgeGroup(1)) * commuting_weights[1] +
                   n.property.base_sim.populations.get_group_total(mio::AgeGroup(2)) * commuting_weights[2] +
                   n.property.base_sim.populations.get_group_total(mio::AgeGroup(3)) * commuting_weights[3] +
                   n.property.base_sim.populations.get_group_total(mio::AgeGroup(4)) * commuting_weights[4] +
                   n.property.base_sim.populations.get_group_total(mio::AgeGroup(5)) * commuting_weights[5];
        });
    const ScalarType num_commuters      = std::accumulate(params_graph.edges().begin(), params_graph.edges().end(), 0.0,
                                                     [mobility_data_commuter](double sum, const auto& e) {
                                                         auto start_node = e.start_node_idx;
                                                         auto end_node   = e.end_node_idx;
                                                         return sum + mobility_data_commuter(start_node, end_node);
                                                     });
    const ScalarType avg_commuter_share = num_commuters / total_population;

    if (mio::mpi::is_root()) {
        std::cout << "avg_commuter_share: " << avg_commuter_share << std::endl;
        std::cout << "avg_traveltime: " << avg_traveltime << std::endl;
    }

    // scale all contact matrices
    for (auto& node : params_graph.nodes()) {
        BOOST_OUTCOME_TRY(
            set_contact_matrices(data_dir, node.property.base_sim.parameters, avg_traveltime, avg_commuter_share));
    }
    return mio::success();
}

// reset population in graph
void init_pop_cologne_szenario(mio::Graph<mio::osecirts::Model<double>, mio::MobilityParameters<double>>& graph,
                               const int id_cologne)
{
    std::vector<std::vector<double>> immunity = {{0.04, 0.61, 0.35}, {0.04, 0.61, 0.35},   {0.075, 0.62, 0.305},
                                                 {0.08, 0.62, 0.3},  {0.035, 0.58, 0.385}, {0.01, 0.41, 0.58}};

    for (auto& node : graph.nodes()) {
        for (auto age = mio::AgeGroup(0); age < mio::AgeGroup(6); age++) {
            auto pop_age = 0.0;
            for (auto inf_state = mio::Index<mio::osecirts::InfectionState>(0);
                 inf_state < mio::osecirts::InfectionState::Count; ++inf_state) {
                pop_age += node.property.populations[{age, inf_state}];
                node.property.populations[{age, inf_state}] = 0.0;
            }
            size_t immunity_index = static_cast<size_t>(age);
            node.property.populations[{age, mio::osecirts::InfectionState::SusceptibleNaive}] =
                pop_age * immunity[immunity_index][0];
            node.property.populations[{age, mio::osecirts::InfectionState::SusceptiblePartialImmunity}] =
                pop_age * immunity[immunity_index][1];
            node.property.populations[{age, mio::osecirts::InfectionState::SusceptibleImprovedImmunity}] =
                pop_age * immunity[immunity_index][2];
        }
        if (node.id == id_cologne) {
            // infect p% of population
            ScalarType p = 0.05;
            for (auto age = mio::AgeGroup(0); age < graph.nodes()[0].property.parameters.get_num_groups(); age++) {
                node.property.populations[{mio::AgeGroup(age), mio::osecirts::InfectionState::InfectedSymptomsNaive}] =
                    node.property.populations[{age, mio::osecirts::InfectionState::SusceptibleNaive}] * p;
                node.property.populations[{age, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity}] =
                    node.property.populations[{age, mio::osecirts::InfectionState::SusceptiblePartialImmunity}] * p;
                node.property.populations[{mio::AgeGroup(age),
                                           mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity}] =
                    node.property.populations[{age, mio::osecirts::InfectionState::SusceptibleImprovedImmunity}] * p;
                node.property.populations[{age, mio::osecirts::InfectionState::SusceptibleNaive}] *= (1 - p);
                node.property.populations[{age, mio::osecirts::InfectionState::SusceptiblePartialImmunity}] *= (1 - p);
                node.property.populations[{age, mio::osecirts::InfectionState::SusceptibleImprovedImmunity}] *= (1 - p);
            }
        }
    }
}

/**
 * @brief Sets the graph nodes for counties or districts.
 * Reads the node ids which could refer to districts or counties and the epidemiological
 * data from json files and creates one node for each id. Every node contains a model.
 * @param[in] params Model Parameters that are used for every node.
 * @param[in] start_date Start date for which the data should be read.
 * @param[in] end_data End date for which the data should be read.
 * @param[in] data_dir Directory that contains the data files.
 * @param[in] population_data_path Path to json file containing the population data.
 * @param[in] stay_times_data_path Path to txt file containing the stay times for the considered local entities.
 * @param[in] is_node_for_county Specifies whether the node ids should be county ids (true) or district ids (false).
 * @param[in, out] params_graph Graph whose nodes are set by the function.
 * @param[in] read_func Function that reads input data for german counties and sets Model compartments.
 * @param[in] node_func Function that returns the county ids.
 * @param[in] scaling_factor_inf Factor of confirmed cases to account for undetected cases in each county.
 * @param[in] scaling_factor_icu Factor of ICU cases to account for underreporting.
 * @param[in] tnt_capacity_factor Factor for test and trace capacity.
 * @param[in] num_days Number of days to be simulated; required to load data for vaccinations during the simulation.
 * @param[in] export_time_series If true, reads data for each day of simulation and writes it in the same directory as the input files.
 */
template <class ReadFunction, class NodeIdFunction, typename FP = double>
mio::IOResult<void>
set_nodes(const mio::osecirts::Parameters<FP>& params, mio::Date start_date, mio::Date end_date,
          const fs::path& data_dir, const std::string& population_data_path, const std::string& stay_times_data_path,
          bool is_node_for_county, mio::ExtendedGraph<mio::osecirts::Model<FP>>& params_graph, ReadFunction&& read_func,
          NodeIdFunction&& node_func, const std::vector<FP>& scaling_factor_inf, FP scaling_factor_icu,
          FP tnt_capacity_factor, const std::vector<std::vector<double>> immunity_population, int num_days = 0,
          bool export_time_series = false, bool rki_age_groups = true, bool masks = false, bool ffp2 = false)
{

    BOOST_OUTCOME_TRY(auto&& duration_stay, mio::read_duration_stay(stay_times_data_path));

    BOOST_OUTCOME_TRY(auto&& node_ids, node_func(population_data_path, is_node_for_county, rki_age_groups));
    std::vector<mio::osecirts::Model<FP>> nodes(node_ids.size(),
                                                mio::osecirts::Model<FP>(int(size_t(params.get_num_groups()))));
    for (auto& node : nodes) {
        node.parameters = params;
    }

    BOOST_OUTCOME_TRY(read_func(nodes, start_date, node_ids, scaling_factor_inf, scaling_factor_icu, data_dir.string(),
                                num_days, immunity_population, export_time_series));

    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {

        auto tnt_capacity = nodes[node_idx].populations.get_total() * tnt_capacity_factor;

        //local parameters
        auto& tnt_value = nodes[node_idx].parameters.template get<mio::osecirts::TestAndTraceCapacity<FP>>();
        tnt_value       = mio::UncertainValue(0.5 * (1.2 * tnt_capacity + 0.8 * tnt_capacity));
        tnt_value.set_distribution(mio::ParameterDistributionUniform(0.8 * tnt_capacity, 1.2 * tnt_capacity));

        //holiday periods
        auto id              = int(mio::regions::CountyId(node_ids[node_idx]));
        auto holiday_periods = mio::regions::get_holidays(mio::regions::get_state_id(id), start_date, end_date);
        auto& contacts       = nodes[node_idx].parameters.template get<mio::osecirts::ContactPatterns<FP>>();
        contacts.get_school_holidays() =
            std::vector<std::pair<mio::SimulationTime, mio::SimulationTime>>(holiday_periods.size());
        std::transform(
            holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(), [=](auto& period) {
                return std::make_pair(mio::SimulationTime(mio::get_offset_in_days(period.first, start_date)),
                                      mio::SimulationTime(mio::get_offset_in_days(period.second, start_date)));
            });

        //uncertainty in populations
        for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = mio::Index<typename mio::osecirts::Model<FP>::Compartments>(0);
                 j < mio::osecirts::Model<FP>::Compartments::Count; ++j) {
                auto& compartment_value = nodes[node_idx].populations[{i, j}];
                compartment_value =
                    mio::UncertainValue(0.5 * (1.1 * double(compartment_value) + 0.9 * double(compartment_value)));
                compartment_value.set_distribution(mio::ParameterDistributionUniform(0.9 * double(compartment_value),
                                                                                     1.1 * double(compartment_value)));
            }
        }

        // Add mobility node
        auto mobility_model = nodes[node_idx];
        mobility_model.populations.set_total(0);

        // reduce transmission on contact due to mask obligation in mobility node
        if (masks) {

            const double fact_surgical_mask = 0.1;
            const double fact_ffp2          = 0.001;

            // first age group not able to (properly) wear masks, second age group only partially
            double factor_mask[] = {1, fact_surgical_mask, fact_ffp2, fact_ffp2, fact_ffp2, fact_ffp2};
            if (!ffp2) {
                for (size_t j = 2; j < 6; j++) {
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
                mobility_model.parameters.template get<mio::osecirts::TransmissionProbabilityOnContact<double>>(),
                transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
        }

        // no vaccination in mobility node
        for (int t_idx = 0; t_idx < num_days; ++t_idx) {
            auto t = mio::SimulationDay((size_t)t_idx);
            for (auto j = mio::AgeGroup(0); j < params.get_num_groups(); j++) {
                mobility_model.parameters.template get<mio::osecirts::DailyPartialVaccinations<double>>()[{j, t}] = 0;
                mobility_model.parameters.template get<mio::osecirts::DailyFullVaccinations<double>>()[{j, t}]    = 0;
                mobility_model.parameters.template get<mio::osecirts::DailyBoosterVaccinations<double>>()[{j, t}] = 0;
            }
        }

        params_graph.add_node(node_ids[node_idx], nodes[node_idx], mobility_model,
                              duration_stay((Eigen::Index)node_idx));
    }
    return mio::success();
}

/**
 * @brief Sets the graph edges.
 * Reads the commuting matrices, travel times and paths from data and creates one edge for each pair of nodes.
 * @param[in] travel_times_path Path to txt file containing the travel times between counties.
 * @param[in] mobility_data_path Path to txt file containing the commuting matrices.
 * @param[in] travelpath_path Path to txt file containing the paths between counties.
 * @param[in, out] params_graph Graph whose nodes are set by the function.
 * @param[in] migrating_compartments Compartments that commute.
 * @param[in] contact_locations_size Number of contact locations.
 * @param[in] commuting_weights Vector with a commuting weight for every AgeGroup.
 */
mio::IOResult<void>
set_edges(const std::string& travel_times_dir, const std::string mobility_data_dir, const std::string& travel_path_dir,
          mio::ExtendedGraph<mio::osecirts::Model<double>>& params_graph,
          std::initializer_list<mio::osecirts::InfectionState>& migrating_compartments, size_t contact_locations_size,
          std::vector<ScalarType> commuting_weights = std::vector<ScalarType>{}, ScalarType theshold_edges = 4e-5)
{
    BOOST_OUTCOME_TRY(auto&& mobility_data_commuter, mio::read_mobility_plain(mobility_data_dir));
    BOOST_OUTCOME_TRY(auto&& travel_times, mio::read_mobility_plain(travel_times_dir));
    BOOST_OUTCOME_TRY(auto&& path_mobility, mio::read_path_mobility(travel_path_dir));

    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            auto& populations = params_graph.nodes()[county_idx_i].property.base_sim.populations;

            // mobility coefficients have the same number of components as the contact matrices.
            // so that the same NPIs/dampings can be used for both (e.g. more home office => fewer commuters)
            auto mobility_coeffs = mio::MobilityCoefficientGroup(contact_locations_size, populations.numel());
            auto num_age_groups =
                (size_t)params_graph.nodes()[county_idx_i].property.base_sim.parameters.get_num_groups();
            commuting_weights =
                (commuting_weights.size() == 0 ? std::vector<ScalarType>(num_age_groups, 1.0) : commuting_weights);

            //commuters
            auto working_population = 0.0;
            auto min_commuter_age   = mio::AgeGroup(2);
            auto max_commuter_age   = mio::AgeGroup(4); //this group is partially retired, only partially commutes
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                working_population += populations.get_group_total(age) * commuting_weights[size_t(age)];
            }
            auto commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) / working_population;
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_index = populations.get_flat_index({age, compartment});
                    mobility_coeffs[size_t(ContactLocation::Work)].get_baseline()[coeff_index] =
                        commuter_coeff_ij * commuting_weights[size_t(age)];
                }
            }

            auto path = path_mobility[county_idx_i][county_idx_j];
            if (static_cast<size_t>(path[0]) != county_idx_i ||
                static_cast<size_t>(path[path.size() - 1]) != county_idx_j)
                std::cout << "Wrong Path for edge " << county_idx_i << " " << county_idx_j << "\n";

            //only add edges with mobility above thresholds for performance
            if (commuter_coeff_ij > theshold_edges) {
                params_graph.add_edge(county_idx_i, county_idx_j, std::move(mobility_coeffs),
                                      travel_times(county_idx_i, county_idx_j),
                                      path_mobility[county_idx_i][county_idx_j]);
            }
        }
    }

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
mio::IOResult<mio::ExtendedGraph<mio::osecirts::Model<double>>> get_graph(const mio::Date start_date,
                                                                          const mio::Date end_date, const int num_days,
                                                                          const std::string& data_dir, bool masks,
                                                                          bool ffp2, bool szenario_cologne, bool edges)
{
    mio::unused(szenario_cologne);
    std::string travel_times_dir = mio::path_join(data_dir, "mobility", "travel_times_pathes.txt");
    std::string durations_dir    = mio::path_join(data_dir, "mobility", "activity_duration_work.txt");

    const std::vector<std::vector<double>> immunity_population = {{0.04, 0.04, 0.075, 0.08, 0.035, 0.01},
                                                                  {0.61, 0.61, 0.62, 0.62, 0.58, 0.41},
                                                                  {0.35, 0.35, 0.305, 0.3, 0.385, 0.58}};

    // global parameters
    const int num_age_groups = 6;
    mio::osecirts::Parameters params(num_age_groups);
    params.get<mio::osecirts::StartDay>() = mio::get_day_in_year(start_date);
    auto params_status                    = set_covid_parameters(params);
    auto contacts_status                  = set_contact_matrices(data_dir, params);
    params.get<mio::osecirts::StartDay>() = mio::get_day_in_year(start_date);

    // create graph
    mio::ExtendedGraph<mio::osecirts::Model<double>> params_graph;

    // set nodes
    auto scaling_factor_infected    = std::vector<double>(size_t(params.get_num_groups()), 1.7);
    auto scaling_factor_icu         = 1.0;
    auto tnt_capacity_factor        = 1.43 / 100000.;
    const auto& read_function_nodes = mio::osecirts::read_input_data_county<mio::osecirts::Model<double>>;
    const auto& node_id_function    = mio::get_node_ids;

    auto set_nodes_status = set_nodes<decltype(read_function_nodes), decltype(node_id_function)>(
        params, start_date, end_date, data_dir,
        mio::path_join(data_dir, "pydata", "Germany", "county_current_population.json"), durations_dir, true,
        params_graph, read_function_nodes, node_id_function, scaling_factor_infected, scaling_factor_icu,
        tnt_capacity_factor, immunity_population, num_days, false, true, ffp2, masks);

    if (!set_nodes_status) {
        return set_nodes_status.error();
    }

    // iterate over all nodes and set the contact_matrix for the mobility models
    for (auto& node : params_graph.nodes()) {
        BOOST_OUTCOME_TRY(set_contact_matrices_transport(data_dir, node.property.mobility_sim.parameters));
    }

    // set edges
    if (edges) {
        auto migrating_compartments = {mio::osecirts::InfectionState::SusceptibleNaive,
                                       mio::osecirts::InfectionState::ExposedNaive,
                                       mio::osecirts::InfectionState::InfectedNoSymptomsNaive,
                                       mio::osecirts::InfectionState::InfectedSymptomsNaive,
                                       mio::osecirts::InfectionState::SusceptibleImprovedImmunity,
                                       mio::osecirts::InfectionState::SusceptiblePartialImmunity,
                                       mio::osecirts::InfectionState::ExposedPartialImmunity,
                                       mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity,
                                       mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity,
                                       mio::osecirts::InfectionState::ExposedImprovedImmunity,
                                       mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity,
                                       mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity};
        auto set_edges_status       = set_edges(
            travel_times_dir, mio::path_join(data_dir, "mobility", "commuter_migration_with_locals.txt"),
            mio::path_join(data_dir, "mobility", "wegketten_ohne_komma.txt"), params_graph, migrating_compartments,
            contact_locations.size(), std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.});

        if (!set_edges_status) {
            return set_edges_status.error();
        }

        // if we have edges/mobility, we also need to scale the local contacts
        BOOST_OUTCOME_TRY(scale_contacts_local(
            params_graph, data_dir, mio::path_join(data_dir, "mobility", "commuter_migration_with_locals.txt"),
            std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.}));
    }

    return params_graph;
}
mio::IOResult<void> run(const std::string data_dir, std::string res_dir, const int num_runs, const int num_days)
{
    // mio::set_log_level(mio::LogLevel::critical);
    auto start_date       = mio::Date(2022, 8, 1);
    auto end_date         = mio::Date(2022, 11, 1);
    constexpr bool masks  = true;
    constexpr bool ffp2   = true;
    const bool edges      = true;
    bool szenario_cologne = false;

    // wenn masks false und ffp2 true, dann error ausgeben
    if constexpr (!masks && ffp2) {
        mio::log_error("ffp2 only possible with masks");
    }

    // auto params_graph = get_graph(num_days, masks);
    BOOST_OUTCOME_TRY(auto&& created,
                      get_graph(start_date, end_date, num_days, data_dir, masks, ffp2, szenario_cologne, edges));
    auto params_graph = created;

    res_dir += "/masks_" + std::to_string(masks) + std::string(ffp2 ? "_ffp2" : "") +
               std::string(szenario_cologne ? "_cologne" : "") + std::string(!edges ? "_no_edges" : "");

    if (mio::mpi::is_root())
        std::cout << "res_dir = " << res_dir << "\n";

    // check if res_dir dir exist, otherwise create
    if (!fs::exists(res_dir)) {
        fs::create_directories(res_dir);
    }

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    // mio::ExtendedGraph<mio::osecirts::Simulation<ScalarType, mio::FlowSimulation<ScalarType, mio::osecirts::Model<double>>>
    // parameter study
    auto parameter_study = mio::ParameterStudy<
        mio::osecirts::Simulation<double, mio::FlowSimulation<double, mio::osecirts::Model<double>>>,
        mio::ExtendedGraph<mio::osecirts::Model<double>>,
        mio::ExtendedGraph<
            mio::osecirts::Simulation<double, mio::FlowSimulation<double, mio::osecirts::Model<double>>>>>(
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
            std::vector<mio::TimeSeries<double>> interpolated_result;
            interpolated_result.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                           std::back_inserter(interpolated_result), [](auto& n) {
                               return interpolate_simulation_result(n.property.base_sim.get_result());
                           });

            auto params = std::vector<mio::osecirts::Model<double>>();
            params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
                           [](auto&& node) {
                               return node.property.base_sim.get_model();
                           });

            auto flows = std::vector<mio::TimeSeries<ScalarType>>{};
            flows.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(flows),
                           [](auto&& node) {
                               auto& flow_node         = node.property.base_sim.get_flows();
                               auto interpolated_flows = mio::interpolate_simulation_result(flow_node);
                               return interpolated_flows;
                           });

            // same for the mobility model (flows are enough here)
            auto params_mobility = std::vector<mio::osecirts::Model<double>>();
            params_mobility.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                           std::back_inserter(params_mobility), [](auto&& node) {
                               return node.property.mobility_sim.get_model();
                           });

            auto flows_mobility = std::vector<mio::TimeSeries<ScalarType>>{};
            flows_mobility.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                           std::back_inserter(flows_mobility), [](auto&& node) {
                               auto& flow_node         = node.property.mobility_sim.get_flows();
                               auto interpolated_flows = mio::interpolate_simulation_result(flow_node);
                               return interpolated_flows;
                           });

            std::cout << "Run " << run_idx << " complete." << std::endl;

            return std::make_tuple(interpolated_result, params, flows, params_mobility, flows_mobility);
        });

    if (ensemble.size() > 0) {
        auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_results.reserve(ensemble.size());
        auto ensemble_params = std::vector<std::vector<mio::osecirts::Model<double>>>{};
        ensemble_params.reserve(ensemble.size());
        auto ensemble_flows = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_flows.reserve(ensemble.size());
        auto ensemble_params_mobility = std::vector<std::vector<mio::osecirts::Model<double>>>{};
        ensemble_params_mobility.reserve(ensemble.size());
        auto ensemble_flows_mobility = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_flows_mobility.reserve(ensemble.size());
        for (auto&& run : ensemble) {
            ensemble_results.emplace_back(std::move(std::get<0>(run)));
            ensemble_params.emplace_back(std::move(std::get<1>(run)));
            ensemble_flows.emplace_back(std::move(std::get<2>(run)));
            ensemble_params_mobility.emplace_back(std::move(std::get<3>(run)));
            ensemble_flows_mobility.emplace_back(std::move(std::get<4>(run)));
        }

        BOOST_OUTCOME_TRY(mio::save_results(ensemble_results, ensemble_params, county_ids, res_dir, false));
        auto result_dir_run_flows = res_dir + "/flows";
        if (mio::mpi::is_root()) {
            boost::filesystem::create_directories(result_dir_run_flows);
            printf("Saving Flow results to \"%s\".\n", result_dir_run_flows.c_str());
        }
        BOOST_OUTCOME_TRY(save_results(ensemble_flows, ensemble_params, county_ids, result_dir_run_flows, false));

        auto result_dir_mobility           = res_dir + "/mobility";
        auto result_dir_run_flows_mobility = result_dir_mobility + "/flows";
        if (mio::mpi::is_root()) {
            boost::filesystem::create_directories(result_dir_run_flows_mobility);
            printf("Saving Flow results to \"%s\".\n", result_dir_run_flows_mobility.c_str());
        }
        BOOST_OUTCOME_TRY(save_results(ensemble_flows_mobility, ensemble_params_mobility, county_ids,
                                       result_dir_run_flows_mobility, false));
    }

    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::warn);
    mio::mpi::init();

    std::string data_dir    = "";
    std::string results_dir = "";
    int num_days            = 10;
    int num_runs            = 2;

    if (argc == 5) {
        // Full specification: <data_dir> <results_dir> <num_runs> <num_days>
        data_dir    = argv[1];
        results_dir = argv[2];
        num_runs    = atoi(argv[3]);
        num_days    = atoi(argv[4]);

        if (mio::mpi::is_root()) {
            printf("Reading data from \"%s\", saving results to \"%s\".\n", data_dir.c_str(), results_dir.c_str());
            printf("Number of runs: %d, Number of days: %d.\n", num_runs, num_days);
        }
    }
    else if (argc == 3) {
        // Partial specification: <data_dir> <results_dir>, use defaults for num_runs and num_days
        data_dir    = argv[1];
        results_dir = argv[2];

        if (mio::mpi::is_root()) {
            printf("Reading data from \"%s\", saving results to \"%s\".\n", data_dir.c_str(), results_dir.c_str());
            printf("Using default values - Number of runs: %d, Number of days: %d.\n", num_runs, num_days);
        }
    }
    else {
        if (mio::mpi::is_root()) {
            printf("Usage:\n");
            printf("2022_omicron_late_phase_mobility <data_dir> <results_dir> <num_runs> <num_days>\n");
            printf("\tRun simulation with data from <data_dir>, saving results to <results_dir>.\n");
            printf("\t<num_runs>: Number of simulation runs.\n");
            printf("\t<num_days>: Number of days to simulate.\n");
            printf("2022_omicron_late_phase_mobility <data_dir> <results_dir>\n");
            printf("\tRun simulation with default num_runs=%d and num_days=%d.\n", num_runs, num_days);
        }
        mio::mpi::finalize();
        return 0;
    }

    // Run the simulation
    auto result = run(data_dir, results_dir, num_runs, num_days);
    if (!result) {
        if (mio::mpi::is_root()) {
            printf("%s\n", result.error().formatted_message().c_str());
        }
        mio::mpi::finalize();
        return -1;
    }

    mio::mpi::finalize();
    return 0;
}
