/**
* Simulation application that was used to produce results of the following publication(s):
* W. Koslow et al, 2022: Appropriate relaxation of non-pharmaceutical interventions minimizes 
* the risk of a resurgence in SARS-CoV-2 infections in spite of the Delta variant
* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Wadim Koslow, Martin Kühn
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
#include "memilio/epidemiology/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/date.h"
#include "ode_secirvvs/parameters.h"
#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/parameter_space.h"
#include "boost/filesystem.hpp"
#include <cstdio>
#include <iomanip>

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
mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters& params, int vacc_effectiveness_szenario)
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
    const double t_imm_min = {30}; // heuristic
    const double t_imm_max = {30}; // heuristic

    array_assign_uniform_distribution(params.get<mio::osecirvvs::IncubationTime>(), incubationTime, incubationTime);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SerialInterval>(), serialIntervalMin,
                                      serialIntervalMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ImmunityInterval1>(), t_imm_min, t_imm_max);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ImmunityInterval2>(), t_imm_min, t_imm_max);
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

    mio::unused(vacc_effectiveness_szenario);
    // if (vacc_effectiveness_szenario == 1) {
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
    // }
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

    params.get<mio::osecirvvs::WaningPartialImmunity>()  = 90.0;
    params.get<mio::osecirvvs::WaningImprovedImmunity>() = 90.0;

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
 * Set NPIs.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param params Object that the NPIs will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_npis(mio::Date start_date, mio::Date end_date, mio::osecirvvs::Parameters& params)
{
    auto& contacts = params.get<mio::osecirvvs::ContactPatterns>();

    mio::unused(start_date);
    mio::unused(end_date);
    mio::unused(params);
    mio::unused(contacts);
    // auto& contact_dampings = contacts.get_dampings();

    // if (test) {
    //     params.get_commuter_nondetection() = 0.85;
    // }
    // else {
    //     params.get_commuter_nondetection() = 1.0;
    // }

    // //weights for age groups affected by an NPI
    // auto group_weights_all     = Eigen::VectorXd::Constant(size_t(params.get_num_groups()), 1.0);
    // auto group_weights_seniors = Eigen::VectorXd::NullaryExpr(size_t(params.get_num_groups()), [](auto&& i) {
    //     return i == 5 ? 1.0 : i == 4 ? 0.5 : 0.0; //65-80 only partially
    // });

    // //helper functions that create dampings for specific NPIs
    // auto contacts_at_home = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::Main)),
    //                                 mio::DampingType(int(Intervention::Home)), t, {size_t(ContactLocation::Home)},
    //                                 group_weights_all);
    // };
    // auto school_closure = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::Main)),
    //                                 mio::DampingType(int(Intervention::SchoolClosure)), t,
    //                                 {size_t(ContactLocation::School)}, group_weights_all);
    // };
    // auto home_office = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::Main)),
    //                                 mio::DampingType(int(Intervention::HomeOffice)), t, {size_t(ContactLocation::Work)},
    //                                 group_weights_all);
    // };
    // auto social_events = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::Main)),
    //                                 mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), t,
    //                                 {size_t(ContactLocation::Other)}, group_weights_all);
    // };
    // auto social_events_work = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::Main)),
    //                                 mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), t,
    //                                 {size_t(ContactLocation::Work)}, group_weights_all);
    // };
    // auto physical_distancing_home = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
    //                                 mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
    //                                 {size_t(ContactLocation::Home)}, group_weights_all);
    // };
    // auto physical_distancing_school = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
    //                                 mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
    //                                 {size_t(ContactLocation::School)}, group_weights_all);
    // };
    // auto physical_distancing_work = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
    //                                 mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
    //                                 {size_t(ContactLocation::Work)}, group_weights_all);
    // };
    // auto physical_distancing_other = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
    //                                 mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
    //                                 {size_t(ContactLocation::Other)}, group_weights_all);
    // };
    // auto senior_awareness = [=](auto t, auto min, auto max) {
    //     auto v = mio::UncertainValue();
    //     assign_uniform_distribution(v, min, max);
    //     return mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::SeniorAwareness)),
    //                                 mio::DampingType(int(Intervention::SeniorAwareness)), t,
    //                                 {size_t(ContactLocation::Home), size_t(ContactLocation::Other)},
    //                                 group_weights_seniors);
    // };

    // //OPEN SCENARIO SPRING
    // auto start_year = mio::Date(2022, 1, 1);
    // double narrow   = 0.05;
    // if (start_year < end_date) {
    //     auto static_open_scenario_spring = mio::SimulationTime(mio::get_offset_in_days(start_year, start_date));
    //     contact_dampings.push_back(contacts_at_home(static_open_scenario_spring, 0.0, 0.0));
    //     contact_dampings.push_back(school_closure(static_open_scenario_spring, 0.0, 0.0));
    //     contact_dampings.push_back(home_office(static_open_scenario_spring, 0.0, 0.0));
    //     contact_dampings.push_back(social_events(static_open_scenario_spring, 0.0, 0.0));
    //     contact_dampings.push_back(social_events_work(static_open_scenario_spring, 0.0, 0.0));
    //     contact_dampings.push_back(physical_distancing_home(static_open_scenario_spring, 0.0, 0.0));
    //     contact_dampings.push_back(physical_distancing_school(static_open_scenario_spring, 0.2 + narrow, 0.4 - narrow));
    //     contact_dampings.push_back(physical_distancing_work(static_open_scenario_spring, 0.2 + narrow, 0.4 - narrow));
    //     contact_dampings.push_back(physical_distancing_other(static_open_scenario_spring, 0.2 + narrow, 0.4 - narrow));
    //     contact_dampings.push_back(senior_awareness(static_open_scenario_spring, 0.0, 0.0));
    // }

    // //OPEN SCENARIO
    // int month_open;
    // if (late) {
    //     month_open = 8;
    // }
    // else {
    //     month_open = 7;
    // }
    // double masks_low, masks_high, masks_low_school, masks_high_school, masks_narrow;
    // if (masks) {
    //     masks_low_school  = 0.2;
    //     masks_high_school = 0.4;
    //     masks_low         = 0.2;
    //     masks_high        = 0.4;
    //     masks_narrow      = narrow;
    // }
    // else {

    //     masks_low_school  = 0.0;
    //     masks_high_school = 0.0;
    //     masks_low         = 0.0;
    //     masks_high        = 0.0;
    //     masks_narrow      = 0.0;
    // }
    // auto start_open = mio::Date(2022, month_open, 1);
    // if (start_open < end_date) {
    //     auto start_summer = mio::SimulationTime(mio::get_offset_in_days(start_open, start_date));
    //     contact_dampings.push_back(contacts_at_home(start_summer, 0.0, 0.0));
    //     contact_dampings.push_back(school_closure(start_summer, 0.0, 0.0));
    //     contact_dampings.push_back(home_office(start_summer, 0.0, 0.0));
    //     contact_dampings.push_back(social_events(start_summer, 0.0, 0.0));
    //     contact_dampings.push_back(social_events_work(start_summer, 0.0, 0.0));
    //     contact_dampings.push_back(physical_distancing_home(start_summer, 0.0, 0.0));
    //     contact_dampings.push_back(physical_distancing_school(start_summer, masks_low_school + masks_narrow,
    //                                                           masks_high_school - masks_narrow));
    //     contact_dampings.push_back(
    //         physical_distancing_work(start_summer, masks_low + masks_narrow, masks_high - masks_narrow));
    //     contact_dampings.push_back(
    //         physical_distancing_other(start_summer, masks_low + masks_narrow, masks_high - masks_narrow));
    //     contact_dampings.push_back(senior_awareness(start_summer, 0.0, 0.0));
    // }

    // auto start_autumn = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2022, 10, 1), start_date));
    // contact_dampings.push_back(contacts_at_home(start_autumn, 0.0, 0.0));
    // contact_dampings.push_back(school_closure(start_autumn, 0.3 + narrow, 0.5 - narrow));
    // // contact_dampings.push_back(home_office(start_autumn, 0.3 + narrow, 0.5 - narrow)); // S3F only
    // contact_dampings.push_back(social_events(start_autumn, 0.3 + narrow, 0.5 - narrow));
    // contact_dampings.push_back(social_events_work(start_autumn, 0.0, 0.0));

    // contact_dampings.push_back(home_office(start_autumn, 0.0 + narrow, 0.2 - narrow)); // S2F

    // //contact_dampings.push_back(school_closure(start_autumn, 0.0 + narrow, 0.2 - narrow)); // S1F
    // //contact_dampings.push_back(home_office(start_autumn, 0.0 + narrow, 0.2 - narrow)); // S1F
    // //contact_dampings.push_back(social_events(start_autumn,  0.0 + narrow, 0.2 - narrow)); // S1F

    // narrow = 0.0;
    // //local dynamic NPIs
    // auto& dynamic_npis        = params.get<mio::osecirvvs::DynamicNPIsInfected>();
    // auto dynamic_npi_dampings = std::vector<mio::DampingSampling>();

    // dynamic_npi_dampings.push_back(contacts_at_home(mio::SimulationTime(0), 0.1 + narrow, 0.3 - narrow));
    // dynamic_npi_dampings.push_back(school_closure(mio::SimulationTime(0), 0.2 + narrow,
    //                                               0.4 - narrow)); //0.25 - 0.25 in autumn
    // dynamic_npi_dampings.push_back(home_office(mio::SimulationTime(0), 0.1 + narrow, 0.3 - narrow));
    // dynamic_npi_dampings.push_back(social_events(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    // dynamic_npi_dampings.push_back(social_events_work(mio::SimulationTime(0), 0.0, 0.0));
    // dynamic_npi_dampings.push_back(physical_distancing_home(mio::SimulationTime(0), 0.0, 0.0));
    // dynamic_npi_dampings.push_back(physical_distancing_school(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    // dynamic_npi_dampings.push_back(physical_distancing_work(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    // dynamic_npi_dampings.push_back(physical_distancing_other(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    // dynamic_npi_dampings.push_back(senior_awareness(mio::SimulationTime(0), 0.0, 0.0));

    // auto dynamic_npi_dampings2 = std::vector<mio::DampingSampling>();
    // dynamic_npi_dampings2.push_back(contacts_at_home(mio::SimulationTime(0), 0.5 + narrow, 0.7 - narrow));
    // dynamic_npi_dampings2.push_back(school_closure(mio::SimulationTime(0), 0.4 + narrow,
    //                                                0.6 - narrow)); //0.25 - 0.25 in autumn
    // dynamic_npi_dampings2.push_back(home_office(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    // dynamic_npi_dampings2.push_back(social_events(mio::SimulationTime(0), 0.7 + narrow, 0.9 - narrow));
    // dynamic_npi_dampings2.push_back(social_events_work(mio::SimulationTime(0), 0.0, 0.0));
    // dynamic_npi_dampings2.push_back(physical_distancing_home(mio::SimulationTime(0), 0.0 + narrow, 0.2 - narrow));
    // dynamic_npi_dampings2.push_back(physical_distancing_school(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    // dynamic_npi_dampings2.push_back(physical_distancing_work(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    // dynamic_npi_dampings2.push_back(physical_distancing_other(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    // dynamic_npi_dampings2.push_back(senior_awareness(mio::SimulationTime(0), 0.0, 0.0));

    // dynamic_npis.set_interval(mio::SimulationTime(1.0));
    // dynamic_npis.set_duration(mio::SimulationTime(14.0));
    // dynamic_npis.set_base_value(100'000);
    // dynamic_npis.set_threshold(35.0, dynamic_npi_dampings);
    // dynamic_npis.set_threshold(100.0, dynamic_npi_dampings2);

    // //school holidays (holiday periods are set per node, see set_nodes)
    // auto school_holiday_value = mio::UncertainValue();
    // assign_uniform_distribution(school_holiday_value, 1.0, 1.0);
    // contacts.get_school_holiday_damping() =
    //     mio::DampingSampling(school_holiday_value, mio::DampingLevel(int(InterventionLevel::Holidays)),
    //                          mio::DampingType(int(Intervention::SchoolClosure)), mio::SimulationTime(0.0),
    //                          {size_t(ContactLocation::School)}, group_weights_all);

    return mio::success();
}

/**
 * Set synthetic population data for testing.
 * Same total populaton but different spread of infection in each county.
 * @param counties parameters for each county.
 */
void set_synthetic_population_data(std::vector<mio::osecirvvs::Model>& counties)
{
    for (size_t county_idx = 0; county_idx < counties.size(); ++county_idx) {
        double nb_total_t0 = 10000, nb_exp_t0 = 2, nb_inf_t0 = 0, nb_car_t0 = 0, nb_hosp_t0 = 0, nb_icu_t0 = 0,
               nb_rec_t0 = 0, nb_dead_t0 = 0;

        nb_exp_t0 = (double)((county_idx % 10 + 1) * 3);

        for (mio::AgeGroup i = 0; i < counties[county_idx].parameters.get_num_groups(); i++) {
            counties[county_idx].populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}]            = nb_exp_t0;
            counties[county_idx].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}] = nb_car_t0;
            counties[county_idx].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]   = nb_inf_t0;
            counties[county_idx].populations[{i, mio::osecirvvs::InfectionState::InfectedSevereNaive}]     = nb_hosp_t0;
            counties[county_idx].populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalNaive}]   = nb_icu_t0;
            counties[county_idx].populations[{i, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}] =
                nb_rec_t0;
            counties[county_idx].populations[{i, mio::osecirvvs::InfectionState::DeadNaive}] = nb_dead_t0;
            counties[county_idx].populations.set_difference_from_group_total<mio::AgeGroup>(
                {i, mio::osecirvvs::InfectionState::SusceptibleNaive}, nb_total_t0);
        }
    }
}

/**
 * Adds county nodes to graph.
 * Reads list counties and populations from files in the data directory. 
 * @param params Parameters that are shared between all nodes.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @param params_graph graph object that the nodes will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_nodes(const mio::osecirvvs::Parameters& params, mio::Date start_date, mio::Date end_date,
                              const fs::path& data_dir,
                              mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>& params_graph,
                              int vacc_campaign_szenario)
{
    namespace de = mio::regions::de;

    BOOST_OUTCOME_TRY(county_ids, mio::get_county_ids((data_dir / "pydata" / "Germany").string()));
    std::vector<mio::osecirvvs::Model> counties(county_ids.size(),
                                                mio::osecirvvs::Model((int)size_t(params.get_num_groups())));
    for (auto& county : counties) {
        county.parameters = params;
    }
    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;
    BOOST_OUTCOME_TRY(
        mio::osecirvvs::read_input_data_county(counties, start_date, county_ids, scaling_factor_infected,
                                               scaling_factor_icu, (data_dir / "pydata" / "Germany").string(),
                                               mio::get_offset_in_days(end_date, start_date), vacc_campaign_szenario));
    //set_synthetic_population_data(counties);

    for (size_t county_idx = 0; county_idx < counties.size(); ++county_idx) {

        //local parameters
        auto tnt_capacity = counties[county_idx].populations.get_total() * 1.43 / 100000.;
        assign_uniform_distribution(counties[county_idx].parameters.get<mio::osecirvvs::TestAndTraceCapacity>(),
                                    0.8 * tnt_capacity, 1.2 * tnt_capacity);

        //holiday periods (damping set globally, see set_npis)
        auto holiday_periods =
            de::get_holidays(de::get_state_id(de::CountyId(county_ids[county_idx])), start_date, end_date);
        auto& contacts = counties[county_idx].parameters.get<mio::osecirvvs::ContactPatterns>();
        contacts.get_school_holidays() =
            std::vector<std::pair<mio::SimulationTime, mio::SimulationTime>>(holiday_periods.size());
        std::transform(
            holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(), [=](auto& period) {
                return std::make_pair(mio::SimulationTime(mio::get_offset_in_days(period.first, start_date)),
                                      mio::SimulationTime(mio::get_offset_in_days(period.second, start_date)));
            });

        //uncertainty in populations
        for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = mio::Index<mio::osecirvvs::InfectionState>(0); j < mio::osecirvvs::InfectionState::Count;
                 ++j) {
                auto& compartment_value = counties[county_idx].populations[{i, j}];
                assign_uniform_distribution(compartment_value, 0.9 * double(compartment_value),
                                            1.1 * double(compartment_value));
            }
        }

        params_graph.add_node(county_ids[county_idx], counties[county_idx]);
    }
    return mio::success();
}

/**
 * Adds edges to graph.
 * Edges represent commuting and other migration between counties.
 * Reads migration from files in the data directory.
 * @param data_dir data directory.
 * @param params_graph graph object that the nodes will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_edges(const fs::path& data_dir,
                              mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>& params_graph)
{
    //migration between nodes
    BOOST_OUTCOME_TRY(
        mobility_data_commuter,
        mio::read_mobility_plain((data_dir / "migration" / "commuter_migration_scaled_2020.txt").string()));
    BOOST_OUTCOME_TRY(mobility_data_twitter,
                      mio::read_mobility_plain((data_dir / "migration" / "twitter_scaled_1252.txt").string()));
    if (mobility_data_commuter.rows() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_commuter.cols() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_twitter.rows() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_twitter.cols() != Eigen::Index(params_graph.nodes().size())) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Mobility matrices do not have the correct size. You may need to run "
                            "transformMobilitydata.py from pycode memilio epidata package.");
    }

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
    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            auto& populations = params_graph.nodes()[county_idx_i].property.populations;
            //mobility coefficients have the same number of components as the contact matrices.
            //so that the same NPIs/dampings can be used for both (e.g. more home office => fewer commuters)
            auto mobility_coeffs = mio::MigrationCoefficientGroup(contact_locations.size(), populations.numel());

            //commuters
            auto working_population = 0.0;
            auto min_commuter_age   = mio::AgeGroup(2);
            auto max_commuter_age   = mio::AgeGroup(4); //this group is partially retired, only partially commutes
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                working_population += populations.get_group_total(age) * (age == max_commuter_age ? 0.33 : 1.0);
            }
            auto commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) /
                                     working_population; //data is absolute numbers, we need relative
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_index = populations.get_flat_index({age, compartment});
                    mobility_coeffs[size_t(ContactLocation::Work)].get_baseline()[coeff_index] =
                        commuter_coeff_ij * (age == max_commuter_age ? 0.33 : 1.0);
                }
            }
            //others
            auto total_population = populations.get_total();
            auto twitter_coeff    = mobility_data_twitter(county_idx_i, county_idx_j) /
                                 total_population; //data is absolute numbers, we need relative
            for (auto age = mio::AgeGroup(0); age < populations.size<mio::AgeGroup>(); ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_idx = populations.get_flat_index({age, compartment});
                    mobility_coeffs[size_t(ContactLocation::Other)].get_baseline()[coeff_idx] = twitter_coeff;
                }
            }

            //only add edges with migration above thresholds for performance
            //thresholds are chosen empirically so that more than 99% of migration is covered, approx. 1/3 of the edges
            if (commuter_coeff_ij > 4e-5 || twitter_coeff > 1e-5) {
                params_graph.add_edge(county_idx_i, county_idx_j, std::move(mobility_coeffs));
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
mio::IOResult<mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>>
create_graph(mio::Date start_date, mio::Date end_date, const fs::path& data_dir, int variant_szenario,
             int vacc_campaign_szenario, int vacc_effectiveness_szenario)
{
    const auto summer_date = mio::Date(2022, 7, 1);

    //global parameters
    const int num_age_groups = 6;
    mio::osecirvvs::Parameters params(num_age_groups);
    params.get<mio::osecirvvs::StartDay>()           = mio::get_day_in_year(start_date);
    params.get<mio::osecirvvs::SzenarioNewVariant>() = variant_szenario;
    params.get_end_dynamic_npis()                    = mio::get_offset_in_days(start_date, summer_date);
    BOOST_OUTCOME_TRY(set_covid_parameters(params, vacc_effectiveness_szenario));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));
    BOOST_OUTCOME_TRY(set_npis(start_date, end_date, params));

    //graph of counties with populations and local parameters
    //and mobility between counties
    mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> params_graph;
    BOOST_OUTCOME_TRY(set_nodes(params, start_date, end_date, data_dir, params_graph, vacc_campaign_szenario));
    BOOST_OUTCOME_TRY(set_edges(data_dir, params_graph));

    return mio::success(params_graph);
}

/**
 * Load the input graph for the parameter study that was previously saved.
 * @param save_dir directory where the graph was saved.
 * @returns created graph or any io errors that happen during reading of the files.
 */
mio::IOResult<mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>> load_graph(const fs::path& save_dir)
{
    return mio::read_graph<mio::osecirvvs::Model>(save_dir.string());
}

/**
 * Save the input graph for the parameter study.
 * @param save_dir directory where the graph will be saved.
 * @returns any io errors that happen during writing of the files.
 */
mio::IOResult<void> save_graph(const mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>& params_graph,
                               const fs::path& save_dir)
{
    return mio::write_graph(params_graph, save_dir.string());
}

/**
 * Create an unconnected graph.
 * Can be used to save space on disk when writing parameters if the edges are not required.
 * @param params parameters for each county node.
 * @param county_ids id of each county node.
 * @return graph with county nodes but no edges.
 */
auto make_graph_no_edges(const std::vector<mio::osecirvvs::Model>& params, const std::vector<int>& county_ids)
{
    //make a graph without edges for writing to file
    auto graph = mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>();
    for (auto i = size_t(0); i < county_ids.size(); ++i) {
        graph.add_node(county_ids[i], params[i]);
    }
    return graph;
}

/**
 * Save the result of a single parameter study run.
 * Creates a new subdirectory for this run.
 * @param result result of the simulation run.
 * @param params parameters used for the simulation run.
 * @param county_ids ids of the county nodes.
 * @param result_dir top level directory for all results of the parameter study.
 * @param run_idx index of the run.
 * @return any io errors that occur during writing of the files.
 */
mio::IOResult<void> save_result(const std::vector<mio::TimeSeries<double>>& result,
                                const std::vector<mio::osecirvvs::Model>& params, const std::vector<int>& county_ids,
                                const fs::path& result_dir, size_t run_idx)
{
    auto result_dir_run = result_dir / ("run" + std::to_string(run_idx));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir_run.string()));
    BOOST_OUTCOME_TRY(mio::save_result(result, county_ids, (int)(size_t)params[0].parameters.get_num_groups(),
                                       (result_dir_run / "Result.h5").string()));
    BOOST_OUTCOME_TRY(
        mio::write_graph(make_graph_no_edges(params, county_ids), result_dir_run.string(), mio::IOF_OmitDistributions));
    return mio::success();
}

/**
 * Save the results of a parameter study.
 * Stores different percentiles and sums of the results and parameters. 
 * @param ensemble_results result of each simulation run.
 * @param ensemble_params parameters used for each simulation run.
 * @param county_ids ids of the county nodes.
 * @param result_dir top level directory for all results of the parameter study.
 * @return any io errors that occur during writing of the files.
 */
mio::IOResult<void> save_results(const std::vector<std::vector<mio::TimeSeries<double>>>& ensemble_results,
                                 const std::vector<std::vector<mio::osecirvvs::Model>>& ensemble_params,
                                 const std::vector<int>& county_ids, const fs::path& result_dir)
{
    //save results and sum of results over nodes
    auto ensemble_result_sum = mio::sum_nodes(ensemble_results);
    auto num_groups          = (int)(size_t)ensemble_params[0][0].parameters.get_num_groups();
    for (size_t i = 0; i < ensemble_result_sum.size(); ++i) {
        BOOST_OUTCOME_TRY((mio::save_result(ensemble_result_sum[i], {0}, num_groups,
                                            (result_dir / ("results_run" + std::to_string(i) + "_sum.h5")).string())));
        BOOST_OUTCOME_TRY((mio::save_result(ensemble_results[i], county_ids, num_groups,
                                            (result_dir / ("results_run" + std::to_string(i) + ".h5")).string())));
    }

    //make directories for percentiles
    auto result_dir_p05 = result_dir / "p05";
    auto result_dir_p25 = result_dir / "p25";
    auto result_dir_p50 = result_dir / "p50";
    auto result_dir_p75 = result_dir / "p75";
    auto result_dir_p95 = result_dir / "p95";
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir_p05.string()));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir_p25.string()));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir_p50.string()));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir_p75.string()));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir_p95.string()));

    //save percentiles of results, summed over nodes
    {
        auto ensemble_results_sum_p05 = mio::ensemble_percentile(ensemble_result_sum, 0.05);
        auto ensemble_results_sum_p25 = mio::ensemble_percentile(ensemble_result_sum, 0.25);
        auto ensemble_results_sum_p50 = mio::ensemble_percentile(ensemble_result_sum, 0.50);
        auto ensemble_results_sum_p75 = mio::ensemble_percentile(ensemble_result_sum, 0.75);
        auto ensemble_results_sum_p95 = mio::ensemble_percentile(ensemble_result_sum, 0.95);

        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_sum_p05, {0}, num_groups, (result_dir_p05 / "Results_sum.h5").string()));
        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_sum_p25, {0}, num_groups, (result_dir_p25 / "Results_sum.h5").string()));
        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_sum_p50, {0}, num_groups, (result_dir_p50 / "Results_sum.h5").string()));
        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_sum_p75, {0}, num_groups, (result_dir_p75 / "Results_sum.h5").string()));
        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_sum_p95, {0}, num_groups, (result_dir_p95 / "Results_sum.h5").string()));
    }

    //save percentiles of results
    {
        auto ensemble_results_p05 = mio::ensemble_percentile(ensemble_results, 0.05);
        auto ensemble_results_p25 = mio::ensemble_percentile(ensemble_results, 0.25);
        auto ensemble_results_p50 = mio::ensemble_percentile(ensemble_results, 0.50);
        auto ensemble_results_p75 = mio::ensemble_percentile(ensemble_results, 0.75);
        auto ensemble_results_p95 = mio::ensemble_percentile(ensemble_results, 0.95);

        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_p05, county_ids, num_groups, (result_dir_p05 / "Results.h5").string()));
        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_p25, county_ids, num_groups, (result_dir_p25 / "Results.h5").string()));
        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_p50, county_ids, num_groups, (result_dir_p50 / "Results.h5").string()));
        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_p75, county_ids, num_groups, (result_dir_p75 / "Results.h5").string()));
        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble_results_p95, county_ids, num_groups, (result_dir_p95 / "Results.h5").string()));
    }

    //save percentiles of parameters
    {
        auto ensemble_params_p05 = ensemble_params_percentile(ensemble_params, 0.05);
        auto ensemble_params_p25 = ensemble_params_percentile(ensemble_params, 0.25);
        auto ensemble_params_p50 = ensemble_params_percentile(ensemble_params, 0.50);
        auto ensemble_params_p75 = ensemble_params_percentile(ensemble_params, 0.75);
        auto ensemble_params_p95 = ensemble_params_percentile(ensemble_params, 0.95);

        auto make_graph = [&county_ids](auto&& params) {
            return make_graph_no_edges(params, county_ids);
        };
        BOOST_OUTCOME_TRY(
            mio::write_graph(make_graph(ensemble_params_p05), result_dir_p05.string(), mio::IOF_OmitDistributions));
        BOOST_OUTCOME_TRY(
            mio::write_graph(make_graph(ensemble_params_p25), result_dir_p25.string(), mio::IOF_OmitDistributions));
        BOOST_OUTCOME_TRY(
            mio::write_graph(make_graph(ensemble_params_p50), result_dir_p50.string(), mio::IOF_OmitDistributions));
        BOOST_OUTCOME_TRY(
            mio::write_graph(make_graph(ensemble_params_p75), result_dir_p75.string(), mio::IOF_OmitDistributions));
        BOOST_OUTCOME_TRY(
            mio::write_graph(make_graph(ensemble_params_p95), result_dir_p95.string(), mio::IOF_OmitDistributions));
    }
    return mio::success();
}

/**
 * Different modes for running the parameter study.
 */
enum class RunMode
{
    Load,
    Save,
};

/**
 * Run the parameter study.
 * Load a previously stored graph or create a new one from data. 
 * The graph is the input for the parameter study.
 * A newly created graph is saved and can be reused.
 * @param mode Mode for running the parameter study.
 * @param data_dir data directory. Not used if mode is RunMode::Load.
 * @param save_dir directory where the graph is loaded from if mode is RunMOde::Load or save to if mode is RunMode::Save.
 * @param result_dir directory where all results of the parameter study will be stored.
 * @returns any io error that occurs during reading or writing of files.
 */
mio::IOResult<void> run(RunMode mode, const fs::path& data_dir, const fs::path& save_dir, const fs::path& result_dir,
                        int variant_szenario, int vacc_campaign_szenario, int vacc_effectiveness_szenario)
{
    mio::Date temp_date;
    temp_date = mio::Date(2022, 6, 1);

    const auto start_date   = temp_date;
    const auto num_days_sim = 10.0;
    const auto end_date     = mio::offset_date_by_days(start_date, int(std::ceil(num_days_sim)));
    const auto num_runs     = 1;

    //create or load graph
    mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> params_graph;
    if (mode == RunMode::Save) {
        BOOST_OUTCOME_TRY(created, create_graph(start_date, end_date, data_dir, variant_szenario,
                                                vacc_campaign_szenario, vacc_effectiveness_szenario));
        BOOST_OUTCOME_TRY(save_graph(created, save_dir));
        params_graph = created;
    }
    else {
        BOOST_OUTCOME_TRY(loaded, load_graph(save_dir));
        params_graph = loaded;
    }

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    //run parameter study
    auto parameter_study =
        mio::ParameterStudy<mio::osecirvvs::Simulation<>>{params_graph, 0.0, num_days_sim, 0.5, num_runs};
    auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
    ensemble_results.reserve(size_t(num_runs));
    auto ensemble_params = std::vector<std::vector<mio::osecirvvs::Model>>{};
    ensemble_params.reserve(size_t(num_runs));
    auto save_result_result = mio::IOResult<void>(mio::success());
    auto run_idx            = size_t(0);
    parameter_study.run(
        [&](auto&& graph) {
            return draw_sample(graph);
        },
        [&](auto results_graph) {
            ensemble_results.push_back(mio::interpolate_simulation_result(results_graph));

            ensemble_params.emplace_back();
            ensemble_params.back().reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                           std::back_inserter(ensemble_params.back()), [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            if (save_result_result) {
                save_result_result =
                    save_result(ensemble_results.back(), ensemble_params.back(), county_ids, result_dir, run_idx);
            }
            std::cout << "run " << run_idx << " complete." << std::endl;
            ++run_idx;
        });
    BOOST_OUTCOME_TRY(save_result_result);
    BOOST_OUTCOME_TRY(save_results(ensemble_results, ensemble_params, county_ids, result_dir));

    return mio::success();
}

int main(int argc, char** argv)
{
    //TODO: proper command line interface to set:
    //- number of runs
    //- start and end date (may be incompatible with runmode::load)
    //- seeds
    //- log level
    //- ...

    mio::set_log_level(mio::LogLevel::warn);

    RunMode mode;
    std::string save_dir;
    std::string data_dir;
    std::string result_dir;
    int variant_szenario            = 1;
    int vacc_campaign_szenario      = 1;
    int vacc_effectiveness_szenario = 1;
    if (argc == 6) {
        mode                        = RunMode::Save;
        data_dir                    = argv[1];
        save_dir                    = argv[2];
        result_dir                  = argv[3];
        variant_szenario            = atoi(argv[4]);
        vacc_campaign_szenario      = atoi(argv[5]);
        vacc_effectiveness_szenario = atoi(argv[6]);

        printf("Szenario: %d, Vaccination campaign szenario: %d, Vaccination effectiveness szenario: %d\n",
               variant_szenario, vacc_campaign_szenario, vacc_effectiveness_szenario);

        printf("Reading data from \"%s\", saving graph to \"%s\".\n", data_dir.c_str(), save_dir.c_str());
    }
    else if (argc == 3) {
        mode       = RunMode::Load;
        save_dir   = argv[1];
        result_dir = argv[2];
        data_dir   = "";
        printf("Loading graph from \"%s\".\n", save_dir.c_str());
    }
    else if (argc == 4) {
        mode       = RunMode::Save;
        data_dir   = argv[1];
        save_dir   = argv[2];
        result_dir = argv[3];
        printf("Reading data from \"%s\", saving graph to \"%s\".\n", data_dir.c_str(), save_dir.c_str());
    }
    else {
        printf("Usage:\n");
        printf("paper1 <data_dir> <save_dir> <result_dir>\n");
        printf("\tMake graph with data from <data_dir> and save at <save_dir>, then run the simulation.\n");
        printf("\tStore the results in <result_dir>\n");
        printf("paper1 <load_dir> <result_dir>\n");
        printf("\tLoad graph from <load_dir>, then run the simulation.\n");
        return 0;
    }

    std::string szenario_info;
    if (variant_szenario == 0) {
        szenario_info = "Omikron (BA5) or a very similar variant stays (none Immune escape variant, disease severity "
                        "remains the same as with omicron).";
    }
    else if (variant_szenario == 1) {
        szenario_info = "New variant with Immune escape but same disease severity";
    }
    else if (variant_szenario == 2) {
        szenario_info = "New variant with Immune escape and disease severity as delta";
    }
    else {
        printf("Input for szenario not definied, has to be between 1-3!");
    }

    printf("Szenario : %s\n", szenario_info.c_str());

    result_dir += "_szenario_" + std::to_string(variant_szenario) + "_vacc_camp_" +
                  std::to_string(vacc_campaign_szenario) + "_vacc_eff_" + std::to_string(vacc_effectiveness_szenario);

    boost::filesystem::path dir(result_dir);
    bool created = boost::filesystem::create_directories(dir);

    if (created) {
        mio::log_info("Directory '{:s}' was created.", dir.string());
    }
    printf("Saving results to \"%s\".\n", result_dir.c_str());

    //mio::thread_local_rng().seed(
    //    {114381446, 2427727386, 806223567, 832414962, 4121923627, 1581162203}); //set seeds, e.g., for debugging
    printf("Seeds: ");
    for (auto s : mio::thread_local_rng().get_seeds()) {
        printf("%u, ", s);
    }
    printf("\n");

    auto result = run(mode, data_dir, save_dir, result_dir, variant_szenario, vacc_campaign_szenario,
                      vacc_effectiveness_szenario);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    return 0;
}
