/*
* Copyright (C) 2020-2024 MEmilio
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
#include "memilio/utils/stl_util.h"
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
mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters<double>& params, bool long_time)
{
    //times
    // TimeExposed and TimeInfectedNoSymptoms are calculated as described in
    // Khailaie et al. (https://doi.org/10.1186/s12916-020-01884-4)
    // given SI_min = 3.935, SI_max = 4.6, INC = 5.2
    const double timeExposedMin            = 2.67;
    const double timeExposedMax            = 4.;
    const double timeInfectedNoSymptomsMin = 1.2;
    const double timeInfectedNoSymptomsMax = 2.53;

    const double timeInfectedSymptomsMin[] = {5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465};
    const double timeInfectedSymptomsMax[] = {8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085};
    const double timeInfectedSevereMin[]   = {3.925, 3.925, 4.85, 6.4, 7.2, 9.};
    const double timeInfectedSevereMax[]   = {6.075, 6.075, 7., 8.7, 9.8, 13.};
    const double timeInfectedCriticalMin[] = {4.95, 4.95, 4.86, 14.14, 14.4, 10.};
    const double timeInfectedCriticalMax[] = {8.95, 8.95, 8.86, 20.58, 19.8, 13.2};

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
    double fac_variant                                 = 1.4;
    const double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                          0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

    const double transmissionProbabilityOnContactMax[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                                          0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};
    const double relativeTransmissionNoSymptomsMin     = 0.5;
    const double relativeTransmissionNoSymptomsMax     = 0.5;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    const double riskOfInfectionFromSymptomaticMin    = 0.0;
    const double riskOfInfectionFromSymptomaticMax    = 0.2;
    const double maxRiskOfInfectionFromSymptomaticMin = 0.4;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.5;
    const double recoveredPerInfectedNoSymptomsMin[]  = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15};
    const double recoveredPerInfectedNoSymptomsMax[]  = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25};
    const double severePerInfectedSymptomsMin[]       = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20};
    const double severePerInfectedSymptomsMax[]       = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25};
    const double criticalPerSevereMin[]               = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35};
    const double criticalPerSevereMax[]               = {0.10, 0.10, 0.10, 0.20, 0.35, 0.45};
    const double deathsPerCriticalMin[]               = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5};
    const double deathsPerCriticalMax[]               = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7};

    const double reducExposedPartialImmunityMin                     = 0.75;
    const double reducExposedPartialImmunityMax                     = 0.85;
    const double reducExposedImprovedImmunityMin                    = 0.281;
    const double reducExposedImprovedImmunityMax                    = 0.381;
    const double reducInfectedSymptomsPartialImmunityMin            = 0.6;
    const double reducInfectedSymptomsPartialImmunityMax            = 0.7;
    const double reducInfectedSymptomsImprovedImmunityMin           = 0.193;
    const double reducInfectedSymptomsImprovedImmunityMax           = 0.293;
    const double reducInfectedSevereCriticalDeadPartialImmunityMin  = 0.05;
    const double reducInfectedSevereCriticalDeadPartialImmunityMax  = 0.15;
    const double reducInfectedSevereCriticalDeadImprovedImmunityMin = 0.041;
    const double reducInfectedSevereCriticalDeadImprovedImmunityMax = 0.141;

    double temp_reducTimeInfectedMild;
    if (long_time) {
        temp_reducTimeInfectedMild = 1.0;
    }
    else {
        temp_reducTimeInfectedMild = 0.5;
    }
    const double reducTimeInfectedMild = temp_reducTimeInfectedMild;

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
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirvvs::Parameters<double>& params)
{
    //TODO: io error handling
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);
    }
    params.get<mio::osecirvvs::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

/**
 * Set NPIs.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param params Object that the NPIs will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_npis(mio::Date start_date, mio::Date end_date, mio::osecirvvs::Parameters<double>& params,
                             bool late, bool masks, bool test)
{
    auto& contacts         = params.get<mio::osecirvvs::ContactPatterns<double>>();
    auto& contact_dampings = contacts.get_dampings();

    if (test) {
        params.get_commuter_nondetection() = 0.85;
    }
    else {
        params.get_commuter_nondetection() = 1.0;
    }

    //weights for age groups affected by an NPI
    auto group_weights_all     = Eigen::VectorXd::Constant(size_t(params.get_num_groups()), 1.0);
    auto group_weights_seniors = Eigen::VectorXd::NullaryExpr(size_t(params.get_num_groups()), [](auto&& i) {
        return i == 5 ? 1.0 : i == 4 ? 0.5 : 0.0; //65-80 only partially
    });

    //helper functions that create dampings for specific NPIs
    auto contacts_at_home = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::Main)),
                                            mio::DampingType(int(Intervention::Home)), t,
                                            {size_t(ContactLocation::Home)}, group_weights_all);
    };
    auto school_closure = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::Main)),
                                            mio::DampingType(int(Intervention::SchoolClosure)), t,
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto home_office = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::Main)),
                                            mio::DampingType(int(Intervention::HomeOffice)), t,
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto social_events = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::Main)),
                                            mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), t,
                                            {size_t(ContactLocation::Other)}, group_weights_all);
    };
    auto social_events_work = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::Main)),
                                            mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), t,
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto physical_distancing_home = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                            {size_t(ContactLocation::Home)}, group_weights_all);
    };
    auto physical_distancing_school = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto physical_distancing_work = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto physical_distancing_other = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                            {size_t(ContactLocation::Other)}, group_weights_all);
    };
    auto senior_awareness = [=](auto t, auto min, auto max) {
        auto v = mio::UncertainValue<double>();
        assign_uniform_distribution(v, min, max);
        return mio::DampingSampling<double>(v, mio::DampingLevel(int(InterventionLevel::SeniorAwareness)),
                                            mio::DampingType(int(Intervention::SeniorAwareness)), t,
                                            {size_t(ContactLocation::Home), size_t(ContactLocation::Other)},
                                            group_weights_seniors);
    };

    //OPEN SCENARIO SPRING
    auto start_year = mio::Date(2021, 1, 1);
    double narrow   = 0.05;
    if (start_year < end_date) {
        auto static_open_scenario_spring = mio::SimulationTime(mio::get_offset_in_days(start_year, start_date));
        contact_dampings.push_back(contacts_at_home(static_open_scenario_spring, 0.0, 0.0));
        contact_dampings.push_back(school_closure(static_open_scenario_spring, 0.0, 0.0));
        contact_dampings.push_back(home_office(static_open_scenario_spring, 0.0, 0.0));
        contact_dampings.push_back(social_events(static_open_scenario_spring, 0.0, 0.0));
        contact_dampings.push_back(social_events_work(static_open_scenario_spring, 0.0, 0.0));
        contact_dampings.push_back(physical_distancing_home(static_open_scenario_spring, 0.0, 0.0));
        contact_dampings.push_back(physical_distancing_school(static_open_scenario_spring, 0.2 + narrow, 0.4 - narrow));
        contact_dampings.push_back(physical_distancing_work(static_open_scenario_spring, 0.2 + narrow, 0.4 - narrow));
        contact_dampings.push_back(physical_distancing_other(static_open_scenario_spring, 0.2 + narrow, 0.4 - narrow));
        contact_dampings.push_back(senior_awareness(static_open_scenario_spring, 0.0, 0.0));
    }

    //OPEN SCENARIO
    int month_open;
    if (late) {
        month_open = 8;
    }
    else {
        month_open = 7;
    }
    double masks_low, masks_high, masks_low_school, masks_high_school, masks_narrow;
    if (masks) {
        masks_low_school  = 0.2;
        masks_high_school = 0.4;
        masks_low         = 0.2;
        masks_high        = 0.4;
        masks_narrow      = narrow;
    }
    else {

        masks_low_school  = 0.0;
        masks_high_school = 0.0;
        masks_low         = 0.0;
        masks_high        = 0.0;
        masks_narrow      = 0.0;
    }
    auto start_open = mio::Date(2021, month_open, 1);
    if (start_open < end_date) {
        auto start_summer = mio::SimulationTime(mio::get_offset_in_days(start_open, start_date));
        contact_dampings.push_back(contacts_at_home(start_summer, 0.0, 0.0));
        contact_dampings.push_back(school_closure(start_summer, 0.0, 0.0));
        contact_dampings.push_back(home_office(start_summer, 0.0, 0.0));
        contact_dampings.push_back(social_events(start_summer, 0.0, 0.0));
        contact_dampings.push_back(social_events_work(start_summer, 0.0, 0.0));
        contact_dampings.push_back(physical_distancing_home(start_summer, 0.0, 0.0));
        contact_dampings.push_back(physical_distancing_school(start_summer, masks_low_school + masks_narrow,
                                                              masks_high_school - masks_narrow));
        contact_dampings.push_back(
            physical_distancing_work(start_summer, masks_low + masks_narrow, masks_high - masks_narrow));
        contact_dampings.push_back(
            physical_distancing_other(start_summer, masks_low + masks_narrow, masks_high - masks_narrow));
        contact_dampings.push_back(senior_awareness(start_summer, 0.0, 0.0));
    }

    auto start_autumn = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2021, 10, 1), start_date));
    contact_dampings.push_back(contacts_at_home(start_autumn, 0.0, 0.0));
    contact_dampings.push_back(school_closure(start_autumn, 0.3 + narrow, 0.5 - narrow));
    // contact_dampings.push_back(home_office(start_autumn, 0.3 + narrow, 0.5 - narrow)); // S3F only
    contact_dampings.push_back(social_events(start_autumn, 0.3 + narrow, 0.5 - narrow));
    contact_dampings.push_back(social_events_work(start_autumn, 0.0, 0.0));

    contact_dampings.push_back(home_office(start_autumn, 0.0 + narrow, 0.2 - narrow)); // S2F

    //contact_dampings.push_back(school_closure(start_autumn, 0.0 + narrow, 0.2 - narrow)); // S1F
    //contact_dampings.push_back(home_office(start_autumn, 0.0 + narrow, 0.2 - narrow)); // S1F
    //contact_dampings.push_back(social_events(start_autumn,  0.0 + narrow, 0.2 - narrow)); // S1F

    narrow = 0.0;
    //local dynamic NPIs
    params.get<mio::osecirvvs::DynamicNPIsImplementationDelay<double>>() = 7;
    auto& dynamic_npis        = params.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms<double>>();
    auto dynamic_npi_dampings = std::vector<mio::DampingSampling<double>>();

    dynamic_npi_dampings.push_back(contacts_at_home(mio::SimulationTime(0), 0.1 + narrow, 0.3 - narrow));
    dynamic_npi_dampings.push_back(school_closure(mio::SimulationTime(0), 0.2 + narrow,
                                                  0.4 - narrow)); //0.25 - 0.25 in autumn
    dynamic_npi_dampings.push_back(home_office(mio::SimulationTime(0), 0.1 + narrow, 0.3 - narrow));
    dynamic_npi_dampings.push_back(social_events(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    dynamic_npi_dampings.push_back(social_events_work(mio::SimulationTime(0), 0.0, 0.0));
    dynamic_npi_dampings.push_back(physical_distancing_home(mio::SimulationTime(0), 0.0, 0.0));
    dynamic_npi_dampings.push_back(physical_distancing_school(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    dynamic_npi_dampings.push_back(physical_distancing_work(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    dynamic_npi_dampings.push_back(physical_distancing_other(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    dynamic_npi_dampings.push_back(senior_awareness(mio::SimulationTime(0), 0.0, 0.0));

    auto dynamic_npi_dampings2 = std::vector<mio::DampingSampling<double>>();
    dynamic_npi_dampings2.push_back(contacts_at_home(mio::SimulationTime(0), 0.5 + narrow, 0.7 - narrow));
    dynamic_npi_dampings2.push_back(school_closure(mio::SimulationTime(0), 0.4 + narrow,
                                                   0.6 - narrow)); //0.25 - 0.25 in autumn
    dynamic_npi_dampings2.push_back(home_office(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    dynamic_npi_dampings2.push_back(social_events(mio::SimulationTime(0), 0.7 + narrow, 0.9 - narrow));
    dynamic_npi_dampings2.push_back(social_events_work(mio::SimulationTime(0), 0.0, 0.0));
    dynamic_npi_dampings2.push_back(physical_distancing_home(mio::SimulationTime(0), 0.0 + narrow, 0.2 - narrow));
    dynamic_npi_dampings2.push_back(physical_distancing_school(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    dynamic_npi_dampings2.push_back(physical_distancing_work(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    dynamic_npi_dampings2.push_back(physical_distancing_other(mio::SimulationTime(0), 0.2 + narrow, 0.4 - narrow));
    dynamic_npi_dampings2.push_back(senior_awareness(mio::SimulationTime(0), 0.0, 0.0));

    dynamic_npis.set_interval(mio::SimulationTime(1.0));
    dynamic_npis.set_duration(mio::SimulationTime(14.0));
    dynamic_npis.set_base_value(100'000);
    dynamic_npis.set_threshold(35.0, dynamic_npi_dampings);
    dynamic_npis.set_threshold(100.0, dynamic_npi_dampings2);

    //school holidays (holiday periods are set per node, see set_nodes)
    auto school_holiday_value = mio::UncertainValue<double>();
    assign_uniform_distribution(school_holiday_value, 1.0, 1.0);
    contacts.get_school_holiday_damping() =
        mio::DampingSampling<double>(school_holiday_value, mio::DampingLevel(int(InterventionLevel::Holidays)),
                                     mio::DampingType(int(Intervention::SchoolClosure)), mio::SimulationTime(0.0),
                                     {size_t(ContactLocation::School)}, group_weights_all);

    return mio::success();
}

/**
 * Set synthetic population data for testing.
 * Only sets immune-naive part of the population. The remaining part is zero.
 * Same total populaton but different spread of infection in each county.
 * @param counties parameters for each county.
 */
void set_synthetic_population_data(std::vector<mio::osecirvvs::Model<double>>& counties)
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
 * Create the input graph for the parameter study.
 * Reads files from the data directory.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @returns created graph or any io errors that happen during reading of the files.
 */
mio::IOResult<mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>>>
get_graph(mio::Date start_date, mio::Date end_date, const fs::path& data_dir, bool late, bool masks, bool test,
          bool long_time)
{
    const auto summer_date = late ? mio::Date(2021, 8, 1) : mio::Date(2021, 7, 1);

    // global parameters
    const int num_age_groups = 6;
    mio::osecirvvs::Parameters<double> params(num_age_groups);
    params.get<mio::osecirvvs::StartDay>() = mio::get_day_in_year(start_date);
    params.get_end_dynamic_npis()          = mio::get_offset_in_days(start_date, summer_date);
    BOOST_OUTCOME_TRY(set_covid_parameters(params, long_time));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));
    BOOST_OUTCOME_TRY(set_npis(start_date, end_date, params, late, masks, test));

    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;
    auto tnt_capacity_factor     = 1.43 / 100000.;
    auto mobile_compartments     = {mio::osecirvvs::InfectionState::SusceptibleNaive,
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

    // graph of counties with populations and local parameters
    // and mobility between counties
    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph;
    const auto& read_function_nodes = mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model<double>>;
    const auto& read_function_edges = mio::read_mobility_plain;
    const auto& node_id_function    = mio::get_node_ids;

    const auto& set_node_function =
        mio::set_nodes<mio::osecirvvs::TestAndTraceCapacity<double>, mio::osecirvvs::ContactPatterns<double>,
                       mio::osecirvvs::Model<double>, mio::MobilityParameters<double>,
                       mio::osecirvvs::Parameters<double>, decltype(read_function_nodes), decltype(node_id_function)>;
    const auto& set_edge_function =
        mio::set_edges<ContactLocation, mio::osecirvvs::Model<double>, mio::MobilityParameters<double>,
                       mio::MobilityCoefficientGroup, mio::osecirvvs::InfectionState, decltype(read_function_edges)>;
    BOOST_OUTCOME_TRY(set_node_function(
        params, start_date, end_date, data_dir,
        mio::path_join((data_dir / "pydata" / "Germany").string(), "county_current_population.json"), true,
        params_graph, read_function_nodes, node_id_function, scaling_factor_infected, scaling_factor_icu,
        tnt_capacity_factor, mio::get_offset_in_days(end_date, start_date), false, true));
    BOOST_OUTCOME_TRY(set_edge_function(data_dir, params_graph, mobile_compartments, contact_locations.size(),
                                        read_function_edges, std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.}));

    return mio::success(params_graph);
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
 * @param save_single_runs [Default: true] Defines if single run results are written to the disk. 
 * @returns any io error that occurs during reading or writing of files.
 */
mio::IOResult<void> run(RunMode mode, const fs::path& data_dir, const fs::path& save_dir, const fs::path& result_dir,
                        bool save_single_runs, bool late, bool masks, bool test, bool high, bool long_time, bool future)
{
    mio::Date temp_date;
    if (future) {
        temp_date = mio::Date(2021, 10, 15);
    }
    else {
        temp_date = mio::Date(2021, 6, 6);
    }
    const auto start_date   = temp_date;
    const auto num_days_sim = 90.0;
    const auto end_date     = mio::offset_date_by_days(start_date, int(std::ceil(num_days_sim)));
    const auto num_runs     = 500;

    //create or load graph
    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph;
    if (mode == RunMode::Save) {
        BOOST_OUTCOME_TRY(auto&& created, get_graph(start_date, end_date, data_dir, late, masks, test, long_time));
        BOOST_OUTCOME_TRY(write_graph(created, save_dir.string()));
        params_graph = created;
    }
    else {
        BOOST_OUTCOME_TRY(auto&& loaded, mio::read_graph<double, mio::osecirvvs::Model<double>>(save_dir.string()));
        params_graph = loaded;
    }

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    //run parameter study
    auto parameter_study =
        mio::ParameterStudy<mio::osecirvvs::Simulation<>>{params_graph, 0.0, num_days_sim, 0.5, num_runs};

    // parameter_study.get_rng().seed(
    //    {114381446, 2427727386, 806223567, 832414962, 4121923627, 1581162203}); //set seeds, e.g., for debugging
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
            return draw_sample(graph, high);
        },
        [&](auto results_graph, auto&& run_idx) {
            auto interpolated_result = mio::interpolate_simulation_result(results_graph);
            auto params              = std::vector<mio::osecirvvs::Model<double>>();
            params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
                           [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            if (save_single_run_result && save_single_runs) {
                save_single_run_result =
                    save_result_with_params(interpolated_result, params, county_ids, result_dir, run_idx);
            }
            std::cout << "run " << run_idx << " complete." << std::endl;
            return std::make_pair(interpolated_result, params);
        });

    if (ensemble.size() > 0) {
        auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_results.reserve(ensemble.size());
        auto ensemble_params = std::vector<std::vector<mio::osecirvvs::Model<double>>>{};
        ensemble_params.reserve(ensemble.size());
        for (auto&& run : ensemble) {
            ensemble_results.emplace_back(std::move(run.first));
            ensemble_params.emplace_back(std::move(run.second));
        }

        BOOST_OUTCOME_TRY(save_single_run_result);
        BOOST_OUTCOME_TRY(save_results(ensemble_results, ensemble_params, county_ids, result_dir, save_single_runs));
    }

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
    mio::mpi::init();

    bool late      = false;
    bool masks     = true;
    bool test      = true;
    bool high      = false;
    bool long_time = false;
    bool future    = false;

    RunMode mode;
    std::string save_dir;
    std::string data_dir;
    std::string result_dir;
    bool save_single_runs = true;
    if (argc == 10) {
        mode       = RunMode::Save;
        data_dir   = argv[1];
        save_dir   = argv[2];
        result_dir = argv[3];
        if (atoi(argv[4]) == 0) {
            save_single_runs = false;
        }
        if (atoi(argv[5]) == 1) {
            high = true;
        }
        else {
            high = false;
        }
        if (atoi(argv[6]) == 1) {
            late = true;
        }
        else {
            late = false;
        }
        if (atoi(argv[7]) == 1) {
            masks = true;
            test  = true;
        }
        else {
            masks = false;
            test  = false;
        }
        if (atoi(argv[8]) == 1) {
            long_time = true;
        }
        else {
            long_time = false;
        }
        if (atoi(argv[9]) == 1) {
            future = true;
        }
        else {
            future = false;
        }
        if (mio::mpi::is_root()) {
            printf("Options: masks set to: %d, late set to: %d, high set to: %d, long set to: %d, future set to: %d\n",
                   (int)masks, (int)late, (int)high, (int)long_time, (int)future);
            printf("Reading data from \"%s\", saving graph to \"%s\".\n", data_dir.c_str(), save_dir.c_str());
            printf("Exporting single run results and parameters: %d.\n", (int)save_single_runs);
        }
    }
    else if (argc == 4) {
        mode       = RunMode::Load;
        save_dir   = argv[1];
        result_dir = argv[2];
        data_dir   = "";
        if (mio::mpi::is_root()) {
            printf("Loading graph from \"%s\".\n", save_dir.c_str());
            printf("Exporting single run results and parameters: %d.\n", (int)save_single_runs);
        }
    }
    else if (argc == 5) {
        mode       = RunMode::Save;
        data_dir   = argv[1];
        save_dir   = argv[2];
        result_dir = argv[3];
        if (atoi(argv[4]) == 0) {
            save_single_runs = false;
        }
        if (mio::mpi::is_root()) {
            printf("Options: masks set to: %d, late set to: %d, high set to: %d, long set to: %d, future set to: %d\n",
                   (int)masks, (int)late, (int)high, (int)long_time, (int)future);
            printf("Reading data from \"%s\", saving graph to \"%s\".\n", data_dir.c_str(), save_dir.c_str());
            printf("Exporting single run results and parameters: %d.\n", (int)save_single_runs);
        }
    }
    else {
        if (mio::mpi::is_root()) {
            printf("Usage:\n");
            printf("2021_vaccination_delta <data_dir> <save_dir> <result_dir> <high> <late> <masks> <long> <future>\n");
            printf("\tMake graph with data from <data_dir> and save at <save_dir>, then run the simulation.\n");
            printf("\tStore the results in <result_dir>\n");
            printf("\t<high> <late> <masks> <long> <future> are either 0 or 1 to define a particular scenario\n");
            printf("2021_vaccination_delta <load_dir> <result_dir>\n");
            printf("\tLoad graph from <load_dir>, then run the simulation.\n");
        }
        mio::mpi::finalize();
        return 0;
    }

    if (future) {
        result_dir += "_future";
    }
    if (long_time) {
        result_dir += "_long";
    }
    if (high) {
        result_dir += "_high";
    }
    if (late) {
        result_dir += "_late";
    }
    if (masks) {
        result_dir += "_mask";
    }
    if (test) {
        result_dir += "_test";
    }
    if (mio::mpi::is_root()) {
        boost::filesystem::path dir(result_dir);
        bool created = boost::filesystem::create_directories(dir);

        if (created) {
            mio::log_info("Directory '{:s}' was created.", dir.string());
        }
        printf("Saving results to \"%s\".\n", result_dir.c_str());
    }

    auto result =
        run(mode, data_dir, save_dir, result_dir, save_single_runs, late, masks, test, high, long_time, future);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        mio::mpi::finalize();
        return -1;
    }
    mio::mpi::finalize();
    return 0;
}
