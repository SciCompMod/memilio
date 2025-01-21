/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Lena Ploetzke
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

#include "memilio/config.h"

#include "lct_secir/parameters_io.h"
#include "lct_secir/parameters.h"
#include "lct_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/date.h"
#include "memilio/utils/time_series.h"
#include "test_data_dir.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include <matchers.h>
#include "json/value.h"

#include <gtest/gtest.h>
#include <string>

// Define synthetic rki data vector without age groups for testing purposes.
std::vector<mio::ConfirmedCasesNoAgeEntry> get_synthetic_rki_data_noage()
{
    Json::Value js(Json::arrayValue);
    std::vector<Json::Value> dates = {"2020-05-26", "2020-05-27", "2020-05-28", "2020-05-29",
                                      "2020-05-30", "2020-05-31", "2020-06-01", "2020-06-02",
                                      "2020-06-03", "2020-06-04", "2020-06-05"};

    for (int day = 0; day < 11; day++) {
        js[day]["Date"]      = dates[day];
        js[day]["Confirmed"] = 100. + day;
        js[day]["Deaths"]    = 2.;
        js[day]["Recovered"] = 0.;
    }
    return mio::deserialize_confirmed_cases_noage(js).value();
}

// Define synthetic rki data vector with age groups for testing purposes.
std::vector<mio::ConfirmedCasesDataEntry> get_synthetic_rki_data_age()
{
    const int num_agegroups = 6;
    Json::Value js(Json::arrayValue);
    std::vector<Json::Value> dates           = {"2020-05-26", "2020-05-27", "2020-05-28", "2020-05-29",
                                                "2020-05-30", "2020-05-31", "2020-06-01", "2020-06-02",
                                                "2020-06-03", "2020-06-04", "2020-06-05"};
    std::vector<Json::Value> age_group_names = {"A00-A04", "A05-A14", "A15-A34", "A35-A59", "A60-A79", "A80+"};
    for (int day = 0; day < 11; day++) {
        for (int age = 0; age < num_agegroups; age++) {
            js[num_agegroups * day + age]["Age_RKI"]   = age_group_names[age];
            js[num_agegroups * day + age]["ID_County"] = 1001;
            js[num_agegroups * day + age]["Date"]      = dates[day];
            js[num_agegroups * day + age]["Confirmed"] = 100. + day * age;
            js[num_agegroups * day + age]["Deaths"]    = age;
            js[num_agegroups * day + age]["Recovered"] = 0.;
        }
    }
    return mio::deserialize_confirmed_cases_data(js).value();
}

// Check that initialization based on synthetic RKI data match a calculated result without age resolution.
TEST(TestLCTParametersIo, ReadPopulationDataRKI)
{
    // Define start date and the total population used for the initialization.
    std::vector<ScalarType> total_population(1, 1000.);
    auto start_date = mio::Date(2020, 6, 1);

    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, 2, 3, 2, 2, 1, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
    model.parameters.get<mio::lsecir::TimeExposed>()[0]            = 2.3;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 1.3;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 2.4;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 1.8;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 1.0;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.2;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.1;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.3;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.2;

    // Calculate initial value vector for subcompartments with RKI data.
    auto read_result =
        mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
            get_synthetic_rki_data_noage(), model.populations, model.parameters, start_date, total_population,
            std::vector<ScalarType>(1, 1.));
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    // Result to compare the simulation result with.
    Eigen::VectorX<ScalarType> compare((Eigen::Index)LctState::Count);
    // Calculate result using that the number of new confirmed cases at each day is one in the synthetic case data
    // and additionally the stay times, the number of subcompartments, and the transition probabilities.
    compare << 889.5, 1. / (2. * 0.8) * 2.3, 1. / (2. * 0.8) * 2.3, 1. / (3. * 0.8) * 1.3, 1. / (3. * 0.8) * 1.3,
        1. / (3. * 0.8) * 1.3, 2.4 / 2., 2.4 / 2., 1.8 / 2. * 0.1, 1.8 / 2. * 0.1, 0.1 * 0.3, 101.39, 2.;

    for (size_t i = 0; i < LctState::Count; i++) {
        EXPECT_NEAR(model.get_initial_values()[i], compare[i], 1e-4) << "at subcompartment number " << i;
    }
}

// Check that initialization based on synthetic RKI data match a calculated result with RKI age groups.
TEST(TestLCTParametersIo, ReadPopulationDataRKIAgeres)
{
    // Define start date and the total population used for the initialization.
    const size_t num_agegroups = 6;
    std::vector<ScalarType> total_population(num_agegroups, 1000.);
    auto start_date = mio::Date(2020, 6, 1);

    using InfState  = mio::lsecir::InfectionState;
    using LctState1 = mio::LctInfectionState<InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using LctState2 = mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model     = mio::lsecir::Model<LctState1, LctState2, LctState1, LctState2, LctState1, LctState2>;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t i = 0; i < num_agegroups; i++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[i]            = 2.3;
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[i] = 1.3;
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[i]   = 2.4;
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[i]     = 1.8;
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[i]   = 1.0;

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[i] = 0.2;
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[i]      = 0.1;
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[i]              = 0.3;
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[i]              = 0.2;
    }

    // Calculate initial value vector for subcompartments with RKI data.
    auto read_result =
        mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesDataEntry>(
            get_synthetic_rki_data_age(), model.populations, model.parameters, start_date, total_population,
            std::vector<ScalarType>(num_agegroups, 1.));
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    // Calculate result using that the number of new confirmed cases at each day is equal to the index of the age group
    // in the synthetic case data and additionally the stay times, the number of subcompartments,
    // and the transition probabilities.
    size_t num_populations = (size_t)InfState::Count * num_agegroups;
    Eigen::VectorX<ScalarType> compare(num_populations);
    for (size_t age = 0; age < num_agegroups; age++) {
        compare[(size_t)InfState::Count * age + (size_t)InfState::Exposed] =
            1. / (1. - model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[age]) *
            model.parameters.get<mio::lsecir::TimeExposed>()[age] * (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::InfectedNoSymptoms] =
            1. / (1. - model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[age]) *
            model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[age] * (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::InfectedSymptoms] =
            model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[age] * (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::InfectedSevere] =
            model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[age] *
            model.parameters.get<mio::lsecir::TimeInfectedSevere>()[age] * (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::InfectedCritical] =
            model.parameters.get<mio::lsecir::CriticalPerSevere>()[age] *
            model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[age] *
            model.parameters.get<mio::lsecir::TimeInfectedCritical>()[age] * (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::Dead] = (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::Recovered] =
            100. + (ScalarType)age * 5. -
            compare.segment((size_t)InfState::Count * age + (size_t)InfState::InfectedSymptoms, 3).sum();
        compare[(size_t)InfState::Count * age + (size_t)InfState::Susceptible] =
            total_population[age] -
            compare
                .segment((size_t)InfState::Count * age + (size_t)InfState::Susceptible + 1, (size_t)InfState::Count - 1)
                .sum();
    }

    mio::TimeSeries<ScalarType> result(model.populations.get_num_compartments());
    // Define TimeSeries as input for the function.
    result.add_time_point(0, model.get_initial_values());
    mio::TimeSeries<ScalarType> population = model.calculate_compartments(result);

    for (size_t i = 0; i < num_populations; i++) {
        EXPECT_NEAR(population.get_value(0)[i], compare[i], 1e-6) << "at subcompartment number " << i;
    }
}

// Check some cases where computation of initial values based on RKI data should fail.
TEST(TestLCTParametersIo, ReadPopulationDataRKIFailure)
{
    std::vector<ScalarType> total_population(1, 1000.);
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, 2, 3, 2, 2, 1, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Define parameters.
    model.parameters.get<mio::lsecir::TimeExposed>()[0]            = 2.3;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 1.3;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 2.4;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 1.8;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 1.0;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.2;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.1;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.3;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.2;

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Case with empty RKI data vector.
    Json::Value js(Json::arrayValue);
    auto start_date = mio::Date(2020, 6, 9);
    auto read_result =
        mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
            mio::deserialize_confirmed_cases_noage(js).value(), model.populations, model.parameters, start_date,
            total_population, std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidFileFormat));

    // Case where start_date is later than maximal provided date in file.
    read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
        get_synthetic_rki_data_noage(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::OutOfRange));

    // Case where not all required dates are provided.
    start_date  = mio::Date(2020, 6, 3);
    read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
        get_synthetic_rki_data_noage(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::OutOfRange));

    // Case with input values that lead to negative entries.
    start_date  = mio::Date(2020, 6, 1);
    read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
        get_synthetic_rki_data_noage(), model.populations, model.parameters, start_date, std::vector<ScalarType>(1, 0.),
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidValue));

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    // Valid case.
    read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
        get_synthetic_rki_data_noage(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(1, 1.));
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
}

// Check some cases where computation of initial values with age groups based on RKI data should fail.
TEST(TestLCTParametersIo, ReadPopulationDataRKIFailureAgeres)
{
    const size_t num_agegroups = 6;
    std::vector<ScalarType> total_population(num_agegroups, 1e6);

    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using Model    = mio::lsecir::Model<LctState, LctState, LctState, LctState, LctState, LctState>;
    Model model;

    // Define parameters.
    for (size_t i = 0; i < num_agegroups; i++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[i]            = 2.3;
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[i] = 1.3;
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[i]   = 2.4;
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[i]     = 1.8;
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[i]   = 1.0;

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[i] = 0.2;
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[i]      = 0.1;
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[i]              = 0.3;
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[i]              = 0.2;
    }

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Case with empty RKI data vector.
    auto start_date = mio::Date(2020, 6, 9);
    Json::Value js(Json::arrayValue);
    auto read_result =
        mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesDataEntry>(
            mio::deserialize_confirmed_cases_data(js).value(), model.populations, model.parameters, start_date,
            total_population, std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidFileFormat));

    // Case where start_date is later than maximal provided date in file.
    start_date  = mio::Date(2021, 1, 1);
    read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesDataEntry>(
        get_synthetic_rki_data_age(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::OutOfRange));

    // Case where not all required dates are provided.
    start_date  = mio::Date(2020, 6, 3);
    read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesDataEntry>(
        get_synthetic_rki_data_age(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::OutOfRange));

    // Case with input values that lead to negative entries.
    start_date  = mio::Date(2020, 6, 1);
    read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesDataEntry>(
        get_synthetic_rki_data_age(), model.populations, model.parameters, start_date,
        std::vector<ScalarType>(num_agegroups, 0.), std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidValue));

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    // Valid case.
    read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesDataEntry>(
        get_synthetic_rki_data_age(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
}
