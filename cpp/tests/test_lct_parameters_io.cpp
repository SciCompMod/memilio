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
#include "memilio/utils/logging.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "test_data_dir.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "matchers.h"

#include "json/value.h"
#include <gtest/gtest.h>
#include <string>

// Define synthetic rki data vector without age groups for testing purposes.
std::vector<mio::ConfirmedCasesNoAgeEntry> get_synthetic_rki_data_noage()
{
    Json::Value js(Json::arrayValue);
    std::vector<Json::Value> dates = {"2020-05-26", "2020-05-27", "2020-05-28", "2020-05-29",
                                      "2020-05-30", "2020-05-31", "2020-06-01", "2020-06-02",
                                      "2020-06-03", "2020-06-04", "2020-06-05", "2020-06-06"};

    for (int day = 0; day < 12; day++) {
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
                                                "2020-06-03", "2020-06-04", "2020-06-05", "2020-06-06"};
    std::vector<Json::Value> age_group_names = {"A00-A04", "A05-A14", "A15-A34", "A35-A59", "A60-A79", "A80+"};
    for (int day = 0; day < 12; day++) {
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

// Define synthetic rki data vector with age groups for testing purposes.
std::vector<mio::DiviEntry> get_synthetic_divi_data()
{
    Json::Value js(Json::arrayValue);
    std::vector<Json::Value> dates = {"2020-05-26", "2020-05-27", "2020-05-28", "2020-05-29",
                                      "2020-05-30", "2020-05-31", "2020-06-01"};
    for (int day = 0; day < 7; day++) {
        js[day]["ICU"]  = 0;
        js[day]["Date"] = dates[day];
    }
    js[6]["ICU"] = 50;
    return mio::deserialize_divi_data(js).value();
}

// Check that initialization based on synthetic RKI data match a calculated result without age resolution.
TEST(TestLCTParametersIo, ReadPopulationDataRKI)
{
    // Define start date and the total population used for the initialization.
    std::vector<ScalarType> total_population(1, 1000.);
    auto start_date = mio::Date(2020, 6, 1);

    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 1, 1, 1>;
    using Model    = mio::lsecir::Model<ScalarType, LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
    model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[0]            = 2.3;
    model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0] = 1.3;
    model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0]   = 2.4;
    model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0]     = 1.8;
    model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0]   = 1.0;

    model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0] = 0.2;
    model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0]      = 0.1;
    model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0]              = 0.3;
    model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0]              = 0.2;

    // Calculate initial value vector for subcompartments with RKI data.
    auto read_result =
        mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
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
    using LctState1 = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using LctState2 = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model     = mio::lsecir::Model<ScalarType, LctState1, LctState2, LctState1, LctState2, LctState1, LctState2>;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t i = 0; i < num_agegroups; i++) {
        model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[i]            = 2.3;
        model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[i] = 1.3;
        model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[i]   = 2.4;
        model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[i]     = 1.8;
        model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[i]   = 1.0;

        model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i] = 0.2;
        model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[i]      = 0.1;
        model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[i]              = 0.3;
        model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[i]              = 0.2;
    }

    // Calculate initial value vector for subcompartments with RKI data.
    auto read_result =
        mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
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
            1. / (1. - model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[age]) *
            model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[age] * (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::InfectedNoSymptoms] =
            1. / (1. - model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[age]) *
            model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[age] * (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::InfectedSymptoms] =
            model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[age] * (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::InfectedSevere] =
            model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[age] *
            model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[age] * (ScalarType)age;
        compare[(size_t)InfState::Count * age + (size_t)InfState::InfectedCritical] =
            model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[age] *
            model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[age] *
            model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[age] * (ScalarType)age;
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

// Check that the scaling using DIVI data is working as expected.
TEST(TestLCTParametersIo, CheckScalingDIVI)
{
    // Define start date and the total population used for the initialization.
    const size_t num_agegroups = 6;
    std::vector<ScalarType> total_population(num_agegroups, 1000.);
    auto start_date = mio::Date(2020, 6, 1);

    using InfState  = mio::lsecir::InfectionState;
    using LctState1 = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using LctState2 = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model     = mio::lsecir::Model<ScalarType, LctState1, LctState2, LctState1, LctState2, LctState1, LctState2>;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t i = 0; i < num_agegroups; i++) {
        model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[i]            = 2.3;
        model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[i] = 1.3;
        model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[i]   = 2.4;
        model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[i]     = 1.8;
        model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[i]   = 1.0;

        model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i] = 0.2;
        model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[i]      = 0.1;
        model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[i]              = 0.3;
        model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[i]              = 0.2;
    }

    // Check that the function get_icu_from_divi_data to get the DIVI data works as expected.
    EXPECT_NEAR(mio::lsecir::details::get_icu_from_divi_data(get_synthetic_divi_data(), start_date).value(), 50, 1e-6);
    // Check additionally for another date.
    EXPECT_NEAR(mio::lsecir::details::get_icu_from_divi_data(get_synthetic_divi_data(), mio::Date(2020, 5, 31)).value(),
                0, 1e-6);

    // Check that the function set_initial_values_from_reported_data WITHOUT DIVI data returns a different number of
    // total InfectedCritical cases than the reported divi data. Otherwise the test case makes no sense.
    auto read_result =
        mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
            get_synthetic_rki_data_age(), model.populations, model.parameters, start_date, total_population,
            std::vector<ScalarType>(num_agegroups, 1.));
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    ScalarType total_InfectedCritical =
        mio::lsecir::details::get_total_InfectedCritical_from_populations<Model::Populations>(model.populations);
    ASSERT_NE(total_InfectedCritical, 50) << "The test does not work because the RKI data already lead to the value "
                                             "from the DIVI data. Please change the synthetic test data.";

    // Check that the function set_initial_values_from_reported_data WITH DIVI data calculates an initial value vector
    // with the correct number of InfectedCritical individuals as defined by the DIVI data.
    read_result = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
        get_synthetic_rki_data_age(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.), get_synthetic_divi_data());
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    total_InfectedCritical =
        mio::lsecir::details::get_total_InfectedCritical_from_populations<Model::Populations>(model.populations);
    EXPECT_NEAR(total_InfectedCritical, 50, 1e-6);

    // Check that the total number of individuals in each age group is not affected.
    // Not tested in for loop as we need constexpr template arguments.
    EXPECT_NEAR(model.populations.get_group_total<0>(), total_population[0], 1e-6);
    EXPECT_NEAR(model.populations.get_group_total<1>(), total_population[1], 1e-6);
    EXPECT_NEAR(model.populations.get_group_total<2>(), total_population[2], 1e-6);
    EXPECT_NEAR(model.populations.get_group_total<3>(), total_population[3], 1e-6);
    EXPECT_NEAR(model.populations.get_group_total<4>(), total_population[4], 1e-6);
    EXPECT_NEAR(model.populations.get_group_total<5>(), total_population[5], 1e-6);

    // Check that we get an error if the date is not part of the input.
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);
    auto read_result_failure =
        mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
            get_synthetic_rki_data_age(), model.populations, model.parameters, mio::Date(2020, 6, 2), total_population,
            std::vector<ScalarType>(num_agegroups, 1.), get_synthetic_divi_data());
    EXPECT_THAT(print_wrap(read_result_failure), IsFailure(mio::StatusCode::OutOfRange));
    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

// Check that the function rescale_to_divi_data() handles each input case well.
TEST(TestLCTParametersIo, CheckRescaleToDIVIDataFunctionCases)
{
    const size_t num_agegroups = 6;
    using InfState             = mio::lsecir::InfectionState;
    using LctState1            = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 5, 1, 1>;
    using LctState2            = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model = mio::lsecir::Model<ScalarType, LctState1, LctState2, LctState1, LctState2, LctState1, LctState2>;
    using Populations = Model::Populations;

    // Initialize a population with zero values.
    Populations pop;
    for (size_t i = 0; i < pop.get_num_compartments(); i++) {
        pop[i] = 0;
    }
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);
    // Check that we get an error if the input for the reported InfectedCritical cases is negative.
    auto status = mio::lsecir::details::rescale_to_divi_data<Populations>(pop, -50, 0);
    EXPECT_THAT(print_wrap(status), IsFailure(mio::StatusCode::InvalidValue));
    // Check that we get an error if the number of Recovered is less than zero after scaling.
    status = mio::lsecir::details::rescale_to_divi_data<Populations>(pop, 50, 0);
    EXPECT_THAT(print_wrap(status), IsFailure(mio::StatusCode::InvalidValue));
    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    // Construct valid case with non-zero values for Recovered and Susceptible compartments.
    // For groups with LctState1.
    for (size_t i = 0; i < 3; i++) {
        size_t first_idx = i * (LctState1::Count + LctState2::Count);
        pop[first_idx + LctState1::get_first_index<InfState::Susceptible>()] = 1000;
        pop[first_idx + LctState1::get_first_index<InfState::Recovered>()]   = 100;
    }
    // For groups with LctState2.
    for (size_t i = 0; i < 3; i++) {
        size_t first_idx = (i + 1) * LctState1::Count + i * LctState2::Count;
        pop[first_idx + LctState2::get_first_index<InfState::Susceptible>()] = 3000;
        pop[first_idx + LctState2::get_first_index<InfState::Recovered>()]   = 300;
    }

    status = mio::lsecir::details::rescale_to_divi_data<Populations>(pop, 50, 0);
    ASSERT_THAT(print_wrap(status), IsSuccess());
    // Check that the reported value of 50 is distributed uniformly as expected.
    // For groups with LctState1.
    for (size_t i = 0; i < 3; i++) {
        size_t first_idx =
            i * (LctState1::Count + LctState2::Count) + LctState1::get_first_index<InfState::InfectedCritical>();
        for (size_t subcomp = 0; subcomp < LctState1::get_num_subcompartments<InfState::InfectedCritical>();
             subcomp++) {
            EXPECT_NEAR(pop[first_idx + subcomp],
                        50. / ((ScalarType)num_agegroups) *
                            (1. / ((ScalarType)LctState1::get_num_subcompartments<InfState::InfectedCritical>())),
                        1e-6);
        }
    }
    // For groups with LctState2.
    for (size_t i = 0; i < 3; i++) {
        size_t first_idx = (i + 1) * LctState1::Count + i * LctState2::Count;
        EXPECT_NEAR(pop[first_idx + LctState2::get_first_index<InfState::InfectedCritical>()],
                    50. / ((ScalarType)num_agegroups), 1e-6);
    }
}

// Check some cases where computation of initial values based on RKI data should fail.
TEST(TestLCTParametersIo, ReadPopulationDataRKIFailure)
{
    std::vector<ScalarType> total_population(1, 1000.);
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 1, 1, 1>;
    using Model    = mio::lsecir::Model<ScalarType, LctState>;
    Model model;

    // Define parameters.
    model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[0]            = 2.3;
    model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0] = 1.3;
    model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0]   = 2.4;
    model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0]     = 1.8;
    model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0]   = 1.0;

    model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0] = 0.2;
    model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0]      = 0.1;
    model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0]              = 0.3;
    model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0]              = 0.2;

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Case with empty RKI data vector.
    Json::Value js(Json::arrayValue);
    auto start_date = mio::Date(2020, 6, 9);
    auto read_result =
        mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
            mio::deserialize_confirmed_cases_noage(js).value(), model.populations, model.parameters, start_date,
            total_population, std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidFileFormat));

    // Case where start_date is later than maximal provided date in file.
    read_result = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
        get_synthetic_rki_data_noage(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::OutOfRange));

    // Case where not all required dates are provided.
    start_date  = mio::Date(2020, 6, 3);
    read_result = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
        get_synthetic_rki_data_noage(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::OutOfRange));

    // Case with input values that lead to negative entries.
    start_date  = mio::Date(2020, 6, 1);
    read_result = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
        get_synthetic_rki_data_noage(), model.populations, model.parameters, start_date, std::vector<ScalarType>(1, 0.),
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidValue));

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    // Valid case.
    read_result = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesNoAgeEntry>(
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
    using LctState = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using Model    = mio::lsecir::Model<ScalarType, LctState, LctState, LctState, LctState, LctState, LctState>;
    Model model;

    // Define parameters.
    for (size_t i = 0; i < num_agegroups; i++) {
        model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[i]            = 2.3;
        model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[i] = 1.3;
        model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[i]   = 2.4;
        model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[i]     = 1.8;
        model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[i]   = 1.0;

        model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i] = 0.2;
        model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[i]      = 0.1;
        model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[i]              = 0.3;
        model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[i]              = 0.2;
    }

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Case with empty RKI data vector.
    auto start_date = mio::Date(2020, 6, 9);
    Json::Value js(Json::arrayValue);
    auto read_result =
        mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
            mio::deserialize_confirmed_cases_data(js).value(), model.populations, model.parameters, start_date,
            total_population, std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidFileFormat));

    // Case where start_date is later than maximal provided date in file.
    start_date  = mio::Date(2021, 1, 1);
    read_result = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
        get_synthetic_rki_data_age(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::OutOfRange));

    // Case where not all required dates are provided.
    start_date  = mio::Date(2020, 6, 3);
    read_result = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
        get_synthetic_rki_data_age(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::OutOfRange));

    // Case with input values that lead to negative entries.
    start_date  = mio::Date(2020, 6, 1);
    read_result = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
        get_synthetic_rki_data_age(), model.populations, model.parameters, start_date,
        std::vector<ScalarType>(num_agegroups, 0.), std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidValue));

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    // Valid case.
    read_result = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
        get_synthetic_rki_data_age(), model.populations, model.parameters, start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
}
