/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "memilio/io/io.h"
#include <matchers.h>

#include <gtest/gtest.h>

// Check that initialization based on synthetic RKI data match previous result.
TEST(TestLCTParametersIo, ReadPopulationDataRKI)
{
    // Define start date and the total population used for the initialization.
    ScalarType total_population = 1000.0;
    auto start_date             = mio::Date(2020, 6, 1);

    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, 2, 3, 2, 2, 1, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
    model.parameters.get<mio::lsecir::TimeExposed>()[0]            = 2.3;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 3.3;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 2.4;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 1.8;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 3.0;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.2;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.1;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.3;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.2;

    // Calculate initial value vector for subcompartments with RKI data.
    auto read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date,
        std::vector<ScalarType>(1, total_population), std::vector<ScalarType>(1, 1.));

    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    // Result of a previous simulation.
    Eigen::VectorXd compare((Eigen::Index)LctState::Count);
    compare << 863.05, 14.30625, 8.53125, 30.1125, 36.1875, 3.8125, 9.88, 3.52, 0.09, 0.25, 0.6888, 27.8712, 1.7;

    for (size_t i = 0; i < LctState::Count; i++) {
        EXPECT_NEAR(model.get_initial_values()[i], compare[i], 1e-4) << "at subcompartment number " << i;
    }
}

// Check that initialization based on synthetic RKI data match previous result with RKI AgeGroups.
TEST(TestLCTParametersIo, ReadPopulationDataRKIAgeres)
{
    // Define start date and the total population used for the initialization.
    const size_t num_agegroups = 6;
    std::vector<ScalarType> total_population(num_agegroups, 1e6);
    auto start_date = mio::Date(2020, 12, 1);

    using InfState  = mio::lsecir::InfectionState;
    using LctState1 = mio::LctInfectionState<InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using LctState2 = mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model     = mio::lsecir::Model<LctState1, LctState2, LctState1, LctState2, LctState1, LctState2>;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t i = 0; i < num_agegroups; i++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[i]            = 2.3;
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[i] = 3.3;
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[i]   = 2.4;
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[i]     = 1.8;
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[i]   = 3.0;

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[i] = 0.2;
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[i]      = 0.1;
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[i]              = 0.3;
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[i]              = 0.2;
    }

    // Calculate initial value vector for subcompartments with RKI data.
    auto read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_age_ma7.json"), start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));

    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    // Result of a previous simulation.
    size_t num_populations = (size_t)InfState::Count * num_agegroups;
    Eigen::VectorXd compare(num_populations);
    compare << 977209.6429, 995.7143, 1372.2143, 793.8, 59.5086, 29.946, 19538.174, 1, 923916.3214, 3781.9821,
        5213.8393, 2969.4571, 226.3171, 114.4997, 63773.8974, 3.6857, 613901.8214, 15635.4821, 21071.8393, 11780.1143,
        882.3657, 447.9634, 336238.6423, 41.7714, 530413.3571, 21446.7321, 28762.9107, 15917.2857, 1185.8914, 596.4754,
        400751.4617, 925.8857, 828383.3571, 8481.8214, 11373.3929, 6158.1429, 451.9086, 224.6966, 138585.452, 6341.2286,
        909054.8214, 6074.875, 8177.3036, 4393.3714, 317.3629, 147.8314, 58419.12, 13415.3143;

    mio::TimeSeries<ScalarType> result(model.populations.get_num_compartments());
    // Define TimeSeries as input for the function.
    result.add_time_point(0, model.get_initial_values());
    mio::TimeSeries<ScalarType> population = model.calculate_compartments(result);

    for (size_t i = 0; i < num_populations; i++) {
        EXPECT_NEAR(population.get_value(0)[i], compare[i], 1e-3) << "at subcompartment number " << i;
    }
}

// Check some cases where computation of initial values for an LCT model based on RKI data should fail.
TEST(TestLCTParametersIo, ReadPopulationDataRKIFailure)
{
    std::vector<ScalarType> total_population(1, 1000.);
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, 2, 3, 2, 2, 1, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
    model.parameters.get<mio::lsecir::TimeExposed>()[0]            = 2.3;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 3.3;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 2.4;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 1.8;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 3.0;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.2;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.08;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.15;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.22;

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Case with empty RKI data file.
    auto start_date  = mio::Date(2020, 6, 9);
    auto read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "test_empty_file.json"), start_date, total_population,
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidFileFormat));

    // Case where start_date is later than maximal provided date in file.
    auto read_result1 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date, total_population,
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result1), IsFailure(mio::StatusCode::OutOfRange));

    // Case where not all required dates are provided.
    start_date        = mio::Date(2020, 6, 6);
    auto read_result2 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date, total_population,
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result2), IsFailure(mio::StatusCode::OutOfRange));

    // Case with input values that are not fitting the data.
    start_date        = mio::Date(2020, 6, 1);
    auto read_result3 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date, std::vector<ScalarType>(1, 0.),
        std::vector<ScalarType>(1, 1.));
    EXPECT_THAT(print_wrap(read_result3), IsFailure(mio::StatusCode::InvalidValue));

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    // Valid case.
    auto read_result4 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date, total_population,
        std::vector<ScalarType>(1, 1.));
    ASSERT_THAT(print_wrap(read_result4), IsSuccess());
}

// Check some cases where computation of initial values for an LCT model with age groups based on RKI data should fail.
TEST(TestLCTParametersIo, ReadPopulationDataRKIFailureAgeres)
{
    const size_t num_agegroups = 6;
    std::vector<ScalarType> total_population(num_agegroups, 1e6);
    auto start_date = mio::Date(2020, 12, 1);

    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using Model    = mio::lsecir::Model<LctState, LctState, LctState, LctState, LctState, LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t i = 0; i < num_agegroups; i++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[i]            = 2.3;
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[i] = 3.3;
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[i]   = 2.4;
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[i]     = 1.8;
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[i]   = 3.0;

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[i] = 0.2;
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[i]      = 0.08;
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[i]              = 0.15;
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[i]              = 0.22;
    }

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Case with empty RKI data file.
    auto read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "test_empty_file.json"), start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result), IsFailure(mio::StatusCode::InvalidFileFormat));

    // Case where start_date is later than maximal provided date in file.
    start_date        = mio::Date(2021, 1, 1);
    auto read_result1 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_age_ma7.json"), start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result1), IsFailure(mio::StatusCode::OutOfRange));

    // Case where not all required dates are provided.
    start_date        = mio::Date(2020, 12, 30);
    auto read_result2 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_age_ma7.json"), start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result2), IsFailure(mio::StatusCode::OutOfRange));

    // Case with input values that are not fitting the data.
    start_date        = mio::Date(2020, 12, 1);
    auto read_result3 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_age_ma7.json"), start_date,
        std::vector<ScalarType>(num_agegroups, 0.), std::vector<ScalarType>(num_agegroups, 1.));
    EXPECT_THAT(print_wrap(read_result3), IsFailure(mio::StatusCode::InvalidValue));

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    // Valid case.
    auto read_result4 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_age_ma7.json"), start_date, total_population,
        std::vector<ScalarType>(num_agegroups, 1.));
    ASSERT_THAT(print_wrap(read_result4), IsSuccess());
}
