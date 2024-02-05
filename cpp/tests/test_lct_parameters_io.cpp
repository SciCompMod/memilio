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
#include "test_data_dir.h"
#include "memilio/io/io.h"
#include <matchers.h>

#include <gtest/gtest.h>

// Check that Initialization based on synthetic RKI data match previous result.
TEST(TestLCTParametersIo, ReadPopulationDataRKI)
{
    ScalarType total_population = 1000.0;
    auto start_date             = mio::Date(2020, 6, 1);

    // Define parameters used for simulation and initialization.
    mio::lsecir::Parameters parameters;
    parameters.get<mio::lsecir::TimeExposed>()            = 2.3;
    parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 3.3;
    parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 2.4;
    parameters.get<mio::lsecir::TimeInfectedSevere>()     = 1.8;
    parameters.get<mio::lsecir::TimeInfectedCritical>()   = 3.0;

    parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.2;
    parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.1;
    parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.3;
    parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.2;

    // Define number of subcompartments.
    std::vector<int> vec_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    // Use subcompartments with a soujourn time of approximately one day in each subcompartment.
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed] =
        (int)round(parameters.get<mio::lsecir::TimeExposed>());
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] =
        (int)round(parameters.get<mio::lsecir::TimeInfectedNoSymptoms>());
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] =
        (int)round(parameters.get<mio::lsecir::TimeInfectedSymptoms>());
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere] =
        (int)round(parameters.get<mio::lsecir::TimeInfectedSevere>());
    // Both realistic distributions for times corresponding to InfectedCritical of the IDE model are exponential distributions.
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical] = 1;
    mio::lsecir::InfectionState infectionState(vec_subcompartments);

    // Calculate initial value vector for subcompartments with RKI data.
    auto read_result =
        mio::lsecir::get_initial_data_from_file(mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date,
                                                infectionState, std::move(parameters), total_population, 1.);

    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    auto init_subcompartments = read_result.value();
    // Previous result.
    Eigen::VectorXd compare(infectionState.get_count());
    compare << 863.05, 14.30625, 8.53125, 30.1125, 36.1875, 3.8125, 9.88, 3.52, 0.09, 0.25, 0.6888, 27.8712, 1.7;

    for (int i = 0; i < infectionState.get_count(); i++) {
        EXPECT_NEAR(init_subcompartments[i], compare[i], 1e-4) << "at subcompartment number " << i;
    }
}

// Check some cases where computation of initial values for an LCT model based on RKI data should fail.
TEST(TestLCTParametersIo, ReadPopulationDataRKIFailure)
{
    ScalarType total_population = 1000.0;

    // Define parameters used for simulation and initialization.
    mio::lsecir::Parameters parameters;
    parameters.get<mio::lsecir::TimeExposed>()            = 2.3;
    parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 3.3;
    parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 2.4;
    parameters.get<mio::lsecir::TimeInfectedSevere>()     = 1.8;
    parameters.get<mio::lsecir::TimeInfectedCritical>()   = 3.0;

    parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.2;
    parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.08;
    parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.15;
    parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.22;

    // Define number of subcompartments.
    std::vector<int> vec_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    // Use subcompartments with a soujourn time of approximately one day in each subcompartment.
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed] =
        (int)round(parameters.get<mio::lsecir::TimeExposed>());
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] =
        (int)round(parameters.get<mio::lsecir::TimeInfectedNoSymptoms>());
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] =
        (int)round(parameters.get<mio::lsecir::TimeInfectedSymptoms>());
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere] =
        (int)round(parameters.get<mio::lsecir::TimeInfectedSevere>());
    // Both realistic distributions for times corresponding to InfectedCritical of the IDE model are exponential distributions.
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical] = 1;
    mio::lsecir::InfectionState infectionState(vec_subcompartments);

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Case where start_date is later than maximal provided date in file.
    auto start_date = mio::Date(2020, 6, 9);
    auto read_result1 =
        mio::lsecir::get_initial_data_from_file(mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date,
                                                infectionState, std::move(parameters), total_population, 1.);

    ASSERT_THAT(print_wrap(read_result1), IsFailure(mio::StatusCode::OutOfRange));
    // Case where not all needed dates are provided.
    start_date = mio::Date(2020, 6, 6);
    auto read_result2 =
        mio::lsecir::get_initial_data_from_file(mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date,
                                                infectionState, std::move(parameters), total_population, 1.);

    ASSERT_THAT(print_wrap(read_result2), IsFailure(mio::StatusCode::OutOfRange));

    // Case with empty RKI data file.
    auto read_result3 =
        mio::lsecir::get_initial_data_from_file(mio::path_join(TEST_DATA_DIR, "test_empty_file.json"), start_date,
                                                infectionState, std::move(parameters), total_population, 1.);

    ASSERT_THAT(print_wrap(read_result3), IsFailure(mio::StatusCode::InvalidFileFormat));

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}