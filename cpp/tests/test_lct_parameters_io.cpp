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

// Check that initialization based on synthetic RKI data match previous result.
TEST(TestLCTParametersIo, ReadPopulationDataRKI)
{
    // Define start date and the total population used for the initialization.
    ScalarType total_population = 1000.0;
    auto start_date             = mio::Date(2020, 6, 1);

    // Use rounded stay times for the number of subcompartments except for InfectedCritical.
    // Template parameters must be compile-time constants, which is why the values are not calculated with parameters.
    using Model    = mio::lsecir::Model<2, 3, 2, 2, 1>;
    using LctState = Model::LctState;
    Model model;

    // Define parameters used for simulation and initialization.
    mio::lsecir::Parameters parameters;
    model.parameters.get<mio::lsecir::TimeExposed>()            = 2.3;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 3.3;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 2.4;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()     = 1.8;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()   = 3.0;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.2;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.1;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.3;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.2;

    // Calculate initial value vector for subcompartments with RKI data.
    auto read_result = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date, total_population, 1.);

    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    // Result of a previous simulation.
    Eigen::VectorXd compare((Eigen::Index)LctState::Count);
    compare << 863.05, 14.30625, 8.53125, 30.1125, 36.1875, 3.8125, 9.88, 3.52, 0.09, 0.25, 0.6888, 27.8712, 1.7;

    for (size_t i = 0; i < LctState::Count; i++) {
        EXPECT_NEAR(model.get_initial_values()[i], compare[i], 1e-4) << "at subcompartment number " << i;
    }
}

// Check some cases where computation of initial values for an LCT model based on RKI data should fail.
TEST(TestLCTParametersIo, ReadPopulationDataRKIFailure)
{
    ScalarType total_population = 1000.0;
    using Model                 = mio::lsecir::Model<2, 3, 2, 2, 1>;
    Model model;
    // Define parameters used for simulation and initialization.
    mio::lsecir::Parameters parameters;
    model.parameters.get<mio::lsecir::TimeExposed>()            = 2.3;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 3.3;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 2.4;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()     = 1.8;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()   = 3.0;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.2;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.08;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.15;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.22;

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Case where start_date is later than maximal provided date in file.
    auto start_date   = mio::Date(2020, 6, 9);
    auto read_result1 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date, total_population, 1.);

    EXPECT_THAT(print_wrap(read_result1), IsFailure(mio::StatusCode::OutOfRange));

    // Case where not all needed dates are provided.
    start_date        = mio::Date(2020, 6, 6);
    auto read_result2 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date, total_population, 1.);

    EXPECT_THAT(print_wrap(read_result2), IsFailure(mio::StatusCode::OutOfRange));

    // Case with empty RKI data file.
    auto read_result3 = mio::lsecir::set_initial_data_from_confirmed_cases<Model>(
        model, mio::path_join(TEST_DATA_DIR, "test_empty_file.json"), start_date, total_population, 1.);

    EXPECT_THAT(print_wrap(read_result3), IsFailure(mio::StatusCode::InvalidFileFormat));

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}