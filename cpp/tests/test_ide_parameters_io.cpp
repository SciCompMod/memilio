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

#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "ide_secir/parameters_io.h"

#include "memilio/config.h"
#include "memilio/utils/date.h"
#include "memilio/io/io.h"

#include "test_data_dir.h"
#include "load_test_data.h"
#include "matchers.h"
#include <gtest/gtest.h>

// Check that initialization based on synthetic RKI data match previous result.
TEST(TestIDEParametersIo, RKIcompareWithPreviousRun)
{
    // Define start date and the total population used for the initialization.
    ScalarType total_population = 80 * 1e6;
    auto start_date             = mio::Date(2020, 6, 3);
    ScalarType dt               = 0.5;
    // Initialize model.
    // The number of deaths will be overwritten if real data is used for initialization. Therefore, an arbitrary number is used for the number of deaths.
    mio::isecir::Model model(mio::TimeSeries<ScalarType>(0., mio::TimeSeries<ScalarType>::Vector::Constant(
                                                                 (int)mio::isecir::InfectionTransition::Count, 1.)),
                             total_population, -1);

    // Set the model parameters so that if the default values are changed, the test is still valid.
    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib((int)mio::isecir::InfectionTransition::Count,
                                                               delaydistribution);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere].set_parameter(1.7);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constfunc(1.0);
    mio::StateAgeFunctionWrapper prob(constfunc);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);

    // Calculate initialization.
    auto status =
        mio::isecir::set_initial_flows(model, dt, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date);

    ASSERT_THAT(print_wrap(status), IsSuccess());

    // Compare with previous run.
    auto compare = load_test_data_csv<ScalarType>("ide-parameters-io-compare.csv");

    ASSERT_EQ(compare.size(), static_cast<size_t>(model.m_transitions.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(model.m_transitions.get_num_elements()) + 1) << "at row " << i;
        ASSERT_NEAR(compare[i][0], model.m_transitions.get_time(i), 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            ASSERT_NEAR(compare[i][j], model.m_transitions.get_value(i)[j - 1], 1e-7) << " at row " << i;
        }
    }
}

// Check some cases where computation of initial values for an IDE model based on RKI data should fail.
TEST(TestIDEParametersIo, ParametersIoRKIFailure)
{
    // Define start date and the total population used for the initialization.
    ScalarType total_population = 80 * 1e6;
    ScalarType dt               = 0.5;

    // Initialize model.
    // The number of deaths will be overwritten if real data is used for initialization. Therefore, an arbitrary number is used for the number of deaths.
    mio::isecir::Model model(mio::TimeSeries<ScalarType>((int)mio::isecir::InfectionTransition::Count),
                             total_population, -1);

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // --- Case where start_date is later than maximal provided date in file.
    auto start_date = mio::Date(2020, 6, 9);
    auto status =
        mio::isecir::set_initial_flows(model, dt, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date);

    ASSERT_THAT(print_wrap(status), IsFailure(mio::StatusCode::OutOfRange));

    // --- Case where not all needed dates are provided.
    start_date = mio::Date(2020, 6, 7);
    status =
        mio::isecir::set_initial_flows(model, dt, mio::path_join(TEST_DATA_DIR, "cases_all_germany.json"), start_date);

    ASSERT_THAT(print_wrap(status), IsFailure(mio::StatusCode::OutOfRange));

    // --- Case with empty RKI data file.
    status =
        mio::isecir::set_initial_flows(model, dt, mio::path_join(TEST_DATA_DIR, "test_empty_file.json"), start_date);

    ASSERT_THAT(print_wrap(status), IsFailure(mio::StatusCode::InvalidFileFormat));

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}