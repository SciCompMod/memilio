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

#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "ide_secir/parameters_io.h"
#include "ide_secir/simulation.h"

#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/io/epi_data.h"
#include "memilio/utils/date.h"
#include "memilio/utils/logging.h"
#include "memilio/io/io.h"

#include "test_data_dir.h"
#include "load_test_data.h"
#include "matchers.h"
#include <gtest/gtest.h>
#include <vector>

// Check that initialization based on synthetic RKI data matches previous result for an IDE model without age resolution.
TEST(TestIDEParametersIo, RKIcompareWithPreviousRun)
{
    // Define number of age groups.
    size_t num_agegroups = 1;
    // Since we only have one age group, the corresponding index is set to 0.
    mio::AgeGroup group = 0;

    // Define start date and the total population used for the initialization.
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 15 * 1e6);
    auto start_date = mio::Date(2020, 06, 2);
    ScalarType dt   = 0.5;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;
    // Initialize model.
    // The number of deaths will be overwritten if real data is used for initialization. Therefore, an arbitrary number
    // is used for the number of deaths.
    // Initial time series for the flows will be also overridden.
    mio::isecir::Model model(
        mio::TimeSeries<ScalarType>(-1, mio::TimeSeries<ScalarType>::Vector::Constant(
                                            (int)mio::isecir::InfectionTransition::Count * num_agegroups, 1.)),
        total_population, mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), -1.),
        num_agegroups);

    // Set the model parameters so that if the default values are changed, the test is still valid.
    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib((int)mio::isecir::InfectionTransition::Count,
                                                               delaydistribution);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]
        .set_distribution_parameter(1.7);

    model.parameters.get<mio::isecir::TransitionDistributions>()[group] = vec_delaydistrib;

    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;

    model.parameters.get<mio::isecir::TransitionProbabilities>()[group] = vec_prob;

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, num_agegroups);
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_agegroups, num_agegroups, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constfunc(1.0);
    mio::StateAgeFunctionWrapper prob(constfunc);

    model.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[group] = prob;
    model.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[group]   = prob;
    model.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[group]   = prob;

    mio::CustomIndexArray<ScalarType, mio::AgeGroup> scale_confirmed_cases =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 1.);

    // Calculate initialization.
    auto status = mio::isecir::set_initial_flows<mio::ConfirmedCasesNoAgeEntry>(
        model, dt, mio::read_confirmed_cases_noage(mio::path_join(TEST_DATA_DIR, "cases_all_germany.json")).value(),
        start_date, scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsSuccess());

    // Compare with previous run.

    std::vector<ScalarType> deaths = {6.};

    std::vector<ScalarType> total_confirmed_cases = {44.5};

    int Di = model.get_state_flat_index((Eigen::Index)mio::isecir::InfectionState::Dead, group);
    EXPECT_NEAR(model.m_populations.get_value(0)[Di], deaths[size_t(group)], 1e-4);
    EXPECT_NEAR(model.m_total_confirmed_cases[group], total_confirmed_cases[size_t(group)], 1e-4);

    // Compare transitions at last time point with results from a previous run that are given here.
    Eigen::VectorX<ScalarType> compare(num_transitions * num_agegroups);
    compare << 30.00000000, 25.50000000, 0.25000000, 6.50000000, 0.89673307, 1.18750000, 0.80038452, 0.80038452,
        0.29817594, 0.29817594;
    mio::isecir::Simulation sim(model, dt);
    ASSERT_EQ(compare.size(), model.m_transitions.get_last_value().size());

    for (int j = 0; j < compare.size(); j++) {
        ASSERT_NEAR(compare[j], model.m_transitions.get_last_value()[j], 1e-7);
    }
}

// Check that initialization based on synthetic RKI data matches previous result for an IDE model with multiple
// age groups.
TEST(TestIDEParametersIo, RKIcompareWithPreviousRunAgeRes)
{
    // Define number of age groups.
    size_t num_agegroups = 6;
    // Define start date and the total population used for the initialization.
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 15 * 1e6);
    auto start_date = mio::Date(2020, 11, 1);
    ScalarType dt   = 0.5;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;
    // Initialize model.
    // The number of deaths will be overwritten if real data is used for initialization. Therefore, an arbitrary number
    // is used for the number of deaths.
    // Initial time series for the flows will be also overridden.
    mio::isecir::Model model(
        mio::TimeSeries<ScalarType>(-1, mio::TimeSeries<ScalarType>::Vector::Constant(
                                            (int)mio::isecir::InfectionTransition::Count * num_agegroups, 1.)),
        total_population, mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), -1.),
        num_agegroups);

    // Set the model parameters so that if the default values are changed, the test is still valid.
    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib((int)mio::isecir::InfectionTransition::Count,
                                                               delaydistribution);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]
        .set_distribution_parameter(1.7);
    for (mio::AgeGroup group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        model.parameters.get<mio::isecir::TransitionDistributions>()[group] = vec_delaydistrib;
    }
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        model.parameters.get<mio::isecir::TransitionProbabilities>()[group] = vec_prob;
    }
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, num_agegroups);
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_agegroups, num_agegroups, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constfunc(1.0);
    mio::StateAgeFunctionWrapper prob(constfunc);

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        model.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[group] = prob;
        model.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[group]   = prob;
        model.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[group]   = prob;
    }

    mio::CustomIndexArray<ScalarType, mio::AgeGroup> scale_confirmed_cases =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 1.);

    // Calculate initialization.
    auto status = mio::isecir::set_initial_flows<mio::ConfirmedCasesDataEntry>(
        model, dt, mio::read_confirmed_cases_data(mio::path_join(TEST_DATA_DIR, "cases_all_age_ma7.json")).value(),
        start_date, scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsSuccess());

    // Compare with previous run.

    std::vector<ScalarType> deaths = {
        1, 2.471428571455, 26.34999999999, 603.621428571465, 3972.41428571431, 7668.84999999998};

    std::vector<ScalarType> total_confirmed_cases = {10269.2857142857, 29615.8571428571, 185321.571428571,
                                                     215386.428571429, 77163.5714285714, 35588.4285714286};

    for (mio::AgeGroup group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        int Di = model.get_state_flat_index((Eigen::Index)mio::isecir::InfectionState::Dead, group);
        EXPECT_NEAR(model.m_populations.get_value(0)[Di], deaths[size_t(group)], 1e-4);
        EXPECT_NEAR(model.m_total_confirmed_cases[group], total_confirmed_cases[size_t(group)], 1e-4);
    }

    // Compare transitions at last time point with results from a previous run that are given here.
    Eigen::VectorX<ScalarType> compare(num_transitions * num_agegroups);
    compare << 336.428571428600, 328.285714285701, 162.000000000000, 163.071428571425, 80.130989648839, 79.803571428575,
        39.476374533415, 39.476374533415, 19.550404043081, 19.550404043081, 1105.714285714297, 1069.857142857200,
        515.714285714250, 163.071428571425, 80.130989648839, 79.803571428575, 39.476374533415, 39.476374533415,
        19.550404043081, 19.550404043081, 5819.000000000000, 5744.000000000000, 2806.428571428551, 163.071428571425,
        80.130989648839, 79.803571428575, 39.476374533415, 39.476374533415, 19.550404043081, 19.550404043081,
        6685.142857142899, 6572.857142857101, 3200.714285714304, 163.071428571425, 80.130989648839, 79.803571428575,
        39.476374533415, 39.476374533415, 19.550404043081, 19.550404043081, 2376.000000000000, 2342.285714285696,
        1142.571428571406, 163.071428571425, 80.130989648839, 79.803571428575, 39.476374533415, 39.476374533415,
        19.550404043081, 19.550404043081, 966.714285714304, 946.142857142797, 457.214285714301, 163.071428571425,
        80.130989648839, 79.803571428575, 39.476374533415, 39.476374533415, 19.550404043081, 19.550404043081;

    mio::isecir::Simulation sim(model, dt);
    ASSERT_EQ(compare.size(), model.m_transitions.get_last_value().size());

    for (int j = 0; j < compare.size(); j++) {
        ASSERT_NEAR(compare[j], model.m_transitions.get_last_value()[j], 1e-7);
    }
}
// Check some cases where computation of initial values for an IDE model without age resolution based on RKI data
// should fail.
TEST(TestIDEParametersIo, ParametersIoRKIFailure)
{
    // Define number of age groups.
    size_t num_agegroups = 1;
    // Since we only have one age group, the corresponding index is set to 0.
    mio::AgeGroup group = 0;

    // Define start date and the total population used for the initialization.
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 15 * 1e6);
    ScalarType dt = 0.5;

    // Initialize model.
    // The number of deaths will be overwritten if real data is used for initialization. Therefore, an arbitrary number
    // is used for the number of deaths.
    mio::isecir::Model model(
        mio::TimeSeries<ScalarType>((size_t)mio::isecir::InfectionTransition::Count * num_agegroups), total_population,
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), -1.), num_agegroups);

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // --- Case with empty RKI data file.
    auto start_date = mio::Date(2020, 06, 02);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> scale_confirmed_cases =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 1.);
    auto status = mio::isecir::set_initial_flows<mio::ConfirmedCasesNoAgeEntry>(
        model, dt, mio::read_confirmed_cases_noage(mio::path_join(TEST_DATA_DIR, "test_empty_file.json")).value(),
        start_date, scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsFailure(mio::StatusCode::InvalidFileFormat));

    // --- Case where start_date is later than maximal provided date in file.
    start_date = mio::Date(2021, 06, 8);
    std::vector<mio::ConfirmedCasesNoAgeEntry> test_data =
        mio::read_confirmed_cases_noage(mio::path_join(TEST_DATA_DIR, "cases_all_germany.json")).value();
    status = mio::isecir::set_initial_flows<mio::ConfirmedCasesNoAgeEntry>(model, dt, test_data, start_date,
                                                                           scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsFailure(mio::StatusCode::OutOfRange));

    // --- Case where not all needed dates from the future are provided.
    start_date = mio::Date(2020, 06, 07);
    status     = mio::isecir::set_initial_flows<mio::ConfirmedCasesNoAgeEntry>(model, dt, test_data, start_date,
                                                                               scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsFailure(mio::StatusCode::OutOfRange));

    // --- Case where not all needed dates from the past are provided.
    start_date = mio::Date(2020, 05, 24);
    status     = mio::isecir::set_initial_flows<mio::ConfirmedCasesNoAgeEntry>(model, dt, test_data, start_date,
                                                                               scale_confirmed_cases);
    // Check that status is Success as just a warning is logged.
    ASSERT_THAT(print_wrap(status), IsSuccess());
    // Check that the flow InfectedNoSymptomsToInfectedSymptoms has actually been set to 0.
    // As start_date is the first date in the data file, there is some infection data for all times >-1
    // (as we interpolate linearly in between -1 and 0). Therefore, we only expect the flow
    // InfectedNoSymptomsToInfectedSymptoms to be zero for times <=-1.
    int INStISy = model.get_transition_flat_index(
        (Eigen::Index)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms, group);
    for (Eigen::Index i = 0; i < model.m_transitions.get_num_time_points() - 2; i++) {
        EXPECT_EQ(0., model.m_transitions.get_value(i)[INStISy]);
    }

    // Reactivate log output.
    mio::set_log_level(mio::LogLevel::warn);

    // --- Valid case.
    start_date = mio::Date(2020, 06, 02);
    status     = mio::isecir::set_initial_flows<mio::ConfirmedCasesNoAgeEntry>(model, dt, test_data, start_date,
                                                                               scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsSuccess());
}

// Check some cases where computation of initial values for an IDE model with multiple age groups based on RKI data
// should fail.
TEST(TestIDEParametersIo, ParametersIoRKIFailureAgeRes)
{
    // Define number of age groups.
    size_t num_agegroups = 6;
    // Define start date and the total population used for the initialization.
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 15 * 1e6);
    ScalarType dt = 0.5;

    // Initialize model.
    // The number of deaths will be overwritten if real data is used for initialization. Therefore, an arbitrary number
    // is used for the number of deaths.
    mio::isecir::Model model(
        mio::TimeSeries<ScalarType>((size_t)mio::isecir::InfectionTransition::Count * num_agegroups), total_population,
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), -1.), num_agegroups);

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // --- Case with empty RKI data file.
    auto start_date = mio::Date(2020, 11, 01);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> scale_confirmed_cases =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 1.);
    auto status = mio::isecir::set_initial_flows<mio::ConfirmedCasesDataEntry>(
        model, dt, mio::read_confirmed_cases_data(mio::path_join(TEST_DATA_DIR, "test_empty_file.json")).value(),
        start_date, scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsFailure(mio::StatusCode::InvalidFileFormat));

    // --- Case where start_date is later than maximal provided date in file.
    start_date = mio::Date(2021, 01, 05);
    std::vector<mio::ConfirmedCasesDataEntry> test_data =
        mio::read_confirmed_cases_data(mio::path_join(TEST_DATA_DIR, "cases_all_age_ma7.json")).value();
    status = mio::isecir::set_initial_flows<mio::ConfirmedCasesDataEntry>(model, dt, test_data, start_date,
                                                                          scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsFailure(mio::StatusCode::OutOfRange));

    // --- Case where not all needed dates from the future are provided.
    start_date = mio::Date(2020, 12, 31);
    status     = mio::isecir::set_initial_flows<mio::ConfirmedCasesDataEntry>(model, dt, test_data, start_date,
                                                                              scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsFailure(mio::StatusCode::OutOfRange));

    // --- Case where not all needed dates from the past are provided.
    start_date = mio::Date(2020, 10, 1);
    status     = mio::isecir::set_initial_flows<mio::ConfirmedCasesDataEntry>(model, dt, test_data, start_date,
                                                                              scale_confirmed_cases);
    // Check that status is Success as just a warning is logged.
    ASSERT_THAT(print_wrap(status), IsSuccess());
    // Check that the flow InfectedNoSymptomsToInfectedSymptoms has actually been set to 0.
    // As start_date is the first date in the data file, there is some infection data for all times >-1
    // (as we interpolate linearly in between -1 and 0). Therefore, we only expect the flow
    // InfectedNoSymptomsToInfectedSymptoms to be zero for times <=-1.
    for (size_t group = 0; group < num_agegroups; ++group) {
        int INStISy = model.get_transition_flat_index(
            (Eigen::Index)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms, group);
        for (Eigen::Index i = 0; i < model.m_transitions.get_num_time_points() - 2; i++) {
            EXPECT_EQ(0., model.m_transitions.get_value(i)[INStISy]);
        }
    }

    // Reactivate log output.
    mio::set_log_level(mio::LogLevel::warn);

    // --- Valid case.
    start_date = mio::Date(2020, 11, 01);
    status     = mio::isecir::set_initial_flows<mio::ConfirmedCasesDataEntry>(model, dt, test_data, start_date,
                                                                              scale_confirmed_cases);

    ASSERT_THAT(print_wrap(status), IsSuccess());
}
