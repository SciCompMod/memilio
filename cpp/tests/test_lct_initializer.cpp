/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "lct_secir/infection_state.h"
#include "lct_secir/initialization.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/math/eigen.h"

#include <gtest/gtest.h>

TEST(TestInitializer, compareWithPrevious)
{
    ScalarType dt = 0.5;
    using Vec     = mio::TimeSeries<ScalarType>::Vector;

    std::vector<int> SubcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Exposed]            = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 3;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = 3;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 2;
    mio::lsecir::InfectionState InfState(SubcompartmentNumbers);

    // Previous result.
    Eigen::VectorXd compare(InfState.get_count());
    compare << 82810889.00545, 850.70432, 970.04980, 315.32890, 391.51799, 391.39351, 565.45854, 580.79267, 85.97421,
        86.02738, 80.26791, 189.53449, 167.57963, 329757.36512, 9710;

    mio::lsecir::Parameters parameters_lct;
    parameters_lct.get<mio::lsecir::TimeExposed>()                      = 3.1;
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms>()           = 3.1;
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms>()             = 6.1;
    parameters_lct.get<mio::lsecir::TimeInfectedSevere>()               = 11.1;
    parameters_lct.get<mio::lsecir::TimeInfectedCritical>()             = 17.1;
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.01;
    mio::ContactMatrixGroup contact_matrix                              = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0].finalize();
    parameters_lct.get<mio::lsecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 1;
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 1;
    parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.1;
    parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.1;
    parameters_lct.get<mio::lsecir::CriticalPerSevere>()              = 0.1;
    parameters_lct.get<mio::lsecir::DeathsPerCritical>()              = 0.1;

    ScalarType total_confirmed_cases = 341223;
    ScalarType deaths                = 9710;
    ScalarType total_population      = 83155031.0;

    // add time points for initialization of transitions
    mio::TimeSeries<ScalarType> init((int)mio::lsecir::InfectionTransition::Count);
    Vec vec_init((int)mio::lsecir::InfectionTransition::Count);
    vec_init[(int)mio::lsecir::InfectionTransition::SusceptibleToExposed]                 = 25.0;
    vec_init[(int)mio::lsecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
    vec_init[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    vec_init[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
    vec_init[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 1.0;
    vec_init[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
    vec_init[(int)mio::lsecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 1.0;
    vec_init[(int)mio::lsecir::InfectionTransition::InfectedSevereToRecovered]            = 1.0;
    vec_init[(int)mio::lsecir::InfectionTransition::InfectedCriticalToDead]               = 1.0;
    vec_init[(int)mio::lsecir::InfectionTransition::InfectedCriticalToRecovered]          = 1.0;
    // add initial time point to time series
    init.add_time_point(-130, vec_init);
    // add further time points until time 0
    while (init.get_last_time() < 0) {
        vec_init *= 1.01;
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    mio::lsecir::Initializer initializer(std::move(init), InfState, std::move(parameters_lct));
    auto init_compartments = initializer.compute_initializationvector(total_population, deaths, total_confirmed_cases);

    for (int i = 0; i < InfState.get_count(); i++) {
        EXPECT_NEAR(init_compartments[i], compare[i], 1e-4) << "at subcompartment number " << i;
    }
}

TEST(TestInitializer, testConstraints)
{
    // Check constraints of Initializer.
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    ScalarType dt = 0.5;

    std::vector<int> SubcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Exposed]            = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 3;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = 3;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 2;
    mio::lsecir::InfectionState InfState(SubcompartmentNumbers);

    // Check wrong form of initial flows.
    mio::TimeSeries<ScalarType> init_wrong_size((int)mio::lsecir::InfectionTransition::Count - 1);
    Eigen::VectorXd vec_wrong_size = Eigen::VectorXd::Ones((int)mio::lsecir::InfectionTransition::Count - 1);
    init_wrong_size.add_time_point(-10, vec_wrong_size);
    init_wrong_size.add_time_point(-9, vec_wrong_size);

    mio::lsecir::Initializer initializer_init_wrong_size(std::move(init_wrong_size), InfState);

    bool constraint_check = initializer_init_wrong_size.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check if last time of initial flows is not zero.
    mio::TimeSeries<ScalarType> init_wrong((int)mio::lsecir::InfectionTransition::Count);
    Eigen::VectorXd vec_init = Eigen::VectorXd::Ones((int)mio::lsecir::InfectionTransition::Count);
    init_wrong.add_time_point(-10, vec_init);
    init_wrong.add_time_point(-9, vec_init);

    mio::TimeSeries<ScalarType> init_copy(init_wrong);
    mio::lsecir::Initializer initializer_init_wrong_last_time(std::move(init_copy), InfState);

    constraint_check = initializer_init_wrong_last_time.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check if time steps of initial flows are not equidistant.
    init_wrong.add_time_point(init_wrong.get_last_time() + 1., vec_init);
    while (init_wrong.get_last_time() < 0) {
        init_wrong.add_time_point(init_wrong.get_last_time() + dt, vec_init);
    }

    mio::lsecir::Initializer initializer_init_wrong_equidistant(std::move(init_wrong), InfState);

    constraint_check = initializer_init_wrong_equidistant.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check large step size.
    mio::TimeSeries<ScalarType> init_wrong_step((int)mio::lsecir::InfectionTransition::Count);
    init_wrong_step.add_time_point(-10, vec_init);
    init_wrong_step.add_time_point(init_wrong_step.get_last_time() + 2., vec_init);
    while (init_wrong_step.get_last_time() < 0) {
        init_wrong_step.add_time_point(init_wrong_step.get_last_time() + dt, vec_init);
    }

    mio::lsecir::Initializer initializer_init_wrong_step(std::move(init_wrong_step), InfState);

    constraint_check = initializer_init_wrong_step.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check with correct flows.
    mio::TimeSeries<ScalarType> init_right((int)mio::lsecir::InfectionTransition::Count);
    init_right.add_time_point(-10, vec_init);
    while (init_right.get_last_time() < 0) {
        init_right.add_time_point(init_right.get_last_time() + dt, vec_init);
    }

    mio::lsecir::Initializer initializer_right(std::move(init_right), InfState);

    constraint_check = initializer_right.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Check with short time frame. Fitting time frame is tested above.
    mio::TimeSeries<ScalarType> init_short((int)mio::lsecir::InfectionTransition::Count);
    init_short.add_time_point(-1., vec_init);
    while (init_short.get_last_time() < 0) {
        init_short.add_time_point(init_short.get_last_time() + dt, vec_init);
    }

    mio::lsecir::Initializer initializer_init_short(std::move(init_short), InfState);
    ScalarType total_confirmed_cases = 341223;
    ScalarType deaths                = 9710;
    ScalarType total_population      = 83155031.0;
    auto initialconditions =
        initializer_init_short.compute_initializationvector(total_population, deaths, total_confirmed_cases);

    for (int i = 2; i < InfState.get_count() - 2; i++) {
        EXPECT_EQ(-1, initialconditions[i]);
    }

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}
/* test num_time_points < calc_time_index and constraints*/