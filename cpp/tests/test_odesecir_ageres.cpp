/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "load_test_data.h"
#include "ode_secir/model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/math/stepper_wrapper.h"
#include <gtest/gtest.h>

TEST(TestSecir, compareAgeResWithSingleRun)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model model(3);
    mio::AgeGroup nb_groups = model.parameters.get_num_groups();
    double fact             = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::StartDay>(60);
    params.set<mio::osecir::Seasonality>(0.2);
    params.get<mio::osecir::TestAndTraceCapacity>() = 35;

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.8;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 9.5;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 7.1;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i]  = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]    = 0.7;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]    = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]    = 0.25;
        params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[i] = 0.45;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]         = 0.2;
        params.get<mio::osecir::CriticalPerSevere>()[i]                 = 0.25;
        params.get<mio::osecir::DeathsPerCritical>()[i]                 = 0.3;
    }

    params.apply_constraints();

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    auto integrator = std::make_shared<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    mio::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model, integrator);

    auto compare = load_test_data_csv<double>("secihurd-compare.csv");

    ASSERT_EQ(compare.size(), static_cast<size_t>(secihurd.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size() - 1, secihurd.get_num_elements() / (size_t)nb_groups) << "at row " << i;
        ASSERT_NEAR(secihurd.get_time(i), compare[i][0], 1e-10) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            double dummy = 0;
            for (size_t k = 0; k < (size_t)nb_groups; k++) {
                dummy += secihurd.get_value(i)[j - 1 + k * (size_t)mio::osecir::InfectionState::Count];
            }
            EXPECT_NEAR(dummy, compare[i][j], 1e-10) << " at row " << i;
        }
    }
}

TEST(TestSecir, compareAgeResWithSingleRunCashKarp)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model model(3);
    mio::AgeGroup nb_groups = model.parameters.get_num_groups();
    double fact             = 1.0 / (double)(size_t)nb_groups;

    model.parameters.set<mio::osecir::StartDay>(60);
    model.parameters.set<mio::osecir::Seasonality>(0.2);

    auto& params = model.parameters;
    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.8;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 9.5;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 7.1;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]                     = fact * nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}]          = fact * nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]            = fact * nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]              = fact * nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]            = fact * nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]                   = fact * nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]                        = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i]  = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]    = 0.7;
        params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[i] = 0.45;
        params.get<mio::osecir::TestAndTraceCapacity>()                 = 35;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]    = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]    = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]         = 0.2;
        params.get<mio::osecir::CriticalPerSevere>()[i]                 = 0.25;
        params.get<mio::osecir::DeathsPerCritical>()[i]                 = 0.3;
    }

    params.apply_constraints();

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    mio::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model, integrator);

    auto compare = load_test_data_csv<double>("secihurd-compare-cashkarp.csv");

    ASSERT_EQ(compare.size(), static_cast<size_t>(secihurd.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size() - 1, secihurd.get_num_elements() / (size_t)nb_groups) << "at row " << i;
        ASSERT_NEAR(secihurd.get_time(i), compare[i][0], 1e-10) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            double dummy = 0;
            for (size_t k = 0; k < (size_t)nb_groups; k++) {
                dummy += secihurd.get_value(i)[j - 1 + k * (size_t)mio::osecir::InfectionState::Count];
            }
            EXPECT_NEAR(dummy, compare[i][j], 1e-10) << " at row " << i;
        }
    }
}
